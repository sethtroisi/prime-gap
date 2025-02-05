// Copyright 2020 Seth Troisi
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <algorithm>
#include <cassert>
#include <chrono>
#include <clocale>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <functional>
#include <iostream>
#include <map>
#include <mutex>
#include <sstream>
#include <thread>
#include <type_traits>
#include <vector>

#include <gmp.h>
#include <omp.h>
#include <primesieve.hpp>

#include "gap_common.h"
#include "modulo_search.h"

using std::cout;
using std::endl;
using std::map;
using std::mutex;
using std::pair;
using std::vector;
using namespace std::chrono;


/**
 * Two MACROS used to validate results
 * GMP_VALIDATE_FACTORS (validates all factors)
 * GMP_VALIDATE_LARGE_FACTORS (validate large factors)
 *
 * GMP_VALIDATE_LARGE_FACTORS only validates the rarer 60+ bit factors
 */
//#define GMP_VALIDATE_FACTORS
#define GMP_VALIDATE_LARGE_FACTORS

// Compresses composite by 50-80%,
// Might make large prime faster but never makes sense because
// Of increased memory size (and time for count unknows)
#define METHOD2_WHEEL   1

// This probably should be optimized to fit in L2/L3
// Related to sizeof(int) * SIEVE_INTERVAL * WHEEL_MAX
// WHEEL should divide config.d
#define METHOD2_WHEEL_MAX (2*3*5*7)

void set_defaults(struct Config& config);
void prime_gap_search(const struct Config& config);
void prime_gap_parallel(struct Config& config);


int main(int argc, char* argv[]) {
    // Display %'d with commas i.e. 12,345
    setlocale(LC_NUMERIC, "");

    Config config = Args::argparse(argc, argv, Args::Pr::SIEVE);

    if (config.verbose >= 2) {
        printf("\tCompiled with GMP %d.%d.%d\n\n",
            __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);
    }

    // More combined sieve specific validation
    {
        set_defaults(config);

        // Both shouldn't be true from gap_common.
        assert(!(config.save_unknowns && config.testing));

        if (!config.save_unknowns && !config.testing) {
            cout << "Must set --save-unknowns" << endl;
            exit(1);
        }

        if (config.sieve_length < 6 * config.p || config.sieve_length > 22 * config.p) {
            int sl_min = ((config.p * 8 - 1) / 500 + 1) * 500;
            int sl_max = ((config.p * 20 - 1) / 500 + 1) * 500;
            printf("--sieve_length(%d) should be between [%d, %d]\n",
                config.sieve_length, sl_min, sl_max);
            exit(1);
        }

        if (config.valid == 0) {
            Args::show_usage(argv[0], Args::Pr::SIEVE);
            exit(1);
        }

        if (config.max_prime > 500'000'000) {
            float m_per = config.max_prime / ((float) config.minc * config.sieve_length);
            if (m_per < .1 && config.p <= 8000) {
                printf("\tmax_prime(%ldB) is probably too large\n",
                    config.max_prime / 1'000'000'000);
            }
        }

        if (config.save_unknowns) {
            std::string fn = Args::gen_unknown_fn(config, ".txt");
            std::ifstream f(fn);
            if (f.good()) {
                printf("\nOutput file '%s' already exists\n", fn.c_str());
                exit(1);
            }
        }

        if (!config.compression && (config.minc * config.sieve_length) > 30'000'000'000L) {
            printf("\tSetting --rle to prevent very large output file\n");
            config.compression = 1;
        }

        if (config.compression == 2) {
            cout << "compression==2 no longer supported" << endl;
            exit(1);
        }
    }

    // Status lines
    if (config.verbose >= 0) {
        printf("\n");
        printf("Testing m * %u#/%u, m = %'ld + [0, %'ld)\n",
            config.p, config.d, config.mstart, config.minc);
    }

    if (config.verbose >= 2 && config.threads > 1) {
        printf("Running with %d threads\n", config.threads);
    }

#ifdef GMP_VALIDATE_FACTORS
    printf("\tValidating factors with GMP\n");
#endif

    if (config.method1) {
        prime_gap_search(config);
    } else {
        prime_gap_parallel(config);
    }
}


void set_defaults(struct Config& config) {
    if (config.valid == 0) {
        // Don't do anything if argparse didn't work.
        return;
    }

    if (config.d % 4 == 0) {
        // AKA min-merit
        config.sieve_length = config.p * config.min_merit;

        // Start from 1
        config.mstart = 1;

        // Large prime near P to make D unique (chosen semi-randomly)
        config.d /= 4;
        vector<uint32_t> P_primes = get_sieve_primes(config.p);
        uint32_t rand_prime = P_primes[P_primes.size() - 2 - (rand() % 10)];
        uint32_t large_p = config.d > 1 ? config.d : rand_prime;
        assert(is_prime_brute(large_p));

        printf("d optimizer for P = %d# | large prime=%d | SL=%d (%.1f merit)\n",
                config.p, large_p, config.sieve_length, config.min_merit);

        /**
         * Secret value to optimize d
         * 1. Test small primorials to find optimal primorial
         * 2. Multiple by large prime (to make unique)
         * 3. test that ~same expected
         */
        vector<uint32_t> primes = {1,2,3,5,7,11,13,17,19,23};
        for (uint32_t lp : {1u, large_p}) {
            config.d = lp;
            for (uint32_t p : primes) {
                // check if large_p already includes p
                if (p != 1 && config.d % p == 0)
                    continue;

                if (__builtin_umul_overflow(config.d, p, &config.d)) {
                    // overflow
                    break;
                }

                // Try searching all values of m (up to 20,000)
                config.minc = std::min(config.d, 20'000U);
                auto expected = count_K_d(config);
                printf("Optimizing | d = %5d * %2d# | %d remaining, %5.0f avg gap | SL insufficient %.3f%% of time\n",
                    lp, p, std::get<1>(expected), std::get<0>(expected), 100 * std::get<2>(expected));
            }
        }
        exit(0);
    }

    mpz_t K;
    double K_log;
    {
        // Suppress log
        int temp = config.verbose;
        config.verbose = -1;

        int K_digits;
        K_stats(config, K, &K_digits, &K_log);

        config.verbose = temp;
    }

    if (config.sieve_length == 0) {
        // Change that a number near K is prime
        // GIVEN no factor of K or D => no factor of P#
        const double prob_prime_coprime_P = prob_prime_coprime(config);

        // factors of K = P#/D
        vector<uint32_t> K_primes = get_sieve_primes(config.p);
        // Remove any factors of D
        K_primes.erase(
            std::remove_if(K_primes.begin(), K_primes.end(),
               [&](uint32_t p){ return config.d % p == 0; }),
            K_primes.end());

        // K = #P/D
        // only numbers K+i has no factor <= p
        //      => (K+i, i) == (K, i) == 1
        //      => only relatively prime i's
        //
        // factors of d are hard because they depend on m*K
        //      some of these m are worse than others so use worst m

        assert( config.p >= 503 );

        // Search till chance of shorter gap is small.
        {
            // Code below is quite slow with larger values of d.
            assert( config.d <= 30030 );

            uint32_t base = mpz_fdiv_ui(K, config.d);

            // count of (m*K) % d over all m
            vector<uint32_t> count_by_mod_d(config.d, 0);
            {
                for (uint64_t mi = 0; mi < config.minc; mi++) {
                    uint64_t m = config.mstart + mi;
                    if (gcd(m, config.d) == 1) {
                        uint32_t center = ((__int128) m * base) % config.d;
                        uint32_t center_down = (config.d - center) % config.d;

                        // distance heading up
                        count_by_mod_d[ center ] += 1;
                        // distance heading up
                        count_by_mod_d[ center_down ] += 1;
                    }
                }
            }

            // Note: By averaging over counts prob_larger could be improve here.
            map<uint32_t, uint32_t> coprime_by_mod_d;
            for (size_t i = 0; i < config.d; i++) {
                if (count_by_mod_d[i] > 0) {
                    coprime_by_mod_d[i] = 0;
                }
            }

            // Keep increasing SL till prob_gap_shorter < 0.8%
            for (size_t tSL = 1; ; tSL += 1) {
                bool any_divisible = false;
                for (int prime : K_primes) {
                    if ((tSL % prime) == 0) {
                        any_divisible = true;
                        break;
                    }
                }
                // Result will be the same as last.
                if (any_divisible) continue;

                // check if tSL is divisible for all center mods
                for (auto& coprime_counts : coprime_by_mod_d) {
                    const auto center = coprime_counts.first;
                    // Some multiple of d will mark this off (for these centers) don't count it.
                    if (gcd(center + tSL, config.d) == 1) {
                        coprime_counts.second += 1;
                    }
                }

                // Find the smallest number of coprimes
                uint32_t min_coprime = tSL;
                for (auto& coprime_counts : coprime_by_mod_d) {
                    min_coprime = std::min(min_coprime, coprime_counts.second);
                }

                // Assume each coprime is independent
                double prob_gap_shorter = pow(1 - prob_prime_coprime_P, min_coprime);

                // This seems to balance PRP fallback and sieve_size
                if (prob_gap_shorter <= 0.008) {
                    config.sieve_length = tSL;
                    printf("AUTO SET: sieve length: %ld (coprime: %d, prob_gap longer %.2f%%)\n",
                        tSL, min_coprime, 100 * prob_gap_shorter);
                    break;
                }
            }
        }
        assert( config.sieve_length > 100 ); // Something went wrong above.
    }

    if (config.max_prime == 0) {
        // each additional numbers removes unknowns / prime
        // and takes log2(prime / sieve_length) time

        // Not worth improving given method2 CTRL+C handling.
        if (K_log >= 1500) {
            config.max_prime =  100'000'000'000;
        } else {
            config.max_prime =   10'000'000'000;
        }
        if (config.method1) {
            printf("Can't use method1 and not set max_prime");
            exit(1);
        }
        if (config.verbose >= 0) {
            printf("AUTO SET: max_prime (log(K) = ~%.0f): %ld\n",
                K_log, config.max_prime);
            printf("WATCH for 'Estimated 2x faster (CTRL+C to stop sieving)' warning");
        }
    }

    mpz_clear(K);
}


static void insert_range_db(
        const struct Config& config,
        long num_rows,
        float time_sieve) {

    DB db_helper(config.search_db.c_str());
    sqlite3 *db = db_helper.get_db();

    const uint64_t rid = db_helper.config_hash(config);
    char sSQL[300];
    sprintf(sSQL,
        "INSERT INTO range(rid, P,D, m_start,m_inc,"
                          "sieve_length, max_prime,"
                          "min_merit,"
                          "num_m,"
                          "time_sieve)"
         "VALUES(%ld,  %d,%d, %ld,%ld,"
                "%d,%ld, %.3f,"
                "%ld,  %.2f)"
         "ON CONFLICT(rid) DO UPDATE SET time_sieve=%.2f",
            rid,  config.p, config.d, config.mstart, config.minc,
            config.sieve_length, config.max_prime,
            config.min_merit,
            num_rows,
            time_sieve, time_sieve);

    char *zErrMsg = nullptr;
    int rc = sqlite3_exec(db, sSQL, nullptr, nullptr, &zErrMsg);
    if (rc != SQLITE_OK) {
        printf("\nrange INSERT failed %d: %s\n",
            rc, sqlite3_errmsg(db));
        exit(1);
    }
}


// Method1


void save_unknowns_method1(
        std::ofstream &unknown_file,
        const uint64_t m, int unknown_p, int unknown_n,
        const unsigned int SL, const vector<char> composite[]) {

    unknown_file << m << " : -" << unknown_p << " +" << unknown_n << " |";

    for (int d = 0; d <= 1; d++) {
        char prefix = "-+"[d];

        for (size_t i = 1; i <= SL; i++) {
            if (!composite[d][i]) {
                unknown_file << " " << prefix << i;
            }
        }
        if (d == 0) {
            unknown_file << " |";
        }
    }
    unknown_file << "\n";
}


void prime_gap_search(const struct Config& config) {
    cout << "METHOD1 NOT IMPLEMENTED" << endl;
    exit(1);
}


// Method 2

bool g_control_c = false;
void signal_callback_handler(int) {
    if (g_control_c) {
        cout << "Caught 2nd CTRL+C stopping now." << endl;
        exit(2);
    } else {
       cout << "Caught CTRL+C stopping and saving after next interval " << endl;
       g_control_c = true;
    }
}


class method2_stats {
    public:
        method2_stats() {};

        method2_stats(
                int thread_i,
                const struct Config& config,
                size_t valid_ms,
                uint64_t threshold,
                double initial_prob_prime
        ) {
            thread = thread_i;

            start_t = high_resolution_clock::now();
            interval_t = high_resolution_clock::now();
            total_unknowns = config.sieve_length * valid_ms;

            if (threshold <= 100000)
               next_mult = 10000;

            prob_prime = initial_prob_prime;
            current_prob_prime = prob_prime;
        }

        // Some prints only happen if thread == 0
        int thread = 0;
        uint64_t next_print = 0;
        uint64_t next_mult = 100000;

        // global and interval start times
        high_resolution_clock::time_point  start_t;
        high_resolution_clock::time_point  interval_t;

        long total_unknowns = 0;
        long small_prime_factors_interval = 0;
        long large_prime_factors_interval = 0;
        // Sum of above two, mostly handled in method2_increment_print
        long prime_factors = 0;

        size_t pi = 0;
        size_t pi_interval = 0;

        uint64_t m_stops = 0;
        uint64_t m_stops_interval = 0;

        uint64_t validated_factors = 0;

        // prob prime after sieve up to some prime threshold
        double current_prob_prime = 0;

        // Constants (more of a stats storage)
        double prp_time_estimate = std::nan("");
        double prob_prime = 0;
        uint64_t last_prime = 0;

        size_t count_coprime_p = 0;

        // Used as a sentinal in method2_large_primes
        bool finished = 0;
};

void method2_increment_print(
        uint64_t prime,
        size_t valid_ms,
        vector<bool> *composite,
        method2_stats &stats,
        const struct Config& config) {

    /**
     * verification requires count_coprime_to_P#
     * Require that first call (next_print = 0) processes all primes up to P
     */
    if (stats.next_print == 0 && stats.count_coprime_p == 0) {
        assert(prime == config.p);

        if (stats.thread == 0) {
            // Other threads don't print details

            if (config.threads > 1 && config.verbose) {
                printf("\nWARNING stats aren't synchronized when "
                       "running with multiple threads(%d)\n\n", config.threads);
            }

            // This sligtly duplicates work below, but we don't care.
            auto   temp = high_resolution_clock::now();
            for (size_t i = 0; i < valid_ms; i++) {
                stats.count_coprime_p += std::count(composite[i].begin(), composite[i].end(), false);
            }
            double interval_count_time = duration<double>(high_resolution_clock::now() - temp).count();
            if (config.verbose >= 2) {
                printf("\t\t counting unknowns takes ~%.1f seconds\n", interval_count_time);
            }
        }
    }

    while (prime >= stats.next_print && stats.next_print < stats.last_prime) {
        //printf("\t\tmethod2_increment_print %'ld >= %'ld\n", prime, stats.next_print);
        const size_t max_mult = 100'000'000'000L * (config.threads > 2 ? 10L : 1L);


        // 10, 20, 30, 40, 50, 100, 200, 300, 400, 500, 1000 ...
        // 60, 70, 80, 90, 100, 120, 150, 200, 300 billion because intervals are wider.
        size_t extra_multiples = prime > ((config.threads > 4) ? 100'000'000 : 1'000'000);
        // With lots of threads small intervals are very fast
        // and large % of time is spent counting unknowns

        // Next time to increment the interval size
        size_t next_next_mult = (5 + 10 * extra_multiples) * stats.next_mult;
        if (stats.next_mult < max_mult && stats.next_print == next_next_mult) {
            stats.next_mult *= 10;
            stats.next_print = 0;
        }

        // 1,2,3,4,5,6,7,8,9,10, SKIP to 12, SKIP to 15
        stats.next_print += stats.next_mult;
        assert(stats.next_print % stats.next_mult == 0);

        if (stats.next_mult < max_mult) {
            int64_t ratio = stats.next_print / stats.next_mult;
            assert(ratio >= 1 && ratio <= 14);

            if (ratio > 10 && ratio < 12) {  // Skip 11 => 12
                stats.next_print = 12 * stats.next_mult;
            } else if (ratio > 12) {  // Skip 13, 14 => 15
                stats.next_print = 15 * stats.next_mult;
            }
        }

        // Never set next_print beyond last_prime
        stats.next_print = std::min(stats.next_print, stats.last_prime);
    }

    bool is_last = (prime == stats.last_prime) || g_control_c;

    if (config.verbose + is_last >= 1) {
        auto   s_stop_t = high_resolution_clock::now();
        // total time, interval time
        double     secs = duration<double>(s_stop_t - stats.start_t).count();
        double int_secs = duration<double>(s_stop_t - stats.interval_t).count();
        uint32_t SIEVE_LENGTH = config.sieve_length;

        if (stats.thread >= 1) {
            printf("Thread %d\t", stats.thread);
        }

        stats.pi += stats.pi_interval;

        printf("%'-10ld (primes %'ld/%ld)\t(seconds: %.2f/%-.1f | per m: %.3g)",
            prime,
            stats.pi_interval, stats.pi,
            int_secs, secs,
            secs / valid_ms);
        if (int_secs > 240) {
            // Add " @ HH:MM:SS" so that it is easier to predict when the next print will happen
            time_t rawtime = std::time(nullptr);
            struct tm *tm = localtime( &rawtime );
            printf(" @ %d:%02d:%02d", tm->tm_hour, tm->tm_min, tm->tm_sec);
        }
        printf("\n");
        stats.interval_t = s_stop_t;

        int verbose = config.verbose + (2 * is_last) + (prime > 1e9) + (stats.thread == 0);
        if (verbose >= 3) {
            stats.prime_factors += stats.small_prime_factors_interval;
            stats.prime_factors += stats.large_prime_factors_interval;
            stats.m_stops += stats.m_stops_interval;

            printf("\tfactors  %'14ld \t"
                   "(interval: %'ld avg m/large_prime interval: %.1f)\n",
                stats.prime_factors,
                stats.small_prime_factors_interval + stats.large_prime_factors_interval,
                1.0 * stats.m_stops_interval / stats.pi_interval);

            // See THEORY.md
            double prob_prime_after_sieve = stats.prob_prime * log(prime) * exp(GAMMA);
            double delta_sieve_prob = (1/stats.current_prob_prime - 1/prob_prime_after_sieve);
            double skipped_prp = 2 * valid_ms * delta_sieve_prob;

            if (is_last || config.threads <= 1) {
                uint64_t t_total_unknowns = 0;
                for (size_t i = 0; i < valid_ms; i++) {
                    t_total_unknowns += std::count(composite[i].begin(), composite[i].end(), false);
                }
                uint64_t new_composites = stats.total_unknowns - t_total_unknowns;

                // count_coprime_sieve * valid_ms also makes sense but leads to smaller numbers
                printf("\tunknowns %'9ld/%-5ld\t"
                       "(avg/m: %.2f) (composite: %.2f%% +%.3f%% +%'ld)\n",
                    t_total_unknowns, valid_ms,
                    1.0 * t_total_unknowns / valid_ms,
                    100.0 - 100.0 * t_total_unknowns / (SIEVE_LENGTH * valid_ms),
                    100.0 * new_composites / (SIEVE_LENGTH * valid_ms),
                    new_composites);

                if (stats.count_coprime_p && prime > 100000 && prime > config.p) {
                    // verify total unknowns & interval unknowns
                    const double prob_prime_coprime_P = prob_prime_coprime(config);

                    float e_unknowns = stats.count_coprime_p * (prob_prime_coprime_P / prob_prime_after_sieve);

                    float delta_composite_rate = delta_sieve_prob * prob_prime_coprime_P;
                    float e_new_composites = stats.count_coprime_p * delta_composite_rate;

                    float error = 100.0 * fabs(e_unknowns - t_total_unknowns) / e_unknowns;
                    float interval_error = 100.0 * fabs(e_new_composites - new_composites) / e_new_composites;

                    if (config.verbose >= 3 || error > 0.1 ) {
                        printf("\tEstimated %.3g unknowns found %.3g (%.2f%% error)\n",
                            e_unknowns, 1.0f * t_total_unknowns, error);
                    }
                    if (config.verbose >= 3 || interval_error > 0.3 ) {
                        printf("\tEstimated %.3g new composites found %.3g (%.2f%% error)\n",
                            e_new_composites, 1.0f * new_composites, interval_error);
                    }
                }
                stats.total_unknowns = t_total_unknowns;
            }

            stats.current_prob_prime = prob_prime_after_sieve;

            double prp_rate = skipped_prp / (int_secs * config.threads);
            if (config.show_timing) {
                printf("\t~ 2x %.2f PRP/m\t\t"
                       "(~ %4.1f skipped PRP => %.1f PRP/%s)\n",
                    1 / stats.current_prob_prime, skipped_prp,
                    prp_rate,
                    config.threads > 1 ? "thread-seconds" : "seconds");
            }
            if (stats.validated_factors) {
                printf("\tValidated %ld factors\n", stats.validated_factors);
            }

            double run_prp_mult = stats.prp_time_estimate / prp_rate;
            if (run_prp_mult > 0.25 && config.show_timing) {
                printf("\t\tEstimated ~%.1fx faster to just run PRP now (CTRL+C to stop sieving)\n",
                    run_prp_mult);
            }

            printf("\n");

            stats.small_prime_factors_interval = 0;
            stats.large_prime_factors_interval = 0;
        }

        stats.pi_interval = 0;
        stats.m_stops_interval = 0;
    }
}

void validate_factor_m_k_x(
        method2_stats& stats,
        mpz_t &test, const mpz_t &K, int64_t m, uint32_t X,
        uint64_t prime) {
#ifdef GMP_VALIDATE_FACTORS
    stats.validated_factors += 1;
    mpz_mul_ui(test, K, m);
    mpz_add_ui(test, test, X);
    uint64_t mod = mpz_fdiv_ui(test, prime);
    assert(mod == 0);
#endif  // GMP_VALIDATE_FACTORS
}


/**
 * TODO better name: RangeStats, KStats, Helpers, Indexes?
 *
 * Helper arrays
 */
class Cached {
    public:
        // mi such that gcd(m_start + mi, D) = 1
        vector<uint32_t> valid_mi;
        // valid_mi.size()
        size_t valid_ms;

        /**
         * m_reindex[mi] = composite[mii] if coprime
         * -1 if gcd(ms + i, D) > 1
         *
         * This is potentially very large use is_m_coprime and is_m_coprime2310
         * to pre check it's a coprime mi before doing the L3/RAM lookup.
         */
        vector<int32_t> m_reindex;
        // if gcd(ms + mi, D) = 1
        vector<bool> is_m_coprime;
        /**
         * is_m_coprime2310[i] = (i, D') == 1
         * D' has only factors <= 11
         * first 2310 values
         */
        vector<bool> is_m_coprime2310;


        // X which are coprime to K
        vector<uint32_t> coprime_X;
        // reindex composite[m][X] for composite[m_reindex[m]][x_reindex[X]]
        // Special 0'th entry stands for all not coprime
        vector<uint32_t> x_reindex;
        // if [x] is coprime to K
        vector<char> is_offset_coprime;


        // reindex composite[m][i] using (m, wheel) (wheel is 1!, 2!, 3!, or 5!)
        // This could be first indexed by x_reindex,
        // Would reduce size from wheel * (SL+1) to wheel * coprime_i

        // Note: Larger wheel eliminates more numbers but takes more space.
        // 6 (saves 2/3 memory), 30 (saves 11/15 memory)
        uint32_t x_reindex_wheel_size;
        vector<uint16_t> x_reindex_wheel[METHOD2_WHEEL_MAX];
        vector<size_t> x_reindex_wheel_count;


        int32_t K_mod2310;

        /** is_comprime2310[i] = (i % 2) && (i % 3) && (i % 5) && (i % 7) && (i % 11)*/
        vector<char> is_coprime2310;

        /**
         * TODO: benchmark adding is_coprime96577 = 13*17*19*23
         * 1 - 2/3 * 4/5 * 6/7 * 10/11 = 58% chance of sharing a factor of 3, 5, 7, or 11
         * 1 - 12/13 * 16/17 * 18/19 * 22/23 = 21% chance of sharing a factor of 13, 17, 19, or 23
         *
         * Could reduce memory contention (is this the slowness on DDR3?)
         */


    Cached(const struct Config& config, const mpz_t &K) {
        const uint32_t P = config.p;
        const uint32_t D = config.d;

        const uint32_t SL = config.sieve_length;

        const vector<uint32_t> P_primes = get_sieve_primes(P);
        assert( P_primes.back() == P);

        // Allocate temp vectors
        m_reindex.resize(config.minc, -1);
        {
            auto temp = is_coprime_and_valid_m(config);
            is_m_coprime = temp.first;
            valid_mi = temp.second;

            for (uint32_t mii = 0; mii < valid_mi.size(); mii++) {
                uint32_t mi = valid_mi[mii];
                m_reindex[mi] = mii;
            }
        }
        valid_ms = valid_mi.size();

        // Includes 0
        is_offset_coprime.resize(SL+1, 1);
        x_reindex.resize(SL+1, 0);

        // reindex composite[m][i] using (m, wheel) (wheel is 1!,2!,3!,5!)
        // This could be first indexed by x_reindex,
        // Would reduce size from wheel * SL to wheel * coprime_i
#if METHOD2_WHEEL
        // Note: Larger wheel eliminates more numbers but takes more space.
        // 6 seems reasonable for larger numbers  (uses 1/3 memory = 33%)
        // 30 is maybe better for smaller numbers (uses 4/15 memory = 26%)
        // 210 is maybe better (uses 24/105 = 23% memory) but x_reindex_wheel might not fit in memory.
        uint32_t wheel = gcd(D, METHOD2_WHEEL_MAX);
        uint32_t reindex_size = wheel * SL * sizeof(uint32_t);
        if (reindex_size > 7 * 1024 * 1024) {
            if (wheel % 7) {
                wheel /= 7;
            } else if (wheel % 5) {
                wheel /= 5;
            }
        }
        x_reindex_wheel_size = wheel;
        assert(x_reindex_wheel_size <= METHOD2_WHEEL_MAX); // Static size in Caches will break;
#else
        x_reindex_wheel_size = 1;
#endif  // METHOD2_WHEEL

        x_reindex_wheel_count.resize(x_reindex_wheel_size, 0);

        for (uint32_t prime : P_primes) {
            if (D % prime != 0) {
                for (size_t x = 0; x <= SL; x += prime) {
                    is_offset_coprime[x] = 0;
                }
            }
        }

        // Center should be marked composite by every prime.
        assert(is_offset_coprime[0] == 0);
        {
            size_t coprime_count = 0;
            for (size_t X = 0; X <= SL; X++) {
                if (is_offset_coprime[X] > 0) {
                    coprime_X.push_back(X);
                    coprime_count += 1;
                    x_reindex[X] = coprime_count;
                }
            }
            assert(coprime_count == coprime_X.size());
        }

        // Start at m_wheel == 0 so that re_index_m_wheel == 1 (D=1) works.

        for (size_t m_wheel = 0; m_wheel < x_reindex_wheel_size; m_wheel++) {
            if (gcd(x_reindex_wheel_size, m_wheel) > 1) continue;
            x_reindex_wheel[m_wheel].resize(SL+1, 0);

            // m * K % wheel => m_wheel % wheel
            uint32_t mod_center = m_wheel * mpz_fdiv_ui(K, x_reindex_wheel_size);

            size_t coprime_count_wheel = 0;
            for (size_t i = 0; i < SL+1; i++) {
                if (is_offset_coprime[i] > 0) {
                    if (gcd(mod_center + i, x_reindex_wheel_size) == 1) {
                        coprime_count_wheel += 1;
                        x_reindex_wheel[m_wheel][i] = coprime_count_wheel;
                    }
                }
            }
            x_reindex_wheel_count[m_wheel] = coprime_count_wheel;

            #if METHOD2_WHEEL
            {
                size_t x_reindex_limit = std::numeric_limits<
                    std::remove_extent_t<decltype(x_reindex_wheel)>::value_type>::max();

                // Only happens when P very large (100K)
                // Fix by changing x_reindex_wheel type to int32_t
                assert(coprime_count_wheel < x_reindex_limit);  // See comment above.
            }
            #endif
        }

        K_mod2310 = mpz_fdiv_ui(K, 2310);

        is_coprime2310.resize(2*3*5*7*11, 1);
        for (int p : {2, 3, 5, 7, 11})
            for (size_t i = 0; i < is_coprime2310.size(); i += p)
                is_coprime2310[i] = 0;

        is_m_coprime2310.resize(2310, 1);
        for (int p : {2, 3, 5, 7, 11})
            if (config.d % p == 0)
                for (int i = 0; i < 2310; i += p)
                    is_m_coprime2310[i] = 0;

        assert(count(is_coprime2310.begin(), is_coprime2310.end(), 1) == 480);
        /**
         * Because 2310 is a multiple of x_reindex_wheel_size,
         * x_reindex_wheel will always be true if prefiltered by is_coprime2310[n % 310]
         */
        assert(2310 % x_reindex_wheel_size == 0);
    }
};


/**
 * Seperate [start_prime, end_prime] into a series of intervals.
 * interval start and ends have nice round numbers
 * percent controls how quickly intervals can grow.
 * [0, 10], [10, 20], [20, 30], ... [90, 100], [100, 200], [200, 300], [300, 400]
 *
 */
std::vector<std::pair<uint64_t, uint64_t>> split_prime_range_to_intervals(
        uint64_t percent, uint64_t start_prime, uint64_t end_prime) {
    assert(percent == 1 || percent == 100);
    /**
     * first ones need to be small enough to print most mults.
     * later larger enough to reduce overhead in prime iterator
     */
    std::vector<std::pair<uint64_t, uint64_t>> intervals;

    // interval_inc, interval_start
    uint64_t i_inc = 10'000;
    uint64_t i_start = 0;
    while(i_start < end_prime) {
        // Next multiple of ten that keeps interval < percent * start
        while (i_inc * 10 * 100 <= i_start * percent) {
            i_inc *= 10;
        }
        uint64_t i_end = std::min(i_start + i_inc, end_prime);
        if (i_end > start_prime) {
            uint64_t first = std::max(i_start, start_prime) + 1;
            intervals.emplace_back(first, i_end);
        }
        i_start = i_end;
    }
    return intervals;
}

void save_unknowns_method2(
        const struct Config& config,
        const mpz_t &K,
        const Cached &caches,
        const vector<bool> *composite) {

    // ----- Open and Save to Output file
    std::ofstream unknown_file;
    {
        std::string fn = Args::gen_unknown_fn(config, ".txt");
        printf("\nSaving unknowns to '%s'\n", fn.c_str());
        unknown_file.open(fn, std::ios::out);
        assert( unknown_file.is_open() ); // Can't open save_unknowns file
    }

    const uint64_t M_start = config.mstart;
    const uint32_t D = config.d;
    const int32_t SL = config.sieve_length;

    const uint32_t neg_K_mod_d = mpz_cdiv_ui(K, D);
    if (D > 1)
        assert(neg_K_mod_d != 0);

    /**
     * coprime_X is [0, SL] and tells if X is coprime to K = P#
     */
    assert( (signed) caches.coprime_X.size() ==
            std::count(caches.is_offset_coprime.begin(), caches.is_offset_coprime.end(), 1) );

    // factors of D = P#/D
    vector<uint32_t> D_primes = get_sieve_primes(config.p);
    D_primes.erase(
        std::remove_if(D_primes.begin(), D_primes.end(), [&](uint32_t p){ return config.d % p != 0; }),
        D_primes.end());
    assert(D_primes.size() <= 9);  // 23# > 2^32
    assert(D_primes.front() >= 2);

    size_t count_a = 0;
    size_t count_b = caches.valid_mi.size() * caches.coprime_X.size();

    #pragma omp parallel for ordered schedule(dynamic, 1) num_threads(config.threads)
    for (size_t mii = 0; mii < caches.valid_mi.size(); mii++) {
        uint64_t mi = caches.valid_mi[mii];
        uint64_t m = M_start + mi;
        assert(gcd(m, D) == 1);
        assert((signed)mii == caches.m_reindex[mi]);

        const auto &comp = composite[mii];
        const auto &x_reindex_m = caches.x_reindex_wheel[m % caches.x_reindex_wheel_size];
        assert(x_reindex_m.size() == (uint64_t) (SL + 1));

        if (config.compression == 2) {
            cout << "compression==2 not implemented here" << endl;
            exit(1);
        }

        std::stringstream line;
        std::stringstream header;
        header << m << " : ";

        // d = 0 yields "|
        for (int d = 0; d <= 1; d++) {
            if (d == 0) {
                line << "|  ";
                header << "-0 ";
                continue;
            }

            assert(d == 1);
            int64_t found = 0;
            line << "|";

            if (config.compression == 1) {
                // RLE
                line << " ";
                int last = 0;

                for (int x : caches.coprime_X) {
                    if (!comp[x_reindex_m[x]]) {
                        found += 1;

                        int delta = x - last;
                        last = x;

                        // Ascii 48 to 122 are all "safe" -> 75 characters -> 5625
                        // Not quite enough so we use 48 + 128 which includes
                        // non printable characters.
                        // assert(0 <= delta && delta < (128L*128));
                        unsigned char upper = 48 + (delta >> 7);
                        unsigned char lower = 48 + (delta & 127);
                        line << upper << lower;
                    }
                }
            } else {
                for (int x : caches.coprime_X) {
                    if (!comp[x_reindex_m[x]]) {
                        line << " +" << x;
                        found += 1;
                    }
                }
            }
            count_a += found;

            line << "\n";
            header << "+" << found << " ";
        }

        #pragma omp ordered
        {
            // unknown_file format is "<m> : -10 : +9 | -1 -5 -10 ... | 3 7 11 ...\n"
            unknown_file << header.str() << line.str();
        }
    }
    if (config.verbose >= 2) {
        printf("\tsaving %ld/%ld (%.1f%%) %s\n",
                count_a, count_b, 100.0f * count_a / count_b,
                config.compression == 2 ? "as bitarray" : "from sieve");
    }
}




method2_stats method2_small_primes(const Config &config, method2_stats &stats,
                          const mpz_t &K,
                          int thread_i,
                          const Cached &caches,
                          const vector<uint32_t> &valid_mi,
                          const uint64_t SMALL_THRESHOLD,
                          vector<bool> *composite) {

    method2_stats temp_stats(thread_i, config, valid_mi.size(), SMALL_THRESHOLD, stats.prob_prime);
    temp_stats.last_prime = stats.last_prime;

    const uint32_t P = config.p;
    const uint32_t D = config.d;

    const uint32_t SL = config.sieve_length;

    const uint32_t x_reindex_wheel_size = caches.x_reindex_wheel_size;

    assert(P < SMALL_THRESHOLD);
    assert(SMALL_THRESHOLD < (size_t) std::numeric_limits<uint32_t>::max());

    mpz_t test;
    mpz_init(test);

    primesieve::iterator iter;
    uint64_t prime = 0;

    while (prime <= SMALL_THRESHOLD) {
        // Handle primes up to (and 1 past) stats.next_mult
        std::vector<std::pair<uint32_t, uint32_t>> p_and_r;

        size_t stop = std::min((prime == 0) ? P : temp_stats.next_print, SMALL_THRESHOLD);
        // Verify this interval contains non-zero numbers (requires all threads call method2_print)
        assert(stop >= prime);
        for (prime = iter.next_prime(); ; prime = iter.next_prime()) {
            temp_stats.pi_interval += 1;

            if (x_reindex_wheel_size % prime == 0) {
                if (thread_i == 0 && config.verbose >= 2) {
                    printf("\t%ld handled by coprime wheel(%d)\n", prime, x_reindex_wheel_size);
                }
                continue;
            }

            // Others are handled by is_offset_coprime above
            if (D % prime == 0 || prime > P) {
                const uint32_t base_r = mpz_fdiv_ui(K, prime);
                p_and_r.push_back({(uint32_t) prime, base_r});
            }

            if (prime >= stop) break;
        }
        if (!p_and_r.empty() && config.verbose >= 3 && thread_i == 0) {
            printf("\tmethod2_small_primes | %ld primes [%d, %d] stop: %ld\n\n",
                p_and_r.size(), p_and_r.front().first, p_and_r.back().first, prime);
        }

        for (uint32_t mi : valid_mi) {
            int32_t mii = caches.m_reindex[mi];
            assert(mii >= 0);

            uint64_t m = config.mstart + mi;
            const auto &x_reindex_m =
#if METHOD2_WHEEL
                caches.x_reindex_wheel[m % x_reindex_wheel_size];
#else
                caches.x_reindex;
#endif  // METHOD2_WHEEL
            vector<bool> &composite_mii = composite[mii];

            bool centerOdd = ((D & 1) == 0) && (m & 1);
            bool lowIsEven = !centerOdd;

            for (const auto &pr : p_and_r) {
                uint64_t a_prime = pr.first;
                uint64_t base_r = pr.second;
                // For each interval that prints

                // Safe as base_r < prime < (2^32-1)
                uint64_t modulo = (m * base_r) % a_prime;
                // negative modulo
                uint32_t first = modulo > 0 ? a_prime - modulo : 0;

                if (first <= SL) {
                    uint32_t shift = a_prime;
                    if (a_prime > 2) {
                        bool evenFromLow = (first & 1) == 0;
                        bool firstIsEven = lowIsEven == evenFromLow;

#ifdef GMP_VALIDATE_FACTORS
                        validate_factor_m_k_x(temp_stats, test, K, config.mstart + mi,
                                              first, a_prime);
                        assert( (mpz_even_p(test) > 0) == firstIsEven );
                        assert( mpz_odd_p(test) != firstIsEven );
#endif  // GMP_VALIDATE_FACTORS

                        if (firstIsEven) {
                            assert( (first >= SL) || composite_mii[x_reindex_m[first]] );
                            // divisible by 2 move to next multiple (an odd multiple)
                            first += a_prime;
                        }

                        // Don't need to count cross off even multiples.
                        shift *= 2;
                    }

                    for (size_t x = first; x <= SL; x += shift) {
                        /**
                         * NOTE: No synchronization issues
                         * Each thread gets a set of mii so no overlap of bits
                         */
                        uint32_t xii = x_reindex_m[x];
                        if (xii > 0) {
                            composite_mii[xii] = true;
                            temp_stats.small_prime_factors_interval += 1;
                        }
                    }
                }
            }
        }

        // Don't print final partial interval
        if (prime >= temp_stats.next_print) {
            method2_increment_print(prime, valid_mi.size(), composite, temp_stats, config);
        }
    }

    #pragma omp critical
    {
        // Update global counters for a couple variables
        if (thread_i == 0) {
            stats.pi          = temp_stats.pi;
            stats.pi_interval = temp_stats.pi_interval;
            stats.count_coprime_p = temp_stats.count_coprime_p;

            stats.next_print = temp_stats.next_print;
            stats.next_mult = temp_stats.next_mult;
        }

        stats.prime_factors += temp_stats.prime_factors;
        stats.small_prime_factors_interval += temp_stats.small_prime_factors_interval;
        stats.validated_factors += temp_stats.validated_factors;
    }

    mpz_clear(test);
    return temp_stats;
}


void method2_medium_primes(const Config &config, method2_stats &stats,
                           const mpz_t &K,
                           int split_i,
                           const Cached &caches,
                           vector<uint32_t> &coprime_X_thread,
                           const uint64_t prime_start, const uint64_t prime_end,
                           vector<bool> *composite) {
    /**
     * NOTE: Theoretically memory access would be better if composite was transposed
     * composite[m,X] would become composite[X,m]
     * Then reverse the loops (for prime, for X) to (for X, for prime)
     * This would make printing harder (like small primes) and means
     * x_reindex_wheel_size wouldn't help.
     * After transpose composite[x, ?] would fit in cache instead of RAM read.
     * This was tried in commit 03e0f6d1 but failed for a number of reasons
     * 1. Takes a long time to invert composite
     * 2. Have to iterate over all primes/inv_k & recalculate m_start_shift for each X
     * 3. REINDEX_WHEEL doesn't work (would require a different reindex)
     */

    /**
     * Similiar to large primes I want to break this up into a long list of
     * <Prime Range, comprime_X_range> then parallel over that list so all
     * work finishes at the same time.
     *
     * This suffers from the same stats tracking problem as small/large primes
     * Here (and for small_primes) it's possible unknowns can be correctly calculated
     * by having each range calculate count unknowns and summing to the main_stats obj.
     *
     * Only tricky bit is either
     * 1. mutex_mi + locking is needed OR
     *      ~1/2 of all "factors" are found at 10-1000x rate of below so locking
     *      much more needed to avoid collision
     * 2. only one thread can execute per coprime_X_thread_chunk
     *      still have issue at edges
     */

    const uint32_t SIEVE_LENGTH = config.sieve_length;
    assert(prime_end <= (size_t)std::numeric_limits<int64_t>::max());
    // Prime can be larger than int32, prime * SIEVE_LENGTH must not overflow int64
    assert(!__builtin_mul_overflow_p(SIEVE_LENGTH, 4 * prime_end, (int64_t) 0));

    const uint64_t M_start = config.mstart;
    const uint64_t M_inc = config.minc;
    assert(config.minc <= (size_t) std::numeric_limits<int32_t>::max());

#if METHOD2_WHEEL
    const uint32_t x_reindex_wheel_size = caches.x_reindex_wheel_size;
#endif  // METHOD2_WHEEL

    mpz_t test, test2;
    mpz_init(test);
    mpz_init(test2);

    const int K_odd = mpz_odd_p(K);

    primesieve::iterator iter(prime_start);
    uint64_t prime = iter.prev_prime();
    assert(prime <= prime_start);
    prime = iter.next_prime();
    assert(prime >= prime_start);
    assert(prime > SIEVE_LENGTH + 1);
    for (; prime <= prime_end; prime = iter.next_prime()) {
        // Only the first split count primes
        if (split_i == 0) {
            stats.pi_interval += 1;
        }

        const uint64_t base_r = mpz_mod_ui(test, K, prime);
        mpz_set_ui(test2, prime);
        assert(mpz_invert(test, test, test2) > 0);

        const int64_t inv_K = mpz_get_ui(test);
        assert(((__int128) inv_K * base_r) % prime == 1);
        const int64_t neg_inv_K = prime - inv_K;

        // -M_start % p
        const int64_t m_start_shift = (prime - (M_start % prime)) % prime;
        const int64_t mi_0_shift = m_start_shift % prime;

        // Lots of expressive (unoptimized) comments and code removed in 9cf1cf40

        // (mi_0 + X) % 2 == (ms + Sl) % 2
        const bool X_M_parity_check = (M_start + SIEVE_LENGTH) & 1;

        size_t small_factors = 0;
        // Find m*K = X, X in [L, R]
        // NOTE: X is positive [0, SL)
        for (int64_t X : coprime_X_thread) {
            // Safe from overflow as 2 * SL * prime < int64
            int64_t mi_0 = (X * neg_inv_K + mi_0_shift) % prime;

            // Check if X parity == m parity
            if (K_odd && ((X ^ mi_0) & 1) == X_M_parity_check) {
                mi_0 += prime;
            }

            uint64_t shift = prime << K_odd; // (1 + K_odd) * prime;

            // Separate loop when shift > M_inc not significantly faster
            for (uint64_t mi = mi_0; mi < M_inc; mi += shift) {
                uint64_t m = M_start + mi;

                uint32_t m_mod2310 = m % 2310;
                // Filters ~80% or more of m where (m, D) != 1
                if (!caches.is_m_coprime2310[m_mod2310])
                    continue;

                // After initial value this increases by (shift * K_mod2310) % 2310
                uint32_t n_mod2310 = ((caches.K_mod2310 * m_mod2310) + X) % 2310;
                if (!caches.is_coprime2310[n_mod2310])
                    continue;

                if (!caches.is_m_coprime[mi])
                    continue;

                int32_t mii = caches.m_reindex[mi];
                assert(mii >= 0);

                small_factors += 1;

#if METHOD2_WHEEL
                uint32_t xii = caches.x_reindex_wheel[m % x_reindex_wheel_size][X];
#else
                uint32_t xii = caches.x_reindex[X];
#endif  // METHOD2_WHEEL

                assert(xii > 0);

                /**
                 * Note: Risk of race condition / contention is reduced by having
                 * each thread handling a different ranche of xii
                 */
                composite[mii][xii] = true;

#ifdef GMP_VALIDATE_FACTORS
                validate_factor_m_k_x(stats, test, K, M_start + mi, X, prime);
                assert( mpz_odd_p(test) );
#endif  // GMP_VALIDATE_FACTORS
            }
        }
        stats.small_prime_factors_interval += small_factors;

        // Should this be moved out of the loop?
        if (split_i == 0 && prime >= stats.next_print && prime != stats.last_prime) {
            // Print counters & stats.
            method2_increment_print(prime, caches.valid_ms, composite, stats, config);
        }
    }
    mpz_clear(test);
    mpz_clear(test2);
}


void method2_large_primes(Config &config, method2_stats &stats,
                          const mpz_t &K,
                          int THREADS,
                          const uint8_t METHOD2_MUTEX_SHIFT,
                          const Cached &caches,
                          const uint64_t MEDIUM_THRESHOLD,
                          vector<bool> *composite) {
    // TODO|XXX: fix for int64 (and in module_search)
    const uint64_t M_start = config.mstart;
    const uint32_t M_inc = config.minc;

    const uint32_t SIEVE_LENGTH = config.sieve_length;
    const uint32_t SL = SIEVE_LENGTH;

    const uint64_t LAST_PRIME = stats.last_prime;

    const uint32_t x_reindex_wheel_size = caches.x_reindex_wheel_size;

    const bool use_lock = config.threads > 1;

    /**
     * This is a LOT of mutexes (given 1-16 threads)
     * Each is 320 bits! so 27 million mi => 1 GB of mutexes.
     * Each mutex is only locked by 1 thread so sharing with a pool
     * reduce size of the vector without substantially impacting contention.
     */
    assert(METHOD2_MUTEX_SHIFT >= 0 && METHOD2_MUTEX_SHIFT < 10);
    vector<mutex> mutex_mi((caches.valid_ms >> METHOD2_MUTEX_SHIFT) + 1);

#if defined GMP_VALIDATE_LARGE_FACTORS && !defined GMP_VALIDATE_FACTORS
    const uint64_t M_end = M_start + M_inc;
    const uint64_t LARGE_PRIME_THRESHOLD = (1LL << 55) / M_end;
    if (LARGE_PRIME_THRESHOLD < LAST_PRIME && config.verbose >= 1) {
        printf("validating factors from primes > %ld\n", LARGE_PRIME_THRESHOLD);
    }
#endif

    const auto intervals = split_prime_range_to_intervals(1, MEDIUM_THRESHOLD, LAST_PRIME);
    if (config.verbose >= 2) {
        printf("\tmethod2_large_primes %ld intervals\n\n", intervals.size());
    }
    if (intervals.size() && intervals.size() < (size_t) THREADS) {
        printf("\n\nCan you ping the thread with your -u <filename>? something is wonky\n");
        printf("\t-u '%s'\n\n", config.unknown_filename.c_str());
    }


    // Setup CTRL+C catcher
    signal(SIGINT, signal_callback_handler);


    uint64_t last_processed = 0;
    map<uint64_t, method2_stats> stats_to_process;

    /**
     * TODO: benchmark firstprivate(cache)
     * "If a SHARED variable in a parallel region is read by the threads executing the region,
     * but not written to by any of the threads, then specify that variable to be FIRSTPRIVATE
     * instead of SHARED. This avoids accessing the variable by dereferencing a pointer, and
     * avoids cache conflicts." -- https://docs.oracle.com/cd/E19059-01/stud.10/819-0501/7_tuning.html
     *
     * But then
     *
     * https://stackoverflow.com/questions/7865555/openmp-shared-vs-firstprivate-performancewise
     * "False sharing does not occur when only reading the variable, that's my understanding on
     * modern processors at least."
     */

    // Better as "for (const auto &interval : intervals) ..." but not supportted by gcc8
    #pragma omp parallel for schedule(dynamic, 1) num_threads(THREADS)
    for (size_t ii = 0; ii < intervals.size(); ii++) {
        const auto &interval = intervals[ii];

        uint64_t first = interval.first;
        uint64_t end = interval.second;

        if (g_control_c) {
            // Can't break in openMP loop, but this has same effect.
            continue;
        }

        mpz_t test;
        mpz_init(test);

        if (config.verbose >= 3) {
            printf("\tmethod2_large_primes(%d) [%'ld, %'ld]\n",
                    omp_get_thread_num(), first, end);
        }

        // Store sentinal for now
        method2_stats test_stats;
        #pragma omp critical
        stats_to_process[end] = test_stats;

        primesieve::iterator it(first - 1);
        for (uint64_t prime = it.next_prime(); prime <= end; prime = it.next_prime()) {
            test_stats.pi_interval += 1;

            // Big improvement over surround_prime is reusing this for each m.
            const uint64_t base_r = mpz_fdiv_ui(K, prime);


            // temp_mod/x from modulo_search_euclid_all_small is faster and helps avoid overflow
            modulo_search_euclid_all_large_next_only(M_start, M_inc, SL, prime, base_r, [&](
                        uint32_t mi, uint64_t temp_mod) {
                assert (mi < M_inc);

                test_stats.m_stops_interval += 1;

                /**
                 * TODO check this
                 * x = (-m * K) % prime
                 *     Computed as
                 * x = (base_r * (M_start + mi)) % prime;
                 */
                assert( temp_mod <= 2 * SL );
                int32_t x = SL - temp_mod;
                if (x < 0) return;

                // Filters ~80% or more of m where (m, D) != 1
                uint64_t m = M_start + mi;
                uint32_t m_mod2310 = m % 2310;
                if (!caches.is_m_coprime2310[m_mod2310])
                    return;

                /**
                 * Check if (m * K + x) has any small factors
                 * Filters 77.2%, Leaves 22.8%
                 */
                uint32_t n_mod2310 = ((caches.K_mod2310 * m_mod2310) + x) % 2310;
                if (!caches.is_coprime2310[n_mod2310])
                    return;

                // Filters ~75% more (coprime to P#) (requires ~10-200kb)
                if (!caches.is_offset_coprime[x])
                    return;

                /**
                 * Would Filters ~60-80%, only ~10-20% after is_m_coprime2310
                 * Requires a M_inc/8 (a few MB) cache so more likely to spill to L3/RAM
                 */
                if (!caches.is_m_coprime[mi])
                    return;

                // ~99% of factors have been skipped at this point (1 in 80-90 will mark a factor)

#if METHOD2_WHEEL
                uint32_t xii = caches.x_reindex_wheel[m % x_reindex_wheel_size][x];
#else
                uint32_t xii = caches.x_reindex[x];
#endif
                assert(xii > 0); // something wrong with n_mod2310 calculation?

                // if coprime with K, try to toggle off factor.
                test_stats.large_prime_factors_interval += 1;

                int32_t mii = caches.m_reindex[mi];
                assert( mii >= 0 );

                if (use_lock) {
                    mutex_mi[mii >> METHOD2_MUTEX_SHIFT].lock();
                    composite[mii][xii] = true;
                    mutex_mi[mii >> METHOD2_MUTEX_SHIFT].unlock();
                } else {
                    composite[mii][xii] = true;
                }

                #ifdef GMP_VALIDATE_FACTORS
                validate_factor_m_k_x(stats, test, K, m, x, prime);
                #elif defined GMP_VALIDATE_LARGE_FACTORS
                if (prime > LARGE_PRIME_THRESHOLD)
                    validate_factor_m_k_x(stats, test, K, m, x, prime);
                #endif
            });
        }
        mpz_clear(test);


        // Normally this is inside the loop but not anymore
        #pragma omp critical
        {
            {
                test_stats.finished = true;
                stats_to_process[end] = test_stats;
                if (config.verbose >= 3) {
                    size_t queued = 0;
                    for(auto&& kv : stats_to_process) queued += kv.second.pi_interval > 0;
                    printf("\tmethod2_large_primes(%d) finished [%'ld, %'ld] (%ld in queue)\n",
                            omp_get_thread_num(), first, end, queued);
                }
            }

            // Walk through other stats (in increasing order) adding if valid
            while (!stats_to_process.empty()) {
                // std::map keys are sorted
                auto min_end = stats_to_process.begin()->first;
                method2_stats temp = stats_to_process[min_end];

                if (config.verbose >= 3) {
                    printf("\tStats[%ld] -> %ld with %ld (%ld vs %ld)\n",
                            stats_to_process.size(),
                            min_end, temp.pi_interval,
                            last_processed, stats.next_print);
                }
                if (temp.finished) {
                    // Stats ready, can process them
                    stats_to_process.erase(min_end);

                    stats.pi_interval += temp.pi_interval;
                    stats.m_stops_interval += temp.m_stops_interval;
                    stats.large_prime_factors_interval += temp.large_prime_factors_interval;

                    assert(last_processed < min_end);
                    assert(last_processed < stats.next_print || stats.next_print == 0);
                    if (min_end >= stats.next_print) {
                        // Print counters & stats.
                        method2_increment_print(min_end, caches.valid_ms, composite, stats, config);
                    }
                    assert(min_end <= stats.next_print);
                    last_processed = min_end;

                } else {
                    break;
                }
            }
        }
        if (config.verbose >= 3) {
            printf("\tmethod2_large_primes(%d) done (%ld)\n",
                omp_get_thread_num(), MEDIUM_THRESHOLD);
        }
    }
    if (config.verbose >= 2) {
        printf("\tmethod2_large_primes done\n");
    }

    // if is_last would truncate .max_prime by 1 million
    if (g_control_c && (last_processed != LAST_PRIME)) {
        uint64_t prime = last_processed;

        // NOTE: the resulting files were sieved by 1 extra prime
        // they will differ from --max_prime=X in a few entries

        if (prime < 1'000'000) {
            cout << "Exit(2) from CTRL+C @ prime=" << prime << endl;
            exit(2);
        }

        cout << "Breaking loop from CTRL+C @ prime=" << prime << endl;
        // reset unknown_filename if cached;
        config.unknown_filename = "";
        config.max_prime = prime - (prime % 1'000'000);
    }
}

// Would be nice to pass const but CTRL+C handler changes max_prime
void prime_gap_parallel(struct Config& config) {
    // Method2
    const uint64_t M_start = config.mstart;
    const uint32_t M_inc = config.minc;

    const uint32_t P = config.p;

    const uint32_t SIEVE_LENGTH = config.sieve_length;
    const uint32_t SL = SIEVE_LENGTH;

    const uint64_t MAX_PRIME = config.max_prime;

    const size_t THREADS = config.threads;

    uint64_t LAST_PRIME = [&] {
        mpz_t test;
        mpz_init(test);

        mpz_set_ui(test, MAX_PRIME);
        mpz_prevprime(test, test);

        uint64_t temp = mpz_get_ui(test);
        mpz_clear(test);
        return temp;
    }();

    assert( LAST_PRIME <= MAX_PRIME && LAST_PRIME + 500 > MAX_PRIME);

    // ----- Generate primes for P
    const vector<uint32_t> P_primes = get_sieve_primes(P);
    assert( P_primes.back() == P);

    // ----- Sieve stats & Merit Stuff
    mpz_t K;
    const double K_log = prob_prime_and_stats(config, K);
    const double N_log = K_log + log(config.mstart);
    const double prob_prime = 1 / N_log - 1 / (N_log * N_log);


    // ----- Allocate memory

    // Various pre-calculated arrays of is_coprime arrays
    const Cached caches(config, K);
    const size_t valid_ms = caches.valid_ms;
    const uint32_t x_reindex_wheel_size = caches.x_reindex_wheel_size;

    const size_t count_coprime_sieve = caches.coprime_X.size();
    assert( count_coprime_sieve % 2 == 0 );

    const auto THRESHOLDS =
        calculate_thresholds_method2(config, count_coprime_sieve, valid_ms);
    const uint64_t SMALL_THRESHOLD = THRESHOLDS.first;
    const uint64_t MEDIUM_THRESHOLD = THRESHOLDS.second;
    if (config.verbose >= 1) {
        printf("sieve_length:  2x %'d\n", config.sieve_length);
        printf("max_prime:        %'ld\n", config.max_prime);
        printf("small_threshold:  %'ld\n", SMALL_THRESHOLD);
        printf("middle_threshold: %'ld\n", MEDIUM_THRESHOLD);
    }

    // SMALL_THRESHOLD must handle all primes that can mark off two items in SIEVE_LENGTH.
    assert( SMALL_THRESHOLD >= SIEVE_LENGTH );
    assert( MEDIUM_THRESHOLD >= SMALL_THRESHOLD );
    assert( MEDIUM_THRESHOLD <= config.max_prime );

    // this controls how many mii (1 << SHIFT) share a mutex
    const uint8_t METHOD2_MUTEX_SHIFT = (THREADS * 1'000'000ul) < valid_ms ? 8 : 1;


    // ----- Timing
    if (config.verbose >= 2) {
        printf("\n");
    }
    // Prints estimate of PRP/s
    const double prp_time_est = prp_time_estimate_composite(
            N_log, config.verbose * config.show_timing);

    // Detailed timing info about different stages
    combined_sieve_method2_time_estimate(config, K, valid_ms, prp_time_est);

    /**
     * Much space is saved via a reindexing scheme
     * composite[mi][x] (0 <= mi < M_inc, 0 <= x <= SL) is re-indexed to
     *      composite[m_reindex[mi]][x_reindex[x]]
     * m_reindex[mi] with (D, M + mi) > 0 are mapped to -1 (and must be handled by code)
     * x_reindex[x]  with (K, x) > 0 are mapped to 0 (and that bit is ignored)
     * x_reindex_wheel[x] same as x_reindex[x]
     */

    // <char> is faster (0-5%?) than <bool>, but uses 8x the memory.
    // Need to change here and in `save_unknowns_method2` signature.
    vector<bool> *composite = new vector<bool>[valid_ms];
    {
        int align_print = 0;
        /**
         * Per m_inc
         *      4 bytes in m_reindex
         *      1 bit   in is_m_coprime
         * Per valid_ms
         *      4  bytes in caches.valid_mi
         *      40 bytes for vector<bool> instance
         *      40>>8 bytes for mutex instance     // in stage2
         *      count_coprime_sieve + 1 bits
         */

        size_t MB = 8 * 1024 * 1024;
        size_t overhead_bits = M_inc * (8 * sizeof(uint32_t) + 1) +
                               valid_ms * 8 * (sizeof(uint32_t) + sizeof(composite[0]));
        if (THREADS > 1) {
            overhead_bits += valid_ms * 8 * sizeof(mutex) >> METHOD2_MUTEX_SHIFT;
        }

        // Per valid_ms
        size_t guess = overhead_bits + valid_ms * (count_coprime_sieve + 1);
        if (config.verbose >= 1) {
            // Using strings instead of printf so sizes can be aligned.
            std::string s_coprime_m = "coprime m    " +
                std::to_string(valid_ms) + "/" + std::to_string(M_inc) + " ";
            std::string s_coprime_i = "coprime i    " +
                std::to_string(count_coprime_sieve / 2) + "/" + std::to_string(SIEVE_LENGTH);
            align_print = s_coprime_m.size();

            printf("%*s", align_print + (int) s_coprime_i.size(), "");
            printf("  ~%'ld MB overhead\n", overhead_bits / MB);
            printf("%s%s, ~%'ld MB\n", s_coprime_m.c_str(), s_coprime_i.c_str(), guess / MB);
        }

        if (x_reindex_wheel_size > 1) {
            // Update guess with first wheel count for OOM prevention check
            size_t guess_avg_count_coprime = caches.x_reindex_wheel_count[1];
            guess = overhead_bits + valid_ms * (guess_avg_count_coprime + 1);
        }

        // Try to prevent OOM, check composite < 10GB allocation,
        if (guess > (size_t) config.max_mem * 1024 * MB) {
            if (config.verbose >= 1 && x_reindex_wheel_size > 1) {
                size_t allocated = guess - overhead_bits;
                printf("%*s", align_print, "");
                printf("coprime wheel %ld/%d, ~%'ldMB\n",
                    allocated / (2 * valid_ms), SIEVE_LENGTH,
                    allocated / 8 / 1024 / 1024);
            }
            printf("\ncombined_sieve expects to use %'ld MB which is greater than %d GB limit\n",
                    guess / MB, config.max_mem);
            printf("\nAdd `--max-mem %ld` to skip this warning\n", (guess / 1024 / MB) + 1);
            std::this_thread::sleep_for(std::chrono::seconds(10));
            exit(1);
        }

        size_t allocated = 0;
        for (size_t i = 0; i < valid_ms; i++) {
            int m_wheel = (M_start + caches.valid_mi[i]) % x_reindex_wheel_size;
            assert(gcd(m_wheel, x_reindex_wheel_size) == 1);

            // +1 reserves extra 0th entry for x_reindex[x] = 0
            composite[i].resize(caches.x_reindex_wheel_count[m_wheel] + 1, false);
            composite[i][0] = true;
            allocated += composite[i].size();
        }
        if (config.verbose >= 1 && x_reindex_wheel_size > 1) {
            printf("%*s", align_print, "");
            printf("coprime wheel %ld/%d, ~%'ld MB\n",
                allocated / (2 * valid_ms), SIEVE_LENGTH,
                allocated / 8 / 1024 / 1024);
        }

        if (config.verbose >= 1) {
            align_print += 1;  // avoid unused warning
            printf("\n");
        }
    }

    // Used for various stats
    method2_stats stats(/* thread */ 0, config, valid_ms, SMALL_THRESHOLD, prob_prime);
    stats.last_prime = LAST_PRIME;
    stats.prp_time_estimate = prp_time_est;

    if (1) { // Small Primes
        /**
         * NOTE: For primes <= SMALL_THRESHOLD, handle per m (with better memory locality)
         * This also avoids need to synchronize access to composite
         * It does make printing slightly harder (see awkward inner loop)
         */

        vector<uint32_t> valid_mi_split[THREADS];
        for (size_t i = 0; i < valid_ms; i++) {
            valid_mi_split[i * THREADS / valid_ms].push_back(caches.valid_mi[i]);
        }

        #pragma omp parallel for num_threads(THREADS)
        for (size_t thread_i = 0; thread_i < THREADS; thread_i++) {
            // Helps keep printing in order
            std::this_thread::sleep_for(milliseconds(50 * thread_i));

            if (config.verbose + (THREADS > 1) >= 3) {
                printf("\tThread %ld method2_small_primes(%'ld/%'ld)\n",
                        thread_i, valid_mi_split[thread_i].size(), valid_ms);
            }

            auto temp_stats = method2_small_primes(
                config, stats, K, thread_i,
                caches,
                valid_mi_split[thread_i],
                SMALL_THRESHOLD, composite);

            #pragma omp critical
            { // Check stats in sync
                stats.interval_t = temp_stats.interval_t;
                stats.total_unknowns = temp_stats.total_unknowns;
                stats.current_prob_prime = temp_stats.current_prob_prime;
            }
        }
    } else {
        // Have to do this to make method2_increment_print happy
        method2_increment_print(config.p, valid_ms, composite, stats, config);

        // XXX: write a random fraction of composite false (and same below).
    }


    if (1) { // Medium Primes
        /**
         * Old (5e782547) parallelization was:
         * each THREADS gets an equal share of work (coprime_X_split)
         *      handles all mediums primes for those comprime_X.
         * Pros:
         *      no locking of composite needed (minus trying to round to multiples of 8)
         * Cons:
         *      stats aren't easier to compute
         *      THREADS = CORES + 1, and last THREAD might take 2x as long
         *          (this really happens that longer threads take 1+ hour longer)
         *
         *
         * New parallelization involves breaking up many extra coprime_X_split
         *      and also breaking up primes into ranges (like large_primes)
         * Pros:
         *      distributes evenly
         *      starts partially works
         * Cons:
         *      more code
         *      still can't do counting in stats (because of wheel)
         *      slowest thread probably unlocks then reclaims least complete range.
         */

        const auto intervals = split_prime_range_to_intervals(100, SMALL_THRESHOLD, MEDIUM_THRESHOLD);
        if (config.verbose >= 1) {
            printf("\tmethod2_medium_primes %ld intervals x %ld THREADS\n\n",
                    intervals.size(), THREADS);
        }

        /**
         * Note: Each extra split means another call to mpz_invert for all primes.
         * This is fast, ~5 minutes for 300M inverts (limit = 7e9) so ~1% overheard
         */
        const size_t NUM_SPLITS = THREADS + 3;
        vector<uint32_t> coprime_X_split[NUM_SPLITS];
        size_t coprime_X_count = caches.coprime_X.size();

        /**
         * To reduce synchronization issues round coprime_X_split to a
         * multiples of 8 so vector<bool>[i] doesn't overlaps between threads
         *
         * XXX: With wheel reindexing this doesn't work!
         * only breaks when coprime_X_split[{0,1,2}] and coprime_X_split[{n-2,n-1,n}]
         * both modify the same composite[m][reindexed] at the sametime
         * could lock/unlock when using boundy elements but would be lots of code.
         */
        {
            size_t xi_start = 0;
            for (size_t t = 0; t < NUM_SPLITS; t++) {
                bool is_final = (t+1) >= NUM_SPLITS;
                size_t xi_end = is_final ?
                    coprime_X_count : ((t+1) * coprime_X_count / NUM_SPLITS);

                // round end down to multiples of 8;
                if (!is_final) {
                    xi_end -= (xi_end % 8);
                }
                assert( xi_start % 8 == 0);
                assert( (xi_end % 8 == 0) || is_final );

                for (; xi_start < xi_end; xi_start++) {
                    coprime_X_split[t].push_back(caches.coprime_X[xi_start]);
                }
            }

            if (config.verbose + (THREADS > 1) >= 2) {
                for (size_t i = 0; i < NUM_SPLITS; i++) {
                    const auto split = coprime_X_split[i];
                    printf("\tSplit %ld method2_medium_primes(%'ld/%'ld) [%d, %d]\n",
                            i, split.size(), coprime_X_count, split.front(), split.back());
                }
            }
        }

        // Can only work on one coprime_X_split interval at a time.
        mutex state_lock;
        vector<int>    split_locked(NUM_SPLITS, 0);
        vector<size_t> split_progress(NUM_SPLITS, 0);

        #pragma omp parallel for num_threads(THREADS)
        for (size_t work_i = 0; work_i < intervals.size() * NUM_SPLITS; work_i++) {
            // Find least finished split (not currently being processed)
            size_t min_p = intervals.size(), max_p = 0, min_index = 0;
            {
                state_lock.lock();
                for (size_t i = 0; i < NUM_SPLITS; i++) {
                    if (split_locked[i] == 0) {
                        max_p = std::max(max_p, split_progress[i]);
                        if (split_progress[i] < min_p) {
                            min_p = split_progress[i];
                            min_index = i;
                        }
                    }
                }

                // There should be some item to work on
                assert(min_p < intervals.size());

                // lock work item (min_p, min_index)
                split_locked[min_index] = 1;

                if ((max_p - min_p) > 2 && config.verbose >= 1) {
                    printf("\tsplit %ld(%ld) falling behind(%ld) @ work item %ld\n",
                        min_index, min_p, max_p, work_i);
                }
                if (config.verbose >= 3) {
                    printf("\twork item %ld, %ld(%ld) | p: [%ld, %ld] X: [%d, %d]\n",
                        work_i, min_index, min_p,
                        intervals[min_p].first, intervals[min_p].second,
                        coprime_X_split[min_index].front(), coprime_X_split[min_index].back());

                }

                state_lock.unlock();
            }

            // Do work item
            method2_medium_primes(config, stats, K,
                                  min_index,
                                  caches,
                                  coprime_X_split[min_index],
                                  intervals[min_p].first, intervals[min_p].second,
                                  composite);

            // Update progress and lock
            {
                state_lock.lock();

                split_progress[min_index]++;
                split_locked[min_index] = 0;

                bool is_final = split_progress[min_index] == intervals.size();
                if (is_final && config.verbose >= 1) {
                    time_t rawtime = std::time(nullptr);
                    struct tm *tm = localtime( &rawtime );
                    printf("\tmethod2_medium_primes(split %ld) done @ %d:%02d:%02d\n",
                        min_index, tm->tm_hour, tm->tm_min, tm->tm_sec);
                }

                state_lock.unlock();
            }
        }

        if (MEDIUM_THRESHOLD >= stats.last_prime) {
            // Handle final print (if no large_prime will be run)
            method2_increment_print(
                stats.last_prime, valid_ms, composite, stats, config);
        }
        if (config.verbose >= 1) {
            printf("\tmethod2_medium_primes all done\n");
        }
    }


    if (1) { // Large Primes
        // THREADS are handled inside function.
        method2_large_primes(
            config, stats,
            K,
            THREADS, METHOD2_MUTEX_SHIFT,
            caches,
            MEDIUM_THRESHOLD,
            composite);
    }


    { // Verify some computation.
        // Likely zeroed in the last interval, but needed if no printing
        stats.m_stops += stats.m_stops_interval;

        // See Merten's Third Theorem
        float expected_m_stops = (log(log(LAST_PRIME)) - log(log(MEDIUM_THRESHOLD))) * SL * M_inc;
        float error_percent = 100.0 * fabs(expected_m_stops - stats.m_stops) / expected_m_stops;
        if (config.verbose >= 3 || error_percent > 0.1 ) {
            printf("Estimated modulo searches (m/prime) error %.2f%%,\t%ld vs expected %.0f\n",
                error_percent, stats.m_stops, expected_m_stops);
        }
    }

    if (config.save_unknowns) {
        auto s_save_t = high_resolution_clock::now();

        save_unknowns_method2(config, K, caches, composite);

        auto s_stop_t = high_resolution_clock::now();
        double   secs = duration<double>(s_stop_t - stats.start_t).count();

        if (config.verbose >= 2) {
            printf("Saving unknowns took %.1f seconds\n",
                    duration<double>(s_stop_t - s_save_t).count());
        }

        insert_range_db(config, valid_ms, THREADS * secs);
    }

    if (config.verbose >= 2) {
        printf("combined sieve Done!\n");
    }

    delete[] composite;

    mpz_clear(K);
}
