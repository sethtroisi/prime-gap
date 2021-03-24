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

// Tweaking this doesn't seem to method1 much.
// method2 is more sensitive and set it's own.
#define SMALL_PRIME_LIMIT_METHOD1       400'000

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

    Config config = Args::argparse(argc, argv);

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
            //exit(1); TODO REENABLE
        }

        if (config.valid == 0) {
            Args::show_usage(argv[0]);
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

        if (!config.compression && (config.minc * config.sieve_length) > 100'000'000'000L) {
            printf("\tSetting --bitcompress to prevent very large output file\n");
            config.compression = 2;
        }
        if (!config.compression && (config.minc * config.sieve_length) > 30'000'000'000L) {
            printf("\tSetting --rle to prevent very large output file\n");
            config.compression = 1;
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
                        uint32_t center = (m * base) % config.d;
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
    //const uint64_t P = config.p;
    const uint64_t D = config.d;
    const uint64_t M_start = config.mstart;
    const uint64_t M_inc = config.minc;

    const unsigned int SIEVE_LENGTH = config.sieve_length;
    const unsigned int SL = SIEVE_LENGTH;

    const uint64_t MAX_PRIME = config.max_prime;

    mpz_t test;
    mpz_init(test);

    if (config.verbose >= 2) {
        printf("\n");
        printf("sieve_length: 2x %'d\n", config.sieve_length);
        printf("max_prime:       %'ld\n", MAX_PRIME);
        printf("\n");
    }

    // ----- Generate primes under SMALL_PRIME_LIMIT_METHOD1
    vector<uint32_t> small_primes;
    primesieve::generate_primes(SMALL_PRIME_LIMIT_METHOD1, &small_primes);

    // ----- Merit / Sieve stats
    mpz_t K;
    prob_prime_and_stats(config, K);


    // ----- Sieve stats
    const size_t SMALL_PRIME_PI = small_primes.size();
    {
        // deals with all primes that can mark off two items in SIEVE_LENGTH.
        assert( SMALL_PRIME_LIMIT_METHOD1 > 2 * SIEVE_LENGTH );
        if (config.verbose >= 1) {
            printf("\tUsing %'ld primes for SMALL_PRIME_LIMIT(%'d)\n\n",
                SMALL_PRIME_PI, SMALL_PRIME_LIMIT_METHOD1);
        }
        assert( small_primes[SMALL_PRIME_PI-1] < SMALL_PRIME_LIMIT_METHOD1);
        assert( small_primes[SMALL_PRIME_PI-1] + 200 > SMALL_PRIME_LIMIT_METHOD1);
    }

    const auto  s_setup_t = high_resolution_clock::now();

    // ----- Allocate memory for a handful of utility functions.

    // Remainders of (p#/d) mod prime
    typedef pair<uint64_t,uint64_t> p_and_r;
    vector<p_and_r> prime_and_remainder;
    prime_and_remainder.reserve(SMALL_PRIME_PI);

    // Big improvement over surround_prime is avoiding checking each large prime.
    // vector<m, vector<pi>> for large primes that only rarely divide a sieve
    int s_large_primes_rem = 0;

    double expected_primes_per = 0;

    // To save space, only save remainder for primes that divide ANY m in range.
    // This helps with memory usage when MAX_PRIME >> SL * MINC;
    auto *large_prime_queue = new vector<p_and_r>[M_inc];
    {
        size_t pr_pi = 0;
        if (config.verbose >= 0) {
            printf("\tCalculating first m each prime divides\n");
        }

        // large_prime_queue size can be approximated by
        // https://en.wikipedia.org/wiki/Meissel–Mertens_constant

        // Print "."s during, equal in length to 'Calculating ...'
        size_t print_dots = 38;

        const size_t expected_primes = primepi_estimate(MAX_PRIME);

        long first_m_sum = 0;

        if (config.verbose >= 0) {
            cout << "\t";
        }
        size_t pi = 0;

        primesieve::iterator it;
        for (uint64_t prime = it.next_prime(); prime <= MAX_PRIME; prime = it.next_prime()) {
            pi += 1;
            if (config.verbose >= 0 && (pi * print_dots) % expected_primes < print_dots) {
                cout << "." << std::flush;
            }

            // Big improvement over surround_prime is reusing this for each m.
            const uint64_t base_r = mpz_fdiv_ui(K, prime);

            if (prime <= SMALL_PRIME_LIMIT_METHOD1) {
                prime_and_remainder.emplace_back(prime, base_r);
                pr_pi += 1;
                continue;
            }

            expected_primes_per += (2.0 * SL + 1) / prime;

            // solve base_r * (M + mi) + (SL - 1)) % prime < 2 * SL
            //   0 <= (base_r * M + SL - 1) + base_r * mi < 2 * SL mod prime
            //
            // let shift = (base_r * M + SL - 1) % prime
            //   0 <= shift + base_r * mi < 2 * SL mod prime
            // add (prime - shift) to all three
            //
            //  (prime - shift) <= base_r * mi < (prime - shift) + 2 * SL mod prime
            uint64_t mi = modulo_search_euclid_gcd(
                    M_start, D, M_inc, SL, prime, base_r);

            // signals mi > M_inc
            if (mi == M_inc) continue;

            assert (mi < M_inc);

            // (M_start + mi) * last_prime < int64 (checked in argparse)
            uint64_t first = (base_r * (M_start + mi) + SL) % prime;
            assert( first <= 2*SL );

            //assert ( gcd(M + mi, D) == 1 );

            large_prime_queue[mi].emplace_back(prime, base_r);
            pr_pi += 1;

            s_large_primes_rem += 1;
            first_m_sum += mi;
        }
        if (config.verbose >= 0) {
            cout << endl;
        }

        assert(prime_and_remainder.size() == small_primes.size());
        if (config.verbose >= 1) {
            printf("\tSum of m1: %ld\n", first_m_sum);
            if (pi == expected_primes) {
                printf("\tPrimePi(%ld) = %ld\n", MAX_PRIME, pi);
            } else {
                printf("\tPrimePi(%ld) = %ld guessed %ld\n", MAX_PRIME, pi, expected_primes);
            }

            printf("\t%ld primes not needed (%.1f%%)\n",
                (pi - SMALL_PRIME_PI) - pr_pi,
                100 - (100.0 * pr_pi / (pi - SMALL_PRIME_PI)));

            double mertens3 = log(log(MAX_PRIME)) - log(log(SMALL_PRIME_LIMIT_METHOD1));
            double theory_count = (2 * SL + 1) * mertens3;
            printf("\texpected large primes/m: %.1f (theoretical: %.1f)\n",
                expected_primes_per, theory_count);
        }
    }
    if (config.verbose >= 0) {
        auto  s_stop_t = high_resolution_clock::now();
        double   secs = duration<double>(s_stop_t - s_setup_t).count();
        printf("\n\tSetup took %.1f seconds\n", secs);
    }


    // ----- Open and Save to Output file
    std::ofstream unknown_file;
    if (config.save_unknowns) {
        std::string fn = Args::gen_unknown_fn(config, ".txt");
        printf("\nSaving to '%s'\n", fn.c_str());
        unknown_file.open(fn, std::ios::out);
        assert( unknown_file.is_open() ); // Can't open save_unknowns file
    }


    // ----- Main sieve loop.

    vector<char> composite[2] = {
        vector<char>(SIEVE_LENGTH+1, 0),
        vector<char>(SIEVE_LENGTH+1, 0)
    };
    assert( composite[0].size() == SIEVE_LENGTH+1 );
    assert( composite[1].size() == SIEVE_LENGTH+1 );

    // Used for various stats
    long  s_tests = 0;
    auto  s_start_t = high_resolution_clock::now();
    long  s_total_unknown = 0;
    long  s_t_unk_prev = 0;
    long  s_t_unk_next = 0;
    long  s_large_primes_tested = 0;

    uint64_t last_mi = M_inc - 1;
    for (; last_mi > 0 && gcd(M_start + last_mi, D) > 1; last_mi -= 1);
    assert(last_mi >= 0 && last_mi < M_inc);
    assert(gcd(M_start + last_mi, D) == 1);

    for (uint64_t mi = 0; mi < M_inc; mi++) {
        const uint64_t m = M_start + mi;
        if (gcd(m, D) > 1) {
            assert( large_prime_queue[mi].empty() );
            continue;
        }

        // Reset sieve array to unknown.
        std::fill_n(composite[0].begin(), SIEVE_LENGTH+1, 0);
        std::fill_n(composite[1].begin(), SIEVE_LENGTH+1, 0);
        // center is always composite.
        composite[0][0] = composite[1][0] = 1;

        // For small primes that we don't do trick things with.
        for (const auto& pr : prime_and_remainder) {
            const uint64_t modulo = (pr.second * m) % pr.first;
            //            const auto& [prime, remainder] = prime_and_remainder[pi];
            //            const uint64_t modulo = (remainder * m) % prime;

            for (size_t x = modulo; x <= SIEVE_LENGTH; x += pr.first) {
                composite[0][x] = true;
            }

            // Not technically correct but fine to skip modulo == 0
            int first_negative = pr.first - modulo;
            assert(first_negative >= 0);
            for (size_t x = first_negative; x <= SIEVE_LENGTH; x += pr.first) {
                composite[1][x] = true;
            }
        }

        // Maybe useful for some stats later.
        // int unknown_small_l = std::count(composite[0].begin(), composite[0].end(), false);
        // int unknown_small_u = std::count(composite[1].begin(), composite[1].end(), false);

        for (const auto& pr : large_prime_queue[mi]) {
            s_large_primes_tested += 1;
            s_large_primes_rem -= 1;

            const auto& prime = pr.first;
            const auto& remainder = pr.second;

            // Large prime should divide some number in SIEVE for this m
            // When done find next mi where prime divides a number in SIEVE.
            const uint64_t modulo = (remainder * m) % prime;

            #ifdef GMP_VALIDATE_FACTORS
                mpz_mul_ui(test, K, m);
                assert(modulo == mpz_fdiv_ui(test, prime));
            #endif  // GMP_VALIDATE_FACTORS

            if (modulo <= SIEVE_LENGTH) {
                // Just past a multiple
                composite[0][modulo] = true;
            } else {
                // Don't have to deal with 0 case anymore.
                int64_t first_positive = prime - modulo;
                assert(first_positive <= SIEVE_LENGTH);  // Bad next m!
                // Just before a multiple
                composite[1][first_positive] = true;
            }

            // Find next mi where primes divides part of SIEVE
            {
                uint64_t start = mi + 1;
                uint64_t next_mi = start + modulo_search_euclid_gcd(
                        M_start + start, D, M_inc - start, SL, prime, remainder);
                if (next_mi == M_inc) continue;

                // (M_start + mi) * prime < int64 (checked in argparse)
                uint64_t mult = (remainder * (M_start + next_mi) + SL) % prime;
                assert(mult < (2 * SL + 1));

                //assert ( gcd(M_start + next_mi, D) == 1 );

                large_prime_queue[next_mi].push_back(pr);
                s_large_primes_rem += 1;
            }
        }
        large_prime_queue[mi].clear();
        large_prime_queue[mi].shrink_to_fit();

        s_tests += 1;

        int unknown_p = std::count(composite[0].begin(), composite[0].end(), false);
        int unknown_n = std::count(composite[1].begin(), composite[1].end(), false);
        s_total_unknown += unknown_p + unknown_n;
        s_t_unk_prev += unknown_p;
        s_t_unk_next += unknown_n;

        // Save unknowns
        if (config.save_unknowns) {
            save_unknowns_method1(
                unknown_file,
                m, unknown_p, unknown_n,
                SL, composite
            );
        }

        bool is_last = (mi == last_mi);

        if ((config.verbose + is_last >= 1) &&
                ((s_tests == 1 || s_tests == 10 || s_tests == 100 || s_tests == 500 || s_tests == 1000) ||
                 (s_tests % 5000 == 0) || is_last) ) {
            auto s_stop_t = high_resolution_clock::now();
            double   secs = duration<double>(s_stop_t - s_start_t).count();
            double t_secs = duration<double>(s_stop_t - s_setup_t).count();

            printf("\t%ld %4d <- unknowns -> %-4d\n",
                    m, unknown_p, unknown_n);

            if (config.verbose + is_last >= 1) {
                // Stats!
                printf("\t    intervals %-10ld (%.2f/sec, with setup per m: %.2g)  %.0f seconds elapsed\n",
                        s_tests, s_tests / secs, t_secs / s_tests, secs);
                printf("\t    unknowns  %-10ld (avg: %.2f), %.2f%% composite  %.2f <- %% -> %.2f%%\n",
                        s_total_unknown, s_total_unknown / ((double) s_tests),
                        100.0 * (1 - s_total_unknown / ((2.0 * SIEVE_LENGTH + 1) * s_tests)),
                        100.0 * s_t_unk_prev / s_total_unknown,
                        100.0 * s_t_unk_next / s_total_unknown);
                printf("\t    large prime remaining: %d (avg/test: %ld)\n",
                        s_large_primes_rem, s_large_primes_tested / s_tests);
            }
        }
    }

    {
        double primes_per_m = s_large_primes_tested / s_tests;
        double error_percent = 100.0 * fabs(expected_primes_per - primes_per_m) /
            expected_primes_per;
        if (config.verbose >= 2 || error_percent > 0.5 ) {
            printf("\n");
            printf("Estimated primes/m error %.2f%%,\t%.1f vs expected %.1f\n",
                error_percent, primes_per_m, expected_primes_per);
        }
    }

    if (config.save_unknowns) {
        auto s_stop_t = high_resolution_clock::now();
        double   secs = duration<double>(s_stop_t - s_setup_t).count();
        insert_range_db(config, s_tests, secs);
    }

    // Should be cleaning up after self.
    for(uint32_t mi = 0; mi < M_inc; mi++)  {
        assert( large_prime_queue[mi].empty() );
    }

    // ----- cleanup

    delete[] large_prime_queue;
    mpz_clear(K);
    mpz_clear(test);
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
                uint32_t thread_i,
                const struct Config& config,
                size_t valid_ms,
                uint64_t threshold,
                double initial_prob_prime
        ) {
            thread = thread_i;

            start_t = high_resolution_clock::now();
            interval_t = high_resolution_clock::now();
            total_unknowns = (2 * config.sieve_length + 1) * valid_ms;

            if (threshold <= 100000)
               next_mult = 10000;

            prob_prime = initial_prob_prime;
            current_prob_prime = prob_prime;
        }

        // Some prints only happen if thread == 0
        int32_t thread = 0;
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

        // TODO Figure out how to set this.
        size_t count_coprime_p = 0;
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
                printf("WARNING stats aren't synchronized when "
                       "running with multiple threads(%d)\n", config.threads);
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
        uint32_t SIEVE_INTERVAL = 2 * config.sieve_length + 1;

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
                    100.0 - 100.0 * t_total_unknowns / (SIEVE_INTERVAL * valid_ms),
                    100.0 * new_composites / (SIEVE_INTERVAL * valid_ms),
                    new_composites);

                if (prime > 100000 && prime > config.p) {
                    // verify total unknowns & interval unknowns
                    const double prob_prime_coprime_P = prob_prime_coprime(config);

                    float e_unknowns = stats.count_coprime_p * (prob_prime_coprime_P / prob_prime_after_sieve);

                    float delta_composite_rate = delta_sieve_prob * prob_prime_coprime_P;
                    float e_new_composites = stats.count_coprime_p * delta_composite_rate;

                    float error = 100.0 * fabs(e_unknowns - t_total_unknowns) / e_unknowns;
                    float interval_error = 100.0 * fabs(e_new_composites - new_composites) / e_new_composites;

                    // TODO: store max error and print at the end.
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
                    skipped_prp / int_secs,
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
        uint64_t prime, uint32_t SL) {
#ifdef GMP_VALIDATE_FACTORS
    stats.validated_factors += 1;
    mpz_mul_ui(test, K, m);
    mpz_sub_ui(test, test, SL);
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
        vector<int32_t> valid_mi;
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
        // if (i, D) == 1 (first 2310 values)
        // XXX: Can I just do is_m_coprime[i] and that will live in cache?
        vector<bool> is_m_coprime2310;


        // X which are coprime to K (X has SIEVE_LENGTH added so x is positive)
        vector<uint32_t> coprime_X;
        // reindex composite[m][X] for composite[m_reindex[m]][x_reindex[X]]
        // Special 0'th entry stands for all not coprime
        vector<uint32_t> x_reindex;
        // if [x+SL] is coprime to K (x has SL added to make value always positive)
        vector<char> is_offset_coprime;


        // reindex composite[m][i] using (m, wheel) (wheel is 1!, 2!, 3!, or 5!)
        // This could be first indexed by x_reindex,
        // Would reduce size from wheel * (2*SL+1) to wheel * coprime_i

        // Note: Larger wheel eliminates more numbers but takes more space.
        // 6 (saves 2/3 memory), 30 (saves 11/15 memory)
        uint32_t x_reindex_m_wheel;
        vector<uint32_t> x_reindex_wheel[METHOD2_WHEEL_MAX];
        vector<size_t> x_reindex_wheel_count;


        int32_t K_mod210;
        int32_t neg_SL_mod210;

        /** is_comprime210[i] = (i % 2) &&(i % 3) && (i % 5) && (i % 7) */
        vector<char> is_coprime210;
};



Cached setup_caches(const struct Config& config,
                    const mpz_t &K) {
    const uint64_t M_start = config.mstart;
    const uint64_t M_inc = config.minc;

    const uint32_t P = config.p;
    const uint32_t D = config.d;

    const uint32_t SL = config.sieve_length;
    const uint32_t SIEVE_INTERVAL = 2 * SL + 1;

    const vector<uint32_t> P_primes = get_sieve_primes(P);
    assert( P_primes.back() == P);

    // Allocate temp vectors
    vector<int32_t> valid_mi;
    vector<int32_t> m_reindex(M_inc, -1);
    vector<bool> is_m_coprime(M_inc, 1);
    {
        for (uint32_t p : P_primes) {
            if (D % p == 0) {
                // mark off any m = m_start + mi that shares factor with d
                uint64_t first = (p - (M_start % p)) % p;
                assert((M_start + first) % p == 0);
                for (uint64_t mi = first; mi < M_inc; mi += p) {
                    is_m_coprime[mi] = 0;
                }
            }
        }

        for (uint32_t mi = 0; mi < M_inc; mi++) {
            if (is_m_coprime[mi]) {
                assert(gcd(M_start + mi, D) == 1);
                m_reindex[mi] = valid_mi.size();
                valid_mi.push_back(mi);
            }
        }
    }

    // if [X] is coprime to K
    vector<char> is_offset_coprime(SIEVE_INTERVAL, 1);
    // reindex composite[m][X] for composite[m_reindex[m]][x_reindex[X]]
    vector<uint32_t> x_reindex(SIEVE_INTERVAL, 0);
    // which X are coprime to K (X has SIEVE_LENGTH added so x is positive)
    vector<uint32_t> coprime_X;

    // reindex composite[m][i] using (m, wheel) (wheel is 1!,2!,3!,5!)
    // This could be first indexed by x_reindex,
    // Would reduce size from wheel * (2*SL+1) to wheel * coprime_i
#if METHOD2_WHEEL
    // Note: Larger wheel eliminates more numbers but takes more space.
    // 6 seems reasonable for larger numbers  (uses 1/3 memory = 33%)
    // 30 is maybe better for smaller numbers (uses 4/15 memory = 26%)
    // 210 is maybe better (uses 24/105 = 23% memory) but x_reindex_wheel might not fit in memory.
    uint32_t wheel = gcd(D, METHOD2_WHEEL_MAX);
    uint32_t reindex_size = wheel * SIEVE_INTERVAL * sizeof(uint32_t);
    if (reindex_size > 7 * 1024 * 1024) {
        if (wheel % 7) {
            wheel /= 7;
        } else if (wheel % 5) {
            wheel /= 5;
        }
    }
    const uint32_t x_reindex_m_wheel = wheel;
    assert(x_reindex_m_wheel <= METHOD2_WHEEL_MAX); // Static size in Caches will break;
    // XXX dynamic use uint16_t when less than 65536 coprimes?
    vector<uint32_t> x_reindex_wheel[x_reindex_m_wheel];
#else
    const uint32_t x_reindex_m_wheel = 1;
    vector<uint32_t> x_reindex_wheel[1];
#endif  // METHOD2_WHEEL

    vector<size_t> x_reindex_wheel_count(x_reindex_m_wheel, 0);

    for (uint32_t prime : P_primes) {
        if (D % prime != 0) {
            uint32_t first = SL % prime;

            assert( 0 <= first && first < prime );
            assert( (SL - first) % prime == 0 );

            for (size_t x = first; x < SIEVE_INTERVAL; x += prime) {
                is_offset_coprime[x] = 0;
            }
        }
    }
    // Center should be marked composite by every prime.
    assert(is_offset_coprime[SL] == 0);
    {
        size_t coprime_count = 0;
        for (size_t X = 0; X < SIEVE_INTERVAL; X++) {
            if (is_offset_coprime[X] > 0) {
                coprime_X.push_back(X);
                coprime_count += 1;
                x_reindex[X] = coprime_count;
            }
        }
        assert(coprime_count == coprime_X.size());
    }

    // Start at m_wheel == 0 so that re_index_m_wheel == 1 (D=1) works.
    uint32_t mod_neg_SL = x_reindex_m_wheel - (SL % x_reindex_m_wheel);

    for (size_t m_wheel = 0; m_wheel < x_reindex_m_wheel; m_wheel++) {
        if (gcd(x_reindex_m_wheel, m_wheel) > 1) continue;
        x_reindex_wheel[m_wheel].resize(SIEVE_INTERVAL, 0);

        // (m * K - SL) % wheel => (m_wheel - SL) % wheel
        uint32_t mod_center = m_wheel * mpz_fdiv_ui(K, x_reindex_m_wheel);
        uint32_t mod_bottom = (mod_center + mod_neg_SL) % x_reindex_m_wheel;

        size_t coprime_count_wheel = 0;
        for (size_t i = 0; i < SIEVE_INTERVAL; i++) {
            if (is_offset_coprime[i] > 0) {
                if (gcd(mod_bottom + i, x_reindex_m_wheel) == 1) {
                    coprime_count_wheel += 1;
                    x_reindex_wheel[m_wheel][i] = coprime_count_wheel;
                }
            }
        }
        x_reindex_wheel_count[m_wheel] = coprime_count_wheel;
    }
    // TODO print largest x_reindex_wheel_count for uint16_t work.

    int32_t K_mod210 = mpz_fdiv_ui(K, 210);
    int32_t neg_SL_mod210 = (210 - (SL % 210)) % 210;

    vector<char> is_coprime210(210, 1);
    for (int p : {2, 3, 5, 7})
        for (int i = 0; i < 210; i += p)
            is_coprime210[i] = 0;

    vector<bool> is_m_coprime2310(2310, 1);
    for (int p : {2, 3, 5, 7, 11})
        if (config.d % p == 0)
            for (int i = 0; i < 2310; i += p)
                is_m_coprime2310[i] = 0;

    assert(count(is_coprime210.begin(), is_coprime210.end(), 1) == 48);
    assert(210 % x_reindex_m_wheel == 0);   // 210 is a multiple of wheel

    Cached cache;

    cache.valid_mi = valid_mi;
    cache.valid_ms = valid_mi.size();

    cache.m_reindex = m_reindex;
    cache.is_m_coprime = is_m_coprime;
    cache.is_m_coprime2310 = is_m_coprime2310;

    cache.coprime_X = coprime_X;
    cache.x_reindex = x_reindex;
    cache.is_offset_coprime = is_offset_coprime;

    cache.x_reindex_m_wheel = x_reindex_m_wheel;
    cache.x_reindex_wheel_count = x_reindex_wheel_count;

    // XXX: because of caches.x_reindex_wheel not being declared here this gets
    // ugly. Some better C++ programer could comment on what to do.
    for(size_t i = 0; i < x_reindex_m_wheel; i++) {
        cache.x_reindex_wheel[i].swap(x_reindex_wheel[i]);
    }

    cache.K_mod210 = K_mod210;
    cache.neg_SL_mod210 = neg_SL_mod210;

    cache.is_coprime210 = is_coprime210;
    return cache;
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
    const uint32_t SIEVE_INTERVAL = 2 * SL + 1;

    const uint32_t neg_K_mod_d = mpz_cdiv_ui(K, D);
    if (D > 1)
        assert(neg_K_mod_d != 0);

    /**
     * coprime_X is [0, 2*SL+1]
     * coprime_prev is [SL-1 ... 0]
     * coprime_next is [SL+1 ... 2 * SL]
     */
    //
    vector<int32_t> coprime_prev;
    vector<int32_t> coprime_next;
    for (int32_t x : caches.coprime_X) {

        // Need to head outwards from SL this is easiest way
        if (x > SL) {
            uint32_t dist = x - SL;
            coprime_prev.push_back(SL - dist);
            coprime_next.push_back(SL + dist);
        }
    }
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
    size_t count_b = 0;

    #pragma omp parallel for ordered schedule(dynamic, 1) num_threads(config.threads)
    for (size_t mii = 0; mii < caches.valid_mi.size(); mii++) {
        uint64_t mi = caches.valid_mi[mii];
        uint64_t m = M_start + mi;
        assert(gcd(m, D) == 1);
        assert((signed)mii == caches.m_reindex[mi]);

        const auto &comp = composite[mii];
        const auto &x_reindex_m = caches.x_reindex_wheel[m % caches.x_reindex_m_wheel];
        assert(x_reindex_m.size() == (size_t) 2 * SL + 1);

        if (config.compression == 2) {
            vector<char> is_offset_fully_coprime(caches.is_offset_coprime);
            for (uint32_t d : D_primes) {
                // First multiple = -(m * K - SL) % d = (m * -K + SL) % d
                // needs +SL for is_offset_fully_coprime which cancels out
                uint64_t first = (m * neg_K_mod_d + SL) % d;
                for (uint64_t mult = first; mult < SIEVE_INTERVAL; mult += d) {
                    assert(comp[x_reindex_m[mult]] == 1);
                    is_offset_fully_coprime[mult] = 0;
                }
            }

            std::stringstream line;

            /**
             * b generally has most bits set it's possible, but any value is possible
             * especially at low sieve depths. So it's nice to avoid `\n` (dec 10) and
             * other ascii control characters (null, space, ...)
             * Could use base64 (6 bits / byte) or ascii85 (6.4 bits / byte), but I decided
             * to just use 1ABCDEFG (7 bits / byte).
             */

            size_t unknowns = 0;
            size_t bytes_written = 0;
            int bit = 0; // index into byte
            unsigned char b = 1 << 7; // current byte

            // XXX: could potentially improve performance by 2x by using wheel of size two
            for (int32_t x : caches.coprime_X) {
                if (is_offset_fully_coprime[x] == 1) {
                    unsigned char is_composite = comp[x_reindex_m[x]];
                    b |= is_composite << bit;
                    unknowns += !is_composite;

                    bit++;
                    if (bit == 7) { // Every 7 items
                        bytes_written++;
                        bit = 0;

                        line << b;
                        b = 1 << 7;  // reset
                    }
                }
            }
            if (bit != 0) {
                bytes_written++;
                line << b;
            }
            count_a += bytes_written;
            count_b += caches.coprime_X.size();

            #pragma omp ordered
            {
                // unknown_file format is "<m> : 19 <bitcount> || <rawbytes> ...\n"
                unknown_file << m << " : " << unknowns << " " << bytes_written << " || " << line.str() << "\n";
            }

            continue;
        }

        std::stringstream line;
        std::stringstream header;
        header << m << " : ";

        for (int d = 0; d <= 1; d++) {
            int64_t found = 0;
            int32_t sign = d == 0 ? -1 : 1;
            line << "|";

            if (config.compression == 1) {
                // RLE
                line << " ";
                int last = SL;

                for (int32_t x : (d == 0 ? coprime_prev : coprime_next)) {
                    if (!comp[x_reindex_m[x]]) {
                        found += 1;

                        int delta = sign * (x - last);
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
                char prefix = "-+"[d];
                for (int32_t x : (d == 0 ? coprime_prev : coprime_next)) {
                    if (!comp[x_reindex_m[x]]) {
                        line << " " << prefix << sign * (x - SL);
                        found += 1;
                    }
                }
            }
            count_a += found;
            count_b += coprime_prev.size();

            line << " \n"[d];
            char prefix = "-+"[d];
            header << prefix << found << " ";
        }

        #pragma omp ordered
        {
            // unknown_file format is "<m> : -10 : +9 | -1 -5 -10 ... | 3 7 11 ...\n"
            unknown_file << header.str() << line.str();
        }
    }
    if (config.verbose >= 2) {
        // TODO compression=2 message

        printf("\tsaving %ld/%ld (%.1f%%) from sieve array\n",
                count_a, count_b, 100.0f * count_a / count_b);
    }
}




method2_stats method2_small_primes(const Config &config, method2_stats &stats,
                          const mpz_t &K,
                          uint32_t thread_i,
                          const Cached &caches,
                          const vector<int32_t> &valid_mi,
                          const uint64_t SMALL_THRESHOLD,
                          vector<bool> *composite) {

    method2_stats temp_stats(thread_i, config, valid_mi.size(), SMALL_THRESHOLD, stats.prob_prime);
    temp_stats.last_prime = stats.last_prime;

    const uint32_t P = config.p;
    const uint32_t D = config.d;

    const uint32_t SIEVE_LENGTH = config.sieve_length;
    const uint32_t SIEVE_INTERVAL = 2 * SIEVE_LENGTH + 1;

    const uint32_t x_reindex_m_wheel = caches.x_reindex_m_wheel;

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

            if (x_reindex_m_wheel % prime == 0) {
                if (thread_i == 0 && config.verbose >= 2) {
                    printf("\t%ld handled by coprime wheel(%d)\n", prime, x_reindex_m_wheel);
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
            const std::vector<uint32_t> &x_reindex_m =
#if METHOD2_WHEEL
                caches.x_reindex_wheel[m % x_reindex_m_wheel];
#else
                caches.x_reindex;
#endif  // METHOD2_WHEEL
            vector<bool> &composite_mii = composite[mii];

            bool centerOdd = ((D & 1) == 0) && (m & 1);
            bool lowIsEven = centerOdd == (SIEVE_LENGTH & 1);

            for (const auto &pr : p_and_r) {
                uint64_t a_prime = pr.first;
                uint64_t base_r = pr.second;
                // For each interval that prints

                // Safe as base_r < prime < (2^32-1)
                uint64_t modulo = (base_r * m) % a_prime;

                // flip = (m * K - SL) % a_prime
                uint32_t flip = modulo + a_prime - ((SIEVE_LENGTH + 1) % a_prime);
                if (flip >= a_prime) flip -= a_prime;

                uint32_t first = a_prime - flip - 1;
                assert(first < a_prime );

                if (first < SIEVE_INTERVAL) {
                    uint32_t shift = a_prime;
                    if (a_prime > 2) {
                        bool evenFromLow = (first & 1) == 0;
                        bool firstIsEven = lowIsEven == evenFromLow;

#ifdef GMP_VALIDATE_FACTORS
                        validate_factor_m_k_x(temp_stats, test, K, config.mstart + mi,
                                              first, a_prime, SIEVE_LENGTH);
                        assert( (mpz_even_p(test) > 0) == firstIsEven );
                        assert( mpz_odd_p(test) != firstIsEven );
#endif  // GMP_VALIDATE_FACTORS

                        if (firstIsEven) {
                            assert( (first >= SIEVE_INTERVAL) || composite_mii[x_reindex_m[first]] );

                            // divisible by 2 move to next multiple (an odd multiple)
                            first += a_prime;
                        }

                        // Don't need to count cross off even multiples.
                        shift *= 2;
                    }

                    for (size_t x = first; x < SIEVE_INTERVAL; x += shift) {
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
                           uint32_t thread_i,
                           const Cached &caches,
                           vector<uint32_t> &coprime_X, // TODO rename to coprime_X_thread
                           const uint64_t SMALL_THRESHOLD, const uint64_t MEDIUM_THRESHOLD,
                           vector<bool> *composite) {
    /**
     * NOTE: Theoretically memory access would be better if composite was transposed
     * composite[m,X] would become composite[X,m]
     * Then reverse the loops (for prime, for X) to (for X, for prime)
     * This would make printing harder (like small primes) and means
     * x_reindex_m_wheel wouldn't help.
     * After transpose composite[x, ?] would fit in cache instead of RAM read.
     * This was tried in commit 03e0f6d1 but failed for a number of reasons
     * 1. Takes a long time to invert composite
     * 2. Have to iterate over all primes/inv_k & recalculate m_start_shift for each X
     * 3. REINDEX_WHEEL doesn't work (would require a different reindex)
     */

    const uint32_t SIEVE_LENGTH = config.sieve_length;
    assert(MEDIUM_THRESHOLD <= (size_t)std::numeric_limits<int64_t>::max());
    // Prime can be larger than int32, prime * SIEVE_LENGTH must not overflow int64
    assert(!__builtin_mul_overflow_p(SIEVE_LENGTH, 4 * MEDIUM_THRESHOLD, (int64_t) 0));

    const uint64_t M_start = config.mstart;
    const uint64_t M_inc = config.minc;
    assert(config.minc <= (size_t) std::numeric_limits<int32_t>::max());

#if METHOD2_WHEEL
    const uint32_t x_reindex_m_wheel = caches.x_reindex_m_wheel;
#endif  // METHOD2_WHEEL

    mpz_t test, test2;
    mpz_init(test);
    mpz_init(test2);

    const int K_odd = mpz_odd_p(K);

    primesieve::iterator iter(SMALL_THRESHOLD+1);
    uint64_t prime = iter.prev_prime();
    assert(prime <= SMALL_THRESHOLD);
    prime = iter.next_prime();
    assert(prime > SMALL_THRESHOLD);
    assert(prime > (2 * SIEVE_LENGTH + 1));
    for (; prime <= MEDIUM_THRESHOLD; prime = iter.next_prime()) {
        // Only the first thread should
        if (thread_i == 0) {
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

        // -SIEVE_LENGTH * inv_K % prime
        const int64_t SL_shift = (SIEVE_LENGTH * inv_K) % prime;
        const int64_t mi_0_shift = (m_start_shift + SL_shift) % prime;

        const bool M_parity = M_start & 1;

        /* Unoptimized expressive code

        - // (m + M_start) * K = (X - SIEVE_LENGTH)
        - // m = (-X*K^-1 - M_start) % p = (X * -(K^-1) - M_start) % p
          int64_t mi_0 = ((prime - dist) * inv_K + m_start_shift) % prime;
        - assert( (base_r * (mi_0 + M_start) + dist) % prime == 0 );

        - int64_t mi_0 = ((prime - dist) * inv_K + m_start_shift) % prime;
        + int64_t mi_0 = (-X * inv_K + SIEVE_LENGTH * inv_K + m_start_shift)

        // Check if X parity == m parity
        - // (X & 1) == X_odd_test <-> ((X + SIEVE_LENGTH) % 2 == 1)
        - const bool X_odd_test = (SIEVE_LENGTH & 1) == 0;
        - const bool M_X_parity = (M_start & 1) ^ (SIEVE_LENGTH & 1);
        - int32_t dist = X - SIEVE_LENGTH;
        - if (((dist ^ mi_0) & 1) == M_X_parity)
        + const bool M_parity = M_start & 1;
        + if (((X ^ mi_0) & 1) == (M_start & 1))

        - / **
        -  * When K is odd
        -  * m parity (even/odd) must not match X parity
        -  * odd m * odd K + odd X   -> even
        -  * even m * odd K + even X -> even
        -  * /
        - size_t m_odd = (M_start + mi_0) & 1;
        - if (((X & 1) == X_odd_test) == m_odd)
        + if (((X^mi_0) & M_start_X_odd_not_same_parity))

            stats.small_prime_factors_interval += 1;
        -    / *
        -    // Doesn't seem to help, slightly tested.
        -    if (D_mod3 && ((dist + K_mod3 * m) % 3 == 0))
        -        continue;
        -    if (D_mod5 && ((dist + K_mod5 * m) % 5 == 0))
        -        continue;
        -    if (D_mod7 && ((dist + K_mod7 * m) % 7 == 0))
        -        continue;
        -    // * /
        */

        size_t small_factors = 0;
        // Find m*K = X, X in [L, R]
        for (int64_t X : coprime_X) {
            // Safe from overflow as SL * prime < int64
            int64_t mi_0 = (X * neg_inv_K + mi_0_shift) % prime;

            uint64_t shift = prime << K_odd; // (1 + K_odd) * prime;
            // Check if X parity == m parity
            if (K_odd && ((X ^ mi_0) & 1) == M_parity) {
                mi_0 += prime;
            }

            // Separate loop when shift > M_inc not significantly faster
            for (uint64_t mi = mi_0; mi < M_inc; mi += shift) {
                uint64_t m = M_start + mi;
                uint32_t m_mod2310 = m % 2310;
                if (!caches.is_m_coprime2310[m_mod2310])
                    continue;

                uint32_t n_mod210 = ((caches.K_mod210 * m_mod2310) + X + caches.neg_SL_mod210) % 210;
                if (!caches.is_coprime210[n_mod210])
                    continue;

                // Wide memory read
                if (!caches.is_m_coprime[mi])
                    continue;

                int32_t mii = caches.m_reindex[mi];
                assert(mii >= 0);

                small_factors += 1;

#if METHOD2_WHEEL
                uint32_t xii = caches.x_reindex_wheel[m % x_reindex_m_wheel][X];
#else
                uint32_t xii = caches.x_reindex[X];
#endif  // METHOD2_WHEEL

                assert(xii > 0);
                composite[mii][xii] = true;

#ifdef GMP_VALIDATE_FACTORS
                validate_factor_m_k_x(stats, test, K, M_start + mi, X, prime, SIEVE_LENGTH);
                assert( mpz_odd_p(test) );
#endif  // GMP_VALIDATE_FACTORS
            }
        }
        stats.small_prime_factors_interval += small_factors;

        if (thread_i == 0 && prime >= stats.next_print && prime != stats.last_prime) {
            // Print counters & stats.
            method2_increment_print(prime, caches.valid_ms, composite, stats, config);
        }
    }
    if (config.verbose >= 2) {
        time_t rawtime = std::time(nullptr);
        struct tm *tm = localtime( &rawtime );
        printf("\tmethod2_medium_primes(thread %d) done @ %d:%02d:%02d\n",
            thread_i, tm->tm_hour, tm->tm_min, tm->tm_sec);
    }


    mpz_clear(test);
    mpz_clear(test2);
}


void method2_large_primes(Config &config, method2_stats &stats,
                          const mpz_t &K,
                          uint32_t THREADS,
                          const Cached &caches,
                          const uint64_t MEDIUM_THRESHOLD,
                          vector<bool> *composite) {
    const uint32_t M_start = config.mstart;
    const uint32_t M_inc = config.minc;

    const uint32_t SIEVE_LENGTH = config.sieve_length;
    const uint32_t SL = SIEVE_LENGTH;

    const uint64_t LAST_PRIME = stats.last_prime;

    const uint32_t x_reindex_m_wheel = caches.x_reindex_m_wheel;

    // This is a LOT of mutexes (given 1-16 threads) might be better to >> 8
    // If that could keep in memory...
    vector<mutex> mutex_mi(caches.valid_ms);

    /**
     * 7# = 2*3*5*7 = 210
     * (1/2) * (2/3) * (4/5) * (6/7) = 23%
     * 77% of numbers have a small factor
     * Helps avoid wide reads to x_reindex_wheel, m_not_coprime ...
     */
    const int32_t K_mod210 = mpz_fdiv_ui(K, 210);
    const int32_t neg_SL_mod210 = (210 - (SL % 210)) % 210;

    const vector<char> is_coprime210 = []() {
        vector<char> temp(210, 1);
        for (int p : {2, 3, 5, 7})
            for (int i = 0; i < 210; i += p)
                temp[i] = 0;
        return temp;
    }();
    assert(count(is_coprime210.begin(), is_coprime210.end(), 1) == 48);
    assert(210 % x_reindex_m_wheel == 0);   // wheel divides 210

    const vector<bool> is_m_coprime2310 = [&]() {
        vector<bool> temp(2310, 1);
        for (int p : {2, 3, 5, 7})
            if (config.d % p == 0)
                for (int i = 0; i < 2310; i += p)
                    temp[i] = 0;
        return temp;
    }();


    // Setup CTRL+C catcher
    signal(SIGINT, signal_callback_handler);


    /**
     * Use un-even intervals
     * first ones need to be small enough to print most mults.
     * later larger to reduce overhead in prime iterator
     */
    std::vector<std::pair<uint64_t, uint64_t>> intervals;
    {
        uint64_t start = 0;
        while(start < LAST_PRIME) {
            const uint64_t interval_inc =
                start >= 2'000'000'000'000 ? 10'000'000'000 :
                    (start >= 200'000'000'000 ? 1'000'000'000 :
                    (start >= 20'000'000'000 ? 100'000'000 :
                    (start >= 2'000'000'000 ? 10'000'000 :
                    (start >= 100'000'000 ? 1'000'000 : 100'000))));
            uint64_t end = std::min(start + interval_inc, LAST_PRIME);

            if (end > MEDIUM_THRESHOLD) {
                uint64_t first = std::max(start, MEDIUM_THRESHOLD) + 1;
                intervals.emplace_back(first, end);
            }
            start = end;
        }

        if (config.verbose >= 2) {
            printf("\tmethod2_large_primes %ld intervals\n", intervals.size());
        }
    }

    const bool use_lock = config.threads > 1;

    uint64_t last_processed = 0;
    map<uint64_t, method2_stats> stats_to_process;

    // Better as "for (const auto &interval : intervals) ..." but not supportted by gcc8
    #pragma omp parallel for schedule(dynamic, 1) num_threads(THREADS)
    for (size_t ii = 0; ii < intervals.size(); ii++) {
        const auto &interval = intervals[ii];

        uint64_t first = interval.first;
        uint64_t end = interval.second;
        // Else it's possible test_stats.pi_interval == 0 below
        assert(first + 1000 < end);


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
            modulo_search_euclid_all_large(M_start, M_inc, SL, prime, base_r, [&](
                        uint32_t mi, uint64_t temp_mod) {
                assert (mi < M_inc);

                test_stats.m_stops_interval += 1;

                /**
                 * x = (SL - m * K) % prime
                 *     Computed as
                 * x =  2*SL - ((SL + m*K) % prime)
                 *     =  SL - m * K
                 *     Requires prime > 2*SL
                 * x = (base_r * (M_start + mi) + SL) % prime;
                 */
                assert( temp_mod <= 2*SL );
                int32_t x = 2*SL - temp_mod;

                // Filters ~80% more (coprime to D)
                uint32_t m = (M_start + mi);
                uint32_t m_mod2310 = m % 2310;
                if (!is_m_coprime2310[m_mod2310])
                    return;

                /**
                 * Check if (m * K + x) has any small factors
                 * Filters 77.2%, Leaves 22.8%
                 *
                 * want positive mod so correct x - SL to x + ...
                 *
                 * XXX: Benchmark
                 * 2x int32 modulo should be faster than 1x int64 modulo
                 * Could save one addition by shifting is_coprime210 table by neg_SL_mod210
                 */
                uint32_t n_mod210 = ((K_mod210 * m_mod2310) + x + neg_SL_mod210) % 210;
                if (!is_coprime210[n_mod210])
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


#if METHOD2_WHEEL
                /**
                 * Filters ~80% (assuming 3 and 5 divide D)
                 * Requires a wide memory read (SL * 8 ~ 1/2 to 3 MB)
                 */
                uint32_t xii = caches.x_reindex_wheel[m % x_reindex_m_wheel][x];
#else
                uint32_t xii = caches.x_reindex[x];
#endif
                assert(xii > 0); // something wrong with n_mod210 calculation?

                // if coprime with K, try to toggle off factor.
                test_stats.large_prime_factors_interval += 1;

                int32_t mii = caches.m_reindex[mi];
                assert( mii >= 0 );

                if (use_lock) {
                    mutex_mi[mii].lock();
                    composite[mii][xii] = true;
                    mutex_mi[mii].unlock();
                } else {
                    composite[mii][xii] = true;
                }

                #ifdef GMP_VALIDATE_FACTORS
                validate_factor_m_k_x(stats, test, K, m, x, prime, SIEVE_LENGTH);
                #elif defined GMP_VALIDATE_LARGE_FACTORS
                if (prime > LARGE_PRIME_THRESHOLD)
                    validate_factor_m_k_x(stats, test, K, m, x, prime, SIEVE_LENGTH);
                #endif
            });
        }
        mpz_clear(test);


        // Normally this is inside the loop but not anymore
        #pragma omp critical
        {
            assert(test_stats.pi_interval > 0);  // Would cause infinite loop below.
            stats_to_process[end] = test_stats;
            if (config.verbose >= 3) {
                size_t queued = 0;
                for(auto&& kv : stats_to_process) queued += kv.second.pi_interval > 0;
                printf("\tmethod2_large_primes(%d) finished [%'ld, %'ld] (%ld in queue)\n",
                        omp_get_thread_num(), first, end, queued);
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
                if (temp.pi_interval > 0) {
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
    const uint32_t M_start = config.mstart;
    const uint32_t M_inc = config.minc;

    const uint32_t P = config.p;

    const uint32_t SIEVE_LENGTH = config.sieve_length;
    const uint32_t SL = SIEVE_LENGTH;
    // SIEVE_INTERVAL includes endpoints [-SL ... K ... SL]
    uint32_t SIEVE_INTERVAL = 2 * SIEVE_LENGTH + 1;

    const uint64_t MAX_PRIME = config.max_prime;

    mpz_t test, test2;
    mpz_init(test);
    mpz_init(test2);

    mpz_set_ui(test, MAX_PRIME);
    mpz_prevprime(test, test);

    uint64_t LAST_PRIME = mpz_get_ui(test);
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
    const Cached caches = setup_caches(config, K);
    const size_t valid_ms = caches.valid_ms;
    const uint32_t x_reindex_m_wheel = caches.x_reindex_m_wheel;

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

    // SMALL_THRESHOLD must handle all primes that can mark off two items in SIEVE_INTERVAL.
    assert( SMALL_THRESHOLD >= SIEVE_INTERVAL );
    assert( MEDIUM_THRESHOLD >= SMALL_THRESHOLD );
    assert( MEDIUM_THRESHOLD <= config.max_prime );

#if defined GMP_VALIDATE_LARGE_FACTORS && !defined GMP_VALIDATE_FACTORS
    // No overflow from gap_common.cpp checks
    const uint32_t M_end = M_start + M_inc;
    const uint64_t LARGE_PRIME_THRESHOLD = (1LL << 55) / M_end;
    if (LARGE_PRIME_THRESHOLD < LAST_PRIME && config.verbose >= 1) {
        printf("validating factors from primes > %ld\n", LARGE_PRIME_THRESHOLD);
    }
#endif

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
     * composite[mi][x] (0 <= mi < M_inc, -SL <= x <= SL) is re-indexed to
     *      composite[m_reindex[mi]][x_reindex[SL + x]]
     * m_reindex[mi] with (D, M + mi) > 0 are mapped to -1 (and must be handled by code)
     * x_reindex[x]  with (K, x) > 0 are mapped to 0 (and that bit is ignored)
     * x_reindex_wheel[x] same as x_reindex[x]
     */

    // <char> is faster (0-5%?) than <bool>, but uses 8x the memory.
    // Need to change here and in `save_unknowns_method2` signature.
    vector<bool> *composite = new vector<bool>[valid_ms];
    {
        int align_print = 0;
        size_t MB = 8 * 1024 * 1024;
        size_t guess = valid_ms * (count_coprime_sieve + 1) / MB;
        if (config.verbose >= 1) {
            align_print = printf("coprime m    %ld/%d,  ", valid_ms, M_inc);
            printf("coprime i     %ld/%d, ~%'ldMB\n",
                count_coprime_sieve / 2, SIEVE_LENGTH, guess);
        }

        if (x_reindex_m_wheel > 1) {
            // Update guess with first wheel count for OOM prevention check
            guess = valid_ms * (caches.x_reindex_wheel_count[1] + 1) / MB;
        }

        // Try to prevent OOM, check composite < 10GB allocation,
        // combined_sieve uses ~5-20% extra space for x_reindex_wheel + extra
        if (guess > (size_t) config.max_mem * 1024) {
            printf("combined_sieve expects to use %'ld MB which is greater than %d GB limit\n",
                    guess, config.max_mem);
            printf("\tAdd `--max-mem %ld` to skip this warning\n", (guess / 1024) + 1);
            exit(1);
        }

        size_t allocated = 0;
        for (size_t i = 0; i < valid_ms; i++) {
            int m_wheel = (M_start + caches.valid_mi[i]) % x_reindex_m_wheel;
            assert(gcd(m_wheel, x_reindex_m_wheel) == 1);

            // +1 reserves extra 0th entry for x_reindex[x] = 0
            composite[i].resize(caches.x_reindex_wheel_count[m_wheel] + 1, false);
            composite[i][0] = true;
            allocated += composite[i].size();
        }
        if (config.verbose >= 1 && x_reindex_m_wheel > 1) {
            printf("%*s", align_print, "");
            printf("coprime wheel %ld/%d, ~%'ldMB\n",
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

    const size_t THREADS = config.threads;

    if (1) { // Small Primes
        /**
         * NOTE: For primes <= SMALL_THRESHOLD, handle per m (with better memory locality)
         * This also avoids need to synchronize access to composite
         * It does make printing slightly harder (see awkward inner loop)
         */

        vector<int32_t> valid_mi_split[THREADS];
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
    }


    if (1) { // Medium Primes
        vector<uint32_t> coprime_X_split[THREADS];
        size_t coprime_X_count = caches.coprime_X.size();
        size_t per_thread = coprime_X_count / THREADS;
        per_thread -= per_thread % 8;

        /**
         * To avoid possible synchronization issues round coprime_X_split to a
         * multiples of 8 so vector<bool>[i]  never overlaps between threads
         */
        for (size_t t = 0; t < THREADS; t++) {
            bool is_final = (t+1) >= THREADS;
            size_t xi_start = t * per_thread;
            size_t xi_end = is_final ? coprime_X_count : xi_start + per_thread;
            assert( xi_start % 8 == 0);
            assert( (xi_end - xi_start) % 8 == 0 || is_final );
            for (size_t xi = xi_start; xi < xi_end; xi++) {
                coprime_X_split[t].push_back(caches.coprime_X[xi]);
            }
        }

        /**
         * Using this parallelization has two side effects
         * 1. If THREADS is CORES+1 would wait 2x as long
         * 2. Most threads finish before stats thread
         *    Then sit idle waiting for final thread to finish.
         */
        #pragma omp parallel for num_threads(THREADS)
        for (size_t thread_i = 0; thread_i < THREADS; thread_i++) {
            // Helps keep printing in order
            std::this_thread::sleep_for(milliseconds(50 * thread_i));

            if (config.verbose + (THREADS > 1) >= 3) {
                printf("\tThread %ld method2_medium_primes(%'ld/%'ld)\n",
                        thread_i, coprime_X_split[thread_i].size(), coprime_X_count);
            }

            method2_medium_primes(config, stats, K,
                                  thread_i,
                                  caches,
                                  coprime_X_split[thread_i],
                                  SMALL_THRESHOLD, MEDIUM_THRESHOLD,
                                  composite);
        }
        if (MEDIUM_THRESHOLD >= stats.last_prime) {
            // Handle final print (if no large_prime will be run)
            method2_increment_print(
                stats.last_prime, valid_ms, composite, stats, config);
        }
        if (config.verbose >= 2) {
            printf("\tmethod2_medium_primes done\n");
        }
    }


    if (1) { // Large Primes
        // THREADS are handled inside function.
        method2_large_primes(
            config, stats,
            K,
            THREADS,
            caches,
            MEDIUM_THRESHOLD,
            composite);
    }


    { // Verify some computation.
        // Likely zeroed in the last interval, but needed if no printing
        stats.m_stops += stats.m_stops_interval;

        // TODO move into print method

        // See Merten's Third Theorem
        float expected_m_stops = (log(log(LAST_PRIME)) - log(log(MEDIUM_THRESHOLD))) * 2*SL * M_inc;
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
    mpz_clear(test);
    mpz_clear(test2);
}
