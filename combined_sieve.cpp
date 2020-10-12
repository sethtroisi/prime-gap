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
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include <gmp.h>

#include "gap_common.h"
#include "modulo_search.h"

using std::cout;
using std::endl;
using std::pair;
using std::map;
using std::vector;
using namespace std::chrono;



// Validate each composite factor (slow)
#ifndef GMP_VALIDATE_FACTORS
    #define GMP_VALIDATE_FACTORS 0
#endif


// Tweaking this doesn't seem to method1 much.
// method2 is more sensative and set it's own.
#define SMALL_PRIME_LIMIT_METHOD1       400'000

void set_defaults(struct Config& config);
void prime_gap_search(const struct Config config);
void prime_gap_parallel(const struct Config config);


int main(int argc, char* argv[]) {
    Config config = Args::argparse(argc, argv);

    if (config.verbose >= 2) {
        printf("\tCompiled with GMP %d.%d.%d\n\n",
            __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);
    }

    set_defaults(config);

    if (config.save_unknowns == 0) {
        cout << "Must set --save-unknowns" << endl;
        exit(1);
    }

    if (config.valid == 0) {
        Args::show_usage(argv[0]);
        exit(1);
    }

    setlocale(LC_NUMERIC, "");
    if (config.verbose >= 0) {
        printf("\n");
        printf("Testing m * %u#/%u, m = %ld + [0, %'ld)\n",
            config.p, config.d, config.mstart, config.minc);
    }
    setlocale(LC_NUMERIC, "C");

    if ((config.max_prime > 2'000'000'000) &&
            (config.max_prime / config.minc > 100'000) &&
            (config.p <= 10000)) {
        printf("\tmax_prime(%ldB) is probably too large\n",
            config.max_prime / 1'000'000'000);
    }

    if (config.save_unknowns) {
        std::string fn = Args::gen_unknown_fn(config, ".txt");
        std::ifstream f(fn);
        if (f.good()) {
            printf("\nOutput file '%s' already exists\n", fn.c_str());
            exit(1);
        }
    }

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
        double N_log = K_log + log(config.mstart);
        double prob_prime_coprime_P = 1 / N_log - 1 / (N_log * N_log);

        // factors of K = P#/D
        vector<uint32_t> K_primes = get_sieve_primes(config.p);
        {
            // Adjust for prob_prime for no primes <= P
            for (auto prime : K_primes) {
                prob_prime_coprime_P /= (1 - 1.0 / prime);
            }

            // Remove any factors of D
            K_primes.erase(
                std::remove_if(K_primes.begin(), K_primes.end(),
                   [&](uint32_t p){ return config.d % p == 0; }));
        }

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
            config.max_prime =  10'000'000'000;
        } else {
            config.max_prime =   1'000'000'000;
        }

        if (config.verbose >= 0) {
            printf("AUTO SET: max_prime (log(K) = ~%.0f): %ld\n",
                K_log, config.max_prime);
        }
    }

    mpz_clear(K);
}


void insert_range_db(
        const struct Config config,
        long num_rows,
        float time_sieve) {

    DB db_helper(config.search_db.c_str());
    sqlite3 *db = db_helper.get_db();

    const uint64_t rid = db_helper.config_hash(config);
    char sSQL[300];
    sprintf(sSQL,
        "INSERT INTO range(rid, P, D, m_start, m_inc,"
                          "sieve_length, max_prime,"
                          "min_merit,"
                          "num_m, num_remaining,"
                          "time_sieve)"
         "VALUES(%ld,  %d,%d, %ld,%ld,"
                " %d,%ld, %.3f,"
                "%ld,%ld,  %.2f)"
         "ON CONFLICT(rid) DO UPDATE SET time_sieve=%.2f",
            rid,  config.p, config.d, config.mstart, config.minc,
            config.sieve_length, config.max_prime,
            config.min_merit,
            num_rows, num_rows,
            time_sieve, time_sieve);

    char *zErrMsg = 0;
    int rc = sqlite3_exec(db, sSQL, NULL, NULL, &zErrMsg);
    if (rc != SQLITE_OK) {
        printf("\nrange INSERT failed %d: %s\n",
            rc, sqlite3_errmsg(db));
        exit(1);
    }
}

// Method1


void save_unknowns_method1(
        std::ofstream &unknown_file,
        const uint64_t mi, int unknown_l, int unknown_u,
        const unsigned int SL, const vector<char> composite[]) {

    unknown_file << mi;
    unknown_file << " : -" << unknown_l << " +" << unknown_u << " |";

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


void prime_gap_search(const struct Config config) {
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
        setlocale(LC_NUMERIC, "");
        printf("\n");
        printf("sieve_length: 2x %'d\n", config.sieve_length);
        printf("max_prime:       %'ld\n", MAX_PRIME);
        printf("\n");
        setlocale(LC_NUMERIC, "C");
    }

    // ----- Generate primes under SMALL_PRIME_LIMIT_METHOD1
    vector<uint32_t> small_primes = get_sieve_primes(SMALL_PRIME_LIMIT_METHOD1);

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

    // To save space, only save remainder for primes that divide ANY m in range.
    // This helps with memory usage when MAX_PRIME >> SL * MINC;
    std::vector<p_and_r> *large_prime_queue = new vector<p_and_r>[M_inc];
    {
        size_t pr_pi = 0;
        if (config.verbose >= 0) {
            printf("\tCalculating first m each prime divides\n");
        }

        // large_prime_queue size can be approximated by
        // https://en.wikipedia.org/wiki/Meissel–Mertens_constant

        // Print "."s during, equal in length to 'Calculat...'
        size_t print_dots = 38;

        // Lookup primepi for common max_prime values.
        const size_t expected_primes = common_primepi.count(MAX_PRIME) ?
            common_primepi[MAX_PRIME] :
            1.04 * MAX_PRIME / log(MAX_PRIME);

        long first_m_sum = 0;
        double expected_large_primes = 0;

        if (config.verbose >= 0) {
            cout << "\t";
        }
        size_t pi = 0;
        get_sieve_primes_segmented_lambda(MAX_PRIME, [&](const uint64_t prime) {
            pi += 1;
            if (config.verbose >= 0 && (pi * print_dots) % expected_primes < print_dots) {
                cout << "." << std::flush;
            }

            // Big improvement over surround_prime is reusing this for each m.
            const uint64_t base_r = mpz_fdiv_ui(K, prime);

            if (prime <= SMALL_PRIME_LIMIT_METHOD1) {
                prime_and_remainder.emplace_back(prime, base_r);
                pr_pi += 1;
                return true;
            }

            expected_large_primes += (2.0 * SL + 1) / prime;

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
            if (mi == M_inc) return true;

            assert (mi < M_inc);

            // (M_start + mi) * last_prime < int64 (checked in argparse)
            uint64_t first = (base_r * (M_start + mi) + SL) % prime;
            assert( first <= 2*SL );

            //assert ( gcd(M + mi, D) == 1 );

            large_prime_queue[mi].emplace_back(prime, base_r);
            pr_pi += 1;

            s_large_primes_rem += 1;
            first_m_sum += mi;

            return true;
        });
        if (config.verbose >= 0) {
            cout << endl;
        }

        assert(prime_and_remainder.size() == small_primes.size());
        if (config.verbose >= 1) {
            printf("\tSum of m1: %ld\n", first_m_sum);
            setlocale(LC_NUMERIC, "");
            if (pi == expected_primes) {
                printf("\tPrimePi(%ld) = %ld\n", MAX_PRIME, pi);
            } else {
                printf("\tPrimePi(%ld) = %ld guessed %ld\n", MAX_PRIME, pi, expected_primes);
            }

            printf("\t%ld primes not needed (%.1f%%)\n",
                (pi - SMALL_PRIME_PI) - pr_pi,
                100 - (100.0 * pr_pi / (pi - SMALL_PRIME_PI)));

            float mertens3 = log(log(MAX_PRIME)) - log(log(SMALL_PRIME_LIMIT_METHOD1));
            printf("\texpected large primes/m: %.1f (theoretical: %.1f)\n",
                expected_large_primes, (2 * SL + 1) * mertens3);
            setlocale(LC_NUMERIC, "C");
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
    long  s_t_unk_low = 0;
    long  s_t_unk_hgh = 0;
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

            for (size_t d = modulo; d <= SIEVE_LENGTH; d += pr.first) {
                composite[0][d] = true;
            }

            // Not technically correct but fine to skip modulo == 0
            int first_negative = pr.first - modulo;
            for (size_t d = first_negative; d <= SIEVE_LENGTH; d += pr.first) {
                composite[1][d] = true;
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

            if (GMP_VALIDATE_FACTORS) {
                mpz_mul_ui(test, K, m);
                assert(modulo == mpz_fdiv_ui(test, prime));
            }

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

        // 2-3% of runtime, could be optimized into save_unknowns loop..
        int unknown_l = std::count(composite[0].begin(), composite[0].end(), false);
        int unknown_u = std::count(composite[1].begin(), composite[1].end(), false);
        s_total_unknown += unknown_l + unknown_u;
        s_t_unk_low += unknown_l;
        s_t_unk_hgh += unknown_u;

        // Save unknowns
        if (config.save_unknowns) {
            save_unknowns_method1(
                unknown_file,
                mi, unknown_l, unknown_u,
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
                    m, unknown_l, unknown_u);

            if (config.verbose + is_last >= 1) {
                // Stats!
                printf("\t    intervals %-10ld (%.2f/sec, with setup per m: %.2g)  %.0f seconds elapsed\n",
                        s_tests, s_tests / secs, t_secs / s_tests, secs);
                printf("\t    unknowns  %-10ld (avg: %.2f), %.2f%% composite  %.2f <- %% -> %.2f%%\n",
                        s_total_unknown, s_total_unknown / ((double) s_tests),
                        100.0 * (1 - s_total_unknown / ((2.0 * SIEVE_LENGTH + 1) * s_tests)),
                        100.0 * s_t_unk_low / s_total_unknown,
                        100.0 * s_t_unk_hgh / s_total_unknown);
                printf("\t    large prime remaining: %d (avg/test: %ld)\n",
                        s_large_primes_rem, s_large_primes_tested / s_tests);
            }
        }
    }

    {
        auto s_stop_t = high_resolution_clock::now();
        double   secs = duration<double>(s_stop_t - s_setup_t).count();
        insert_range_db(config, s_tests, secs);
    }

    // Should be cleaning up after self.
    for(uint32_t mi = 0; mi < M_inc; mi++)  {
        assert( large_prime_queue[mi].size() == 0 );
    }

    // ----- cleanup

    delete[] large_prime_queue;
    mpz_clear(K);
    mpz_clear(test);
}


// Method 2


void save_unknowns_method2(
        const struct Config config,
        const vector<int32_t> &valid_mi,
        const vector<int32_t> &m_reindex,
        const vector<uint32_t> &i_reindex,
        const vector<bool> *composite) {

    // ----- Open and Save to Output file
    std::ofstream unknown_file;
    {
        std::string fn = Args::gen_unknown_fn(config, ".txt");
        printf("\nSaving unknowns to '%s'\n", fn.c_str());
        unknown_file.open(fn, std::ios::out);
        assert( unknown_file.is_open() ); // Can't open save_unknowns file
    }

    const uint32_t M_start = config.mstart;
    const uint32_t D = config.d;
    const uint32_t SL = config.sieve_length;

    const size_t count_coprime_sieve = *std::max_element(i_reindex.begin(), i_reindex.end());

    for (uint64_t mi : valid_mi) {
        assert(gcd(M_start + mi, D) == 1);
        int32_t mii = m_reindex[mi];
        assert( mii >= 0 );

        const auto& comp = composite[mii];

        // composite[0] isn't a real entry.
        const size_t size_side = count_coprime_sieve / 2;
        auto real_begin = comp.begin() + 1;
        size_t unknown_l = std::count(real_begin, real_begin + size_side, false);
        size_t unknown_u = std::count(real_begin + size_side, comp.end(), false);

        unknown_file << mi << " : -" << unknown_l << " +" << unknown_u << " |";
        for (int d = 0; d <= 1; d++) {
            char prefix = "-+"[d];
            size_t found = 0;

            for (size_t i = 1; i <= SL; i++) {
                int a = SL + (2*d - 1) * i;
                if (!comp[i_reindex[a]]) {
                    unknown_file << " " << prefix << i;
                    found += 1;
                }
            }
            if (d == 0) {
                unknown_file << " |";
                assert( found == unknown_l );
            } else {
                assert( found == unknown_u );
            }
        }
        unknown_file << "\n";
    }
}


bool g_control_c = false;
void signal_callback_handler(int signum) {
    if (g_control_c) {
        cout << "Caught 2nd CTRL+C stopping now." << endl;
        exit(2);
    } else {
       cout << "Caught CTRL+C stopping and saving after next interval " << endl;
       g_control_c = true;
    }
}


// Would be nice to pass const but CTRL+C handler changes max_prime
void prime_gap_parallel(struct Config config) {
    // Method2
    const uint32_t M_start = config.mstart;
    const uint32_t M_inc = config.minc;
    const uint32_t P = config.p;
    const uint32_t D = config.d;

    const uint32_t SIEVE_LENGTH = config.sieve_length;
    const uint32_t SL = SIEVE_LENGTH;

    const uint64_t MAX_PRIME = config.max_prime;

    mpz_t test;
    mpz_init(test);

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

    // Also prints estimate if verbose >= 2
    const double prp_time_est = prp_time_estimate_composite(N_log, config.verbose);

    // ----- Allocate memory
    // SIEVE_INTERVAL includes endpoints [-SL ... K ... SL]
    uint32_t SIEVE_INTERVAL = 2 * SIEVE_LENGTH + 1;

    vector<int32_t> valid_mi;
    vector<int32_t> m_reindex(M_inc, -1);
    size_t valid_ms = 0;
    {
        for (uint32_t mi = 0; mi < M_inc; mi++) {
            if (gcd(M_start + mi, D) == 1) {
                m_reindex[mi] = valid_ms;
                valid_mi.push_back(mi);
                valid_ms++;
            }
        }
    }
    assert(valid_ms == valid_mi.size());

    // which [i] are coprime to K
    vector<char> coprime_composite(SIEVE_INTERVAL+1, 1);
    // reindex composite[m][i]
    vector<uint32_t> i_reindex(SIEVE_INTERVAL+1, 0);
    {
        for (uint32_t prime : P_primes) {
            if (D % prime != 0) {
                uint32_t first = SIEVE_LENGTH % prime;
                if (GMP_VALIDATE_FACTORS) {
                    mpz_set(test, K);
                    mpz_sub_ui(test, test, SIEVE_LENGTH);
                    mpz_add_ui(test, test, first);
                    assert( 0 == mpz_fdiv_ui(test, prime) );
                }

                assert( 0 <= first && first < prime );
                assert( (SIEVE_LENGTH - first) % prime == 0 );

                for (size_t d = first; d <= SIEVE_INTERVAL; d += prime) {
                    coprime_composite[d] = 0;
                }
            }
        }
        // Center should be marked composite by every prime.
        assert(coprime_composite[SL] == 0);

        size_t count = 1;
        for (size_t i = 0; i < SIEVE_INTERVAL; i++) {
            if (coprime_composite[i] > 0) {
                i_reindex[i] = count++;
            }
        }
    }
    const size_t count_coprime_sieve = *std::max_element(i_reindex.begin(), i_reindex.end());
    assert( count_coprime_sieve % 2 == 0 );

    // Rough math is MULT  vs  log2(MULT) * (M_inc/valid_ms)
    const float SMALL_MULT = std::max(8.0, log(8) * M_inc / valid_ms);
    const uint64_t SMALL_THRESHOLD = SMALL_MULT * 2 * SIEVE_LENGTH;

    // SMALL_THRESHOLD mult deal with all primes that can mark off two items in SIEVE_LENGTH.
    assert( SMALL_THRESHOLD > (2 * SIEVE_LENGTH + 1)  );

    if (config.verbose >= 1) {
        setlocale(LC_NUMERIC, "");
        printf("sieve_length: 2x %'d\n", config.sieve_length);
        printf("max_prime:       %'ld   small_threshold:  %'ld (%.1f x SL)\n",
            config.max_prime, SMALL_THRESHOLD, 2 * SMALL_MULT);
        //printf("last prime :  %'ld\n", LAST_PRIME);
        setlocale(LC_NUMERIC, "C");
    }

    /**
     * Much space is saved via a reindexing scheme
     * composite[mi][x] (0 <= mi < M_inc, -SL <= x <= SL) is reindexed to
     *      composite[m_reindex[mi]][i_reindex[SL + x]]
     * m_reindex[mi] with (D, M + mi) > 0 are mapped to -1 (and must be handled by code)
     * i_reindex[x]  with (K, x) > 0 are mapped to 0 (and that bit is ignored)
     */

    // <char> is faster (0-5%?) than <bool>, but uses 8x the memory.
    // Need to change here and in `save_unknowns_method2` signature.
    vector<bool> *composite = new vector<bool>[valid_ms];
    {
        if (config.verbose >= 1) {
            printf("coprime m    %ld/%d,  coprime i %ld/%d,  ~%'ldMB\n",
                valid_ms, M_inc, count_coprime_sieve / 2, SIEVE_LENGTH,
                valid_ms * count_coprime_sieve / 8 / 1024 / 1024);
            printf("\n");
        }
        for (size_t i = 0; i < valid_ms; i++) {
            // Improve this setup (segfaults at ~0.5GB required)
            // +1 reserves extra 0th entry for i_reindex[x] = 0
            composite[i].resize(count_coprime_sieve + 1, false);
            composite[i][0] = true;
        };
    }

    // Used for various stats
    auto  s_start_t = high_resolution_clock::now();
    auto  s_interval_t = high_resolution_clock::now();
    long  s_total_unknowns = SIEVE_INTERVAL * valid_ms;
    long  s_prime_factors = 0;
    long  s_small_prime_factors_interval = 0;
    long  s_large_prime_factors_interval = 0;
    long  m_stops_interval = 0;

    uint64_t  s_next_print = 0;
    uint64_t  next_mult = SMALL_THRESHOLD <= 10000 ? 10000 : 100000;
    double s_prp_needed = 1 / prob_prime;

    // Setup CTRL+C catcher
    signal(SIGINT, signal_callback_handler);

    // Note: Handling small primes here had better localized memory access
    // But wasn't worth the extra code IMHO.

    const int K_mod3 = mpz_fdiv_ui(K, 3); // K % 3
    const int K_mod5 = mpz_fdiv_ui(K, 5); // K % 5
    const int K_mod7 = mpz_fdiv_ui(K, 7); // K % 7
    const int D_mod2 = D % 2 == 0;
    const int D_mod3 = D % 3 == 0;
    const int D_mod5 = D % 5 == 0;
    const int D_mod7 = D % 7 == 0;

    size_t pi = 0;
    size_t pi_interval = 0;
    get_sieve_primes_segmented_lambda(MAX_PRIME, [&](const uint64_t prime) {
        pi_interval += 1;
        // Big improvement over surround_prime is reusing this for each m.
        const uint64_t base_r = mpz_fdiv_ui(K, prime);

        if (prime <= SMALL_THRESHOLD) {
            // Handled by coprime_composite above
            if (D % prime != 0 && prime <= P) {
                return true;
            }

            for (uint32_t mi : valid_mi) {
                int32_t mii = m_reindex[mi];
                assert(mii >= 0);

                uint64_t m = M_start + mi;
                uint64_t modulo = (base_r * m) % prime;

                // flip = (m * K - SL) % prime
                uint32_t flip = modulo + prime - ((SIEVE_LENGTH+1) % prime);
                if (flip >= prime) flip -= prime;

                uint32_t first = prime - flip - 1;
                assert( first < prime );

                uint32_t shift = prime;
                if (prime > 2) {
                    bool centerOdd = ((D & 1) == 0) && (m & 1);
                    bool lowIsEven = centerOdd == (SIEVE_LENGTH & 1);
                    bool evenFromLow = (first & 1) == 0;
                    bool firstIsEven = lowIsEven == evenFromLow;

                    if (GMP_VALIDATE_FACTORS) {
                        mpz_mul_ui(test, K, M_start + mi);
                        mpz_sub_ui(test, test, SIEVE_LENGTH);
                        assert( (mpz_even_p(test) > 0) == lowIsEven );

                        mpz_add_ui(test, test, first);
                        assert( (mpz_even_p(test) > 0) == firstIsEven );

                        assert( 0 == mpz_fdiv_ui(test, prime) );
                    }

                    if (firstIsEven) {
                        // divisible by 2 move to next multiple (an odd multiple)

                        assert( (first >= SIEVE_INTERVAL) || composite[mii][i_reindex[first]] );
                        first += prime;
                    }

                    // Don't need to count cross off even multiples.
                    shift = 2*prime;
                }

                for (size_t d = first; d < SIEVE_INTERVAL; d += shift) {
                    composite[mii][i_reindex[d]] = true;
                    s_small_prime_factors_interval += 1;
                }
            }
        } else {
            modulo_search_euclid_all_small(M_start, M_inc, SL, prime, base_r, [&](
                        const uint32_t mi, uint64_t first) {
                assert (mi < M_inc);

                m_stops_interval += 1;

                // With D even, (ms + mi) must be odd (or share a factor of 2)
                // Helps avoid wide memory read
                if (((D & 1) == 0) && ((M_start & 1) == (mi & 1)))
                    return;

                int32_t mii = m_reindex[mi];
                if (mii < 0)
                    return;

                // Returning first from modulo_search_euclid_all_small is
                // slightly faster on benchmarks, and slightly faster here

                // first = (SL - m * K) % prime
                //     Computed as
                // first =  2*SL - ((SL + m*K) % prime)
                //       =  SL - m * K
                //     Requires prime > 2*SL
                //uint64_t first = (base_r * (M_start + mi) + SL) % prime;
                assert( first <= 2*SL );
                first = 2*SL - first;

                if (GMP_VALIDATE_FACTORS) {
                    mpz_mul_ui(test, K, M_start + mi);
                    mpz_sub_ui(test, test, SIEVE_LENGTH);
                    mpz_add_ui(test, test, first);
                    uint64_t mod = mpz_fdiv_ui(test, prime);
                    assert( mod == 0 );
                }

                int64_t dist = first - SIEVE_LENGTH;
                uint32_t m = M_start + mi;
                if (D_mod2 && (dist & 1))
                    return;
                if (D_mod3 && ((dist + K_mod3 * m) % 3 == 0))
                    return;
                if (D_mod5 && ((dist + K_mod5 * m) % 5 == 0))
                    return;
                if (D_mod7 && ((dist + K_mod7 * m) % 7 == 0))
                    return;

                if (!coprime_composite[first]) {
                    return;
                }

                composite[mii][i_reindex[first]] = true;
                s_large_prime_factors_interval += 1;
            });
        }

        if (prime >= s_next_print) {
            if (prime >= s_next_print) {
                if (s_next_print == 5 * next_mult) {
                    next_mult = 10 * next_mult;
                    s_next_print = 0;
                }
                s_next_print += next_mult;
                s_next_print = std::min(s_next_print, LAST_PRIME);
            }

            auto   s_stop_t = high_resolution_clock::now();
            // total time, interval time
            double     secs = duration<double>(s_stop_t - s_start_t).count();
            double int_secs = duration<double>(s_stop_t - s_interval_t).count();

            // See THEORY.md
            double prob_prime_after_sieve = prob_prime * log(prime) * exp(GAMMA);

            uint64_t t_total_unknowns = 0;
            for (size_t i = 0; i < valid_ms; i++) {
                t_total_unknowns += std::count(composite[i].begin(), composite[i].end(), false);
            }
            uint64_t new_composites = s_total_unknowns - t_total_unknowns;
            double skipped_prp = 2 * valid_ms * (s_prp_needed - 1/prob_prime_after_sieve);

            pi += pi_interval;
            s_prime_factors += s_small_prime_factors_interval;
            s_prime_factors += s_large_prime_factors_interval;

            bool is_last = (prime == LAST_PRIME) || g_control_c;

            setlocale(LC_NUMERIC, "");
            if (config.verbose + is_last >= 1) {
                printf("%'-10ld (primes %'ld/%'ld)\t(seconds: %.2f/%-.1f | per m: %.2g)\n",
                    prime,
                    pi_interval, pi,
                    int_secs, secs,
                    secs / valid_ms);
            }

            if ((config.verbose + 2*is_last + (prime > 1e9)) >= 2) {
                printf("\tfactors  %'14ld \t"
                       "(interval: %'ld avg m/large_prime interval: %.1f)\n",
                    s_prime_factors,
                    s_small_prime_factors_interval + s_large_prime_factors_interval,
                    1.0 * m_stops_interval / pi_interval);
                // count_coprime_sieve * valid_ms also makes sense but leads to smaller numbers
                printf("\tunknowns %'9ld/%-5ld\t"
                       "(avg/m: %.2f) (composite: %.2f%% +%.3f%% +%'ld)\n",
                    t_total_unknowns, valid_ms,
                    1.0 * t_total_unknowns / valid_ms,
                    100.0 - 100.0 * t_total_unknowns / (SIEVE_INTERVAL * valid_ms),
                    100.0 * new_composites / (SIEVE_INTERVAL * valid_ms),
                    new_composites);

                printf("\t~ 2x %.2f PRP/m\t\t"
                       "(~ %4.1f skipped PRP => %.1f PRP/seconds)\n",
                    1 / prob_prime_after_sieve,
                    skipped_prp,
                    skipped_prp / int_secs);

                double run_prp_mult = int_secs / (prp_time_est * skipped_prp);
                if (run_prp_mult > 4) {
                    printf("\t\tEstimated ~%.1fx faster to just run PRP now (CTRL+C to stop sieving)\n",
                        run_prp_mult);
                }

                printf("\n");

                s_total_unknowns = t_total_unknowns;
                s_interval_t = s_stop_t;
                s_prp_needed = 1 / prob_prime_after_sieve;

                s_small_prime_factors_interval = 0;
                s_large_prime_factors_interval = 0;
                m_stops_interval = 0;
                pi_interval = 0;
            }
            setlocale(LC_NUMERIC, "C");

            // if is_last would truncate .max_prime by 1 million
            if (g_control_c && (prime != LAST_PRIME)) {
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

                return false;
            }

            #ifdef SAVE_INCREMENTS
            if (config.save_unknowns && prime > 1e8 && prime != LAST_PRIME) {
                // reset unknown_filename if cached;
                config.unknown_filename = "";
                size_t old = config.max_prime;
                config.max_prime = prime - (prime % 1'000'000);
                save_unknowns_method2(
                    config,
                    valid_mi, m_reindex, i_reindex,
                    composite);
                config.max_prime = old;
            }
            #endif // SAVE_INCREMENTS
        }

        return true;
    });

    // TODO assert s_prime_factors is close to (2 * SL * valid_m) * (log(log(SL)) + M_c)

    if (config.save_unknowns) {
        save_unknowns_method2(
            config,
            valid_mi, m_reindex, i_reindex,
            composite);

        auto s_stop_t = high_resolution_clock::now();
        double   secs = duration<double>(s_stop_t - s_start_t).count();
        insert_range_db(config, valid_mi.size(), secs);
    }

    delete[] composite;

    mpz_clear(K);
    mpz_clear(test);
}

