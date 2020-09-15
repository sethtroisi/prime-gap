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

// Tweaking this doesn't seem to method1 much.
// method2 seems very sensative and is controlled elsewhere.
#define SIEVE_SMALL       400'000

void set_defaults(struct Config& config);
void prime_gap_search(const struct Config config);
void prime_gap_parallel(const struct Config config);


int main(int argc, char* argv[]) {
    Config config = argparse(argc, argv);

    if (config.verbose >= 2) {
        printf("\tCompiled with GMP %d.%d.%d\n\n",
            __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);
    }

    set_defaults(config);

    if (config.save_unknowns == 0) {
        cout << "Must set --save-unknowns" << endl;
        return 1;
    }

    if (config.valid == 0) {
        show_usage(argv[0]);
        return 1;
    }

    setlocale(LC_NUMERIC, "");
    if (config.verbose >= 0) {
        printf("\n");
        printf("Testing m * %u#/%u, m = %ld + [0, %'ld)\n",
            config.p, config.d, config.mstart, config.minc);
    }
    setlocale(LC_NUMERIC, "C");

    if ((config.sieve_range > 2'000'000'000) &&
            (config.sieve_range / config.minc > 100'000) &&
            (config.p <= 10000)) {
        printf("\tsieve_range(%ldB) is probably too large\n",
            config.sieve_range / 1'000'000'000);
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
        double prob_prime_coprime = 1 / (K_log + log(config.mstart));

        // factors of K = P#/D
        vector<uint32_t> K_primes = get_sieve_primes(config.p);
        {
            for (size_t pi = K_primes.size()-1; ; pi--) {
                prob_prime_coprime /= 1 - 1.0 / K_primes[pi];
                if (config.d % K_primes[pi] == 0) {
                    K_primes.erase(K_primes.begin() + pi);
                }
                if (pi == 0) break;
            }
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
                    // Some factor of d will mark this off (for these centers) don't count it.
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
                double prob_gap_shorter = pow(1 - prob_prime_coprime, min_coprime);

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

    if (config.sieve_range == 0) {
        // each additional numbers removes unknowns / prime
        // and takes log2(prime / sieve_length) time

        // TODO improve this.
        // Potentially based on:
        //  sieve_length
        //  min_merit
        if (K_log >= 1500) {
            // Largest supported number right now
            config.sieve_range = 4'000'000'000;
            // 2020-02-09 tuning notes
            //  P 1627 SL=8192
            //      log(t) ~= 1600, ~= 4 primes in SL
            //       100M => 620/s | 293 unknowns
            //       200M => 549/s | 282 unknowns
            //       400M => 480/s | 273 unknowns
            //          81.6 PRP / test => 1.70s/test
            //       800M => 440/s | 264 unknowns
            //          78.6 PRP / test => 1.75s/test
            //      1600M => 440/s | 254 unknowns
            //          76.2 PRP / test => 1.78s/test
        } else {
            config.sieve_range =   1'000'000'000;
        }

        if (config.verbose >= 0) {
            printf("AUTO SET: sieve range (log(K) = ~%.0f): %ld\n",
                K_log, config.sieve_range);
        }
    }

    mpz_clear(K);
}


double prob_prime_and_stats(
        const struct Config& config, mpz_t &K, const vector<uint32_t> &primes) {

    int K_digits;
    double K_log;
    K_stats(config, K, &K_digits, &K_log);

    double prob_prime = 1 / (K_log + log(config.mstart));

    // From Mertens' 3rd theorem
    double unknowns_after_sieve = 1 / (log(config.sieve_range) * exp(GAMMA));
    double prob_prime_after_sieve = prob_prime / unknowns_after_sieve;

    size_t count_coprime = 0;
    double prob_prime_coprime = 0;
    double prob_gap_hypothetical = prop_gap_larger(
        config, prob_prime, &prob_prime_coprime, &count_coprime);

    if (config.verbose >= 2) {
        printf("\n");
        printf("\tavg %.0f left from %.3f%% of %u interval with %ldM sieve\n",
                count_coprime * (unknowns_after_sieve / prob_prime_coprime),
                100 * unknowns_after_sieve,
                config.sieve_length, config.sieve_range/1'000'000);
        printf("\t%.3f%% of %d digit numbers are prime\n",
                100 * prob_prime, K_digits);
        printf("\t%.3f%% of tests should be prime (%.1fx speedup)\n",
                100 * prob_prime_after_sieve, 1 / unknowns_after_sieve);
        printf("\t~2x%.1f = %.1f PRP tests per m\n",
                1 / prob_prime_after_sieve, 2 / prob_prime_after_sieve);
        printf("\tsieve_length=%d is insufficient ~~%.3f%% of time\n",
                config.sieve_length, 100 * prob_gap_hypothetical);
        printf("\n");
    }

    return prob_prime;
}


void prime_gap_search(const struct Config config) {
    const uint64_t M_start = config.mstart;
    const uint64_t M_inc = config.minc;
    //const uint64_t P = config.p;
    const uint64_t D = config.d;

    const unsigned int SIEVE_LENGTH = config.sieve_length;
    const unsigned int SL = SIEVE_LENGTH;

    const uint64_t SIEVE_RANGE = config.sieve_range;

    mpz_t test;
    mpz_init(test);

    if (config.verbose >= 2) {
        setlocale(LC_NUMERIC, "");
        printf("\n");
        printf("sieve_length: 2x %'d\n", config.sieve_length);
        printf("sieve_range:  %'ld\n", config.sieve_range);
        printf("\n");
        setlocale(LC_NUMERIC, "C");
    }

    // ----- Generate primes under SIEVE_RANGE.
    vector<uint32_t> small_primes = get_sieve_primes(SIEVE_SMALL);

    // ----- Merit Stuff
    mpz_t K;
    prob_prime_and_stats(config, K, small_primes);


    // ----- Sieve stats
    const size_t SIEVE_SMALL_PRIME_PI = small_primes.size();
    {
        // SIEVE_SMALL deals with all primes that can mark off two items in SIEVE_LENGTH.
        assert( SIEVE_SMALL > 2 * SIEVE_LENGTH );
        if (config.verbose >= 1) {
            printf("\tUsing %'ld primes for SIEVE_SMALL(%'d)\n\n",
                SIEVE_SMALL_PRIME_PI, SIEVE_SMALL);
        }
        assert( small_primes[SIEVE_SMALL_PRIME_PI-1] < SIEVE_SMALL);
        assert( small_primes[SIEVE_SMALL_PRIME_PI-1] + 200 > SIEVE_SMALL);
    }

    auto  s_setup_t = high_resolution_clock::now();

    // ----- Allocate memory for a handful of utility functions.

    // Remainders of (p#/d) mod prime
    typedef pair<uint64_t,uint64_t> p_and_r;
    vector<p_and_r> prime_and_remainder;

    // Big improvement over surround_prime is avoiding checking each large prime.
    // vector<m, vector<pi>> for large primes that only rarely divide a sieve
    int s_large_primes_rem = 0;

    // To save space, only save remainder for primes that divide ANY m in range.
    // This helps with memory usage when SIEVE_RANGE >> SL * MINC;
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

        size_t expected_primes = 1.01 * SIEVE_RANGE / log(SIEVE_RANGE);

        long first_m_sum = 0;
        double expected_large_primes = 0;

        if (config.verbose >= 0) {
            cout << "\t";
        }
        size_t pi = 0;
        get_sieve_primes_segmented_lambda(SIEVE_RANGE, [&](const uint64_t prime) {
            pi += 1;
            if (config.verbose >= 0 && (pi * print_dots) % expected_primes < print_dots) {
                cout << "." << std::flush;
            }

            // Big improvement over surround_prime is reusing this for each m.
            const uint64_t base_r = mpz_fdiv_ui(K, prime);

            if (prime <= SIEVE_SMALL) {
                prime_and_remainder.emplace_back(prime, base_r);
                pr_pi += 1;
                return true;
            }

            expected_large_primes += (2.0 * SL - 1) / prime;

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

            __int128 mult = (__int128) base_r * (M_start + mi) + (SL - 1);
            assert( mult % prime < (2*SL-1) );

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
            printf("\tPrimePi(%ld) = %ld guessed %ld\n", SIEVE_RANGE, pi, expected_primes);

            printf("\t%ld primes not needed (%.1f%%)\n",
                (pi - SIEVE_SMALL_PRIME_PI) - pr_pi,
                100 - (100.0 * pr_pi / (pi - SIEVE_SMALL_PRIME_PI)));
            printf("\texpected large primes/m: %.1f (theoretical: %.1f)\n",
                expected_large_primes,
                (2 * SL - 1) * (log(log(SIEVE_RANGE)) - log(log(SIEVE_SMALL))));
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
        std::string fn = gen_unknown_fn(config, ".txt");
        printf("\nSaving unknowns to '%s'\n", fn.c_str());
        unknown_file.open(fn, std::ios::out);
        assert( unknown_file.is_open() ); // Can't open save_unknowns file
    }


    // ----- Main sieve loop.
    if (config.verbose >= 1) {
        cout << "\nStarting m=" << M_start << "\n" << endl;
    }

    // vector<bool> uses bit indexing which is ~5% slower.
    vector<char> composite[2] = {
        vector<char>(SIEVE_LENGTH, 0),
        vector<char>(SIEVE_LENGTH, 0)
    };
    assert( composite[0].size() == SIEVE_LENGTH );
    assert( composite[1].size() == SIEVE_LENGTH );

    // Used for various stats
    long  s_tests = 0;
    auto  s_start_t = high_resolution_clock::now();
    long  s_total_unknown = 0;
    long  s_t_unk_low = 0;
    long  s_t_unk_hgh = 0;
    long  s_large_primes_tested = 0;

    uint64_t last_mi = M_inc - 1;
    for (; last_mi > 0 && gcd(M_start + last_mi, D) > 1; last_mi -= 1);
    assert(last_mi > 0 && last_mi < M_inc);

    for (uint64_t mi = 0; mi < M_inc; mi++) {
        uint64_t m = M_start + mi;
        if (gcd(m, D) > 1) {
            assert( large_prime_queue[mi].empty() );
            continue;
        }

        // Reset sieve array to unknown.
        std::fill_n(composite[0].begin(), SIEVE_LENGTH, 0);
        std::fill_n(composite[1].begin(), SIEVE_LENGTH, 0);
        // center is always composite.
        composite[0][0] = composite[1][0] = 1;

        // For small primes that we don't do trick things with.
        for (const auto& pr : prime_and_remainder) {
            const uint64_t modulo = (pr.second * m) % pr.first;
            //            const auto& [prime, remainder] = prime_and_remainder[pi];
            //            const uint64_t modulo = (remainder * m) % prime;

            for (size_t d = modulo; d < SIEVE_LENGTH; d += pr.first) {
                composite[0][d] = true;
            }

            // Not technically correct but fine to skip modulo == 0
            int first_negative = pr.first - modulo;
            for (size_t d = first_negative; d < SIEVE_LENGTH; d += pr.first) {
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

            if (0) {
                mpz_mul_ui(test, K, m);
                uint64_t mod = mpz_fdiv_ui(test, prime);
                assert( mod == modulo );
            }

            if (modulo < SIEVE_LENGTH) {
                // Just past a multiple
                composite[0][modulo] = true;
            } else {
                // Don't have to deal with 0 case anymore.
                int64_t first_positive = prime - modulo;
                assert( first_positive < SIEVE_LENGTH); // Bad next m!
                // Just before a multiple
                composite[1][first_positive] = true;
            }

            // Find next mi where primes divides part of SIEVE
            {
                uint64_t start = mi + 1;
                uint64_t next_mi = start + modulo_search_euclid_gcd(
                        M_start + start, D, M_inc - start, SL, prime, remainder);
                if (next_mi == M_inc) continue;

                //assert( (remainder * (M + next_mi) + (SL - 1)) % prime < (2*SL-1) );
                __int128 mult = (__int128) remainder * (M_start + next_mi) + (SL - 1);
                assert ( mult % prime <= (2 * SL - 1) );

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
            unknown_file << mi;
            unknown_file << " : -" << unknown_l << " +" << unknown_u << " |";

            for (int d = 0; d <= 1; d++) {
                char prefix = "-+"[d];

                for (size_t i = 1; i < SL; i++) {
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
                        100.0 * (1 - s_total_unknown / (2.0 * (SIEVE_LENGTH - 1) * s_tests)),
                        100.0 * s_t_unk_low / s_total_unknown,
                        100.0 * s_t_unk_hgh / s_total_unknown);
                printf("\t    large prime remaining: %d (avg/test: %ld)\n",
                        s_large_primes_rem, s_large_primes_tested / s_tests);
            }
        }
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


// Would be nice to pass const but CTRL+C handler changes sieve_range
void prime_gap_parallel(struct Config config) {
    // Method2
    const uint32_t M_start = config.mstart;
    const uint32_t M_inc = config.minc;
    const uint32_t P = config.p;
    const uint32_t D = config.d;

    const uint32_t SIEVE_LENGTH = config.sieve_length;
    const uint32_t SL = SIEVE_LENGTH;

    const uint64_t SIEVE_RANGE = config.sieve_range;

    // TODO: Would be nice to have a prev prime.
    uint64_t LAST_PRIME = SIEVE_RANGE;
    {
        bool prime;
        do {
            prime = true;
            if (LAST_PRIME % 2 == 0)
                LAST_PRIME--;

            for (uint64_t p = 3; p*p <= LAST_PRIME; p += 2) {
                if (LAST_PRIME % p == 0) {
                    prime = false;
                    LAST_PRIME--;
                    break;
                }
            }
        } while (prime == false);
    }

    // TODO increase this when valid_ms is very small
    // TODO seems not well calibrated right now, lower to 40x?
    const uint64_t SMALL_THRESHOLD = 50 * SIEVE_LENGTH;

    // SMALL_THRESHOLD deals with all primes that can mark off two items in SIEVE_LENGTH.
    assert( SMALL_THRESHOLD > 2 * SIEVE_LENGTH );

    mpz_t test;
    mpz_init(test);

    // ----- Generate primes for P
    const vector<uint32_t> P_primes = get_sieve_primes(P);
    assert( P_primes.back() == P);

    // ----- Sieve stats & Merit Stuff
    mpz_t K;
    double prob_prime = prob_prime_and_stats(config, K, P_primes);

    // ----- Allocate memory
    uint32_t SIEVE_INTERVAL = 2 * SIEVE_LENGTH - 1;

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
    vector<char> coprime_composite(SIEVE_INTERVAL, 1);
    // reindex composite[m][i]
    vector<uint32_t> i_reindex(SIEVE_INTERVAL, 0);
    {
        coprime_composite[0] = 0;
        for (uint32_t prime : P_primes) {
            if (D % prime != 0) {
                int first = (SIEVE_LENGTH-1) % prime;
                assert ( ((SIEVE_LENGTH-1) - first) % prime == 0 );
                for (size_t d = first; d < SIEVE_INTERVAL; d += prime) {
                    coprime_composite[d] = 0;
                }
            }
        }

        size_t count = 1;
        for (size_t i = 0; i < SIEVE_INTERVAL; i++) {
            if (coprime_composite[i] > 0) {
                i_reindex[i] = count++;
            }
        }
    }
    const size_t count_coprime_sieve = *std::max_element(i_reindex.begin(), i_reindex.end());

    if (config.verbose >= 1) {
        setlocale(LC_NUMERIC, "");
        printf("sieve_length: 2x %'d\n", config.sieve_length);
        printf("sieve_range: %'ld   small_threshold:  %'ld\n", config.sieve_range, SMALL_THRESHOLD);
        //printf("last prime :  %'ld\n", LAST_PRIME);
        setlocale(LC_NUMERIC, "C");
    }

    // <bool> is slower than <char>, but uses 1/8th the memory.
    vector<bool> composite[valid_ms];
    {
        for (size_t i = 0; i < valid_ms; i++) {
            // Improve this setup.
            composite[i].resize(count_coprime_sieve + 1, false);
        };
        if (config.verbose >= 1) {
            printf("coprime m    %ld/%d,  coprime i %ld/%d,  ~%'ldMB\n",
                valid_ms, M_inc, count_coprime_sieve / 2, SIEVE_LENGTH,
                valid_ms * count_coprime_sieve / 8 / 1024 / 1024);
            printf("\n");
        }
    }

    // Used for various stats
    auto  s_start_t = high_resolution_clock::now();
    auto  s_interval_t = high_resolution_clock::now();
    long  s_total_unknowns = SIEVE_INTERVAL * valid_ms;
    long  s_prime_factors = 0;
    long  s_small_prime_factors_interval = 0;
    long  s_large_prime_factors_interval = 0;
    uint64_t  s_next_print = 0;
    uint64_t  next_mult = SMALL_THRESHOLD <= 10000 ? 10000 : 100000;
    double s_prp_needed = 1 / prob_prime;

    // Setup CTRL+C catcher
    signal(SIGINT, signal_callback_handler);

    // Note: Handling small primes here had better localized memory access
    // But wasn't worth the extra code IMHO.

    size_t pi = 0;
    size_t pi_interval = 0;
    get_sieve_primes_segmented_lambda(SIEVE_RANGE, [&](const uint64_t prime) {
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

                uint32_t flip = modulo + prime - (SIEVE_LENGTH % prime);
                if (flip >= prime) flip -= prime;
                uint32_t first = prime - flip - 1;

                uint32_t shift = prime;
                if (prime > 2) {
                    //bool lowIsEven = ((SIEVE_LENGTH % 2) == 0) ^ ((D % 2) == 1) || ((m % 2) == 0);
                    //bool evenFrom = (first % 2) == 0;

                    bool lowIsEven = (D & 1) || ((m & 1) == 0);
                    bool evenFromCenter = (first & 1) == (SIEVE_LENGTH & 1);

                    if (lowIsEven ^ evenFromCenter) {
                        // Same parity (divisible by 2) move to next prime
                        assert( (first >= SIEVE_INTERVAL) ||
                                (composite[mii][i_reindex[first]] > 0) );
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
            modulo_search_euclid_all_small(M_start, M_inc, SL, prime, base_r, [&](const uint32_t mi) {
                int32_t mii = m_reindex[mi];
                if (mii < 0) {
                    return;
                }

                // TODO validate no overflow
                // TODO return first from lambda?
                uint64_t first = (base_r * (M_start + mi) + (SL-1)) % prime;
                assert( first < SIEVE_INTERVAL );
                first = SIEVE_INTERVAL - first - 1;
                assert( 0 <= first && first < SIEVE_INTERVAL );

                if (!coprime_composite[first]) {
                    return;
                }

                if (0) {
                    mpz_mul_ui(test, K, M_start + mi);
                    mpz_sub_ui(test, test, SIEVE_LENGTH-1);
                    mpz_add_ui(test, test, first);
                    uint64_t mod = mpz_fdiv_ui(test, prime);
                    assert( mod == 0 );
                }

                composite[mii][i_reindex[first]] = true;
                s_large_prime_factors_interval += 1;
            });
        }

        if (prime >= s_next_print) {
            if (prime >= s_next_print) {
                if (s_next_print == 5 * next_mult) {
                    next_mult = 2 * s_next_print;
                    s_next_print = 0;
                }
                s_next_print += next_mult;
                s_next_print = std::min(s_next_print, LAST_PRIME);
            }

            auto   s_stop_t = high_resolution_clock::now();
            // total time, interval time
            double     secs = duration<double>(s_stop_t - s_start_t).count();
            double int_secs = duration<double>(s_stop_t - s_interval_t).count();

            double unknowns_after_sieve = 1 / (log(prime) * exp(GAMMA));
            double prob_prime_after_sieve = prob_prime / unknowns_after_sieve;

            uint64_t t_total_unknowns = 0;
            for (size_t i = 0; i < valid_ms; i++) {
                t_total_unknowns += std::count(composite[i].begin(), composite[i].end(), false);
            }
            uint64_t saved_prp = s_total_unknowns - t_total_unknowns;
            double skipped_prp = valid_ms * (s_prp_needed - 1/prob_prime_after_sieve);

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
                printf("\tfactors  %'9ld \t\t(interval: %'ld, avg m/large_prime interval: %.1f)\n",
                    s_prime_factors,
                    s_small_prime_factors_interval + s_large_prime_factors_interval,
                    1.0 * s_large_prime_factors_interval / pi_interval);
                printf("\tunknowns %'9ld/%-5ld\t(avg/m: %.2f) (composite: %.2f%% +%.3f%%)\n",
                    t_total_unknowns, valid_ms,
                    1.0 * t_total_unknowns / valid_ms,
                    100.0 - 100.0 * t_total_unknowns / (SIEVE_INTERVAL * valid_ms),
                    100.0 * saved_prp / (SIEVE_INTERVAL * valid_ms));

                printf("\t~ 2x %.2f PRP/m\t\t(%ld new composites ~ %4.1f skipped PRP => %.1f PRP/seconds)\n\n",
                    1 / prob_prime_after_sieve,
                    saved_prp,
                    skipped_prp,
                    1.0 * skipped_prp / int_secs);
            }
            setlocale(LC_NUMERIC, "C");

            if (config.verbose + is_last >= 1) {
                s_total_unknowns = t_total_unknowns;
                s_interval_t = s_stop_t;
                s_prp_needed = 1 / prob_prime_after_sieve;

                s_small_prime_factors_interval = 0;
                s_large_prime_factors_interval = 0;
                pi_interval = 0;
            }

            // if is_last would truncate .sieve_range by 1 million
            if (g_control_c && (prime != LAST_PRIME)) {
                // NOTE: the resulting files were sieved by 1 extra prime
                // they will differ from --sieve_range=X in a few entries

                if (prime < 1'000'000) {
                    cout << "Exit(2) from CTRL+C @ prime=" << prime << endl;
                    exit(2);
                }

                cout << "Breaking loop from CTRL+C @ prime=" << prime << endl;
                config.sieve_range = prime - (prime % 1'000'000);

                return false;
            }
        }

        return true;
    });

    // TODO assert s_prime_factors is close to (2 * SL * valid_m) * (log(log(SL)) + M_c)

    if (config.save_unknowns) {
        // ----- Open and Save to Output file
        std::ofstream unknown_file;
        {
            std::string fn = gen_unknown_fn(config, ".txt");
            printf("\nSaving unknowns to '%s'\n", fn.c_str());
            unknown_file.open(fn, std::ios::out);
            assert( unknown_file.is_open() ); // Can't open save_unknowns file
        }

        for (uint64_t mi : valid_mi) {
            assert(gcd(M_start + mi, D) == 1);
            int32_t mii = m_reindex[mi];
            assert( mii >= 0 );

            const auto& comp = composite[mii];

            const size_t size_side = count_coprime_sieve / 2;
            size_t unknown_l = std::count(comp.begin(), comp.begin() + size_side, false);
            size_t unknown_u = std::count(comp.begin() + size_side, comp.end(), false);

            unknown_file << mi << " : -" << unknown_l << " +" << unknown_u << " |";
            for (int d = 0; d <= 1; d++) {
                char prefix = "-+"[d];

                for (size_t i = 1; i < SL; i++) {
                    int a = (SL-1) + (2*d -1) * i;
                    if (!comp[i_reindex[a]]) {
                        unknown_file << " " << prefix << i;
                    }
                }
                if (d == 0) {
                    unknown_file << " |";
                }
            }
            unknown_file << "\n";
        }
    }

    mpz_clear(K);
    mpz_clear(test);
}

