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
#include <cmath>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <gmp.h>

#include "gap_common.h"

using std::cout;
using std::endl;
using std::vector;
using namespace std::chrono;


void prime_gap_test(const struct Config config);


int main(int argc, char* argv[]) {
    Config config = argparse(argc, argv);
    if (config.valid == 0) {
        show_usage(argv[0]);
        return 1;
    }

    if (config.verbose >= 2) {
        printf("\tCompiled with GMP %d.%d.%d\n\n",
            __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);
    }

    if (config.sieve_length == 0) {
        cout << "Must set sieve-length for " << argv[0] << endl;
        show_usage(argv[0]);
        return 1;
    }

    if (config.sieve_range == 0) {
        cout << "Must set sieve-length for " << argv[0] << endl;
        show_usage(argv[0]);
        return 1;
    }

    printf("\n");
    printf("Testing m * %d#/%d, m = %ld + [0, %ld)\n",
        config.p, config.d, config.mstart, config.minc);

    printf("\n");
    printf("sieve_length: 2x%d\n", config.sieve_length);
    printf("sieve_range:  %ld\n", config.sieve_range);
    printf("\n");

    printf("run_prp:  %d\n", config.run_prp);
    printf("\n");

    prime_gap_test(config);
}


void prime_gap_test(const struct Config config) {
    const uint64_t M_start = config.mstart;
    const uint64_t M_inc = config.minc;
    const uint64_t P = config.p;
    const uint64_t D = config.d;
    const float min_merit = config.min_merit;

    const unsigned int SIEVE_LENGTH = config.sieve_length;
    const unsigned int SL = SIEVE_LENGTH;

    // ----- Merit Stuff
    mpz_t K;
    mpz_init(K);
    mpz_primorial_ui(K, P);
    assert( 0 == mpz_tdiv_q_ui(K, K, D) );

    int K_digits = mpz_sizeinbase(K, 10);
    float K_log;
    {
        long exp;
        double mantis = mpz_get_d_2exp(&exp, K);
        K_log = log(mantis) + log(2) * exp;
        float m_log = log(M_start);
        int K_bits   = mpz_sizeinbase(K, 2);

        printf("K = %d bits, %d digits, log(K) = %.2f\n",
            K_bits, K_digits, K_log);
        printf("Min Gap ~= %d (for merit > %.1f)\n\n",
            (int) (min_merit * (K_log + m_log)), min_merit);
    }

    // ----- Open unknown input file
    std::ifstream unknown_file;
    {
        std::string fn = gen_unknown_fn(config, ".txt");
        printf("\tReading unknowns from '%s'\n", fn.c_str());
        unknown_file.open(fn, std::ios::in);
        assert( unknown_file.is_open() ); // Can't open save_unknowns file
        assert( unknown_file.good() );    // Can't open save_unknowns file
    }

    // used in next_prime
    assert( P <= 80000 );
    vector<uint32_t> primes = get_sieve_primes(80000);

    // ----- Allocate memory for a handful of utility functions.

    // Remainders of (p#/d) mod prime
    int *remainder   = (int*) malloc(sizeof(int) * primes.size());
    {
        for (size_t pi = 0; pi < primes.size(); pi++) {
            const long prime = primes[pi];

            // Big improvement over surround_prime is reusing this for each m.
            long mod = mpz_fdiv_ui(K, prime);
            assert( 0 <= mod && mod < prime );
            remainder[pi] = mod;
        }
    }

    // ----- Sieve stats
    {
        assert( config.sieve_range >= 1e6 );
        //From Mertens' 3rd theorem
        double unknowns_after_sieve = 1 / (log(config.sieve_range) * exp(GAMMA));

        double prob_prime = 1 / (K_log + log(M_start));
        double prob_prime_coprime = 1;
        double prob_prime_after_sieve = prob_prime / unknowns_after_sieve;

        for (size_t pi = 0; primes[pi] <= P; pi++) {
            if (D % primes[pi] != 0) {
                prob_prime_coprime *= (1 - 1.0/primes[pi]);
            }
        }

        int count_coprime = SL-1;
        for (size_t i = 1; i < SL; i++) {
            for (int prime : primes) {
                if ((unsigned) prime > P) break;
                if ((i % prime) == 0 && (D % prime) != 0) {
                    count_coprime -= 1;
                    break;
                }
            }
        }
        double chance_coprime_composite = 1 - prob_prime / prob_prime_coprime;
        double prob_gap_shorter_hypothetical = pow(chance_coprime_composite, count_coprime);

        // count_coprime already includes some parts of unknown_after_sieve
        printf("\t%.3f%% of sieve should be unknown (%ldM) ~= %.0f\n",
            100 * unknowns_after_sieve,
            config.sieve_range/1'000'000,
            count_coprime * (unknowns_after_sieve / prob_prime_coprime));
        printf("\t%.3f%% of %d digit numbers are prime\n",
            100 * prob_prime, K_digits);
        printf("\t%.3f%% of tests should be prime (%.1fx speedup)\n",
            100 * prob_prime_after_sieve, 1 / unknowns_after_sieve);
        printf("\t~2x%.1f=%.1f PRP tests per m\n",
            1 / prob_prime_after_sieve, 2 / prob_prime_after_sieve);
        printf("\tsieve_length=%d is insufficient ~%.2f%% of time\n",
            SIEVE_LENGTH, 100 * prob_gap_shorter_hypothetical);
        cout << endl;
    }


    // ----- Main sieve loop.
    cout << "\nStarting m=" << M_start << "\n" << endl;

    // vector<bool> uses bit indexing which is ~5% slower.
    vector<char> composite[2] = {
        vector<char>(SIEVE_LENGTH, 0),
        vector<char>(SIEVE_LENGTH, 0)
    };
    assert( composite[0].size() == SIEVE_LENGTH );
    assert( composite[1].size() == SIEVE_LENGTH );

    // Used for various stats
    auto  s_start_t = high_resolution_clock::now();
    uint32_t  s_tests     = 0;
    uint64_t  s_total_unknown = 0;
    uint64_t  s_t_unk_low = 0;
    uint64_t  s_t_unk_hgh = 0;
    uint64_t  s_total_prp_tests = 0;
    uint64_t  s_gap_out_of_sieve_prev = 0;
    uint64_t  s_gap_out_of_sieve_next = 0;
    float     s_best_merit_interval = 0;
    uint64_t  s_best_merit_interval_m = 0;

    for (uint32_t mi = 0; mi < M_inc; mi++) {
        long m = M_start + mi;
        if (gcd(m, D) > 1) {
            continue;
        }

        // Reset sieve array to all composite.
        std::fill_n(composite[0].begin(), SIEVE_LENGTH, 1);
        std::fill_n(composite[1].begin(), SIEVE_LENGTH, 1);

        int unknown_l;
        int unknown_u;
        // Read a line from the file
        {
            uint32_t mtest;
            unknown_file >> mtest;
            if (mtest != mi ) {
                cout << "Mismatched mi " << mtest << " vs " << mi << endl;
            }
            std::string delim;
            unknown_file >> delim;
            assert( delim == ":" );

            unknown_file >> unknown_l;
            unknown_l *= -1;
            unknown_file >> unknown_u;

            unknown_file >> delim;
            assert( delim == "|" );

            int c;
            for (int k = 0; k < unknown_l; k++) {
                unknown_file >> c;
                composite[0][-c] = 0;
            }
            unknown_file >> delim;
            assert( delim == "|" );

            for (int k = 0; k < unknown_u; k++) {
                unknown_file >> c;
                composite[1][c] = 0;
            }
        }

        int unknown_l_test = std::count(composite[0].begin(), composite[0].end(), false);
        int unknown_u_test = std::count(composite[1].begin(), composite[1].end(), false);
        assert( unknown_l == unknown_l_test );
        assert( unknown_u == unknown_u_test );

        s_total_unknown += unknown_l + unknown_u;
        s_t_unk_low += unknown_l;
        s_t_unk_hgh += unknown_u;

        // TODO break out to function, also count tests.
        int prev_p_i = 0;
        int next_p_i = 0;
        if (config.run_prp) {
            mpz_t center, ptest;
            mpz_init(center); mpz_init(ptest);
            mpz_mul_ui(center, K, m);

            for (unsigned int i = 1; (next_p_i == 0 || prev_p_i == 0) && i < SL; i++) {
                if (prev_p_i == 0 && !composite[0][i]) {
                    s_total_prp_tests += 1;

                    mpz_sub_ui(ptest, center, i);
                    if (mpz_probab_prime_p(ptest, 25)) {
                        prev_p_i = i;
                    }
                }
                if (next_p_i == 0 && !composite[1][i]) {
                    s_total_prp_tests += 1;

                    mpz_add_ui(ptest, center, i);
                    if (mpz_probab_prime_p(ptest, 25)) {
                        next_p_i = i;
                    }
                }
            }

            if (next_p_i == 0) {
                s_gap_out_of_sieve_next += 1;
                // Using fallback to slower gmp routine
                //cout << "\tfalling back to mpz_nextprime" << endl;
                mpz_add_ui(ptest, center, SIEVE_LENGTH - 1);
                mpz_nextprime(ptest, ptest);
                mpz_sub(ptest, ptest, center);
                next_p_i = mpz_get_ui(ptest);
            }

            if (prev_p_i == 0) {
                s_gap_out_of_sieve_prev += 1;
                /*
                // REALLY UGLY FALLBACK
                cout << "\tUGLY prevprime hack" << endl;
                mpz_sub_ui(ptest, center, 2*SIEVE_LENGTH-1);
                mpz_nextprime(ptest, ptest);
                if (mpz_cmp(ptest, center) > 1) {
                    cout << m << "What!" << endl;
                    exit(1);
                }

                while (mpz_cmp(ptest, center) < 0) {
                    // save distance
                    mpz_sub(center, center, ptest);
                    prev_p_i = mpz_get_ui(center);
                    mpz_add(center, center, ptest);
                    mpz_nextprime(ptest, ptest);
                }
                // */

                // /*
                // Medium ugly fallback.
                for (int i = SIEVE_LENGTH; ; i++) {
                    bool composite = false;
                    for (size_t pi = 0; pi < primes.size(); pi++) {
                        const long prime = primes[pi];
                        long modulo = (remainder[pi] * m) % prime;
                        if (i % prime == modulo) {
                            composite = true;
                            break;
                        }
                    }
                    if (!composite) {
                        mpz_sub_ui(ptest, center, i);
                        if (mpz_probab_prime_p(ptest, 25)) {
                            prev_p_i = i;
                            break;
                        }
                    }
                }
                // */
            }

            int gap = next_p_i + prev_p_i;
            float merit = gap / (K_log + log(m));

            if (merit > min_merit)  {
                // TODO write to file.
                printf("%d  %.4f  %ld * %ld#/%ld -%d to +%d\n",
                    gap, merit, m, P, D, prev_p_i, next_p_i);
            }
            if (merit > s_best_merit_interval) {
                s_best_merit_interval = merit;
                s_best_merit_interval_m = m;
            }

            mpz_clear(center); mpz_clear(ptest);
        }

        s_tests += 1;
        if ( false || (s_tests == 1 || s_tests == 10 || s_tests == 100 ||
                               s_tests == 500 || s_tests == 1000) ||
              (s_tests % 5000 == 0) ||
              (s_tests == M_inc) ) {
            auto s_stop_t = high_resolution_clock::now();
            double   secs = duration<double>(s_stop_t - s_start_t).count();

            printf("\t%ld %4d <- unknowns -> %-4d\t%4d <- gap -> %-4d\n",
                m,
                unknown_l, unknown_u,
                prev_p_i, next_p_i);

            // Stats!
            if (config.verbose > 0) {
                printf("\t    tests     %-10d (%.2f/sec)  %.0f seconds elapsed\n",
                    s_tests, s_tests / secs, secs);
                printf("\t    unknowns  %-10ld (avg: %.2f), %.2f%% composite  %.2f%% <- %% -> %.2f%%\n",
                    s_total_unknown, s_total_unknown / ((double) s_tests),
                    100.0 * (1 - s_total_unknown / (2.0 * (SIEVE_LENGTH - 1) * s_tests)),
                    100.0 * s_t_unk_low / s_total_unknown,
                    100.0 * s_t_unk_hgh / s_total_unknown);
                if (config.run_prp) {
                    printf("\t    prp tests %-10ld (avg: %.2f) (%.1f tests/sec)\n",
                        s_total_prp_tests,
                        s_total_prp_tests / (float) s_tests,
                        s_total_prp_tests / secs);
                    printf("\t    fallback prev_gap %ld (%.1f%%), next_gap %ld (%.1f%%)\n",
                        s_gap_out_of_sieve_prev, 100.0 * s_gap_out_of_sieve_prev / s_tests,
                        s_gap_out_of_sieve_next, 100.0 * s_gap_out_of_sieve_next / s_tests);
                    printf("\t    best merit this interval: %.3f (at m=%ld)\n",
                        s_best_merit_interval, s_best_merit_interval_m);
                }

                s_best_merit_interval = 0;
                s_best_merit_interval_m = -1;
            }
        }
    }

    // ----- cleanup

    mpz_clear(K);
}

