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
        printf("Compiled with GMP %d.%d.%d\n",
            __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);
    }

    if( !has_prev_prime_gmp() ) {
        cout << "See Notes in README.md for instructions on using dev GMPlib" << endl;
        return 1;
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

    setlocale(LC_NUMERIC, "");
    if (config.verbose >= 0) {
        printf("\n");
        printf("Testing m * %d#/%d, m = %ld + [0, %'ld)\n",
            config.p, config.d, config.mstart, config.minc);
    }

    if (config.verbose >= 2) {
        printf("\n");
        printf("sieve_length: 2x %'d\n", config.sieve_length);
        printf("sieve_range:  %'ld\n", config.sieve_range);
        printf("\n");

        printf("run_prp:  %d\n", config.run_prp);
        printf("\n");
    }
    setlocale(LC_NUMERIC, "C");

    prime_gap_test(config);
}


void load_and_verify_unknowns(
        const uint64_t mi,
        const unsigned int SIEVE_LENGTH,
        std::ifstream &unknown_file,
        vector<int32_t> (&unknowns)[2]) {
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
            assert( 1 <= -c && -c < (int) SIEVE_LENGTH );
            unknowns[0].push_back(-c);
        }
        unknown_file >> delim;
        assert( delim == "|" );

        for (int k = 0; k < unknown_u; k++) {
            unknown_file >> c;
            assert( 1 <= c && c < (int) SIEVE_LENGTH );
            unknowns[1].push_back(c);
        }
    }

    assert( unknown_l >= 0 && unknown_u >= 0 );
    assert( (size_t) unknown_l == unknowns[0].size() );
    assert( (size_t) unknown_u == unknowns[1].size() );

    assert( is_sorted(unknowns[0].begin(), unknowns[0].begin()) );
    assert( is_sorted(unknowns[1].begin(), unknowns[1].begin()) );
}

void test_interval(
        const uint64_t m, const mpz_t &K, const size_t SIEVE_LENGTH,
        size_t &s_total_prp_tests,
        size_t &s_gap_out_of_sieve_prev, size_t &s_gap_out_of_sieve_next,
        vector<int32_t> (&unknowns)[2],
        int &next_p, int &prev_p) {

    mpz_t center, ptest;
    mpz_init(center); mpz_init(ptest);
    mpz_mul_ui(center, K, m);

    for (auto low : unknowns[0]) {
        s_total_prp_tests += 1;

        mpz_sub_ui(ptest, center, low);
        if (mpz_probab_prime_p(ptest, 25)) {
            prev_p = low;
            break;
        }
    }
    for (auto high : unknowns[1]) {
        s_total_prp_tests += 1;

        mpz_add_ui(ptest, center, high);
        if (mpz_probab_prime_p(ptest, 25)) {
            next_p = high;
            break;
        }
    }

    if (prev_p == 0) {
        s_gap_out_of_sieve_prev += 1;

        mpz_sub_ui(ptest, center, SIEVE_LENGTH - 1);
        mpz_prevprime(ptest, ptest);
        mpz_sub(ptest, center, ptest);
        prev_p = mpz_get_ui(ptest);
    }

    if (next_p == 0) {
        s_gap_out_of_sieve_next += 1;

        mpz_add_ui(ptest, center, SIEVE_LENGTH - 1);
        mpz_nextprime(ptest, ptest);
        mpz_sub(ptest, ptest, center);
        next_p = mpz_get_ui(ptest);
    }

    mpz_clear(center); mpz_clear(ptest);
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
    int K_digits;
    double K_log;
    K_stats(config, K, &K_digits, &K_log);
    {
        float m_log = log(M_start);
        if (config.verbose >= 1) {
            printf("Min Gap ~= %d (for merit > %.1f)\n",
                (int) (min_merit * (K_log + m_log)), min_merit);
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

        vector<uint32_t> primes = get_sieve_primes(P);
        for (auto prime : primes) {
            assert( prime <= P );
            if (gcd(D, prime) == 1) {
                prob_prime_coprime *= (1 - 1.0/prime);
            }
        }

        // TODO gcd_ui(K, i)
        size_t count_coprime = SL-1;
        for (size_t i = 1; i < SL; i++) {
            for (auto prime : primes) {
                assert( prime <= P );

                if ((i % prime) == 0 && (D % prime) != 0) {
                    count_coprime -= 1;
                    break;
                }
            }
        }

        double chance_coprime_composite = 1 - prob_prime / prob_prime_coprime;
        double prob_gap_shorter_hypothetical = pow(chance_coprime_composite, count_coprime);

        if (config.verbose >= 2) {
            // count_coprime already includes some parts of unknown_after_sieve
            printf("\n");
            printf("\t%.3f%% of sieve should be unknown (%ldM) ~= %.0f\n",
                100 * unknowns_after_sieve,
                config.sieve_range/1'000'000,
                count_coprime * (unknowns_after_sieve / prob_prime_coprime));
            printf("\t%.3f%% of %d digit numbers are prime\n",
                100 * prob_prime, K_digits);
            printf("\t%.3f%% of tests should be prime (%.1fx speedup)\n",
                100 * prob_prime_after_sieve, 1 / unknowns_after_sieve);
            printf("\t~2x%.1f ~ %.1f PRP tests per m\n",
                1 / prob_prime_after_sieve, 2 / prob_prime_after_sieve);
            printf("\tsieve_length=%d is insufficient ~%.2f%% of time\n",
                SIEVE_LENGTH, 100 * prob_gap_shorter_hypothetical);
        }
    }

    // ----- Open unknown input file
    std::ifstream unknown_file;
    {
        std::string fn = gen_unknown_fn(config, ".txt");
        if (config.verbose >= 1) {
            printf("\nReading unknowns from '%s'\n", fn.c_str());
        }
        unknown_file.open(fn, std::ios::in);
        assert( unknown_file.is_open() ); // Can't open save_unknowns file
        assert( unknown_file.good() );    // Can't open save_unknowns file
    }

    uint64_t valid_ms = 0;
    for (uint64_t mi = 0; mi < M_inc; mi++) {
        if (gcd(M_start + mi, D) == 1) {
            valid_ms++;
        }
    }
    assert(valid_ms > 0);

    uint64_t first_mi = 0;
    for (; first_mi > 0 && gcd(M_start + first_mi, D) > 1; first_mi++);
    assert(first_mi >= 0 && first_mi < M_inc);

    uint64_t last_mi = M_inc - 1;
    for (; last_mi > 0 && gcd(M_start + last_mi, D) > 1; last_mi--);
    assert(last_mi > 0 && last_mi < M_inc);

    // ----- Main sieve loop.
    if (config.verbose >= 1) {
        printf("\nStarting m=%ld\ttesting %ld m (%ld to %ld)\n\n",
            M_start, valid_ms, first_mi, last_mi);
    }

    // Used for various stats
    auto  s_start_t = high_resolution_clock::now();
    uint32_t  s_tests     = 0;
    size_t    s_total_unknown = 0;
    size_t    s_t_unk_low = 0;
    size_t    s_t_unk_hgh = 0;
    size_t    s_total_prp_tests = 0;
    size_t    s_gap_out_of_sieve_prev = 0;
    size_t    s_gap_out_of_sieve_next = 0;
    float     s_best_merit_interval = 0;
    size_t    s_best_merit_interval_m = 0;

    for (uint32_t mi = 0; mi < M_inc; mi++) {
        long m = M_start + mi;
        if (gcd(m, D) > 1) {
            continue;
        }

        vector<int32_t> unknowns[2];

        load_and_verify_unknowns(
            mi, SIEVE_LENGTH, unknown_file, unknowns);

        size_t unknown_l = unknowns[0].size();
        size_t unknown_u = unknowns[1].size();

        s_total_unknown += unknown_l + unknown_u;
        s_t_unk_low += unknown_l;
        s_t_unk_hgh += unknown_u;

        // TODO break out to function, also count tests.
        int prev_p = 0;
        int next_p = 0;
        if (config.run_prp) {
            test_interval(
                m, K, SIEVE_LENGTH,
                s_total_prp_tests,
                s_gap_out_of_sieve_prev, s_gap_out_of_sieve_next,
                unknowns, prev_p, next_p);

            assert( prev_p > 0 && next_p > 0 );

            int gap = next_p + prev_p;
            float merit = gap / (K_log + log(m));

            if (merit > min_merit)  {
                // TODO write to file.
                printf("%-5d %.4f  %ld * %ld#/%ld -%d to +%d\n",
                    gap, merit, m, P, D, prev_p, next_p);
            }
            if (merit > s_best_merit_interval) {
                s_best_merit_interval = merit;
                s_best_merit_interval_m = m;
            }

        }

        s_tests += 1;
        bool is_last = (mi == last_mi);
        if ( is_last || (s_tests == 1 || s_tests == 10 || s_tests == 100 ||
                         s_tests == 500 || s_tests == 1000)
                     || (s_tests % 5000 == 0)) {
            auto s_stop_t = high_resolution_clock::now();
            double   secs = duration<double>(s_stop_t - s_start_t).count();

            if ((config.verbose + is_last) >= 1) {
                printf("\tm=%ld %4ld <- unknowns -> %-4ld\t%4d <- gap -> %-4d\n",
                    m,
                    unknown_l, unknown_u,
                    prev_p, next_p);

                // Stats!
                if (s_tests > secs) {
                    printf("\t    tests     %-10d (%.2f/sec)  %.0f seconds elapsed\n",
                        s_tests, s_tests / secs, secs);
                } else {
                    printf("\t    tests     %-10d (%.2f secs/test)  %.0f seconds elapsed\n",
                        s_tests, secs / s_tests, secs);
                }

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

                    if (config.verbose >= 2) {
                        if (s_gap_out_of_sieve_prev + s_gap_out_of_sieve_next > 0) {
                            printf("\t    fallback prev_gap %ld (%.1f%%), next_gap %ld (%.1f%%)\n",
                                s_gap_out_of_sieve_prev, 100.0 * s_gap_out_of_sieve_prev / s_tests,
                                s_gap_out_of_sieve_next, 100.0 * s_gap_out_of_sieve_next / s_tests);
                        }
                        printf("\t    best merit this interval: %.3f (at m=%ld)\n",
                            s_best_merit_interval, s_best_merit_interval_m);
                    }
                }

                s_best_merit_interval = 0;
                s_best_merit_interval_m = -1;
            }
        }
    }

    // ----- cleanup

    mpz_clear(K);
}

