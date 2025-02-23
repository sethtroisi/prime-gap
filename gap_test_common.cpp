// Copyright 2021 Seth Troisi
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

#include "gap_test_common.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <sstream>

#include "gap_common.h"

using std::cout;
using std::endl;
using std::vector;


void StatsCounters::process_results(
        const Config &config,
        long m, bool is_last,
        size_t unknown_l, size_t unknown_u,
        int prev_p, int next_p,
        int p_tests, int n_tests,
        float merit) {

    s_tests += 1;

    s_total_unknown += unknown_l + unknown_u;
    s_t_unk_prev += unknown_l;
    s_t_unk_next += unknown_u;

    // TODO break out s_p_tests, s_n_tests;
    s_total_prp_tests += p_tests + n_tests;

    // TODO s_side_skips

    s_total_merit += merit;
    s_total_prev_p += std::max(0, prev_p);
    s_total_next_p += std::max(0, next_p);

    if (merit > s_best_merit_interval) {
        s_best_merit_interval = merit;
        s_best_merit_interval_m = m;
    }

    if (possible_print_stats(config, m, is_last, unknown_l, unknown_u, prev_p, next_p)) {
        s_best_merit_interval = 0;
        s_best_merit_interval_m = 0;
    }
}

bool StatsCounters::possible_print_stats(
        const Config &config,
        long m, bool is_last,
        size_t unknown_l, size_t unknown_u,
        int prev_p, int next_p) const {

    // truncate to a nearby multiple of 10000 (avoid making zero)
    size_t print_interval = 1800 * s_tests_per_second;
    if (print_interval > 10000)
        print_interval -= (print_interval % 10000);

    // if s_tests = {1,3,5} * 10 ^ x
    bool is_power_print = (s_tests == 1);
    for (uint64_t p = 100; p <= s_tests; p *= 10) {
        is_power_print |= (s_tests == p) || (s_tests == 3*p) || (s_tests == 5*p);
    }

    if ( is_last || is_power_print || (print_interval > 0 && s_tests % print_interval == 0) ) {
        auto s_stop_t = std::chrono::high_resolution_clock::now();
        double   secs = std::chrono::duration<double>(s_stop_t - s_start_t).count();
        s_tests_per_second = s_tests / secs;

        if ((config.verbose + is_last) >= 1) {
            printf("\tm=%ld %4ld <- unknowns -> %-4ld\t%4d <- gap -> %-4d\n",
                m,
                unknown_l, unknown_u,
                prev_p, next_p);

            // Stats!
            if (s_tests > secs) {
                printf("\t    tests     %-10lu (%.2f/sec)  %.0f seconds elapsed\n",
                    s_tests, s_tests / secs, secs);
            } else {
                printf("\t    tests     %-10lu (%.2f secs/test)  %.0f seconds elapsed\n",
                    s_tests, secs / s_tests, secs);
            }

            printf("\t    unknowns  %-10ld (avg: %.2f), %.2f%% composite  %.2f%% <- %% -> %.2f%%\n",
                s_total_unknown, s_total_unknown / ((double) s_tests),
                100.0 * (1 - s_t_unk_next / ((double) config.sieve_length * s_tests)),
                100.0 * s_t_unk_prev / s_total_unknown,
                100.0 * s_t_unk_next / s_total_unknown);
            printf("\t    prp tests %-10ld (avg: %.2f) (%.1f tests/sec)\n",
                s_total_prp_tests,
                s_total_prp_tests / (float) s_tests,
                s_total_prp_tests / secs);

            if (config.verbose >= 2) {
                // Suppress for now, this is almost exactly 100 - extra sieves prev_gap %
                /*
                if (s_skips_after_one_side) {
                    printf("\t    only next_prime %ld (%.2f%%)\n",
                        s_skips_after_one_side, 100.0 * s_skips_after_one_side / s_tests);
                }
                */
                if (s_gap_out_of_sieve_prev + s_gap_out_of_sieve_next > 0) {
                    printf("\t    extra sieves prev_gap %ld (%.2f%%), next_gap %ld (%.2f%%)\n",
                        s_gap_out_of_sieve_prev, 100.0 * s_gap_out_of_sieve_prev / s_tests,
                        s_gap_out_of_sieve_next, 100.0 * s_gap_out_of_sieve_next / s_tests);
                }

                printf("\t    best merit this interval: %.3f (at m=%ld)\n",
                    s_best_merit_interval, s_best_merit_interval_m);

                /*
                printf("\t    avg merit: %.3f, avg gap: %.1f (%lu + %lu)\n",
                    s_total_merit / s_tests,
                    ((float) s_total_prev_p + s_total_next_p) / s_tests,
                    s_total_prev_p / s_tests,
                    s_total_next_p / s_tests
                );
                */
            }
            return true;
        }
    }
    return false;
}



void test_interval_cpu(
        const uint64_t m, const mpz_t &K, const size_t SIEVE_LENGTH,
        size_t &s_total_prp_tests,
        size_t &s_gap_out_of_sieve_prev, size_t &s_gap_out_of_sieve_next,
        vector<int32_t> (&unknowns)[2],
        int &prev_p, int &next_p) {

    mpz_t center, prime_test;
    mpz_init(center); mpz_init(prime_test);
    mpz_mul_ui(center, K, m);

    for (auto low : unknowns[0]) {
        s_total_prp_tests += 1;
        assert(low > 0);

        mpz_sub_ui(prime_test, center, low);
        if (mpz_probab_prime_p(prime_test, 25)) {
            prev_p = low;
            break;
        }
    }
    for (auto high : unknowns[1]) {
        s_total_prp_tests += 1;
        assert(high > 0);

        mpz_add_ui(prime_test, center, high);
        if (mpz_probab_prime_p(prime_test, 25)) {
            next_p = high;
            break;
        }
    }

    if (prev_p == 0) {
        s_gap_out_of_sieve_prev += 1;

        // If only one-sided then do the thing
        mpz_sub_ui(prime_test, center, unknowns[0].size() > 2 ? SIEVE_LENGTH : 1);
        mpz_prevprime(prime_test, prime_test);
        mpz_sub(prime_test, center, prime_test);
        prev_p = mpz_get_ui(prime_test);
    }

    if (next_p == 0) {
        s_gap_out_of_sieve_next += 1;

        // Double checks SL which is fine.
        mpz_add_ui(prime_test, center, SIEVE_LENGTH);
        mpz_nextprime(prime_test, prime_test);
        mpz_sub(prime_test, prime_test, center);
        next_p = mpz_get_ui(prime_test);
    }

    mpz_clear(center); mpz_clear(prime_test);
}

void sieve_interval_cpu(
        const uint64_t m, const mpz_t &K,
        const vector<std::pair<uint32_t, uint32_t>> p_and_r,
        const int32_t sieve_start, const int32_t sieve_length,
        vector<int32_t> &unknowns) {

    // To make math easy we always start at the bottom of the sieve and head positive
    // This is natural for next_prime and returns the numbers backwards for a prev_prime

    char composite[sieve_length];
    std::fill_n(composite, sieve_length, 0);

    // Handle 2 because it's weird
    {
        bool start_mod_2 = mpz_odd_p(K) ^ (m & 1) ^ (sieve_start & 1);
        for (int32_t i = 1 - start_mod_2; i < sieve_length; i += 2) {
            composite[i] = 1;
        }
    }

    int64_t negative_m = -m;

    // TODO add an assert that m*p doesn't overflow
    for (const auto& [p, r] : p_and_r) {
        // (-(m * K + sieve_start) % r)
        int32_t center_mod = ((int64_t) r * negative_m - sieve_start) % p;
        if (center_mod < 0) {
            center_mod += p;
        }
        if (center_mod < sieve_length) {
            if (r >= (uint32_t) sieve_length) {
                composite[center_mod] = 1;
            } else {
                // Could do by 2*p but not worth the extra code IMO
                for (int32_t i = center_mod; i < sieve_length; i += p) {
                    composite[i] = 1;
                }
            }
        }
    }

    for (size_t i = 0; i < (uint32_t) sieve_length; i++) {
        if (composite[i] == 0)
            unknowns.push_back(i);
    }
}

