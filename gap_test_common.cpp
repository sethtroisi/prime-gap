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

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "gap_common.h"

using std::cout;
using std::endl;
using std::vector;


void load_and_verify_unknowns(
        const int compression,
        const uint64_t m,
        const int SIEVE_LENGTH,
        std::ifstream &unknown_file,
        vector<int32_t> (&unknowns)[2]) {
    int unknown_l = 0;
    int unknown_u = 0;

    {
        uint32_t m_test = 0;
        unknown_file >> m_test;
        if (m_test != m ) {
            cout << "Mismatched m " << m_test << " vs " << m << endl;
            assert (m_test == m);
        }
        std::string delim;
        unknown_file >> delim;
        assert( delim == ":" );

        unknown_file >> unknown_l;
        unknown_l *= -1;
        unknown_file >> unknown_u;

        unknown_file >> delim;
        assert( delim == "|" );

        if (compression) {
            for (int group = 0; group <= 1; group++) {
                unsigned char c1 = 0, c2 = 0;
                int accum = 0;
                int unknown_count = group == 0 ? unknown_l : unknown_u;
                for (int k = 0; k < unknown_count; k++) {
                    unknown_file >> c1;
                    unknown_file >> c2;
                    int delta = 128 * (c1 - 48) + (c2 - 48);
                    accum += delta;
                    assert( 1 <= accum && accum <= SIEVE_LENGTH );
                    unknowns[group].push_back(accum);
                }
                if (group == 0) {
                    unknown_file >> delim;
                    assert( delim == "|" );
                }
            }
        } else {
            int c = 0;
            for (int k = 0; k < unknown_l; k++) {
                unknown_file >> c;
                assert( 1 <= -c && -c <= SIEVE_LENGTH );
                unknowns[0].push_back(-c);
            }
            unknown_file >> delim;
            assert( delim == "|" );

            for (int k = 0; k < unknown_u; k++) {
                unknown_file >> c;
                assert( 1 <= c && c <= SIEVE_LENGTH );
                unknowns[1].push_back(c);
            }
        }
    }

    assert( unknown_l >= 0 && unknown_u >= 0 );
    assert( (size_t) unknown_l == unknowns[0].size() );
    assert( (size_t) unknown_u == unknowns[1].size() );

    assert( is_sorted(unknowns[0].begin(), unknowns[0].begin()) );
    assert( is_sorted(unknowns[1].begin(), unknowns[1].begin()) );
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

        // Double checks SL which is fine.
        mpz_sub_ui(prime_test, center, SIEVE_LENGTH);
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
