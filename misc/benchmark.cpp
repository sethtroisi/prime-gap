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

#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>
#include <set>
#include <vector>

#include <gmp.h>

#include "modulo_search.h"

using std::cout;
using std::endl;
using std::set;
using std::vector;
using namespace std::chrono;


#define DELTA_SINCE(start_time) duration<double>( \
    high_resolution_clock::now() - start_time).count();

uint32_t modulo_search_one_op(uint32_t p, uint32_t A, uint32_t L, uint32_t R) {
    return (p % A) + L - R;
}

// uint32_t modulo_search_brute(uint32_t p, uint32_t A, uint32_t L, uint32_t R)

// uint32_t modulo_search_euclid_small(uint32_t p, uint32_t a, uint32_t l, uint32_t r)

// uint64_t modulo_search_euclid(uint64_t p, uint64_t a, uint64_t l, uint64_t r)

// uint64_t modulo_search_euclid_gcd(
//        uint64_t M, uint64_t D, uint64_t max_m, uint64_t SL,
//        uint64_t prime, uint64_t base_r)

// void modulo_search_euclid_all(
//        uint64_t M, uint64_t max_m, uint64_t SL,
//        uint64_t prime, uint64_t base_r,
//        std::function<void (uint64_t)> lambda)


// Method 1 use modulo_search_euclid_gcd (_gcd avoids returns where ((M + mi), D) which are handled by a different D)
/*
            uint64_t mi = modulo_search_euclid_gcd(
                    M_start, D, M_inc, SL, prime, base_r);
*/

// Method 2 uses modulo_search_euclid_all (using lambda return for valid mi with calling many times)
/*
          modulo_search_euclid_all(M_start, M_inc, SL, prime, base_r, [&](const uint64_t mi) {
                uint64_t first = (base_r * (M_start + mi) + (SL-1)) % prime;
                assert( first < SIEVE_INTERVAL );
                first = SIEVE_INTERVAL - first - 1;
                assert( 0 <= first && first < SIEVE_INTERVAL );
          });
*/



/**
 * Generate count primes which are between [2 ** (bits-1), 2 ** bits)
 * Store into save
 */
void generate_primes(int bits, size_t count, vector<uint64_t> &save) {
    mpz_t p;

    uint64_t start = (1UL << (bits-1)) + (1UL << (bits-3));
    uint64_t start_mod = 1000;
    while (start_mod * 100 < start) start_mod *= 10;

    start -= (start % start_mod);
    assert(start >= (1UL << (bits-1)) );

    mpz_init_set_ui(p, start);

    for (size_t i = 1; i <= count; i++) {
        mpz_nextprime(p, p);
        save.push_back(mpz_get_ui(p));
        if (0) {
            if (i <= 3 || i + 2 > count) {
                printf("\t\t%2d, %7ld  %ld\n", bits, i, mpz_get_ui(p));
            }
        }
    }
    assert( (1Ul << (bits-1)) <= save.front() );
    assert( save.back() < (1UL << bits) );

    mpz_clear(p);
}

void generate_PALR(
        int bits, size_t count, size_t S,
        vector<uint64_t> &primes,
        vector<uint64_t> &A,
        vector<uint64_t> &L,
        vector<uint64_t> &R) {

    generate_primes(bits, count, primes);

    std::mt19937 mt_rand(S);
    for (uint64_t p : primes) {
        assert( p > S );

        // A = random number 1 to P
        // L = random number 1 to P - S
        // R = L + S

        uint64_t a = 0;
        while (a == 0) a = mt_rand() % p;

        uint64_t l = 0;
        while (l == 0) l = mt_rand() % (p - S);

        A.push_back(a);
        L.push_back(l);
        R.push_back(l + S);
    }
}

// Create a type for pointer to modulo_search uint32 signature
typedef uint32_t(*modulo_search_uint32_sig)(uint32_t p, uint32_t a, uint32_t l, uint32_t r);


void benchmark_method_small(
        const char* benchmark_row, const char* ref_name,
        int bits, size_t count,
        const vector<uint64_t> &primes,
        const vector<uint64_t> &A,
        const vector<uint64_t> &L,
        const vector<uint64_t> &R,
        modulo_search_uint32_sig ref_func) {

    assert( count <= primes.size() );
    assert( (primes.size() == A.size()) && (A.size() == L.size()) && (A.size() == R.size()) );

    auto t_start = high_resolution_clock::now();

    size_t found = 0;
    for (size_t i = 0; i < primes.size(); i++) {
        uint32_t m = ref_func(primes[i], A[i], L[i], R[i]);

        uint64_t t = (m * A[i]) % primes[i];
        if ((L[i] <= t) && (t <= R[i])) {
            found++;
        }
        assert( ref_func == modulo_search_one_op ||
                ((L[i] <= t) && (t <= R[i])) );
    }

    auto time = DELTA_SINCE(t_start);

    printf(benchmark_row,
        bits, count, ref_name,
        found, found, time, time * 1e9 / count);
}

void benchmark_method_large(
        const char* benchmark_row, const char* ref_name,
        int bits, size_t count,
        size_t SL, size_t max_m,
        const vector<uint64_t> &primes,
        const vector<uint64_t> &A,
        const vector<uint64_t> &L,
        const vector<uint64_t> &R,
        int method) {

    assert( count <= primes.size() );
    assert( (primes.size() == A.size()) && (A.size() == L.size()) && (A.size() == R.size()) );

    auto t_start = high_resolution_clock::now();

    size_t found = 0, found2 = 0;
    for (size_t i = 0; i < primes.size(); i++) {
        const uint64_t M = 123;
        uint64_t p = primes[i];
        uint64_t base_r = A[i];

        uint64_t m;
        if (method == 0) {
            uint64_t l = (A[i] + SL - 1) % p;
            if (l <= 2*SL-2) {
                m = 0;
            } else {
                l = p - l;
                uint64_t r = l + 2*SL-2;
                m = modulo_search_euclid(p, A[i], l, r);
            }

            uint64_t m1 = modulo_search_euclid_gcd(1, 1, m+1, SL, p, base_r);
            uint64_t m2 = modulo_search_euclid_gcd2(1, 1, m+1, SL, p, base_r);

            assert( m == m1 );
            assert( m == m2 );
            found++;

        } else if (method == 1) {
            m = modulo_search_euclid(p, A[i], L[i], R[i]);
            found++;
            uint64_t t = ((__int128) m * A[i]) % p;
            assert( (L[i] <= t) && (t <= R[i]) );

        } else if (method == 2) {
            // M = 1, D = 1, max_m = 1e9
            // D = 1 mean only zero/one modulo_search.
            m = modulo_search_euclid_gcd(M, 1, max_m, SL, p, base_r);
            found++;
            if (m == max_m) continue;

            uint64_t t = ((__int128) base_r * (M + m)) % p;
            assert( (t < SL) || (t + SL) > p );

        } else if (method == 3) {
            // M = 1, D = 1, max_m = 1e9
            // D = 1 mean only zero/one modulo_search.
            m = modulo_search_euclid_gcd2(M, 1, max_m, SL, p, base_r);
            found++;
            if (m == max_m) continue;

            uint64_t t = ((__int128) base_r * (M + m)) % p;
            assert( (t < SL) || (t + SL) > p );

        } else {
            uint64_t previous = found2;

            modulo_search_euclid_all(
                M, max_m, SL, p, base_r,
                [&](const uint64_t mi) {
                    found2++;
                    uint64_t t = ((__int128) base_r * (M + mi)) % p;
                    assert( (t < SL) || (t + SL) > p );
                }
            );
            // Did we find any m for this prime?
            found += (found2 > previous);

            // one extra
            found2++;

        }
    }

    auto time = DELTA_SINCE(t_start);

    found2 = found2 == 0 ? found : found2;

    printf(benchmark_row,
        bits, count, ref_name,
        found, found2, time, time * 1e9 / found2);
}


void benchmark(int bits, size_t count) {
    auto t_setup = high_resolution_clock::now();

    // TODO: describe this somewhere
    size_t SL = 100000;
    size_t S = 2 * SL - 1;  // R - L = S

    vector<uint64_t> primes, A, L, R;
    generate_PALR(bits, count, S, primes, A, L, R);
    assert( primes.size() == count );

    if (0) {
        double secs = DELTA_SINCE(t_setup);
        printf("\tBenchmark %ld x %d bits, setup %.3f seconds (%ld to %ld)\n",
            primes.size(), bits,
            secs,
            primes.front(), primes.back());
        printf("\t\t%ld to %ld\n", primes.front(), primes.back());
    }

    // Header
    cout << endl;
    const auto benchmark_row = "\t| %d x %7ld | modulo_search_%-15s | %-8ld | %-8ld | %-7.4f | %7.0f |\n";
    printf("\t| bits x count | method_name%18s | found    | total    | time(s) | ns/iter |\n", "");

    if (bits < 32) {
        // Per Method benchmark
        benchmark_method_small(
            benchmark_row, "one_op", bits, count,
            primes, A, L, R, modulo_search_one_op);

        if (bits != 32) {
            benchmark_method_small(
                benchmark_row, "brute", bits, count,
                primes, A, L, R, modulo_search_brute);
        }

        benchmark_method_small(
            benchmark_row, "euclid_small", bits, count,
            primes, A, L, R, modulo_search_euclid_small);
    }

//    for (size_t max_m : {1'000, 100'000}) {
    for (size_t max_m : {1'000}) {
//        cout << "max_m: " << max_m << endl;

        benchmark_method_large(
            benchmark_row, "verify", bits, count,
            SL, max_m, primes, A, L, R, 0);

        benchmark_method_large(
            benchmark_row, "euclid", bits, count,
            SL, max_m, primes, A, L, R, 1);

        benchmark_method_large(
            benchmark_row, "euclid_gcd", bits, count,
            SL, max_m, primes, A, L, R, 2);

        benchmark_method_large(
            benchmark_row, "euclid_gcd2", bits, count,
            SL, max_m, primes, A, L, R, 3);

        // Set max_m so ~1 + 1 m per p
        size_t all_max_m = std::min(max_m, std::max((size_t) 10, 2 * primes.front() / S));
        benchmark_method_large(
            benchmark_row, "euclid_all", bits, count,
            SL, all_max_m, primes, A, L, R, 4);
    }
}

int main(int argc, char **argv) {
    set<int> benchmark_sizes = {25, 30, 31, 32, 35, 40, 45};
    //set<int> benchmark_sizes = {25, 30, 31};
    //set<int> benchmark_sizes = {35, 40, 45};
    //set<int> benchmark_sizes = {55};

    // Input validation
    assert(argc >= 2);

    uint64_t count_primes = atol(argv[1]);
    assert( (count_primes > 0) && (count_primes < 1000001) );
    for (int i = 2; i < argc; i++) {
        int t = atoi(argv[i]);
        assert((10 < t) && (t <= 62));
        benchmark_sizes.insert(t);
    }

    // Status output
    cout << "Using " << count_primes << " primes during benchmark" << endl;

    // Benchmarking
    cout << endl;
    cout << "Starting benchmarking" << endl;

    for (int bits : benchmark_sizes)
        benchmark(bits, count_primes);

    return 0;
}
