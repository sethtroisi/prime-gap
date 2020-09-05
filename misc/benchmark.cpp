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
#include <utility>
#include <vector>

#include <gmp.h>

#include "modulo_search.h"

using std::cout;
using std::endl;
using std::pair;
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


// Method 1 seems to use modulo_search_euclid_gcd (to avoid recalling)
/*
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
            if (mi == M_inc) return;

            assert(mi < M_inc);

            __int128 mult = (__int128) base_r * (M_start + mi) + (SL - 1);
            assert( mult % prime < (2*SL-1) );
*/
// Method 2 seems to use modulo_search_euclid_all (to get all
/*
          modulo_search_euclid_all(M_start, M_inc, SL, prime, base_r, [&](const uint64_t mi) {
                uint64_t first = (base_r * (M_start + mi) + (SL-1)) % prime;
                assert( first < SIEVE_INTERVAL );
                first = SIEVE_INTERVAL - first - 1;
                assert( 0 <= first && first < SIEVE_INTERVAL );
          });
*/




void generate_primes(int bits, uint64_t count, vector<uint64_t> &save) {
    mpz_t p;

    uint64_t start = (1UL << bits) + (1UL << (bits -2));
    uint64_t start_mod = 1000;
    while (start_mod * 100 < start) start_mod *= 10;

    start -= (start % start_mod);
    assert(start >= (1UL << bits) );

    mpz_init_set_ui(p, start);

    for (size_t i = 1; i <= count; i++) {
        mpz_nextprime(p, p);
        save.push_back(mpz_get_ui(p));
        if (i <= 3 || i + 2 > count) {
            printf("\t\t%2d, %7ld  %ld\n", bits, i, mpz_get_ui(p));
        }
    }
    cout << endl;

    assert( (1Ul << bits) <= save.front() );
    assert( save.back() < (2UL << bits) );

    mpz_clear(p);
}


// Create a type for pointer to modulo_search uint32 signature
typedef uint32_t(*modulo_search_uint32_sig)(uint32_t p, uint32_t a, uint32_t l, uint32_t r);

// uint32_t modulo_search_one_op(uint32_t p, uint32_t A, uint32_t L, uint32_t R) {

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

pair<size_t, double>
benchmark_method_small(
        size_t to_test,
        const vector<uint64_t> &primes,
        const vector<uint64_t> &A,
        const vector<uint64_t> &L,
        const vector<uint64_t> &R,
        modulo_search_uint32_sig ref_func) {

    assert( to_test <= primes.size() );
    assert( (primes.size() == A.size()) &&
            (A.size() == L.size()) &&
            (A.size() == R.size()) );

    auto t_start = high_resolution_clock::now();

    size_t verified = 0;
    for (size_t i = 0; i < primes.size(); i++) {
        uint32_t m = ref_func(primes[i], A[i], L[i], R[i]);

        uint64_t t = (m * A[i]) % primes[i];
        if ((L[i] <= t) && (t <= R[i])) {
            verified++;
        }
    }

    auto time = DELTA_SINCE(t_start);

    return {verified, time};
}


void benchmark(size_t S, int bits, const vector<uint64_t> &primes) {
    auto t_setup = high_resolution_clock::now();

    vector<uint64_t> A, L, R;

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

    if (0) {
        double secs = DELTA_SINCE(t_setup);
        fprintf(stderr, "\tBenchmark %ld x %d bits, setup %.3f seconds (%ld to %ld)\n",
            primes.size(), bits,
            secs,
            primes.front(), primes.back());
    }

    size_t count = primes.size();

    size_t brute_count = count;
    if (bits >= 30) {
        // Brute is exponetial, lower the count based on bit size.
        brute_count = std::max(10L, (long) (1e9 * 4 * 10 / pow(2.0, bits)));
    }

    // Header
    cout << endl;
    const auto benchmark_row = "\t| %d x %7ld | modulo_search_%-15s | %-8ld | %-7.4f |\n";
    printf("\t| bits x count | method_name%18s | verified | time(s) | us/iter |\n", "");

    if (bits <= 32) {
        // Per Method benchmark
        auto v_m1 = benchmark_method_small(count, primes, A, L, R, modulo_search_one_op);
        printf(benchmark_row, bits, count, "one_op", v_m1.first, v_m1.second);

        auto v_m2 = benchmark_method_small(count, primes, A, L, R, modulo_search_brute);
        printf(benchmark_row, bits, brute_count, "brute", v_m2.first, v_m2.second);

        auto v_m3 = benchmark_method_small(count, primes, A, L, R, modulo_search_euclid_small);
        printf(benchmark_row, bits, count, "euclid_small", v_m3.first, v_m3.second);
    }

// uint32_t modulo_search_euclid_small(uint32_t p, uint32_t a, uint32_t l, uint32_t r)
}

int main(int argc, char **argv) {
    vector<uint64_t> primes[64]; // indexed by number of bits (in start)
    set<int> benchmark_sizes = {20, 30, 40};

    size_t S = 100000; // R - L

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
    for (int bits : benchmark_sizes)
        cout << "\tBenchmarking modulo_search @ " << bits << " bits with S = " << S << endl;
    cout << endl;

    cout << "Using " << count_primes << " primes during benchmark" << endl;

    // Generating primes
    auto t_find_primes = high_resolution_clock::now();

    for (int bits : benchmark_sizes)
        generate_primes(bits, count_primes, primes[bits]);

    double secs = DELTA_SINCE(t_find_primes);
    printf("Found %ld x %ld Primes in %.1f seconds\n",
        benchmark_sizes.size(), count_primes, secs);

    // Benchmarking
    cout << endl;
    cout << "Starting benchmarking" << endl;
    cout << endl;

    for (int bits : benchmark_sizes)
        benchmark(S, bits, primes[bits]);

    return 0;
}
