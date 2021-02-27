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

#include <cassert>
#include <cmath>
#include <iostream>
#include <random>

#include <primesieve.hpp>
#include <benchmark/benchmark.h>

#include "../modulo_search.h"

using std::vector;


/**
 * Generate count primes which are between [2 ** (bits-1), 2 ** bits)
 * Store into save
 */
vector<uint64_t> generate_primes(int bits, size_t count) {
    // Start the primes at a nice round number.
    vector<uint64_t> primes;
    primes.reserve(count);

    uint64_t start = (1UL << (bits-1)) + (1UL << (bits-3));
    uint64_t start_mod = 1000;
    while (start_mod * 100 < start) start_mod *= 10;

    start -= (start % start_mod);
    assert(start >= (1UL << (bits-1)) );

    uint64_t max_prime = (1UL << bits) - 1;

    primesieve::iterator it(start);
    for(uint64_t prime = it.next_prime(); prime <= max_prime; prime = it.next_prime()) {
        primes.push_back(prime);
        if (primes.size() == count) break;
    }

    assert( (1UL << (bits-1)) <= primes.front() );
    assert( primes.back() < (1UL << bits) );

    return primes;
}


void generate_PALR(
        int bits, size_t count, size_t S,
        vector<uint64_t> &primes,
        vector<uint64_t> &A,
        vector<uint64_t> &L,
        vector<uint64_t> &R) {

    auto prime_tmp = generate_primes(bits, count);
    primes.swap(prime_tmp);

    A.reserve(count);
    L.reserve(count);
    R.reserve(count);

    std::mt19937_64 mt_rand(bits + count);
    for (uint64_t p : primes) {
        // Only need 2 * S, but assert plenty of space
        assert( p > 4*S );

        // A = random number 1 to P
        // L = random number 1 to P - S
        // R = L + S

        uint64_t a = 0;
        while (a == 0) a = mt_rand() % p;

        // TODO: for i <= 10 generate using L = SL, R = -SL
        // and find A[i] such that m=i is a solution.
        // Tricky way is to find modular_inverse (i, p)

        uint64_t l = 0;
        while (l <= S) l = mt_rand() % (p - S);

        assert(a < p);
        // avoid chance of overflow in l + S
        assert(l < (p - S));

        A.push_back(a);
        L.push_back(l);
        R.push_back(l + S);
    }

    if (primes.size() < count) {
        assert(bits <= 25);
        // Ran out of primes for small bit count
        primes.reserve(count);

        for (size_t i = 0; primes.size() < count; i++) {
            primes.push_back(primes[i]);
            A.push_back(A[i]);
            L.push_back(L[i]);
            R.push_back(R[i]);
        }
    }
}


static void BM_module_search_euclid(benchmark::State& state) {
    assert(state.max_iterations <= 50'000'000);
    size_t count = state.max_iterations;
    assert(count > 0);

    // state doesn't change while we iterate
    int bits = state.range(0);
    size_t SL = state.range(1);
    vector<uint64_t> primes, A, L, R;
    generate_PALR(bits, count, 2 * SL + 1, primes, A, L, R);

    assert( (primes.size() == count) );
    assert( (primes.size() == A.size()) &&
            (primes.size() == L.size()) &&
            (primes.size() == R.size()) );

    size_t i = 0;
    for (auto _ : state) {
        uint64_t m = modulo_search_euclid(primes[i], A[i], L[i], R[i]);
        benchmark::DoNotOptimize(m);

        i++;
    }
}


static void BM_module_search_euclid_stack(benchmark::State& state) {
    assert(state.max_iterations <= 50'000'000);
    size_t count = state.max_iterations;
    assert(count > 0);

    int bits = state.range(0);
    size_t SL = state.range(1);
    vector<uint64_t> primes, A, L, R;
    generate_PALR(bits, count, 2 * SL + 1, primes, A, L, R);

    assert( (primes.size() == count) );
    assert( (primes.size() == A.size()) &&
            (primes.size() == L.size()) &&
            (primes.size() == R.size()) );

    size_t i = 0;
    for (auto _ : state) {
        uint64_t m = modulo_search_euclid_stack(primes[i], A[i], L[i], R[i]);
        benchmark::DoNotOptimize(m);

        i++;
    }
}

static void BM_module_search_euclid_verify_both(benchmark::State& state) {
    assert(state.max_iterations <= 50'000'000);
    size_t count = state.max_iterations;
    assert(count > 0);

    int bits = state.range(0);
    size_t SL = state.range(1);
    vector<uint64_t> primes, A, L, R;
    generate_PALR(bits, count, 2 * SL + 1, primes, A, L, R);

    assert( (primes.size() == count) );
    assert( (primes.size() == A.size()) &&
            (primes.size() == L.size()) &&
            (primes.size() == R.size()) );

    size_t i = 0;
    for (auto _ : state) {
        assert(bits == state.range(0));
        assert(SL == (size_t) state.range(1));

        uint64_t p = primes[i];
        uint64_t m = modulo_search_euclid_stack(p, A[i], L[i], R[i]);
        uint64_t m2 = modulo_search_euclid(p, A[i], L[i], R[i]);
        assert(m == m2);

        uint64_t t = ((__int128) m * A[i]) % p;
        assert( (L[i] <= t) && (t <= R[i]) );

        i++;
    }
}


static void LargeBitArguments(benchmark::internal::Benchmark* benchmark) {
    benchmark
        // {Number of bits, SL}
        ->Args({25, 15'000})
        ->Args({25, 100'000})
    //    ->Args({30, 15'000})
    //    ->Args({31, 15'000})
    //    ->Args({32, 15'000})
        ->Args({35, 15'000})
    //    ->Args({40, 15'000})
        ->Args({45, 15'000})
    //    ->Args({55, 15'000})
        ->Args({60, 15'000})
        ->Args({60, 100'000});
}


BENCHMARK(BM_module_search_euclid)->Apply(LargeBitArguments);
BENCHMARK(BM_module_search_euclid_stack)->Apply(LargeBitArguments);
BENCHMARK(BM_module_search_euclid_verify_both)->Apply(LargeBitArguments);


BENCHMARK_MAIN();
