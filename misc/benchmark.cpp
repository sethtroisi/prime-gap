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
#include <cstring>
#include <iostream>
#include <random>
#include <set>

#include <gmp.h>
#include <primesieve.hpp>

#include "../modulo_search.h"

using std::cout;
using std::endl;
using std::set;
using std::vector;
using namespace std::chrono;



#define DELTA_SINCE(start_time) duration<double>( \
    high_resolution_clock::now() - start_time).count();

const uint32_t M_START = 1000;


#include <thread>
#include <x86intrin.h>

// global frequency
double freq;

double approx_Hz(unsigned sleeptime) {
    auto t_start = high_resolution_clock::now();
    uint64_t elapsed_cycles = __rdtsc();
    std::this_thread::sleep_for(std::chrono::milliseconds(sleeptime));
    elapsed_cycles = __rdtsc() - elapsed_cycles;
    auto time = DELTA_SINCE(t_start);
    auto time2 = DELTA_SINCE(t_start);

    return elapsed_cycles / (time - 2 * (time2 - time));
}


/**
 * Generate count primes which are between [2 ** (bits-1), 2 ** bits)
 * Store into save
 */
void generate_primes(int bits, size_t count, vector<uint64_t> &save) {
    // Start the primes at a nice round number.
    uint64_t start = (1UL << (bits-1)) + (1UL << (bits-3));
    uint64_t start_mod = 1000;
    while (start_mod * 100 < start) start_mod *= 10;

    start -= (start % start_mod);
    assert(start >= (1UL << (bits-1)) );

    uint64_t max_prime = (1UL << bits) - 1;

    /*
    mpz_t p;
    mpz_init_set_ui(p, start);
    for (size_t i = 1; i <= count; i++) {
        mpz_nextprime(p, p);
        uint64_t t = mpz_get_ui(p);
        if (t > max_prime) {
            break;
        }
        save.push_back(t);
        if (0) {
            if (i <= 3 || i + 2 > count) {
                printf("\t\t%2d, %7ld  %ld\n", bits, i, t);
            }
        }
    }
    mpz_clear(p);
    */

    primesieve::iterator it(start);
    for(uint64_t prime = it.next_prime(); prime <= max_prime; prime = it.next_prime()) {
        save.push_back(prime);
        if (save.size() == count) break;
    }

    assert( (1UL << (bits-1)) <= save.front() );
    assert( save.back() < (1UL << bits) );
}


void generate_PALR(
        int bits, size_t count, size_t S,
        vector<uint64_t> &primes,
        vector<uint64_t> &A,
        vector<uint64_t> &L,
        vector<uint64_t> &R) {

    generate_primes(bits, count, primes);

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
}


void benchmark_primorial_modulo(
        const char* benchmark_row, const char* filter,
        unsigned int P, int bits, size_t count,
        const vector<uint64_t> &primes) {

    mpz_t K;
    mpz_init(K);
    mpz_primorial_ui(K, P);

    auto log2 = mpz_sizeinbase(K, 2);
    char ref_name[50];
    sprintf(ref_name, "%d# mod <%d bit>p", P, bits);

    if (strstr(ref_name, filter)) {
        auto t_start = high_resolution_clock::now();

        uint64_t z = 0;
        for (auto p : primes) {
            z += mpz_fdiv_ui(K, p);
        }

        auto time = DELTA_SINCE(t_start);

        printf(benchmark_row,
            log2, count, ref_name,
            primes.size(), z % 10000,
            time, time * 1e9 / primes.size(),
            time * freq / primes.size() / ((log2 - 1) / 64 + 1));
    }

    mpz_clear(K);
}


uint32_t modulo_search_single_mod(uint32_t p, uint32_t A, uint32_t L, uint32_t R) {
    return (p % A) + L - R;
}


// Create a type for pointer to modulo_search uint32 signature
typedef uint32_t(*modulo_search_uint32_sig)(uint32_t p, uint32_t a, uint32_t l, uint32_t r);

void benchmark_method_small(
        const char* benchmark_row, const char* ref_name, const char* filter,
        int bits, size_t count,
        const vector<uint64_t> &primes,
        const vector<uint64_t> &A,
        const vector<uint64_t> &L,
        const vector<uint64_t> &R,
        modulo_search_uint32_sig ref_func) {
    if (strstr(ref_name, filter) == NULL) return;

    //assert( count <= primes.size() );
    assert( (primes.size() == A.size()) && (A.size() == L.size()) && (A.size() == R.size()) );

    auto t_start = high_resolution_clock::now();

    size_t found = 0;
    for (size_t i = 0; i < primes.size(); i++) {
        uint32_t m = ref_func(primes[i], A[i], L[i], R[i]);
        uint64_t t = (m * A[i]) % primes[i];
        if ((L[i] <= t) && (t <= R[i])) {
            found++;
        }
        assert( ref_func == modulo_search_single_mod ||
                ((L[i] <= t) && (t <= R[i])) );
    }

    auto time = DELTA_SINCE(t_start);

    printf(benchmark_row,
        bits, count, ref_name,
        found, count, time,
        time * 1e9 / count, time * freq / count);
}

/* Used to double check modulo_search_euclid at one point
uint64_t modulo_search_brute_large(uint64_t p, uint32_t max_m, uint64_t A, uint64_t L, uint64_t R) {

    uint64_t goal = R - L;
    uint64_t temp = p - L;
    for (uint32_t i = 0; ; i++) {
        if (temp <= goal) {
            return i;
        }
        temp += A;
        if (temp >= p) temp -= p;
    }
}
*/

void benchmark_method_large(
        const char* benchmark_row, const char* ref_name, const char* filter,
        int bits, size_t count,
        size_t SL, size_t max_m,
        const vector<uint64_t> &primes,
        const vector<uint64_t> &A,
        const vector<uint64_t> &L,
        const vector<uint64_t> &R,
        int method) {
    if (strstr(ref_name, filter) == NULL) return;

    assert( count <= primes.size() );
    assert( (primes.size() == A.size()) && (A.size() == L.size()) && (A.size() == R.size()) );

    auto t_start = high_resolution_clock::now();

    size_t found = 0, found2 = 0;
    for (size_t i = 0; i < primes.size(); i++) {
        uint64_t p = primes[i];
        uint64_t base_r = A[i];

        uint64_t m;
        if (method == 0) {
            // Can't use L, R as euclid_gcd doesn't accept those.
            uint64_t l = (A[i] + SL) % p;
            if (l <= 2*SL) {
                m = 0;
            } else {
                l = p - l;
                uint64_t r = l + 2*SL;
                m = modulo_search_euclid(p, A[i], l, r);
            }
            // Only verify interesting solutions here
            if (m == 0 || m >= max_m)
                continue;

            uint64_t m1 = modulo_search_euclid_gcd(1, 1, max_m, SL, p, base_r);
            uint64_t m2 = modulo_search_euclid_gcd2(1, 1, max_m, SL, p, base_r);

            if (m != m1 || m != m2) {
                printf("Mismatch | max_m: %ld | m: %ld, m1: %ld, m2: %ld\n",
                    max_m, m, m1, m2);
            }

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
            m = modulo_search_euclid_gcd(M_START, 1, max_m, SL, p, base_r);
            found++;
            if (m == max_m) continue;
            found2++;

            uint64_t t = ((__int128) base_r * (M_START + m)) % p;
            assert( (t <= SL) || (t + SL) >= p );

        } else if (method == 3) {
            // D = 1 mean only zero/one modulo_search.
            m = modulo_search_euclid_gcd2(M_START, 1, max_m, SL, p, base_r);
            found++;
            if (m == max_m) continue;
            found2++;

            uint64_t t = ((__int128) base_r * (M_START + m)) % p;
            assert( (t <= SL) || (t + SL) >= p );

        } else if (method == 4) {
            uint64_t previous = found2;
            modulo_search_euclid_all_small(
                M_START, max_m, SL, p, base_r,
                [&](uint32_t mi, uint64_t first) {
                    found2++;
                    uint64_t t = (((base_r * M_START) % p) + ((base_r * mi) % p) + SL) % p;
                    assert( first == t );
                    assert( (t <= 2 * SL) );
                }
            );
            // Did we find any m for this prime?
            found += (found2 > previous);
            // to account for the final search
            found2++;

        } else {
            uint64_t previous = found2;

            modulo_search_euclid_all_large(
                M_START, max_m, SL, p, base_r,
                [&](uint32_t mi, uint64_t first) {
                    found2++;
                    __int128 t = base_r;
                    t *= (M_START + mi);
                    t += SL;
                    t %= p;
                    assert( first == t );
                    assert( t <= 2 * SL );
                }
            );
            // Did we find any m for this prime?
            found += (found2 > previous);
            // to account for the final search
            found2++;
        }
    }

    auto time = DELTA_SINCE(t_start);

    found2 = found2 == 0 ? found : found2;

    printf(benchmark_row,
        bits, count, ref_name,
        found, found2,
        time, time * 1e9 / found2, time * freq / found2);
}


void benchmark(int bits, size_t count, const char* filter) {
    /**
     * These control how many quickly solutions are found
     * And how likely large primes are to find solutions
     */
    const size_t SL = bits <= 32 ? 100'000 :
        (bits <= 42 ? 1'000'000 : 100'000'000);
    const size_t S = 2 * SL + 1;  // R - L = S

    vector<uint64_t> primes, A, L, R;
    generate_PALR(bits, count, S, primes, A, L, R);
    assert( (primes.size() == count) || (bits <= 25));

    size_t N = primes.size();

    // Header
    cout << endl;
    const auto benchmark_row = "\t| %5d x %8ld | %-30s | %-8ld | %-8ld | %-7.4f | %7.0f | %-8.1f    |\n";
    printf("\t|  bits x count    | method_name%19s | found    | total    | time(s) | ns/iter | cycles/iter |\n", "");
    auto padding = "-------------------------------------------------";
    printf("\t|%.18s|%.32s|%.10s|%.10s|%.9s|%.9s|%.13s|\n",
        padding, padding, padding, padding, padding, padding, padding);

    if (bits < 32) {
        // Per Method benchmark
        benchmark_method_small(
            benchmark_row, "single_mod_op", filter, bits, N,
            primes, A, L, R, modulo_search_single_mod);

        if (bits != 32) {
            benchmark_method_small(
                benchmark_row, "modulo_search_brute", filter, bits, N / 10,
                primes, A, L, R, modulo_search_brute);
        }

        benchmark_method_small(
            benchmark_row, "modulo_search_euclid_small", filter, bits, N,
            primes, A, L, R, modulo_search_euclid_small);
    }

    size_t max_m_overflow = std::numeric_limits<uint64_t>::max() / primes.back() - 10;

    benchmark_method_large(
        benchmark_row, "modulo_search_verify", filter, bits, N,
        SL, max_m_overflow, primes, A, L, R, 0);

    benchmark_method_large(
        benchmark_row, "modulo_search_euclid", filter, bits, N,
        SL, max_m_overflow, primes, A, L, R, 1);

    benchmark_method_large(
        benchmark_row, "modulo_search_euclid_gcd", filter, bits, N,
        SL, max_m_overflow, primes, A, L, R, 2);

    benchmark_method_large(
        benchmark_row, "modulo_search_euclid_gcd2", filter, bits, N,
        SL, max_m_overflow, primes, A, L, R, 3);

    /**
     * Try to set max_m so it finds 1m per p
     * Limit m so doesn't overflow 64 bits
     */
    size_t two_found_per = std::max(2UL, (2 * primes.front()) / S);
    size_t max_m = std::min(max_m_overflow, two_found_per);
    if (max_m == max_m_overflow) {
        cout << "max_m_overflow:" << max_m_overflow << endl;
    }

    benchmark_method_large(
        benchmark_row, "modulo_search_euclid_all_small", filter, bits, N,
        SL, max_m, primes, A, L, R, 4);

    benchmark_method_large(
        benchmark_row, "modulo_search_euclid_all_large", filter, bits, N,
        SL, max_m, primes, A, L, R, 5);

    if (strstr("# mod < bits>p", filter) != NULL) {
        printf("\n");
        printf("\t|  bits x count   | method_name%19s | found    | total    | time(s) | ns/iter | cycles/limb |\n", "");
        printf("\t|%.17s|%.32s|%.10s|%.10s|%.9s|%.9s|%.13s|\n",
            padding, padding, padding, padding, padding, padding, padding);
    }
    for (auto P : {503, 1009, 1999, 5003, 10007, 20011}) {
        benchmark_primorial_modulo(
            benchmark_row, filter,
            P,
            bits, N, primes);
    }
}

int main(int argc, char **argv) {
    set<int> benchmark_sizes = {25, 30, 31, 32, 35, 40, 45, 55, 60};
    //set<int> benchmark_sizes = {25, 30, 31};
    //set<int> benchmark_sizes = {35, 40, 45};
    //set<int> benchmark_sizes = {55};

    // Input validation
    assert(argc == 2 || argc == 3);

    long count_primes = atol(argv[1]);
    //assert( (count_primes > 0) && (count_primes <= 10'000'000) );

    char empty_filter[] = "";
    char *filter = empty_filter;
    if (argc == 3) {
        filter = argv[2];
    }

    // Status output
    cout << "Using " << count_primes << " primes during benchmark" << endl;

    // Benchmarking
    cout << endl;

    freq = approx_Hz(250);
    printf("Starting benchmarking (%1.2g Hz)\n", freq);

    for (int bits : benchmark_sizes)
        benchmark(bits, count_primes, filter);

    return 0;
}
