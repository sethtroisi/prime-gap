// Copyright 2025 Seth Troisi
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

#pragma once

#include <cassert>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <vector>

#include <cuda.h>
#include <gmp.h>
#include <primesieve.hpp>
#include <boost/dynamic_bitset.hpp>

#include "gap_common.h"
#include "cached.h"

using std::cout;
using std::endl;
using std::vector;
using boost::dynamic_bitset;
using namespace std::chrono;


// TODO figure out what to set here
#define GRID_SIZE 16
#define BLOCK_SIZE 16


// support routines
void cuda_check(cudaError_t status, const char *action=NULL, const char *file=NULL, int32_t line=0) {
  // check for cuda errors

  if(status!=cudaSuccess) {
    printf("CUDA error occurred: %s\n", cudaGetErrorString(status));
    if(action!=NULL)
      printf("While running %s   (file %s, line %d)\n", action, file, line);
    exit(1);
  }
}
#define CUDA_CHECK(action) cuda_check(action, #action, __FILE__, __LINE__)


__global__ void count_set_bits_in_array(char *array, uint64_t bytes) {
    uint32_t set = 0;
    for (size_t i = 0; i < bytes; i++) {
        set += array[i] > 0;
    }
    printf("%u/%lu bytes are non zero\n", set, bytes);
}

/** Called by host executed on device. */
__global__ void method2_medium_primes_kernal(
    /** Caches section **/
    uint16_t *X_reindex,
    char *is_m_coprime2310,
    char *is_coprime2310,
    // Maybe later? is_m_coprime
    int32_t *m_reindex,
    /** End caches section */

    /** config section **/
    uint64_t M_start,
    uint32_t M_inc,
    uint32_t K_mod2310,

    char *composite,

    uint32_t num_primes,
    uint32_t *primes,
    uint32_t *remainders, // r = K mod p
    int32_t *neg_inv_Ks,  // r^1 mod p

    uint32_t num_coprimes,
    uint32_t *coprime_X_thread
) {
    int index = threadIdx.x + blockIdx.x * BLOCK_SIZE;
    uint32_t total_threads = GRID_SIZE * BLOCK_SIZE;

    uint32_t print_mult = 10000;
    uint32_t next_print = print_mult;
    uint32_t next_mult = 5 * print_mult;

    uint32_t small_factors = 0;
    for (size_t i = 0; i < 2310; i++) {
        small_factors += is_m_coprime2310[i];
        if (index == 0) {
            printf("is_m_coprime2310[%lu] = %u | %u\n", i, (uint32_t) is_m_coprime2310[i], small_factors);
        }
    }

    if (small_factors != 1111)
        return;

    for (uint32_t pi = 0; pi < num_primes; pi++) {
        const uint32_t prime = primes[pi];
        const uint32_t base_r = remainders[pi];
        const int32_t neg_inv_K = neg_inv_Ks[pi];

        if (1) {
            // as large as prime^2
            uint64_t t = ((uint64_t) neg_inv_K) * base_r;
            if (index == 0 && prime >= next_print) {
                printf("prime(%u): %u\n", pi, prime);

                if (next_print == next_mult) {
                    print_mult *= 10;
                    next_print = print_mult;
                    next_mult = 5 * print_mult;
                } else {
                    next_print += print_mult;
                }
            }
            assert(t % prime == (prime-1));
            assert(base_r < prime);
            assert(0 < neg_inv_K);
            assert(neg_inv_K < prime);
        }

        // -M_start % p
        int64_t mi_0_shift = prime - (M_start % prime);
        if (mi_0_shift == prime) {
            mi_0_shift = 0;
        }

        const uint8_t M_parity_check = M_start & 1;
        uint32_t shift = prime << 1;

        // Find m*K = X, X in [L, R]
        // NOTE: X is positive [0, SL]
        for (size_t cxti = index; cxti < num_coprimes; cxti += total_threads) {
            int64_t X = coprime_X_thread[cxti];
            // Safe from overflow as (SL * prime + prime) < int64
            int64_t mi_0 = (X * neg_inv_K + mi_0_shift) % prime;
            mi_0 += (((X ^ mi_0) & 1) == M_parity_check) ? prime : 0;

            uint64_t mi = mi_0;
            for (; mi < M_inc; mi += shift) {
                uint64_t m = M_start + mi;
                uint32_t m_mod2310 = m % 2310;

                if (index == 0 && m_mod2310 == 2303) {
                    printf("%lu -> %u -> %u\n", m, m_mod2310, (uint16_t)is_m_coprime2310[m_mod2310]);
                }
                // Filters ~80% or more of m where (m, D) != 1
                if (!is_m_coprime2310[m_mod2310])
                    continue;

                small_factors += 1;

                // After initial value this increases by (shift * K_mod2310) % 2310
                uint32_t n_mod2310 = ((K_mod2310 * m_mod2310) + X) % 2310;
                if (!is_coprime2310[n_mod2310])
                    continue;

                // TODO re-add is_m_coprime test.

                int32_t mii = m_reindex[mi];
                if (mii == 0)
                    continue;

                uint32_t xii = X_reindex[X];
                // TODO composite_line_size -> num_coprimes?
                composite[mii * num_coprimes + xii] = true;
                //small_factors += 1;
            }
        }
    }

    printf(" GPU thread: %d -> %u factors\n", index, small_factors);
}

/**
 * TODO moved somewhere sharde with sieve_small
 * Return a^-1 mod p
 * a^-1 * a mod p = 1
 * assumes that gcd(a, p) = 1
 */
int32_t invert_tmp(int32_t a, int32_t p) {
    // Use extended euclidean algorithm to get
    // a * x + p * y = 1
    // inverse of a is then x
    int x = 1, y = 0;
    int x1 = 0, y1 = 1, a1 = a, p1 = p;
    while (p1) {
        int32_t q = a1 / p1;

        int32_t t = x1;
        x1 = x - q * x1;
        x = t;

        t = y1;
        y1 = y - q * y1;
        y = t;

        t = p1;
        p1 = a1 - q * p1;
        a1 = t;
    }
    return x < 0 ? (p + x) : x;
}

class GPUSieve {
    public:
        cudaStream_t runner_stream;

        // Cache stuff
        uint32_t num_coprimes;
        uint32_t *coprime_X;

        uint16_t *X_reindex;
        char     *is_coprime2310;
        char     *is_m_coprime2310;
        int32_t  *m_reindex;

        // Cached prime stuff
        uint32_t num_primes;
        uint32_t *primes;
        uint32_t *remainders; // r = K mod p
        int32_t *neg_inv_Ks;  // r^1 mod p

        // Maybe later? is_m_coprime
        char *composite;

        uint64_t M_start;
        uint32_t K_mod2310;

        ~GPUSieve() {
            CUDA_CHECK(cudaFree(coprime_X));
            CUDA_CHECK(cudaFree(X_reindex));
            CUDA_CHECK(cudaFree(is_coprime2310));
            CUDA_CHECK(cudaFree(is_m_coprime2310));
            CUDA_CHECK(cudaFree(m_reindex));

            CUDA_CHECK(cudaFree(primes));
            CUDA_CHECK(cudaFree(remainders));
            CUDA_CHECK(cudaFree(neg_inv_Ks));

            CUDA_CHECK(cudaStreamDestroy(runner_stream));
        }

        GPUSieve(
                    const Config &config,
                    const mpz_t &K,
                    const Cached &caches,
                    const uint64_t prime_start, const uint64_t prime_end
        ) {
            auto T0 = high_resolution_clock::now();

            CUDA_CHECK(cudaSetDevice(0));
            CUDA_CHECK(cudaStreamCreate(&runner_stream));

            // TODO: only use this if batch takes > 100ms
            // Reduces GPU_THREAD cpu from 100% while waiting
            //CUDA_CHECK(cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync));
            CUDA_CHECK(cudaSetDeviceFlags(cudaDeviceScheduleSpin));


            { // Compute prime stuff and copy over
                vector<uint32_t> host_primes;
                vector<uint32_t> host_remainders;
                vector<int32_t> host_neg_inv_Ks;

                primesieve::iterator iter(prime_start);
                uint32_t prime = iter.next_prime();
                num_primes = 0;
                for (; prime <= prime_end; prime = iter.next_prime()) {
                    const uint64_t base_r = mpz_fdiv_ui(K, prime);
                    const int32_t inv_K = invert_tmp(base_r, prime);
                    assert( (inv_K * base_r) % prime == 1 );
                    const int64_t neg_inv_K = (prime - inv_K) % prime;
                    assert( (neg_inv_K * base_r) % prime == (prime-1) );

                    host_primes.push_back(prime);
                    host_remainders.push_back(base_r);
                    host_neg_inv_Ks.push_back(neg_inv_K);
                    num_primes += 1;
                }

                const size_t bytes = sizeof(uint32_t) * num_primes;
                CUDA_CHECK(cudaMalloc(&primes, bytes));
                CUDA_CHECK(cudaMemcpy(primes, host_primes.data(), bytes, cudaMemcpyHostToDevice));

                CUDA_CHECK(cudaMalloc(&remainders, bytes));
                CUDA_CHECK(cudaMemcpy(remainders, host_remainders.data(), bytes, cudaMemcpyHostToDevice));

                CUDA_CHECK(cudaMalloc(&neg_inv_Ks, bytes));
                CUDA_CHECK(cudaMemcpy(neg_inv_Ks, host_neg_inv_Ks.data(), bytes, cudaMemcpyHostToDevice));
                // TODO use stream synchronize
                cudaDeviceSynchronize();
            }

            num_coprimes = caches.coprime_X.size();
            const size_t X_bytes = sizeof(uint32_t) * num_coprimes;
            CUDA_CHECK(cudaMalloc(&coprime_X, X_bytes));
            CUDA_CHECK(cudaMemcpy(coprime_X, (void*)caches.coprime_X.data(), X_bytes, cudaMemcpyHostToDevice));

            {
                // From uint32_t -> uint16_t
                vector<uint16_t> X_reindex_tmp;
                for (uint32_t x_i : caches.x_reindex) X_reindex_tmp.push_back(x_i);
                CUDA_CHECK(cudaMalloc(&X_reindex, X_bytes/2));
                CUDA_CHECK(cudaMemcpy(X_reindex, X_reindex_tmp.data(), X_bytes/2, cudaMemcpyHostToDevice));
                count_set_bits_in_array<<<1, 1>>>((char*) X_reindex, X_bytes/2);
            }

            CUDA_CHECK(cudaMalloc(&is_coprime2310, 2310));
            CUDA_CHECK(cudaMemcpy(is_coprime2310, (void*)caches.is_coprime2310.data(), 2310, cudaMemcpyHostToDevice));

            {
                // From bitset to char
                char is_m_coprime2310_tmp[2310];
                for (size_t i = 0; i < 2310; i++) is_m_coprime2310_tmp[i] = caches.is_m_coprime2310[i];
                CUDA_CHECK(cudaMalloc(&is_m_coprime2310, 2310));
                CUDA_CHECK(cudaMemcpy(is_m_coprime2310, is_m_coprime2310_tmp, 2310, cudaMemcpyHostToDevice));
                // TODO use stream synchronize
                cudaDeviceSynchronize();
            }

            const size_t m_reindex_bytes = sizeof(int32_t) * caches.m_reindex.size();
            CUDA_CHECK(cudaMalloc(&m_reindex, m_reindex_bytes));
            CUDA_CHECK(cudaMemcpy(m_reindex, caches.m_reindex.data(), m_reindex_bytes, cudaMemcpyHostToDevice));

            const size_t composite_bytes = sizeof(char) * num_coprimes * caches.valid_ms;
            CUDA_CHECK(cudaMalloc(&composite, composite_bytes));
            CUDA_CHECK(cudaMemset(composite, 0, composite_bytes));

            printf("GPUSieve(): malloced: 4x %d, %lu, %lu, 2x%d, %lu, %lu\n",
                    4*num_primes, X_bytes, X_bytes/2, 2310, m_reindex_bytes, composite_bytes);

            M_start   = config.mstart;
            K_mod2310 = caches.K_mod2310;

            auto T1 = high_resolution_clock::now();
            auto gpu_setup_ms = duration_cast<milliseconds>(T1 - T0).count();
            cudaDeviceSynchronize();
            cout << "GPU setup: " << gpu_setup_ms << " ms" << endl;
        }

        void run_sieve(
            uint64_t M_start, uint32_t M_inc,
            const Cached &caches, dynamic_bitset<uint32_t> &output_composite
        ) {
            assert( M_start == this->M_start );

            return;

            if (0) { // Run GPU Sieve!
                auto T0 = high_resolution_clock::now();
                method2_medium_primes_kernal<<<GRID_SIZE, BLOCK_SIZE>>>(
                    this->X_reindex,
                    this->is_m_coprime2310,
                    this->is_coprime2310,
                    // Maybe later? is_m_coprime
                    this->m_reindex,

                    M_start,
                    M_inc,
                    this->K_mod2310,

                    this->composite,

                    this->num_primes,
                    this->primes,
                    this->remainders,
                    this->neg_inv_Ks,

                    this->num_coprimes,
                    this->coprime_X // TODO pass a portion of this or something IDK
                );

                // TODO streamsynchronize after async memcopy
                cudaDeviceSynchronize();

                auto T1 = high_resolution_clock::now();
                auto kernel_ms = duration_cast<milliseconds>(T1 - T0).count();
                cout << "GPU sieve: " << kernel_ms << " ms" << endl;
            }

            if (0) { // Parse results back to composite.
                auto T0 = high_resolution_clock::now();
                // Segmented by 1M rows
                const size_t valid_m = caches.valid_ms;
                const size_t segment_size = 1'000'000;
                const size_t segment_bytes = sizeof(char) * segment_size * num_coprimes;
                char* host_composite = (char*) malloc(segment_bytes);
                assert( host_composite != NULL );
                //CUDA_CHECK(cudaMallocHost((void**) &host_composite, segment_bytes));
                printf("\thost buffer: %lu\n", segment_bytes);

                char* composite_start = composite;
                for (size_t start_mii = 0; start_mii < valid_m; start_mii += segment_size) {
                    size_t last_mii = std::min(start_mii + segment_size, valid_m);
                    printf("\tCopying over mi [%lu, %lu)\n", start_mii, last_mii);

                    size_t chunk_bytes = (last_mii - start_mii) * num_coprimes * sizeof(char);
                    printf("Hi %lu -> %p, %p, %lu\n", start_mii, host_composite, composite_start, chunk_bytes); // TODO cleanup
                    assert( 0 < chunk_bytes && chunk_bytes <= segment_bytes );
                    assert( chunk_bytes == segment_bytes || last_mii == valid_m );
                    CUDA_CHECK(cudaMemcpy(host_composite, composite_start, chunk_bytes, cudaMemcpyDeviceToHost));
                    composite_start += chunk_bytes;

                    // TODO could do something smart like build up chunks of dynamic bitset and commit them.
                    size_t had_factor = 0;
                    size_t offset = 0;
                    for(size_t mii = start_mii; mii < last_mii; mii++) {
                        uint64_t m = M_start + caches.valid_mi[mii];
                        const auto &x_reindex_m = caches.x_reindex_wheel[m % caches.x_reindex_wheel_size];
                        size_t m_offset = mii * caches.composite_line_size;

                        for (size_t xi = 0; xi < num_coprimes; xi++) {
                            // if [mii][xi] is composite
                            if (host_composite[offset++]) {
                                had_factor += 1;
                                auto X = caches.coprime_X[xi];
                                auto xii = x_reindex_m[X];
                                output_composite[m_offset + xii] = 1;
                            }
                        }
                    }
                    printf("\tmi [%lu, %lu) -> %lu factors\n", start_mii, last_mii, had_factor);
                }

                // TODO free host_composite
                auto T1 = high_resolution_clock::now();
                auto bitfiddling_ms = duration_cast<milliseconds>(T1 - T0).count();
                cout << "GPU copy-back: " << bitfiddling_ms << " ms" << endl;
            }

        }
};
