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

#include <iostream>

#include <cuda.h>
#include <gmp.h>


using std::cout
using std::endl;
using std::vector;


// TODO figure out what to set here
#define GRID_SIZE 16
#define BLOCK_SIZE 4


#define checkCuda(expr) {                          \
  auto status = (expr);                             \
  if (status != cudaSuccess) {                       \
    cerr << "cuda Error on line " << __LINE__ << ": " \
         << cudaGetErrorString(status) << endl;        \
    exit(EXIT_FAILURE);                                 \
  }                                                      \
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

    char *composites,

    uint32_t num_primes,
    uint32_t *primes,
    uint32_t *remainders, // r = K mod p
    int32_t *neg_inv_Ks,     // r^1 mod p
    uint32_t *coprime_X_thread,
) {
    int index = threadIdx.x + blockIdx.x * BLOCK_SIZE;

    // Code assumes K is odd.

    for (uint32_t pi = 0; pi < num_primes; pi++) {
        const uint32_t prime = primes;
        const uint32_t base_r = remainders[pi];
        const int32_t neg_inv_Ks = neg_inv_Ks

        if (1) {
            // as large as prime^2
            uint64_t t = inv_K * base_r
            assert(t % prime == 1);
            assert(base_r < prime);
            assert(0 < neg_inv_Ks);
            assert(neg_inv_Ks < prime);
        }

        // -M_start % p
        int64_t mi_0_shift = prime - (M_start % prime);
        if (m_start_shift == prime) {
            m_start_shift = 0;
        }

        const uint8_t M_parity_check = M_start & 1;
        uint32_t shift = prime << 1;

        // Find m*K = X, X in [L, R]
        // NOTE: X is positive [0, SL]
        for (int64_t X : coprime_X_thread) {
            // Safe from overflow as (SL * prime + prime) < int64
            int64_t mi_0 = (X * neg_inv_K + mi_0_shift) % prime;
            mi_0 += (((X ^ mi_0) & 1) == M_parity_check) ? prime : 0;

            uint64_t mi = mi_0;
            for (; mi < M_inc; mi += shift) {
                uint64_t m = M_start + mi;
                uint32_t m_mod2310 = m % 2310;

                // Filters ~80% or more of m where (m, D) != 1
                if (!is_m_coprime2310[m_mod2310])
                    continue;

                // After initial value this increases by (shift * K_mod2310) % 2310
                uint32_t n_mod2310 = ((K_mod2310 * m_mod2310) + X) % 2310;
                if (!is_coprime2310[n_mod2310])
                    continue;

                // TODO readd is_m_coprime test.

                int32_t mii = m_reindex[mi];
                if (mii == 0)
                    continue;

                uint32_t xii = X_reindex[X];
                composite[mii * caches.composite_line_size + xii] = true;
            }
        }
    }
}

#define CUDA_CHECK(action) cuda_check(action, #action, __FILE__, __LINE__)

class GPUSieve {
    public:
        cudaStream_t runner_stream;

        // Cache stuff
        uint32_t *coprime_X,
        uint16_t *X_reindex;
        char     *is_coprime2310,
        char     *is_m_coprime2310;
        int32_t  *m_reindex,

        // Cached prime stuff
        uint32_t num_primes,
        uint32_t *primes,
        uint32_t *remainders, // r = K mod p
        int32_t *neg_inv_Ks,     // r^1 mod p

        // Maybe later? is_m_coprime
        char *composites

        uint32_t K_mod2310;

        GPUSieve(
            const Config &config,
            const mpz_t &K,
            const Cached &caches,
            const uint64_t prime_start, const uint64_t prime_end,
        );

        ~GPUSieve() {
            CUDA_CHECK(cudaFree(coprime_X));
            CUDA_CHECK(cudaFree(X_reindex));
            CUDA_CHECK(cudaFree(is_coprime2310))
            CUDA_CHECK(cudaFree(is_m_coprime2310))
            CUDA_CHECK(cudaFree(m_reindex))

            CUDA_CHECK(cudaFree(primes))
            CUDA_CHECK(cudaFree(remainders))
            CUDA_CHECK(cudaFree(neg_inv_Ks))

            CUDA_CHECK(cudaStreamDestroy(runner_stream)
        }

        run_sieve(
            uint64_t M_start, uint32_t M_inc,
            const Cached &caches, dynamic_bitset<uint32_t> &composite);
};


/**
 * TODO moved somewhere sharde with sieve_small
 * Return a^-1 mod p
 * a^-1 * a mod p = 1
 * assumes that gcd(a, p) = 1
 */
int32_t invert(int32_t a, int32_t p) {
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


GPUSieve::GPUSieve(
            const Config &config,
            const mpz_t &K,
            const Cached &caches,
            const uint64_t prime_start, const uint64_t prime_end,
) {
    //auto T0 = chrono::high_resolution_clock::now();

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
        for (; prime <= prime_end; prime = iter.next_prime()) {
            const uint64_t base_r = mpz_fdiv_ui(K, prime);
            const int32_t inv_K = invert(base_r, prime);
            assert((inv_K * base_r) % prime == 1);
            const int64_t neg_inv_K = (prime - inv_K) % prime;

            host_primes.push_back(prime);
            host_remainders.push_back(base_r);
            host_neg_inv_Ks.push_back(neg_inv_K);
            num_primes += 1;
        }

        const size_t bytes = sizeof(uint32_t) * num_primes;
        checkCuda(cudaMalloc(&primes, bytes);
        checkCuda(cudaMemcpy(host_primes, primes, bytes, cudaMemcpyHostToDevice));

        checkCuda(cudaMalloc(&remainders, bytes);
        checkCuda(cudaMemcpy(host_remainders, remainders, bytes, cudaMemcpyHostToDevice));

        checkCuda(cudaMalloc(&neg_inv_Ks, bytes);
        checkCuda(cudaMemcpy(host_neg_inv_Ks, neg_inv_Ks, bytes, cudaMemcpyHostToDevice));
    }

        // Cache stuff
        uint32_t *coprime_X,
        uint16_t *X_reindex;
        char     *is_coprime2310,
        char     *is_m_coprime2310;
        int32_t  *m_reindex,

        // Maybe later? is_m_coprime
        char *composites

        uint32_t K_mod2310;



    //auto T1 = chrono::high_resolution_clock::now();
    //auto gpu_setup_ms = chrono::duration_cast<chrono::milliseconds>(T1 - T0).count();
    //cout << "GPU setup " << gpu_setup_ms << " ms" << endl;
}

run_sieve(
    uint64_t M_start, uint32_t M_inc,
    const Cached &caches, dynamic_bitset<uint32_t> &composite
) {

  method2_medium_primes_kernal<<<GRID_SIZE, BLOCK_SIZE>>>(
      d_is_prime, d_div_mods, d_primes, pi_start, pi_end, d_results);

  checkCuda(cudaMemcpy(results, d_results, pi_end, cudaMemcpyDeviceToHost));
  }
