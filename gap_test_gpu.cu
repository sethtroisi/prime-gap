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
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <gmp.h>

#include "gap_common.h"
#include "gap_test_common.h"
#include "miller_rabin.h"

using std::cout;
using std::endl;
using std::vector;
using namespace std::chrono;

#ifdef GPU_BITS
#define BITS GPU_BITS
#else
#define BITS 1024
#endif


#define WINDOW_BITS ((BITS <= 1024) ? 5 : 6)

void prime_gap_test(const struct Config config);


int main(int argc, char* argv[]) {
    Config config = Args::argparse(argc, argv);
    if (config.valid == 0) {
        Args::show_usage(argv[0]);
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
        Args::show_usage(argv[0]);
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
        printf("max_prime:       %'ld\n", config.max_prime);
        printf("\n");
    }
    setlocale(LC_NUMERIC, "C");

    prime_gap_test(config);
}


class GPUBatched {
    public:
        // current index;
        int i;

        // If z[i] should be tested
        vector<bool>  active;
        // Result from GPU
        vector<int>  result;

        // number to check if prime
        vector<mpz_t*> z;

        // index into 'open' (BatchedM)
        vector<int> open_i;

        // if this is p_i or n_i
        vector<int> open_pn;
        // offset
        vector<int> open_off_i;
};

class BatchedM {
    public:
        long m;
        mpz_t center;
        vector<int32_t> unknowns[2];

        bool p_found = false, n_found = false;
        int prev_p = 0, next_p = 0;

        // if this entry needs to be handled manually
        bool overflow = false;

        int p_tests = 0;
        int n_tests = 0;

};


void prime_gap_test(const struct Config config) {
    const uint64_t M_start = config.mstart;
    const uint64_t M_inc = config.minc;
    const uint64_t P = config.p;
    const uint64_t D = config.d;
    const float min_merit = config.min_merit;

    const unsigned int SIEVE_LENGTH = config.sieve_length;

    // ----- Merit / Sieve stats
    mpz_t K;
    double K_log = prob_prime_and_stats(config, K);
    {
        float m_log = log(M_start);
        if (config.verbose >= 1) {
            printf("Min Gap ~= %d (for merit > %.1f)\n",
                (int) (min_merit * (K_log + m_log)), min_merit);
        }
    }

    // ----- Open unknown input file
    int compression;
    std::ifstream unknown_file;
    {
        std::string fn = Args::gen_unknown_fn(config, ".txt");
        if (config.verbose >= 1) {
            printf("\nReading unknowns from '%s'\n", fn.c_str());
        }
        unknown_file.open(fn, std::ios::in);
        assert( unknown_file.is_open() ); // Can't open save_unknowns file
        assert( unknown_file.good() );    // Can't open save_unknowns file

        compression = Args::guess_compression(config, unknown_file);
    }

    // TODO this is kinda slow (and blocks program start) for very large numbers
    uint64_t valid_ms = 0;
    for (uint64_t mi = 0; mi < M_inc; mi++) {
        if (gcd(M_start + mi, D) == 1) {
            valid_ms++;
        }
    }
    assert(valid_ms > 0 && valid_ms <= M_inc);

    uint64_t first_mi = 0;
    for (; first_mi > 0 && gcd(M_start + first_mi, D) > 1; first_mi++);
    assert(first_mi < M_inc);

    uint64_t last_mi = M_inc - 1;
    for (; last_mi > 0 && gcd(M_start + last_mi, D) > 1; last_mi--);
    assert(last_mi > 0 && last_mi < M_inc);

    // ----- Main sieve loop.
    if (config.verbose >= 1) {
        printf("\n%ld tests M_start(%ld) + mi(%ld to %ld)\n\n",
            valid_ms, M_start, first_mi, last_mi);
    }

    // Used for various stats
    StatsCounters stats(high_resolution_clock::now());

    // High throughput, low overheard (wasted PRP) method
    // Load unknowns for {2,4,8,16,32,64} m's
    // Keep track of {p_i, n_i} for each m

    const size_t BATCHED_M = 256;
    const size_t BATCH_GPU = 1024;
    assert( BATCH_GPU % BATCHED_M == 0);
    const size_t CANDIDATES_PER_M = BATCH_GPU / BATCHED_M;

    /**
     * Originally 8 which has highest throughput but only if we have LOTS of instances
     * this helps reduce the number of parallel instances needed
     */
    const int THREADS_PER_INSTANCE = 16;
    const int ROUNDS = 1;

    // Setup test runner
    printf("BITS=%d\tWINDOW_BITS=%d\n", BITS, WINDOW_BITS);
    printf("PRP/BATCH=%ld\tM/BATCH=%ld (candidates/M=%ld)\n",
            BATCH_GPU, BATCHED_M, CANDIDATES_PER_M);
    printf("THREADS/PRP=%d\n", THREADS_PER_INSTANCE);

    {
        size_t N_bits = mpz_sizeinbase(K, 2) + log2(M_start + last_mi);
        for (size_t bits : {512, 1024, 1536, 2048, 3036, 4096}) {
            if (N_bits <= bits) {
                if (bits < BITS) {
                    printf("\nFASTER WITH `make gap_test_gpu BITS=%ld` (may require `make clean`)\n\n", bits);
                }
                break;
            }
        }
        assert( N_bits < BITS ); // See last debug line.
        assert( BITS <= (1 << (2 * WINDOW_BITS)) );
    }

    // XXX: params1024, params2048 with *runner1024, *runner2048 and only new one of them.
    typedef mr_params_t<THREADS_PER_INSTANCE, BITS, WINDOW_BITS> params;
    test_runner_t<params> runner(BATCH_GPU, ROUNDS);

    std::unique_ptr<BatchedM> open[BATCHED_M];
    size_t open_count = 0;

    GPUBatched gpu_batch;
    // XXX: This is ugly hack because you can't create mpz_t vector easily
    mpz_t z_array[BATCH_GPU];

    { // Initialize gpu_batch variables
        gpu_batch.active.resize(BATCH_GPU, 0);
        gpu_batch.result.resize(BATCH_GPU, -1);

        for (size_t i = 0; i < BATCH_GPU; i++) {
            mpz_init(z_array[i]);
            gpu_batch.z.push_back(&z_array[i]);
        }

        gpu_batch.open_i.resize(BATCH_GPU, -1);
        gpu_batch.open_pn.resize(BATCH_GPU, -1);
        gpu_batch.open_off_i.resize(BATCH_GPU, -1);
    }

    mpz_t test_z;
    mpz_init(test_z);

    // Main loop
    uint32_t mi = 0;
    while (mi < M_inc || open_count > 0) {
        // Add a new M if space
        while (open_count < BATCHED_M && mi < M_inc) {
            int m = M_start + mi;
            if (gcd(m, D) > 1) {
                mi++;
                continue;
            }
            BatchedM *test = new BatchedM;
            test->m = m;

            load_and_verify_unknowns(
                compression, M_start + mi, SIEVE_LENGTH, unknown_file, test->unknowns);

            mpz_init(test->center);
            mpz_mul_ui(test->center, K, test->m);

            // First first unassigned slot
            for (size_t i = 0; i < BATCHED_M; i++) {
                if (!open[i]) {
                    open[i].reset(test);
                    open_count++;
                    //cout << "open[" << i << "] = " << open[i]->m << endl;
                    break;
                }
            }
            mi++;
        }

        // Grap some entries from each item in M
        {
            gpu_batch.i = 0;
            // Turn off all entries in gpu_batch
            std::fill_n(gpu_batch.active.begin(), BATCH_GPU, false);

            for (size_t i = 0; i < BATCHED_M; i++) {
                if (!open[i]) {
                    assert(mi >= M_inc); // we are out of M to add to open
                    continue;
                }
                BatchedM &test = *open[i];

                for (size_t j = 0; j < CANDIDATES_PER_M; j++) {
                    assert(! (test.p_found && test.n_found) );

                    int gpu_i = gpu_batch.i;
                    gpu_batch.open_i[gpu_i] = i;
                    gpu_batch.result[gpu_i] = -1;

                    if (!test.p_found) {
                        if ((size_t) test.p_tests < test.unknowns[0].size()) {
                            gpu_batch.open_pn[gpu_i] = 0;
                            mpz_sub_ui(*gpu_batch.z[gpu_i], test.center, test.unknowns[0][test.p_tests]);
                            gpu_batch.open_off_i[gpu_i] = test.p_tests++;
                        } else {
                            // Haven't found previous prime, but run out of unknowns to test
                            test.prev_p = -1;
                            test.overflow = 1; // Indicates prev side has overflowed and should be processed out of band.
                            break;
                        }
                    }

                    if (test.p_found) {
                        if ((size_t)test.n_tests < test.unknowns[1].size()) {
                            assert(!test.n_found);
                            gpu_batch.open_pn[gpu_i] = 1;
                            mpz_add_ui(*gpu_batch.z[gpu_i], test.center, test.unknowns[1][test.n_tests]);
                            gpu_batch.open_off_i[gpu_i] = test.n_tests++;
                        } else {
                            // Haven't found next prime, but run out of unknowns to test
                            test.next_p = -1;
                            test.overflow = 1; // Indicates next side has overflowed and should be processed out of band.
                            break;
                        }
                    }

                    //gmp_printf("batch[%d] = %d,%d = %d | %Zd\n", gpu_i, i, j, test.m, *gpu_batch.z[gpu_i]);
                    gpu_batch.active[gpu_i] = true;
                    gpu_batch.i++;
                }
            }
        }

        // Run this large batch
        {
            runner.run_test(gpu_batch.z, gpu_batch.result);
        }

        // Decode and possible finalize M
        {
            for (size_t i = 0; i < BATCH_GPU; i++) {
                if (!gpu_batch.active[i]) {
                    continue;
                }
                assert (gpu_batch.result[i] == 0 || gpu_batch.result[i] == 1);

                size_t open_i = gpu_batch.open_i[i];
                /*
                if (open[open_i] && open[open_i]->m == 2141) {
                    const BatchedM &test = *open[open_i];
                    bool n_side = gpu_batch.open_pn[i];
                    int offset = test.unknowns[n_side][gpu_batch.open_off_i[i]];
                    char sgn = "+-"[n_side];
                    cout << "m=" << test.m << " " << sgn << offset
                         << " | " << test.p_found << " " << test.n_found
                         << " | " << test.prev_p << " " << test.next_p
                         << " | " << gpu_batch.result[i] << endl;
                }  // */

                if (gpu_batch.result[i]) {
                    BatchedM *test = open[open_i].get();
                    if (test == nullptr) {
                        //cout << "Found two primes for open[" << open_i << "] from batch of " << CANDIDATES_PER_M << endl;
                        continue;
                    }


                    int offset_i = gpu_batch.open_off_i[i];
                    if (gpu_batch.open_pn[i] == 0) {
                        if (test->p_found) {
                            /*
                            cout << "Found two previous primes for m=" << test->m << endl;
                            cout << "\t" << test->prev_p << " vs "
                                 << test->unknowns[0][offset_i] << "(" << offset_i << ")" << endl;
                            */
                            continue;
                        }

                        // prev_prime found
                        assert(test->p_tests > 0 );
                        test->p_found = true;
                        test->prev_p = test->unknowns[0][offset_i];
                    } else {
                        if (test->n_found) {
                            /*
                             *
                            cout << "Found two next primes for m=" << test->m << endl;
                            cout << "\t" << test->next_p << " vs "
                                 << test->unknowns[1][offset_i] << "(" << offset_i << ")" << endl;
                            */
                            continue;
                        }

                        // next_prime found (and done)
                        assert(test->p_found );
                        assert(test->n_tests > 0 );
                        test->n_found = true;
                        //cout << "Setting " << test->m << " as " << test->unknowns[1][gpu_batch.open_off_i[i]] << endl;;
                        test->next_p = test->unknowns[1][offset_i];
                    }
                }
            }
        }

        // Finalize any finished (or overflowed open)
        {
            for (size_t i = 0; i < BATCHED_M; i++) {
                if (!open[i])
                    continue;

                BatchedM &test = *open[i];

                int prev_p = test.prev_p;
                int next_p = test.next_p;

                if (test.overflow) {
                    if (prev_p == -1) {
                        assert(test.p_tests > 0);
                        //cout << "gap_out_of_sieve_prev m=" << test.m << endl;
                        mpz_sub_ui(test_z, test.center, SIEVE_LENGTH);
                        mpz_prevprime(test_z, test_z);
                        mpz_sub(test_z, test.center, test_z);
                        prev_p = test.prev_p = mpz_get_ui(test_z);
                        test.p_found = true;
                        test.overflow = 0;
                        continue;
                    }
                    if (next_p == -1) {
                        assert(test.p_tests > 0);
                        //cout << "gap_out_of_sieve_next m=" << test.m << endl;
                        mpz_add_ui(test_z, test.center, SIEVE_LENGTH);
                        mpz_nextprime(test_z, test_z);
                        mpz_sub(test_z, test_z, test.center);
                        next_p = test.next_p = mpz_get_ui(test_z);
                        test.n_found = true;
                        test.overflow = 0;
                    }
                }

                if (!test.p_found || !test.n_found) {
                    continue;
                }

                assert( prev_p > 0 && next_p > 0 );

                float merit = (next_p + prev_p) / (K_log + log(test.m));
                if (merit > min_merit)  {
                    // TODO: Record finished mi in log file / db.
                    printf("%-5d %.4f  %ld * %ld#/%ld -%d to +%d\n",
                        (next_p + prev_p), merit, test.m, P, D, prev_p, next_p);
                }

                bool is_last = (mi >= M_inc) && open_count == 0;
                stats.process_results(config, test.m, is_last,
                    test.unknowns[0].size(), test.unknowns[1].size(),
                    prev_p, next_p,
                    test.p_tests, test.n_tests, merit);

                open[i].reset();  // Delete open[i]
                open_count--;
            }
        }
    }

    // ----- cleanup
    {
        mpz_clear(K);
        mpz_clear(test_z);
        for (size_t i = 0; i < BATCH_GPU; i++)
            mpz_clear(*gpu_batch.z[i]);
    }
}

