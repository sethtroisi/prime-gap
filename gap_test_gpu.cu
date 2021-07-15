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

// TODO accept from makefile
#define BITS 2048
#define WINDOW_BITS 6

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

    if (BITS > (1 << (2 * WINDOW_BITS))) {
        printf("Not enought WINDOW_BITS(%d) for BITS(%d)\n", WINDOW_BITS, BITS);
        return 1;
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

        int p_i = 0, n_i = 0;
        bool p_found = false, n_found = false;

        int p_tests = 0; // Should equal p_i + 1 with zero overheard
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

    uint64_t valid_ms = 0;
    for (uint64_t mi = 0; mi < M_inc; mi++) {
        if (gcd(M_start + mi, D) == 1) {
            valid_ms++;
        }
    }
    assert(valid_ms > 0);

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

    const size_t BATCHED_M = 32;
    const size_t BATCH_GPU = 128;

    /**
     * Originally 8 which has highest throughput but only if we have LOTS of instances
     * this helps reduce the number of parallel instances needed
     */
    const int THREADS_PER_INSTANCE = 16;
    const int ROUNDS = 1;

    assert( BATCH_GPU % BATCHED_M == 0);
    const size_t INSTANCES_PER_M = BATCH_GPU / BATCHED_M;

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
            open_count++;
            for (size_t i = 0; i < BATCHED_M; i++) {
                if (!open[i]) {
                    open[i].reset(test);
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


                for (size_t j = 0; j < INSTANCES_PER_M; j++) {
                    assert(! (test.p_found && test.n_found) );

                    int gpu_i = gpu_batch.i++; // Increment after we get value (now points at next value)
                    gpu_batch.open_i[gpu_i] = i;
                    gpu_batch.active[gpu_i] = true;
                    gpu_batch.result[gpu_i] = -1;

                    if (!test.p_found) {
                        gpu_batch.open_pn[gpu_i] = 0;
                        mpz_sub_ui(*gpu_batch.z[gpu_i], test.center, test.unknowns[0][test.p_tests]);
                        gpu_batch.open_off_i[gpu_i] = test.p_tests++;
                        if ((size_t) test.p_tests == test.unknowns[0].size()) break;
                        // TODO handle overflow of unknowns[0]
                    } else {
                        assert(!test.n_found);
                        gpu_batch.open_pn[gpu_i] = 1;
                        mpz_add_ui(*gpu_batch.z[gpu_i], test.center, test.unknowns[1][test.n_tests]);
                        gpu_batch.open_off_i[gpu_i] = test.n_tests++;
                        if ((size_t)test.n_tests == test.unknowns[1].size()) break;
                        // TODO handle overflow of unknowns[1]
                    }
        //            gmp_printf("batch[%d] = %d,%d = %d | %Zd\n", gpu_i, i, j, test.m, *gpu_batch.z[gpu_i]);
                }
            }
        }

        // Run this large batch
        {
            typedef mr_params_t<THREADS_PER_INSTANCE, BITS, WINDOW_BITS> params;
            run_test<params>(gpu_batch.z, gpu_batch.result, ROUNDS);
        }

        // Decode and possible finalize M
        {
            for (size_t i = 0; i < BATCH_GPU; i++) {
                if (!gpu_batch.active[i]) {
                    continue;
                }
                assert (gpu_batch.result[i] == 0 || gpu_batch.result[i] == 1);

                /*
                size_t open_i = gpu_batch.open_i[i];
                if (open[open_i]) {
                    const BatchedM &test = *open[open_i];
                    bool n_side = gpu_batch.open_pn[i];
                    int offset = test.unknowns[n_side][gpu_batch.open_off_i[i]];
                    char sgn = "+-"[n_side];
                    cout << "m=" << test.m << " " << sgn << offset << " | " << gpu_batch.result[i] << endl;
                }  // */

                if (gpu_batch.result[i]) {
                    size_t open_i = gpu_batch.open_i[i];
                    if (!open[open_i]) {
                        //cout << "Found two primes for open[" << open_i << "] from batch of " << INSTANCES_PER_M << endl;
                        continue;
                    }

                    if (gpu_batch.open_pn[i] == 0) {
                        if (open[open_i]->p_found) {
                            /*
                            cout << "Found two previous primes for m=" << open[open_i]->m << endl;
                            cout << "\t" << open[open_i]->p_i << " " << gpu_batch.open_off_i[i] << endl;
                            cout << "\t" << open[open_i]->unknowns[0][open[open_i]->p_i] <<
                                    " "  << open[open_i]->unknowns[0][gpu_batch.open_off_i[i]] << endl;
                            */
                            continue;
                        }

                        // prev_prime found
                        assert(open[open_i]->p_tests > 0 );
                        open[open_i]->p_found = true;
                        open[open_i]->p_i = gpu_batch.open_off_i[i];
                    } else {
                        // next_prime found (and done)
                        BatchedM *test = open[open_i].release();
                        open_count--;

                        assert(test->p_found );
                        assert(test->n_tests > 0 );
                        test->n_found = true;
                        test->n_i = gpu_batch.open_off_i[i];

                        int prev_p = test->unknowns[0][test->p_i];
                        int next_p = test->unknowns[1][test->n_i];
                        assert( prev_p > 0 && next_p > 0 );

                        float merit = (next_p + prev_p) / (K_log + log(test->m));
                        if (merit > min_merit)  {
                            // TODO: write to file or database
                            printf("%-5d %.4f  %ld * %ld#/%ld -%d to +%d\n",
                                (next_p + prev_p), merit, test->m, P, D, prev_p, next_p);
                            // TODO: Record finished mi in log file / db.
                        }

                        bool is_last = (mi >= M_inc) && open_count == 0;
                        stats.process_results(config, test->m, is_last,
                            test->unknowns[0].size(), test->unknowns[1].size(),
                            prev_p, next_p,
                            test->p_tests, test->n_tests, merit);
                    }
                }
            }
        }
    }

    // ----- cleanup
    {
        mpz_clear(K);
        for (size_t i = 0; i < BATCH_GPU; i++)
            mpz_clear(*gpu_batch.z[i]);
    }
}

