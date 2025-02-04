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
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <unistd.h>
#include <unordered_map>
#include <vector>

// pthread_setname_np
#include <pthread.h>

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

/**
 * BATCH_GPU is 2^n >= 1024
 * SEQUENTIAL_IN_BATCH = {1,2,4}
 *      1 => 0 overhead
 *      2 => 0.5 extra PRP/m
 *      4 => 1.5 extra PRP/M
 *
 * BATCHED_M is number of M loaded at the same time
 * ------------------
 * Try:
 *  1024,   16384,  1,  8,  1
 *  ???2048,    4096,  2,  8,  1
 *  ???4096,    2048,  4,  16, 1
 *
 */
const size_t BATCH_GPU = 1024; //2*8192;
const size_t SEQUENTIAL_IN_BATCH = 2;
const size_t BATCHED_M = 2 * BATCH_GPU * 120 / 100 / SEQUENTIAL_IN_BATCH;  // 10% extra

/**
 * Originally 8 which has highest throughput but only if we have LOTS of instances
 * this helps reduce the number of parallel instances needed
 */
const int THREADS_PER_INSTANCE = 16;
const int ROUNDS = 1;

//************************************************************************

void prime_gap_test(const struct Config config);


int main(int argc, char* argv[]) {
    Config config = Args::argparse(argc, argv, Args::Pr::TEST_GPU);
    if (config.valid == 0) {
        Args::show_usage(argv[0], Args::Pr::TEST_GPU);
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
        Args::show_usage(argv[0], Args::Pr::TEST_GPU);
        return 1;
    }

    setlocale(LC_NUMERIC, "");
    if (config.verbose >= 0) {
        printf("\n");
        printf("Testing m * %d#/%d, m = %ld + [0, %'ld)\n",
            config.p, config.d, config.mstart, config.minc);
    }

    if (config.mskip > 0) {
        printf("\tskipping m < %'ld\n", config.mskip);
        assert(config.mskip >= config.mstart);
        assert(config.mskip < (config.mstart + config.minc));
    }

    setlocale(LC_NUMERIC, "C");

    // Determine compression
    {
        std::string fn = Args::gen_unknown_fn(config, ".txt");
        std::ifstream unknown_file(fn, std::ios::in);
        assert( unknown_file.is_open() ); // Can't open save_unknowns file
        assert( unknown_file.good() );    // Can't open save_unknowns file
        config.compression = Args::guess_compression(config, unknown_file);
    }

    prime_gap_test(config);
}


class GPUBatch {
    public:
        enum State { EMPTY, READY, RESULT_WRITTEN };
        State state = EMPTY;

        // current index;
        int i;

        // number to check if prime
        vector<mpz_t*> z;
        // XXX: This is an ugly hack because you can't create mpz_t vector easily
        mpz_t *z_array;

        // If z[i] should be tested
        vector<bool>  active;
        // Result from GPU
        vector<int>  result;

        // index into 'processing' (DataM)
        vector<int64_t> data_i;

        vector<int> unknown_i;

        GPUBatch(size_t n) {
            elements = n;

            z_array = (mpz_t *) malloc(n * sizeof(mpz_t));
            for (size_t i = 0; i < n; i++) {
                // MAYBE FIXES MY STALL ISSUE?
                mpz_init(z_array[i]);
                z.push_back(&z_array[i]);
            }

            active.resize(n, 0);
            result.resize(n, -1);

            data_i.resize(n, -1);
            unknown_i.resize(n, -1);
        }

        ~GPUBatch() {
            for (size_t i = 0; i < elements; i++) {
                mpz_clear(z_array[i]);
            }
        }

    private:
        size_t elements;
};

class DataM {
    public:
        /**
         * Elements in READY state can ONLY be modified by load_thread
         * Elements in RUNNING are either part of a GPU batch in overflowed queue
         */
        DataM() {};
        DataM(long m): m(m) {};

        enum State { READY, RUNNING, OVERFLOW_DONE };
        State state = READY;

        long m;
        mpz_t center;
        vector<int32_t> unknowns[2];

        bool p_found = false, n_found = false;
        int prev_p = 0, next_p = 0;

        // if this entry needs to be handled manually
        bool overflow = false;

        size_t p_tests = 0;
        size_t n_tests = 0;
};


/** Shared state between threads */
std::atomic<bool> is_running;

/**
 * Note: Uses a double batched system
 * C++ Thread is preparing batch_a (even more m)
 * While GPU runs batch_b
 */
vector<GPUBatch> batches = {{BATCH_GPU}, {BATCH_GPU}};

std::mutex overflow_mtx;
std::condition_variable overflow_cv;
vector<std::shared_ptr<DataM>> overflowed;

void run_gpu_thread(const struct Config config) {
    pthread_setname_np(pthread_self(), "RUN_GPU_THREAD");

    // XXX: params1024, params2048 with *runner1024, *runner2048 and only new one of them.
    typedef mr_params_t<THREADS_PER_INSTANCE, BITS, WINDOW_BITS> params;
    test_runner_t<params> runner(BATCH_GPU, ROUNDS);

    size_t processed_batches = 0;
    size_t no_batch_count_ms = 0;
    while (is_running) {
        bool no_batch = true;
        for (GPUBatch& batch : batches) {
            if (batch.state == GPUBatch::State::READY) {
                if (batch.i != BATCH_GPU) {
                    size_t test_active = 0;
                    for (size_t gpu_i = 0; gpu_i < BATCH_GPU; gpu_i++) {
                        test_active += batch.active[gpu_i];
                        if (!batch.active[gpu_i]) {
                            // This prevents the GPU from stalling if z was never initalized.
                            mpz_set_ui(*batch.z[gpu_i], 7);
                        }
                    }
                    printf("Partial batch %d/%ld actual: %lu\n", batch.i, BATCH_GPU, test_active);
                }
                // Run batch on GPU and wait for results to be set
                runner.run_test(batch.z, batch.result);
                batch.state = GPUBatch::State::RESULT_WRITTEN;
                no_batch = false;
                processed_batches++;
            }
        }
        if (no_batch) {
            // Waiting doesn't count till 1st batch is ready
            if (config.verbose >= 0 && processed_batches > 0) {
                no_batch_count_ms += 100;
                printf("Waiting on batch%ld => %.1f seconds\n",
                        no_batch_count_ms / 100, no_batch_count_ms / 1000.0);
            }
            usleep(250000); // 250ms
        }
    }

    if (config.verbose >= 1) {
        printf("Processed %'ld batches\n", processed_batches);
    }
}

void run_overflow_thread(const struct Config config) {
    mpz_t prime_test;
    mpz_init(prime_test);

    std::unique_lock<std::mutex> lock(overflow_mtx);

    while (true) {
        overflow_cv.wait(lock, []{ return overflowed.size() || !is_running; });
        if (!is_running) break;

        while (overflowed.size()) {
            DataM& interval = *overflowed.back(); overflowed.pop_back();
            lock.unlock();  // Allow main thread to add more things while we process
            assert (interval.overflow && interval.state == DataM::State::RUNNING);

            // NOTE: Overhead to doing this while GPU waits seems small (<1% of candidates)
            // But is actually A LOT because 40x slower. Becomes ~20-40% overhead quickly.

            if (interval.next_p == -1) {
                assert(interval.n_tests > 0);

                //cout << "gap_out_of_sieve_next m=" << interval.m << endl;
                mpz_add_ui(prime_test, interval.center, config.sieve_length);
                mpz_nextprime(prime_test, prime_test);
                mpz_sub(prime_test, prime_test, interval.center);

                interval.next_p = mpz_get_ui(prime_test);
                //cout << "gap_out_of_sieve_next m=" << interval.m << " -> " << interval.next_p << endl;
                interval.n_found = true;
                interval.overflow = 0;
            }
            if (interval.prev_p == -1) {
                // It took two years to mpz_prevprime into gmp.
                // I'm so proud and excited to get to use it here.

                assert(interval.p_tests == 0);
                //cout << "gap_out_of_sieve_prev m=" << interval.m << endl;
                mpz_prevprime(prime_test, interval.center);
                mpz_sub(prime_test, interval.center, prime_test);

                interval.prev_p = mpz_get_ui(prime_test);
                //cout << "gap_out_of_sieve_prev m=" << interval.m << " -> " << interval.prev_p << endl;
                interval.p_found = true;
                interval.overflow = 0;
            }

            // Mark interval as finished processing
            // NOTE: don't mark as READY or race_condition can happen in load
            interval.state = DataM::State::OVERFLOW_DONE;

            lock.lock(); // Lock so that overflow_cv / unlock waits correctly
        }
    }

    mpz_clear(prime_test);
}

void load_batch_thread(const struct Config config, const size_t QUEUE_SIZE) {
    // TODO ask C++ person if I need to worry about CPU doing cache invalidation with this setup
    // if batch is RESULT_WRITTEN | read result back to DataM processing | update to EMPTY
    // if batch is EMPTY          | load data from DataM processing      | update to READY, unlock GPU thread
    // if all batches EMPTY, wait(thread_sync)

    mpz_t K;
    double K_log;
    std::ifstream unknown_file;

    // Used for various stats
    StatsCounters stats(high_resolution_clock::now());

    std::unordered_map<int64_t, std::shared_ptr<DataM>> processing;

    const uint64_t P = config.p;
    const uint64_t D = config.d;
    const uint64_t M_start = config.mstart;
    const uint64_t M_inc = config.minc;

    const float min_merit = config.min_merit;

    // Print Header info & Open unknown_fn
    {

        // ----- Merit / Sieve stats
        K_log = prob_prime_and_stats(config, K);
        {
            float m_log = log(M_start);
            if (config.verbose >= 1) {
                printf("Min Gap ~= %d (for merit > %.1f)\n",
                    (int) (min_merit * (K_log + m_log)), min_merit);
            }
        }

        // ----- Open unknown input file
        {
            std::string fn = Args::gen_unknown_fn(config, ".txt");
            if (config.verbose >= 1) {
                printf("\nReading unknowns from '%s'\n", fn.c_str());
            }
            unknown_file.open(fn, std::ios::in);
            assert( unknown_file.is_open() ); // Can't open save_unknowns file
            assert( unknown_file.good() );    // Can't open save_unknowns file
        }

        uint64_t first_mi = 0;
        for (; first_mi > 0 && gcd(M_start + first_mi, D) > 1; first_mi++);
        assert(first_mi < M_inc);

        uint64_t last_mi = M_inc - 1;
        for (; last_mi > 0 && gcd(M_start + last_mi, D) > 1; last_mi--);
        assert(last_mi > 0 && last_mi < M_inc);

        // ----- Main sieve loop.
        if (config.verbose >= 1) {
            uint64_t valid_ms = count_num_m(M_start, M_inc, D);
            assert(valid_ms > 0 && valid_ms <= M_inc);

            printf("\n%ld tests M_start(%ld) + mi(%ld to %ld)\n\n",
                valid_ms, M_start, first_mi, last_mi);
        }
    }

    // For compressed lines
    BitArrayHelper helper(config, K);

    // Main loop
    uint64_t mi = 0;
    while (mi < M_inc || !processing.empty()) {
        usleep(500); // 0.5ms
        for (GPUBatch& batch : batches) {
            // If batch is ready to have new data loaded
            if (batch.state == GPUBatch::State::EMPTY) {
                // Add new DataM if free space
                for (; processing.size() < QUEUE_SIZE && mi < M_inc; mi++) {
                    uint64_t m = M_start + mi;
                    if (gcd(m, D) > 1) continue;

                    std::string line;
                    // Loop can be pragma omp parallel if this is placed in critical section
                    std::getline(unknown_file, line);

                    std::istringstream iss_line(line);

                    // Can skip if m < M_RESUME without parsing line here
                    if (m < config.mskip) continue;

                    auto test = std::make_shared<DataM>(m);

                    uint64_t m_parsed = parse_unknown_line(
                        config, helper, m, iss_line, test->unknowns[0], test->unknowns[1]);
                    assert(m_parsed == (uint64_t) m);

                    mpz_init(test->center);
                    mpz_mul_ui(test->center, K, test->m);

                    processing[test->m] = test;
                }

                // Grap some entries from each item in M
                {
                    batch.i = 0;
                    // Turn off all entries in batch
                    std::fill_n(batch.active.begin(), BATCH_GPU, false);
                    // Mark all results as invalid
                    std::fill_n(batch.result.begin(), BATCH_GPU, -1);

                    for (auto& pair : processing) {
                        auto& interval = *pair.second;
                        if (interval.state != DataM::State::READY or interval.overflow) {
                            // Already part of some other batch
                            continue;
                        }

                        for (size_t j = 0; j < SEQUENTIAL_IN_BATCH; j++) {
                            assert(! (interval.p_found && interval.n_found) );

                            int gpu_i = batch.i;  // [GPU] batch index
                            batch.data_i[gpu_i] = interval.m;  // [Data] index for GPU Batch

                            // One sided only runs positive side.
                            assert(!interval.n_found);
                            if (interval.n_tests < interval.unknowns[1].size()) {
                                mpz_add_ui(*batch.z[gpu_i], interval.center, interval.unknowns[1][interval.n_tests]);
                                batch.unknown_i[gpu_i] = interval.n_tests++;
                            } else {
                                // Haven't found next prime, but run out of unknowns to test
                                interval.next_p = -1;
                                interval.overflow = 1; // Indicates next side has overflowed unknowns
                                break;
                            }

                            //gmp_printf("batch[%d] = %d,%d = %d | %Zd\n", gpu_i, i, j, interval.m, *batch.z[gpu_i]);
                            interval.state = DataM::State::RUNNING;
                            batch.active[gpu_i] = true;
                            batch.i++;
                            if (batch.i == BATCH_GPU) break;
                        }
                        if (batch.i == BATCH_GPU) break;
                    }

                    // Every batch should be full unless we are almost done
                    // technically if many overflowed results this could not be true.
                    assert( (mi >= M_inc) || (batch.i == BATCH_GPU) );
                }

                // Mark batch as ready for GPU processing
                batch.state = GPUBatch::State::READY;
            }

            // If PRP result has been written to all entries by GPU
            if (batch.state == GPUBatch::State::RESULT_WRITTEN) {
                // Read results, mark any found primes, and possible finalize m-interval
                {
                    for (size_t i = 0; i < BATCH_GPU; i++) {
                        if (!batch.active[i]) {
                            continue;
                        }
                        // Verify GPU really did write the result
                        assert (batch.result[i] == 0 || batch.result[i] == 1);

                        DataM &interval = *processing.at(batch.data_i[i]);
                        // Mark interval as being ready again
                        interval.state = DataM::State::READY;

                        if (batch.result[i]) {
                            // Found prime in last partial batch of unknowns, no longer overflowed
                            interval.overflow = 0;

                            int offset_i = batch.unknown_i[i];
                            if (interval.n_found) {
                                /*
                                cout << "Found two next primes for m=" << interval.m << endl;
                                cout << "\t" << interval.next_p << " vs "
                                     << interval.unknowns[1][offset_i] << "(" << offset_i << ")" << endl;
                                */
                                continue;
                            }

                            // next_prime found (and done)
                            assert(interval.n_tests > 0 );
                            interval.n_found = true;
                            interval.next_p = interval.unknowns[1][offset_i];
                        }
                    }
                }

                // Finalize any finished (or overflowed) results from processing
                {
                    // Push Out-Of-Sieve gaps to overflow queue and notify that thread
                    {
                        bool pushed_to_overflow = false;
                        for (auto& pair : processing) {
                            auto& interval = *pair.second;
                            if (interval.overflow && interval.state == DataM::State::READY) {
                                if (interval.next_p == -1) {
                                    assert(interval.n_tests > 0);
                                    stats.s_gap_out_of_sieve_next += 1;
                                }

                                if (interval.prev_p == -1) {
                                    assert(interval.n_found);
                                    stats.s_gap_out_of_sieve_prev += 1;
                                }

                                // Push to overflow and wake up that thread
                                interval.state = DataM::State::RUNNING;
                                {
                                    std::unique_lock<std::mutex> lock(overflow_mtx);
                                    // TODO I THINK THERE'S BAD WRITES IF OVERFLOWED CHANGES POINTS AND STUFF
                                    // MAYBE I CAN MAKE SHARED POINTERS OR SOMETHING
                                    // OR I CAN PUSH TO A DIFFERENT QUEUE TBD NOT SURE
                                    overflowed.push_back(pair.second);
                                    pushed_to_overflow = true;
                                }
                            }
                        }
                        // TODO print warning if overflowed.size() is very large
                        if (pushed_to_overflow) {
                            overflow_cv.notify_one();
                        }
                    }

                    {
                        // Update any items finished in overflow as ready to be loaded into batches again
                        for (auto& pair : processing) {
                            auto& interval = *pair.second;
                            if (interval.state == DataM::State::OVERFLOW_DONE) {
                                interval.state = DataM::State::READY;
                            }
                        }
                    }

                    // Ugly code that allows for remove during iteration
                    auto it = processing.begin();
                    while (it != processing.end()) {
                        auto& interval = *it->second;

                        int prev_p = interval.prev_p;
                        int next_p = interval.next_p;

                        // Potentially do Side-Skip if next_p is not very large.
                        // Only consider if next_p just found.
                        if (interval.n_found && interval.prev_p == 0 && !interval.p_found) {
                            // TODO improve this with constant and logging
                            float next_merit = next_p / (K_log + log(interval.m));
                            /**
                             * TODO better math
                             * With Y = 24
                             * 50% of gaps with merit > 24 merit have prev > 12 merit
                             *      only test 1/2^(12-3) = 1/512 gaps
                             * 75% of gaps with merit > 24 merit have prev > 6 merit
                             *      test 1/2^(6-3) = 1/8 gaps
                             */
                            float MIN_MERIT_TO_CONTINUE = min_merit / 2 - 2;

                            if (next_merit < MIN_MERIT_TO_CONTINUE) {
                                stats.s_skips_after_one_side += 1;

                                bool is_last = (mi >= M_inc) && processing.size() == 1;
                                stats.process_results(config, interval.m, is_last,
                                    interval.unknowns[0].size(), interval.unknowns[1].size(),
                                    prev_p, next_p,
                                    interval.p_tests, interval.n_tests, next_merit);

                                mpz_clear(interval.center);
                                it = processing.erase(it);  // Erase this element
                                continue;
                            }

                            //cout << "Queued prev_p for check " << interval.m << endl;

                            // Mark this for overflow
                            // TODO not clear if it's bad that not pushed to overflow here.
                            interval.prev_p = -1;
                            interval.overflow = 1; // Indicates a side has overflowed
                            continue;
                        }

                        if (!interval.p_found || !interval.n_found) {
                            ++it;
                            continue;
                        }
                        assert( prev_p > 0 && next_p > 0 );

                        float merit = (next_p + prev_p) / (K_log + log(interval.m));
                        if (merit > min_merit)  {
                            // TODO: Record finished mi in log file / db.
                            printf("%-5d %.4f  %ld * %ld#/%ld -%d to +%d\n",
                                (next_p + prev_p), merit, interval.m, P, D, prev_p, next_p);
                        }

                        bool is_last = (mi >= M_inc) && processing.size() == 1;
                        stats.process_results(config, interval.m, is_last,
                            interval.unknowns[0].size(), interval.unknowns[1].size(),
                            prev_p, next_p,
                            interval.p_tests, interval.n_tests, merit);

                        mpz_clear(interval.center);
                        it = processing.erase(it);  // Erase this element
                    }
                }

                // Result batch to EMPTY
                batch.state = GPUBatch::State::EMPTY;
            }
        }
    }

    // ----- cleanup
    {
        mpz_clear(K);
    }
}


void prime_gap_test(struct Config config) {
    // Setup test runner
    printf("BITS=%d\tWINDOW_BITS=%d\n", BITS, WINDOW_BITS);
    printf("PRP/BATCH=%ld\tM/BATCH=%ld\n",
            BATCH_GPU, BATCH_GPU/SEQUENTIAL_IN_BATCH);
    printf("THREADS/PRP=%d\n", THREADS_PER_INSTANCE);

    assert( BATCH_GPU == 1024 || BATCH_GPU == 2048 || BATCH_GPU == 4096 ||
            BATCH_GPU == 8192 || BATCH_GPU ==16384 || BATCH_GPU ==32768 );
    assert( SEQUENTIAL_IN_BATCH == 1 || SEQUENTIAL_IN_BATCH == 2 || SEQUENTIAL_IN_BATCH == 4 );

    {
        mpz_t K;
        init_K(config, K);
        size_t N_bits = mpz_sizeinbase(K, 2) + log2(config.mstart + config.minc);
        mpz_clear(K);

        // P# roughly 349, 709, 1063, 1447
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

    is_running = true;

    std::thread load_thread(load_batch_thread, config, BATCHED_M);
    std::thread gpu_thread(run_gpu_thread, config);
    std::thread overflow_thread(run_overflow_thread, config);

    load_thread.join();

    is_running = false;
    overflow_cv.notify_one();  // wake up overflow thread

    gpu_thread.join();
    overflow_thread.join();
}

