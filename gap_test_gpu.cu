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
#include <ranges>
#include <sstream>
#include <string>
#include <thread>
#include <unistd.h>
#include <unordered_map>
#include <vector>

// pthread_setname_np
#include <pthread.h>

#include <gmp.h>
#include <primesieve.hpp>

#include "gap_common.h"
#include "gap_test_common.h"
#include "combined_sieve_small.h"
#include "miller_rabin.h"

using std::cout;
using std::cerr;
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
 *  1025,   16384,  1,  8,  1
 *  ???2048,    4096,  2,  8,  1
 *  ???4096,    2048,  4,  16, 1
 *
 */
const size_t BATCH_GPU = 1024; //2*8192;
const size_t SEQUENTIAL_IN_BATCH = 1;
const size_t BATCHED_M = 2 * BATCH_GPU * 130 / 100 / SEQUENTIAL_IN_BATCH;  // 30% extra for overflow

// 1080Ti - 1024, 2, 16 -> 170-180K!
// A100   - 4096, 4, 8 -> 327K!


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
    // TODO renable these checks
    /*
    {
        set_defaults(config);

        // Both shouldn't be true from gap_common.
        assert(!(config.save_unknowns && config.testing));

        if (!config.save_unknowns && !config.testing) {
            cout << "Must set --save-unknowns" << endl;
            exit(1);
        }

        if (config.sieve_length < 6 * config.p || config.sieve_length > 22 * config.p) {
            int sl_min = ((config.p * 8 - 1) / 500 + 1) * 500;
            int sl_max = ((config.p * 20 - 1) / 500 + 1) * 500;
            printf("--sieve_length(%d) should be between [%d, %d]\n",
                config.sieve_length, sl_min, sl_max);
            exit(1);
        }

        if (config.valid == 0) {
            Args::show_usage(argv[0], Args::Pr::SIEVE);
            exit(1);
        }


        if (config.max_prime > 100'000'000) {
            printf("\tmax_prime(%ldM) is probably too large\n",
                config.max_prime / 1'000'000);
        }

        if (config.save_unknowns) {
            std::string fn = Args::gen_unknown_fn(config, ".txt");
            std::ifstream f(fn);
            if (f.good()) {
                printf("\nOutput file '%s' already exists\n", fn.c_str());
                exit(1);
            }
        }

        if (!config.compression) {
            printf("\tSetting --rle to prevent very large output file\n");
            config.compression = 1;
        }

        if (config.compression != 1) {
            cout << "only --rle compression is supported" << endl;
            exit(1);
        }
    }
    */

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
        printf("Testing m * %d#/%d\n", config.p, config.d);
    }

    if (config.mskip != 0) {
        cout << "Use m_start not mskip" << endl;
        return 1;
    }

    setlocale(LC_NUMERIC, "C");

    // Always use RLE
    config.compression = 1;

    // Delete unknown_fn so it gets recreated by sieve
    config.unknown_filename = "";
    prime_gap_test(config);
}


class GPUBatch {
    public:
        enum State { EMPTY, READY, RESULT_WRITTEN };
        State state = EMPTY;

        // current index;
        int i;

        // true if running out of m (used to suppress "No results ready" warning)
        bool end_of_file;

        // number to check if prime
        vector<mpz_t*> z;
        // XXX: This is an ugly hack because you can't create mpz_t vector easily
        mpz_t *z_array;

        // If z[i] should be tested
        vector<bool>  active;
        // Result from GPU
        vector<int>  result;

        // index into 'processing', value is m
        vector<int64_t> data_m;

        // index into the interval[data_m].unknowns
        vector<int> unknown_i;

        GPUBatch(size_t n) {
            elements = n;

            z_array = (mpz_t *) malloc(n * sizeof(mpz_t));
            for (size_t i = 0; i < n; i++) {
                mpz_init(z_array[i]);
                z.push_back(&z_array[i]);
            }

            active.resize(n, 0);
            result.resize(n, -1);

            data_m.resize(n, -1);
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
        DataM(long m, unsigned short sl): m(m), sieve_length(sl) {};

        const long m;
        const unsigned short sieve_length = 0;

        // READY means that load thread should pull next unknown and run it
        // RUNNING means that next unkown is currently running on GPU
        // OVERFLOWED means that load thread ran out of unknowns to load
        // SIEVING means that cpu sieve thread is currently processing
        // PRIMED means that a prime endpoint (next_p or prev_p has just been found)
        enum State : uint8_t { READY, RUNNING, OVERFLOWED, SIEVING, PRIMED};
        State state = READY;

        enum Side : uint8_t { NEXT_P, PREV_P };
        Side side = NEXT_P;

        mpz_t center;
        int32_t unknown_start[2];
        vector<int32_t> unknowns[2];
        // Index into unknowns
        uint8_t p_index = 0;
        uint8_t n_index = 0;

        bool p_found = false;
        bool n_found = false;

        // Both are recored as positive
        int prev_p = 0;
        int next_p = 0;

        // Number of tests made so far
        size_t p_tests = 0;
        size_t n_tests = 0;
};


class SieveResult {
    public:
        SieveResult(const struct Config config): config(config) {};

        enum State { CONFIGURED, SIEVED, TESTING };
        State state = CONFIGURED;

        struct Config config;

        std::unique_ptr<SieveOutput> sieved;
};


/** Shared state between threads */
std::atomic<bool> is_running;

/**
 * Note: Uses a double batched system
 * C++ Thread is preparing batch_a (even more m)
 * While GPU runs batch_b
 */
vector<GPUBatch> batches = {{BATCH_GPU}, {BATCH_GPU}};

std::mutex sieve_mtx;
vector<std::shared_ptr<SieveResult>> sieveds;

std::mutex interval_mtx;
std::condition_variable overflow_cv;
vector<std::shared_ptr<DataM>> overflowed;


void run_gpu_thread(int verbose) {
    pthread_setname_np(pthread_self(), "RUN_GPU_THREAD");

    // XXX: params1024, params2048 with *runner1024, *runner2048 and only new one of them.
    typedef mr_params_t<THREADS_PER_INSTANCE, BITS, WINDOW_BITS> params;
    test_runner_t<params> runner(BATCH_GPU, ROUNDS);

    size_t processed_batches = 0;
    size_t no_batch_count_ms = 0;
    auto last_finished = high_resolution_clock::now();
    bool is_close_to_end = false;
    while (is_running) {
        bool no_batch = true;
        for (GPUBatch& batch : batches) {
            //printf("BB %d\n", (int)batch.state);
            if (batch.state == GPUBatch::State::READY) {
                if (batch.i != BATCH_GPU) {
                    if (!batch.end_of_file) {
                        printf("Partial batch %d vs %lu\n", batch.i, BATCH_GPU);
                    }
                }
                assert(std::count(batch.active.begin(), batch.active.end(), 1) == batch.i);
                // Run batch on GPU and wait for results to be set
                if (true) {
                    runner.run_test(batch.z, batch.result);
                } else {
                    // Return true for 1/10 results (helps not overflow sieve)
                    for (size_t gpu_i = 0; gpu_i < BATCH_GPU; gpu_i++) {
                        if (batch.active[gpu_i]) {
                            batch.result[gpu_i] = (std::rand() % 10) == 1;
                        }
                    }
                }
                batch.state = GPUBatch::State::RESULT_WRITTEN;
                no_batch = false;
                processed_batches++;
                is_close_to_end = batch.end_of_file;
                last_finished = high_resolution_clock::now();
            }
        }
        if (no_batch) {
            // Waiting doesn't count till 1st batch is ready
            auto now = high_resolution_clock::now();
            auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(
                    now - last_finished).count();
            if (verbose >= 0 && processed_batches > 0 && delta > 50) {
                // When we run out of m's it's normal for batches to be empty.
                if (!is_close_to_end) {
                    printf("No results ready for batch %ld. Total wait %.1f seconds\n",
                            processed_batches, no_batch_count_ms / 1000.0);
                }
            }
            uint32_t sleep_ms = delta < 10 ? 1 : (delta < 100 ? 5 : 20);
            no_batch_count_ms += sleep_ms;
            usleep(sleep_ms * 1000);
        }
    }

    if (verbose >= 1) {
        printf("Processed %'ld batches\n", processed_batches);
    }
}

void run_sieve_thread(void) {
    pthread_setname_np(pthread_self(), "SIEVE_THREAD");

    std::unique_lock<std::mutex> lock(interval_mtx, std::defer_lock);
    while (is_running) {
        std::shared_ptr<SieveResult> to_sieve = nullptr;

        lock.lock();
        // Check if any config to start processing
        for (auto& sieved : sieveds) {
            if (sieved->state == SieveResult::CONFIGURED) {
                to_sieve = sieved;
                break;
            }
        }
        lock.unlock();

        if (to_sieve == nullptr) {
            usleep(1'000'000); // 1,000ms
            continue;
        }

        //auto M_end = to_sieve->config.mstart + to_sieve->config.minc;
        //printf("\tStarting Combined Sieve %ld to %ld\n", to_sieve->config.mstart, M_end);
        auto s_save_t = high_resolution_clock::now();
        // DO MORE FOR FIRST RUN THEN LESS AFTER
        to_sieve->config.threads = 6;
        to_sieve->config.verbose -= 3;
        auto result = prime_gap_parallel(to_sieve->config);
        to_sieve->config.verbose += 3;
        auto s_stop_t = high_resolution_clock::now();
        printf("\tCombined Sieve took %.1f seconds\n",
               duration<double>(s_stop_t - s_save_t).count());

        lock.lock();
        to_sieve->state = SieveResult::SIEVED;
        to_sieve->sieved.swap(result);
        lock.unlock();
    }
}

void run_overflow_thread(const mpz_t &K_in) {
    pthread_setname_np(pthread_self(), "CPU_OVERFLOW_SIEVE_THREAD");
    mpz_t K;
    mpz_init_set(K, K_in);

    // TODO make this a #define or based on config or something
    // TODO how do I get a reference to K
    std::vector<std::pair<uint32_t, uint32_t>> p_and_r;
    primesieve::iterator iter;
    uint64_t prime = iter.next_prime();
    assert (prime == 2);  // we skip 2 which is the oddest prime.
    for (prime = iter.next_prime(); ; prime = iter.next_prime()) {
        const uint32_t base_r = mpz_fdiv_ui(K, prime);
        p_and_r.emplace_back((uint32_t) prime, base_r);
    }

    std::unique_lock<std::mutex> lock(interval_mtx, std::defer_lock);

    size_t tested = 0;
    // TODO where else should # of overflows be computed?

    while (is_running) {
        // Important so that overflow_cv / unlock waits correctly
        if (!lock.owns_lock())
            lock.lock();
        // Lock IS NOT held while waiting.
        overflow_cv.wait(lock, []{ return overflowed.size() || !is_running; });

        while (overflowed.size()) {
            DataM& interval = *overflowed.back(); overflowed.pop_back();
            assert(interval.state == DataM::State::OVERFLOWED);
            interval.state = DataM::State::SIEVING;
            lock.unlock();  // Allow main thread to add more things while we process
            // We are safeish to modify interval (expect state) without the lock.

            // One and only one is overflowed at a time
            assert((interval.next_p == -1) ^ (interval.prev_p == -1));

            int32_t sieve_start = 0;
            unsigned sl = interval.sieve_length;

            if (interval.side == DataM::Side::PREV_P) {
                assert( interval.prev_p == -1 );
                interval.unknown_start[0] -= sl;
                // sieves always start at bottom of the range to make sieve math easier
                sieve_start = interval.unknown_start[0] - sl;
            } else {
                assert( interval.next_p == -1 && interval.n_tests > 0);
                interval.unknown_start[1] += sl;
                sieve_start = interval.unknown_start[1];
            }

            vector<int32_t> unknowns;
            sieve_interval_cpu(interval.m, K, p_and_r, sieve_start, sl, unknowns);

            if (interval.prev_p == -1) {
                interval.unknowns[0].clear();
                for (auto u : unknowns) {
                    interval.unknowns[0].push_back(sieve_start + u);
                }
            } else {
                // clear existing unknowns and replace
                interval.unknowns[1].clear();
                assert(sieve_start < 0);
                for (auto iter = unknowns.rbegin(); iter != unknowns.rend(); iter++) {
                    auto u = *iter;
                    assert(u > 0);
                    interval.unknowns[1].push_back(sieve_start + u);
                }
            }

            lock.lock();
            interval.state = DataM::State::READY;
            tested += 1;
        }
    }

    cout << "\tOverflowed " << tested << " intervals" << endl;
    mpz_clear(K);
}


void batch_run_config(
        std::shared_ptr<const SieveResult> result,
        const size_t QUEUE_SIZE,
        StatsCounters& stats) {
    // So we can modify
    Config config = result->config;
    const auto& output = *result->sieved;

    mpz_t K;
    double K_log;

    std::unordered_map<int64_t, std::shared_ptr<DataM>> processing;

    const uint64_t P = config.p;
    const uint64_t D = config.d;
    const uint64_t M_start = config.mstart;
    const uint64_t M_inc = config.minc;

    const float min_merit = config.min_merit;
    // See THEORY.md!
    // Added const is small preference for doing less prev_p
    const float MIN_MERIT_TO_CONTINUE = 2.6 + std::log2(min_merit * std::log(2) + 1);

    // Print Header info & Open unknown_fn
    {

        // ----- Merit / Sieve stats
        int tmp_verbose = config.verbose;
        if (stats.s_tests > 0) config.verbose = 0;
        K_log = prob_prime_and_stats(config, K);
        if (stats.s_tests > 0) config.verbose = tmp_verbose;
        {
            float m_log = log(M_start);
            if (stats.s_tests == 0 && config.verbose >= 1) {
                printf("Min Gap ~= %d (for merit > %.1f)\n",
                    (int) (min_merit * (K_log + m_log)), min_merit);
                printf("Min Gap to continue ~= %d (for merit = %.1f)\n",
                    (int) (MIN_MERIT_TO_CONTINUE * (K_log + m_log)),
                    MIN_MERIT_TO_CONTINUE);
            }
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

    // Two delta trackers for SieveOutput
    size_t result_i = 0;
    uint64_t result_m = output.m_start;

    // Main loop
    uint64_t mi = 0;
    while (mi < M_inc || !processing.empty()) {
        usleep(100); // 0.1ms

        bool any_pushed_to_overflow = false;
        for (GPUBatch& batch : batches) {
            // If batch is ready to have new data loaded
            if (batch.state == GPUBatch::State::EMPTY) {
                // Add new DataM if free space in processing
                for (; processing.size() < QUEUE_SIZE && mi < M_inc; mi++) {
                    uint64_t m = M_start + mi;
                    if (gcd(m, D) > 1) continue;

                    result_m += std::get<0>(output.m_inc[result_i]);
                    int64_t found = std::get<1>(output.m_inc[result_i]);
                    const auto m_unknown_deltas = output.unknowns[result_i];
                    result_i++;

                    assert(result_m == m);

                    auto test = std::make_shared<DataM>(m, config.sieve_length);

                    test->unknowns[1].reserve(m_unknown_deltas.size());
                    int32_t offset = 0;
                    for (uint16_t delta : m_unknown_deltas) {
                        offset += delta;
                        test->unknowns[1].push_back(offset);
                    }

                    assert(m_unknown_deltas.size() > 2);
                    assert(test->unknowns[1].size() == m_unknown_deltas.size());
                    assert(test->unknowns[1].back() <= config.sieve_length);

                    mpz_init(test->center);
                    mpz_mul_ui(test->center, K, test->m);

                    processing[test->m] = test;
                }

                // Not a full count but (only covers the portion till GPU_BATCH is full).
                size_t currently_overflowed = 0;

                // Grap some entries from each item in M
                {
                    batch.i = 0;
                    // Turn off all entries in batch
                    std::fill_n(batch.active.begin(), BATCH_GPU, false);
                    // Mark all results as invalid
                    std::fill_n(batch.result.begin(), BATCH_GPU, -1);
                    for (size_t i = 0; i < BATCH_GPU; i++) {
                        // This prevents the GPU from stalling in partial batches
                        mpz_set_ui(*batch.z[i], 7);
                    }

                    {
                        std::unique_lock<std::mutex> lock(interval_mtx);
                        for (auto& pair : processing) {
                            auto& interval = *pair.second;
                            currently_overflowed += DataM::State::OVERFLOWED || DataM::State::SIEVING;
                            if (interval.state != DataM::State::READY) {
                                continue;  // Not ready to be part of a batch.
                            }

                            int32_t last_offset = 0;

                            for (size_t j = 0; j < SEQUENTIAL_IN_BATCH; j++) {
                                // interval shouldn't have any ends (as only 1 side is tested in this loop);
                                assert(!interval.n_found);
                                if (interval.n_index == interval.unknowns[1].size()) {
                                    /* Don't push to overflow queue while GPU still has sieve to check.
                                     * leads to race condition with overflow thread & uglier code.
                                     */
                                    if (j == 0) {
                                        interval.next_p = -1;
                                        interval.state = DataM::State::OVERFLOWED;
                                        stats.s_gap_out_of_sieve_next += 1;
                                        overflowed.push_back(pair.second);
                                        any_pushed_to_overflow = true;
                                    }
                                    break;
                                }

                                interval.state = DataM::State::RUNNING;
                                int gpu_i = batch.i;  // [GPU] batch index
                                auto offset = interval.unknowns[1][interval.n_index];
                                //printf("%d = %d,%d = %d*K+%d\n", gpu_i, i, j, interval.m, offset);

                                // entries should grow away from 0 in magnitude
                                assert(abs(last_offset) < abs(offset));
                                last_offset = offset;

                                mpz_add_ui(*batch.z[gpu_i], interval.center, offset);
                                batch.data_m[gpu_i] = interval.m;  // [Data] index for GPU Batch
                                batch.unknown_i[gpu_i] = interval.n_index++;
                                batch.active[gpu_i] = true;

                                batch.i++;
                                if (batch.i == BATCH_GPU) break;
                            }
                            if (batch.i == BATCH_GPU) break;
                        }
                    }

                    // Every batch should be full unless we are almost done
                    // technically if many overflowed results this could not be true.
                    if (!( (mi >= M_inc) || (batch.i == BATCH_GPU) )) {
                        printf("Partial load @ %lu/%lu -> %d/%lu | %lu\n",
                               mi, M_inc, batch.i, BATCH_GPU, currently_overflowed);
                    }
                }

                // TODO maybe handle this differently
                // Mark batch as ready for GPU processing
                batch.end_of_file = (mi == M_inc);
                // Needs to be last so that GPU_thread doesn't read other parts early.
                batch.state = batch.i > 0 ? GPUBatch::State::READY : GPUBatch::State::EMPTY;

            } else if (batch.state == GPUBatch::State::RESULT_WRITTEN) {
                // If GPU PRP result are ready:
                // Read results, mark any found primes, and possible finalize m-interval
                {
                    for (size_t i = 0; i < BATCH_GPU; i++) {
                        if (!batch.active[i]) {
                            continue;
                        }
                        // Verify GPU really did write the result
                        assert (batch.result[i] == 0 || batch.result[i] == 1);

                        // Doesn't lock as none of these will be overflows
                        DataM &interval = *processing.at(batch.data_m[i]);
                        assert( interval.state == DataM::State::RUNNING );

                        int side_i;
                        if (interval.side == DataM::Side::PREV_P) {
                            side_i = 0;
                            interval.p_tests++;
                        } else {
                            side_i = 1;
                            interval.n_tests++;
                        }

                        if (!batch.result[i]) {
                            // Mark interval as being ready again
                            interval.state = DataM::State::READY;
                        } else {
                            // Found prime!

                            // next_prime found (and done)
                            int offset_i = batch.unknown_i[i];

                            int found_offset = interval.unknowns[side_i][offset_i];

                            // Two primes happens fairly regularly with SEQUENTIAL_IN_BATCH > 1
                            // Ignore the further away prime
                            if (side_i == 0) {
                                if (interval.n_found) {
                                    assert( SEQUENTIAL_IN_BATCH > 1 );
                                    assert( interval.next_p < found_offset );
                                    continue;
                                }
                                assert(interval.next_p == 0);
                                interval.n_found = true;
                                interval.next_p = found_offset;
                            } else {
                                if (interval.p_found) {
                                    assert( SEQUENTIAL_IN_BATCH > 1 );
                                    assert( abs(interval.prev_p) < abs(found_offset) );
                                    continue;
                                }
                                assert(interval.prev_p == 0);
                                interval.p_found = true;
                                interval.prev_p = found_offset;
                            }
                            interval.state = DataM::State::PRIMED;
                        }
                    }
                }

                // Result batch to EMPTY
                batch.state = GPUBatch::State::EMPTY;
            }


            if (any_pushed_to_overflow) {
                // CPU sieving thread will start if unlocked and notified
                overflow_cv.notify_one();
                // TODO print warning if overflowed.size() is very large
            }

            {
                // Process & remove finished intervals.
                std::unique_lock<std::mutex> lock(interval_mtx);
                for (auto it = processing.begin(); it != processing.end(); /* increment in loop */) {
                    auto& interval = *it->second;

                    if (interval.state != DataM::State::PRIMED) {
                        continue;
                    }

                    int prev_p = interval.prev_p;
                    int next_p = interval.next_p;

                    assert(interval.n_found);
                    assert(interval.n_tests > 0);
                    assert(next_p > 0);

                    // Potentially do Side-Skip if next_p is not very large (and just found)
                    if (!interval.p_found) {
                        assert( interval.p_tests == 0 );
                        float next_merit = next_p / (K_log + log(interval.m));

                        if (next_merit >= MIN_MERIT_TO_CONTINUE) {
                            // Mark this for prev_p sieve
                            interval.side = DataM::Side::PREV_P;
                            interval.prev_p = -1;
                            interval.state = DataM::State::OVERFLOWED;
                            continue;
                        }

                        assert(prev_p == 0);
                        stats.s_skips_after_one_side += 1;
                    } else {
                        assert(interval.p_found);
                        assert(interval.p_tests > 0);
                        assert(prev_p > 0);
                    }

                    float merit = (next_p + prev_p) / (K_log + log(interval.m));
                    if (merit > min_merit)  {
                        // TODO: Record finished mi in log file / db.
                        printf("%-5d %.4f  %ld * %ld#/%ld -%d to +%d\n",
                            (next_p + prev_p), merit, interval.m, P, D, prev_p, next_p);
                    }

                    bool is_last = (mi >= M_inc) && processing.size() == 1;
                    // TODO this doesn't seem to always be working now
                    stats.process_results(config, interval.m, is_last,
                        interval.unknowns[0].size(), interval.unknowns[1].size(),
                        prev_p, next_p,
                        interval.p_tests, interval.n_tests, merit);

                    mpz_clear(interval.center);
                    it = processing.erase(it);  // Erase this element
                }
            }
        }
    }

    // ----- cleanup
    {
        mpz_clear(K);
    }
}

void coordinator_thread(struct Config global_config) {
    pthread_setname_np(pthread_self(), "MAIN_LOAD_AND_BATCH_THREAD");

    size_t skipped_ms = 0;

    // Used for various stats
    StatsCounters stats(high_resolution_clock::now());

    // TODO make some get_sieved_by_status method.
    std::unique_lock<std::mutex> lock(interval_mtx, std::defer_lock);
    while (is_running) {
        std::shared_ptr<SieveResult> to_test = nullptr;

        lock.lock();
        {
            int needed = 2;
            for (auto& sieved : sieveds) {
                if (sieved->state == SieveResult::CONFIGURED || sieved->state == SieveResult::SIEVED) {
                    needed -= 1;
                }
            }
            while (needed > 0) {
                // Add new config to the queue
                auto this_config = std::make_shared<SieveResult>(global_config);
                printf("\tQueued %'lu to %'lu\n",
                        global_config.mstart,
                        global_config.mstart + global_config.minc);
                sieveds.push_back(this_config);
                needed -= 1;

                global_config.mstart += global_config.minc;
            }
        }

        for (auto& sieved : sieveds) {
            if (sieved->state == SieveResult::SIEVED) {
                if (!to_test || to_test->config.mstart > sieved->config.mstart) {
                    to_test = sieved;
                }
            }
        }
        if (to_test) {
            to_test->state = SieveResult::TESTING;
        }
        lock.unlock();

        if (to_test == nullptr) {
            skipped_ms += 1000;
            if (stats.s_tests > 0) {
                printf("\tNothing to test (%.1f seconds paused)\n", skipped_ms / 1000.0);
            }
            usleep(1'000'000); // 1,000ms
            continue;
        }

        // TODO: Figure out how to avoid partial batches as each batch finishes up.
        // Have a middle person (this thread?) merge them into a new queue or something.
        printf("\tTesting %'lu to %'lu\n",
                to_test->config.mstart,
                to_test->config.mstart + to_test->config.minc);
        batch_run_config(to_test, BATCHED_M, stats);

        lock.lock();
        sieveds.erase(std::remove(sieveds.begin(), sieveds.end(), to_test), sieveds.end());
        lock.unlock();
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

    mpz_t K;
    {
        init_K(config, K);
        size_t N_bits = mpz_sizeinbase(K, 2) + log2(config.mstart);

        // P# roughly 349, 709, 1063, 1447
        for (size_t bits : {512, 1024, 1536, 2048, 3036, 4096}) {
            if (N_bits <= bits) {
                if (bits < BITS) {
                    printf("\nFASTER WITH `make gap_test_gpu BITS=%ld` (may require `make clean`)\n\n", bits);
                }
                break;
            }
        }
        assert( (N_bits + 1) < BITS ); // See last debug line.
        assert( BITS <= (1 << (2 * WINDOW_BITS)) );
    }

    is_running = true;

    std::thread main_thread(coordinator_thread, config);
    std::thread sieve_thread(run_sieve_thread);
    std::thread gpu_thread(run_gpu_thread, config.verbose);
    std::thread overflow_sieve_thread(run_overflow_thread, &K);

    main_thread.join();
    // Tell other threads to quit
    is_running = false;
    cout << "All done. Goodbye!" << endl;

    sieve_thread.join();
    gpu_thread.join();
    overflow_cv.notify_all();  // wake up all overflow thread
    overflow_sieve_thread.join();
    mpz_clear(K);
}

