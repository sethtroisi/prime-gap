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
#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <deque>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <queue>
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
using std::deque;
using std::vector;
using namespace std::chrono;

#ifdef GPU_BITS
const int BITS = GPU_BITS;
#else
const int BITS = 1024;
#endif

const int WINDOW_BITS = (BITS <= 1024) ? 5 : 6;
const int THREADS_PER_INSTANCE = (BITS <= 512) ? 4 : 8;

/**
 * GPU_BATCHES the number of simultanious batches to create & queue.
 * GPU_BATCH_SIZE is 2^n | best is between 4K and 16K.
 */
const size_t GPU_BATCHES = 2;
const size_t GPU_BATCH_SIZE = 4 * 1024;

/********** BENCHMARKING ***********/
// 701#
// GPU    - BATCH TPI -> PRP/second
// 1080Ti - 2K,   8 -> 250K
// 1080Ti - 4K,   8 -> 266K
// 1080Ti - 8K,   8 -> 282K
// 1080Ti - 16K,  8 ->

// Remember to set sm_80
// A100   - 4K,   8 -> 930K
// A100   - 8K,   8 -> 1050K
// A100   - 16K,  8 -> 1085K!

// 347# as high as 140K m/sec
// 1080Ti - 8K, 4  -> 1680K
// 1080Ti - 4K, 4  -> 1690K!
// 1080Ti - 2K, 4  -> 1330K

// A100   - 8K,   4 ->
// A100   - 16K,  4 -> 2210K (40% utilization)
// A100   - 32K,  4 ->

// 257# as high as 340K m/sec!
// 1080Ti - 2K,  4 -> 2500K!
// 1080Ti - 4K,  4 -> 3300K!
// 1080Ti - 8K,  4 -> 3060K!

// A100   - 16K, 4 ->

/********** BENCHMARKING ***********/


// Always use 1.
const int ROUNDS = 1;
// 20-60% extra for overflow is very reasonable.
const size_t QUEUE_SIZE = 140 * GPU_BATCHES * GPU_BATCH_SIZE / 100;

// From 701# I believe.
// 1M -> 2000/second
// 200K -> 4000/second
const size_t CPU_SIEVE_LIMIT = 70'000;

const size_t COMBINED_SIEVE_THREADS = 8;

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

    if (config.sieve_length == 0) {
        cout << "Must set sieve-length for " << argv[0] << endl;
        Args::show_usage(argv[0], Args::Pr::TEST_GPU);
        return 1;
    }

    if (config.sieve_length < 6 * config.p || config.sieve_length > 22 * config.p) {
        int sl_min = ((config.p * 8 - 1) / 500 + 1) * 500;
        int sl_max = ((config.p * 20 - 1) / 500 + 1) * 500;
        printf("--sieve_length(%d) should be between [%d, %d]\n",
            config.sieve_length, sl_min, sl_max);
        return 1;
    }

    if (100'000 < config.max_prime && config.max_prime > 500'000'000) {
        printf("\tmax_prime(%'ld) should be between 100K and 500M\n", config.max_prime);
    }

    if (config.compression != 0) {
        cout << argv[0] << " Doesn't support any compression options." << endl;
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

    // Delete unknown_fn so it gets recreated by sieve
    config.unknown_filename = "";

    prime_gap_test(config);

    cout << "gap_test_gpu ended" << endl;
}


class SieveResult {
    public:
        SieveResult(const struct Config config): config(config) {};

        enum State { CONFIGURED, SIEVED, DONE };
        State state = CONFIGURED;

        struct Config config;

        uint64_t last_m = (uint64_t) -1;
        size_t current_index = 0;
        std::unique_ptr<SieveOutput> result;
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
        int32_t unknown_start[2] = {0, 0};
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


class GPUBatch {
    public:
        enum State : uint8_t { EMPTY, READY, DONE };
        std::atomic<State> state = EMPTY;

        // Used in debugging Batch Timing.
        time_point<high_resolution_clock> fill_start;
        time_point<high_resolution_clock> fill_end;
        time_point<high_resolution_clock> gpu_start;
        time_point<high_resolution_clock> gpu_end;
        time_point<high_resolution_clock> results_start;
        time_point<high_resolution_clock> results_end;

        // current index;
        size_t i;

        // number to check if prime
        vector<mpz_t*> z;
        // XXX: This is an ugly hack because you can't create mpz_t vector easily
        mpz_t *z_array;

        // If z[i] should be tested
        vector<char>  active;
        // Result from GPU
        vector<int>  result;

        // Intervals for each element
        vector<DataM*> intervals; // TODO was data_m

        // index into the interval.unknowns
        vector<int> unknown_i;

        // For signaling
        std::mutex m;
        std::condition_variable cv;

        explicit GPUBatch(size_t n) {
            elements = n;

            z_array = (mpz_t *) malloc(n * sizeof(mpz_t));
            for (size_t i = 0; i < n; i++) {
                mpz_init(z_array[i]);
                z.push_back(&z_array[i]);
            }

            active.resize(n, 0);
            result.resize(n, -1);

            intervals.resize(n, nullptr);
            unknown_i.resize(n, -1);
        }

        ~GPUBatch() {
            cout << "~GPUBatch" << endl;
            for (size_t i = 0; i < elements; i++) {
                mpz_clear(z_array[i]);
            }
            free(z_array);
        }

        GPUBatch(const GPUBatch&) = delete;
        GPUBatch& operator=(const GPUBatch&) = delete;

    private:
        size_t elements;
};


/** Shared state between threads */
std::atomic<bool> is_running;
std::atomic<bool> queue_new_work;

// Don't read from sieveds without holding sieve_mtx
std::mutex sieve_mtx;
vector<std::unique_ptr<SieveResult>> sieveds;

// Don't mutate an interval (especially status) without holding interval_mtx
std::mutex interval_mtx;
std::condition_variable overflow_cv;
// deque (double ended queue) avoids a degenerate case of large gap getting stuck
// if this can't keep up. Try to avoid falling behind, but this is an extra safety.
deque<DataM*> overflowed;

void run_gpu_thread(int verbose, int runner_num, GPUBatch& batch) {
    try {
        // TODO use runner_num in threadname.
        pthread_setname_np(pthread_self(), "RUN_GPU_THREAD");

        // TODO test changing cudaDeviceScheduleBlockingSync to cudaDeviceScheduleYield or cudaDeviceScheduleSpin
        typedef mr_params_t<THREADS_PER_INSTANCE, BITS, WINDOW_BITS> params;
        test_runner_t<params> runner(GPU_BATCH_SIZE, ROUNDS);

        size_t processed_batches = 0;
        std::unique_lock lock(batch.m, std::defer_lock);
        while (is_running) {
            lock.lock();
            batch.cv.wait(lock, [&] { return batch.state == GPUBatch::State::READY ||!is_running; });

            if (!is_running) break;

            assert(batch.state == GPUBatch::State::READY);
            // Active items are all at the front of the batch.
            auto mid = batch.active.begin();
            std::advance(mid, batch.i);
            assert(std::count(batch.active.begin(), mid, 1) == batch.i);
            assert(std::count(mid,   batch.active.end(), 1) == 0);
            batch.gpu_start = high_resolution_clock::now();
            lock.unlock();

            // Run batch on GPU and wait for results to be set
            if (1) {
                //printf("GPU(%d): Starting batch %lu\n", runner_num, processed_batches);
                runner.run_test(batch.i, batch.z, batch.result);
                //printf("GPU(%d): Finished batch %lu\n", runner_num, processed_batches);
            } else {
                // Return true for 1/10 results (helps not overflow sieve)
                for (size_t gpu_i = 0; gpu_i < GPU_BATCH_SIZE; gpu_i++) {
                    if (batch.active[gpu_i]) {
                        batch.result[gpu_i] = (std::rand() % 10) == 1;
                    }
                }
            }

            lock.lock();
            batch.gpu_end = high_resolution_clock::now();
            batch.state = GPUBatch::State::DONE;
            processed_batches += 1;

            // let CPU thread one.
            lock.unlock();
            batch.cv.notify_one();
        }

        if (verbose >= 1) {
            printf("GPU(%d): Processed %'ld batches\n", runner_num, processed_batches);
        }
    } catch (const std::exception &e) {
        cout << "ERROR in run_gpu_thread" << endl;
        cout << e.what() << endl;
        is_running = false;
    }
}

void run_sieve_thread(void) {
    try {
        pthread_setname_np(pthread_self(), "SIEVE_THREAD");

        std::unique_lock<std::mutex> lock(sieve_mtx, std::defer_lock);
        while (queue_new_work) {
            SieveResult *to_sieve = nullptr;

            lock.lock();
            // Check if any config to start processing
            for (auto& sieved : sieveds) {
                if (sieved->state == SieveResult::CONFIGURED) {
                    to_sieve = sieved.get();
                    break;
                }
            }
            // Safe to unlock because CONFIGURED are never changed by anyone else.
            lock.unlock();

            if (to_sieve == nullptr) {
                usleep(1'000'000); // 1,000ms
                continue;
            }

            auto M_end = to_sieve->config.mstart + to_sieve->config.minc;
            auto s_start_t = high_resolution_clock::now();
            to_sieve->config.threads = COMBINED_SIEVE_THREADS;
            to_sieve->config.verbose -= 3;
            auto result = prime_gap_parallel(to_sieve->config);
            to_sieve->config.verbose += 3;
            auto s_stop_t = high_resolution_clock::now();
            printf("\tCombined Sieve (%ldM to %ldM) took %.1f seconds\n",
                   to_sieve->config.mstart / 1'000'000, M_end / 1'000'000,
                   duration<double>(s_stop_t - s_start_t).count());

            lock.lock();
            if (queue_new_work) {
                to_sieve->state = SieveResult::SIEVED;
                to_sieve->last_m = result->m_start;
                to_sieve->result.swap(result);
            } else {
                to_sieve->state = SieveResult::DONE;
            }
            lock.unlock();
        }

        // Delete any CONFIGURED but not started
        lock.lock();
        for (auto& sieved : sieveds) {
            if (sieved->state == SieveResult::CONFIGURED) {
                sieved->state = SieveResult::DONE;
            }
        }
        lock.unlock();
    } catch (const std::exception &e) {
        cout << "ERROR in run_sieve_thread" << endl;
        cout << e.what() << endl;
        is_running = false;
    }
}

void run_overflow_thread(const mpz_t &K_in) {
    try {
        pthread_setname_np(pthread_self(), "CPU_OVERFLOW_SIEVE_THREAD");
        mpz_t K;
        mpz_init_set(K, K_in);

        std::vector<std::pair<uint32_t, uint32_t>> p_and_r;
        primesieve::iterator iter;
        uint64_t prime = iter.next_prime();
        assert (prime == 2);  // we skip 2 which is the oddest prime.
        for (prime = iter.next_prime(); prime < CPU_SIEVE_LIMIT; prime = iter.next_prime()) {
            const uint32_t base_r = mpz_fdiv_ui(K, prime);
            p_and_r.emplace_back((uint32_t) prime, base_r);
        }

        std::unique_lock<std::mutex> lock(interval_mtx, std::defer_lock);
        size_t tested = 0;

        while (is_running) {
            // Important so that overflow_cv / unlock waits correctly
            if (!lock.owns_lock())
                lock.lock();
            // Lock IS NOT held while waiting.
            overflow_cv.wait(lock, []{ return overflowed.size() || !is_running; });

            while (is_running && overflowed.size()) {
                // 16K batch might dump 160-400 per batch.
                if (tested % 10'000 == 0 && overflowed.size() > 1000) {
                    printf("\tCPU Sieve Queue: %lu open, %lu processed\n",
                            overflowed.size(), tested);
                }

                DataM& interval = *overflowed.front(); overflowed.pop_front();
                assert(interval.state == DataM::State::OVERFLOWED);
                interval.state = DataM::State::SIEVING;

                int32_t sieve_start = 0;
                int32_t sl = interval.sieve_length;

                if (interval.side == DataM::Side::PREV_P) {
                    assert( interval.prev_p == 0 );
                    interval.unknown_start[0] -= sl;
                    sieve_start = interval.unknown_start[0];
                } else {
                    assert( interval.next_p == 0 && interval.n_tests > 0);
                    interval.unknown_start[1] += sl;
                    sieve_start = interval.unknown_start[1];
                }

                lock.unlock();  // Allow main thread to add more things while we process
                vector<int32_t> unknowns;
                sieve_interval_cpu(interval.m, K, p_and_r, sieve_start, sl, unknowns);
                lock.lock();

                assert(5 <= unknowns.size() && unknowns.size() <= (((size_t) sl) / 4));
                assert(unknowns.front() >= 0 && unknowns.back() <= sl);

                if (interval.side == DataM::Side::PREV_P) {
                    interval.unknowns[0].clear();
                    interval.unknowns[0].reserve(unknowns.size());
                    interval.p_index = 0;
                    assert(sieve_start + sl <= 0);
                    for (auto iter = unknowns.rbegin(); iter != unknowns.rend(); iter++) {
                        auto u = *iter;
                        assert(u >= 0);
                        assert( sieve_start + u < 0 );
                        interval.unknowns[0].push_back(sieve_start + u);
                    }
                } else {
                    interval.unknowns[1].clear();
                    interval.unknowns[1].reserve(unknowns.size());
                    interval.n_index = 0;
                    for (auto u : unknowns) {
                        interval.unknowns[1].push_back(sieve_start + u);
                    }
                    assert( sieve_start <= interval.unknowns[1].front() );
                    assert( interval.unknowns[1].front() < interval.unknowns[1].back() );
                    assert( interval.unknowns[1].back() < sieve_start + sl );
                }

                interval.state = DataM::State::READY;
                tested += 1;
            }
        }

        cout << "\tOverflowed " << tested << " intervals" << endl;
        mpz_clear(K);
    } catch (const std::exception &e) {
        cout << "ERROR in run_overflow_thread" << endl;
        cout << e.what() << endl;
        is_running = false;
    }
}


size_t add_to_processing(
        const mpz_t &K,
        const uint32_t D,
        std::vector<std::unique_ptr<DataM>> &processing) {
    // Add the n smallest m range from sieveds to processing.

    SieveResult *sieve = nullptr;
    {
        // Only need the lock while reading state
        std::lock_guard lock(sieve_mtx);
        for (auto& test : sieveds) {
            if (test->state == SieveResult::SIEVED) {
                if (sieve == nullptr || test->result->m_start < sieve->result->m_start) {
                    sieve = test.get();
                }
            }
        }
    }

    if (!sieve) {
        // Nothing to add!
        return 0;
    }

    const SieveOutput& result = *(sieve->result);

    int32_t sl = result.sieve_length;
    assert( sl == (int32_t) sieve->config.sieve_length );

    const auto& coprime_X = result.coprime_X;

    size_t added = 0;

    //printf("Adding %lu from m_start: %lu @ index: %lu/%lu\n",
    //    add, result.m_start, sieve->current_index, result.m_inc.size());

    for (auto& interval : processing) {
        if (interval.get() != nullptr) {
            continue;
        }

        // Updated last_m (current_m during this loop)
        // Update current_index (next_index during the loop)

        const auto [m_add, found] = result.m_inc[sieve->current_index];
        const auto m_unknown_deltas = result.unknowns[sieve->current_index];
        assert( (size_t) found == m_unknown_deltas.size() );
        sieve->last_m += m_add;
        sieve->current_index++;

        uint64_t m = sieve->last_m;
        assert( gcd(m, D) == 1 );

        added += 1;
        interval.reset(new DataM(m, sl));
        interval->unknowns[1].reserve(found);

        int32_t offset = 0;
        for (size_t j = 0; j < (unsigned) found; j++) {
            auto delta = m_unknown_deltas[j];
            offset += delta;
            interval->unknowns[1].push_back(coprime_X[offset]);
        }

        assert( (size_t) offset < coprime_X.size() );
        assert( m_unknown_deltas.size() > 2 );
        assert( interval->unknowns[1].size() == (size_t) found );
        assert( interval->unknowns[1].back() <= sl );

        // No longer need a copy in SieveResult.
        sieve->result->unknowns[sieve->current_index - 1].clear();
        // Try to reclaim memory.
        sieve->result->unknowns[sieve->current_index - 1].shrink_to_fit();

        mpz_init(interval->center);
        mpz_mul_ui(interval->center, K, interval->m);

        if (sieve->current_index == result.m_inc.size()) {
            std::lock_guard lock(sieve_mtx);
            sieve->state = SieveResult::DONE;
            break;
        }
    }

    return added;
}


bool process_finished_batch(GPUBatch& batch) {
    bool found_any_primes = false;
    std::lock_guard lock(interval_mtx);
    for (size_t i = 0; i < GPU_BATCH_SIZE; i++) {
        if (!batch.active[i]) {
            continue;
        }
        // Verify GPU really did write the result
        assert (batch.result[i] == 0 || batch.result[i] == 1);

        DataM &interval = *batch.intervals[i];
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

            int offset_i = batch.unknown_i[i];
            int found_offset = interval.unknowns[side_i][offset_i];

            // Shouldn't have found a previous prime.
            if (interval.side == DataM::Side::PREV_P) {
                assert( found_offset < 0 );
                assert( !interval.p_found );
                assert( interval.prev_p == 0 );
                interval.p_found = true;
                interval.prev_p = -found_offset;
            } else {
                assert( !interval.n_found );
                assert(interval.next_p == 0);
                interval.n_found = true;
                interval.next_p = found_offset;
            }
            interval.state = DataM::State::PRIMED;
            found_any_primes = true;
        }
    }

    return found_any_primes;
}


void fill_batch(
        GPUBatch& batch,
        std::vector<std::unique_ptr<DataM>> &processing,
        StatsCounters& stats) {
    assert( batch.state == GPUBatch::State::EMPTY);

    // Not a full count but (only covers the portion till GPU_BATCH_SIZE is full).
    size_t num_overflowed = 0;
    size_t num_tombstoned = 0;
    size_t num_primed = 0;
    size_t num_running = 0;

    // Grap some entries from each item in M
    batch.i = 0;
    // Turn off all entries in batch
    std::fill_n(batch.active.begin(), GPU_BATCH_SIZE, false);
    // Mark all results as invalid
    std::fill_n(batch.result.begin(), GPU_BATCH_SIZE, -1);

    bool any_pushed_to_overflow = false;
    {
        std::lock_guard lock(interval_mtx);
        for (const auto& row : processing) {
            if (!row) {
                num_tombstoned += 1;
                continue;
            }

            DataM& interval = *row;
            if (interval.state != DataM::State::READY) {
                num_overflowed += (
                    interval.state == DataM::State::OVERFLOWED ||
                    interval.state == DataM::State::SIEVING
                );
                num_primed += interval.state == DataM::State::PRIMED;
                num_running += interval.state == DataM::State::RUNNING;
                continue;  // Not ready to be part of a batch.
            }

            int32_t last_offset = 0;

            auto side = interval.side;
            int32_t index = 0;
            int32_t offset = 0;
            if (side == DataM::Side::NEXT_P) {
                assert(!interval.n_found);
                if (interval.n_index == interval.unknowns[1].size()) {
                    interval.state = DataM::State::OVERFLOWED;
                    stats.s_gap_out_of_sieve_next += 1;
                    overflowed.push_back(row.get());
                    any_pushed_to_overflow = true;
                    continue;
                }
                index = interval.n_index++;
                offset = interval.unknowns[1][index];
            } else {
                assert(interval.n_found);
                assert(!interval.p_found);
                if (interval.p_index == interval.unknowns[0].size()) {
                    interval.state = DataM::State::OVERFLOWED;
                    stats.s_gap_out_of_sieve_prev += 1;
                    overflowed.push_back(row.get());
                    any_pushed_to_overflow = true;
                    continue;
                }
                index = interval.p_index++;
                offset = interval.unknowns[0][index];
            }

            interval.state = DataM::State::RUNNING;
            int gpu_i = batch.i;  // [GPU] batch index

            if (offset < 0) {
                assert( side == DataM::Side::PREV_P );
                assert( offset < last_offset );  // Moving away from center.
                mpz_sub_ui(*batch.z[gpu_i], interval.center, -offset);
            } else {
                assert( side == DataM::Side::NEXT_P );
                assert( offset > last_offset );  // Moving away from center.
                mpz_add_ui(*batch.z[gpu_i], interval.center, offset);
            }
            batch.intervals[gpu_i] = row.get();
            batch.unknown_i[gpu_i] = index;
            batch.active[gpu_i] = true;

            batch.i++;
            last_offset = offset;
            if (batch.i == GPU_BATCH_SIZE) break;
        }
    }
    if (any_pushed_to_overflow) {
        // CPU sieving thread will start if unlocked and notified
        overflow_cv.notify_one();
    }

    assert( batch.i <= GPU_BATCH_SIZE);

    // Batches should be full unless lots of overflowed results.
    if (batch.i > 0 && batch.i < GPU_BATCH_SIZE && num_tombstoned < GPU_BATCH_SIZE) {
        printf("Partial load @ %lu -> %lu/%lu | size: %lu | running: %lu, primed: %lu, overflowed: %lu, tombstoned: %lu\n",
            batch.i > 0 ? batch.intervals[0]->m : 0,
            batch.i, GPU_BATCH_SIZE,
            processing.size(),
            num_running, num_primed, num_overflowed, num_tombstoned);
    }
}


void create_gpu_batches(const struct Config og_config) {

    // gap / 2 up to 60 merit
    uint64_t distance_counts[10000] = {};

    try {
        cout << endl;

        // K is initialized in prob_prime_and_stats
        mpz_t K;
        double K_log = prob_prime_and_stats(og_config, K);
        const uint64_t P = og_config.p;
        const uint64_t D = og_config.d;

        const float min_merit = og_config.min_merit;

        // See THEORY.md! Added const is small preference for doing less prev_p.
        const float MIN_MERIT_TO_CONTINUE = 2.6 + std::log2(min_merit * std::log(2) + 1);

        // Print Header info
        if (og_config.verbose >= 1) {
            setlocale(LC_NUMERIC, "");
            // ----- Merit / Sieve stats
            float m_log = log(og_config.minc);
                printf("Min Gap ~= %'d (for merit > %.1f)\n",
                    (int) (min_merit * (K_log + m_log)), min_merit);
                printf("Min Gap to continue ~= %'d (for merit = %.1f)\n",
                    (int) (MIN_MERIT_TO_CONTINUE * (K_log + m_log)),
                    MIN_MERIT_TO_CONTINUE);

            const uint64_t M_start = og_config.mstart;
            const uint64_t M_inc = og_config.minc;
            uint64_t valid_ms = count_num_m(M_start, M_inc, D);
            assert(valid_ms > 0 && valid_ms <= M_inc);
            printf("\nTesting ranges of %'ld ~ %'ld m per range.\n\n", M_inc, valid_ms);
            setlocale(LC_NUMERIC, "C");
        }

        // Create QUEUE_SIZE empty unique_ptrs.
        std::vector<std::unique_ptr<DataM>> processing(QUEUE_SIZE);

        while (is_running) {
            usleep(50'000); // 50ms
            // Wait for first data before intitalizing StatsCounter to get more stable numbers.
            if (add_to_processing(K, D, processing)) {
                break;
            }
        }

        // Used for various stats
        StatsCounters stats(high_resolution_clock::now());

        /* Note: Uses a double batched system
         * C++ Thread is preparing batch_a (even more m), while GPU runs batch_b */
        std::array<GPUBatch, GPU_BATCHES> gpu_batches = {
            GPUBatch(GPU_BATCH_SIZE),
            GPUBatch(GPU_BATCH_SIZE),
            //GPUBatch(GPU_BATCH_SIZE),
        };

        std::thread gpu_threads[GPU_BATCHES];
        for(size_t i = 0; i < GPU_BATCHES; i++) {
            gpu_threads[i] = std::thread(run_gpu_thread, og_config.verbose, i, std::ref(gpu_batches[i]));
        }

        // Silly but that's what life is.
        std::queue<int> open_gpu;

        //auto s_start_t = high_resolution_clock::now();
        //auto s_batch0_t = s_start_t;
        //auto s_batch1_t = s_start_t;

        // Main loop
        while (is_running) {
            /**
             * Try to fill all batches
             * queue on gpu all ready batches
             * wait for result from the 1st batch sent
             */

            // Add new DataM for any tombstoned items.
            add_to_processing(K, D, processing);

            for (size_t i = 0; i < GPU_BATCHES; i++) {
                GPUBatch& batch = gpu_batches[i];

                if (batch.state != GPUBatch::State::EMPTY)
                    continue;

                batch.fill_start = high_resolution_clock::now();
                fill_batch(batch, processing, stats);
                batch.fill_end = high_resolution_clock::now();

                if (batch.i > 0) {
                    batch.state = GPUBatch::State::READY;

                    open_gpu.push(i);

                    // Start the batch
                    batch.cv.notify_one();
                }
            }

            // Wait for the next batch to be done.
            {
                if (open_gpu.empty()) {
                    // Why would this happen, empty and nothing to queue?
                    usleep(50'000); // 50ms
                } else {
                    int i = open_gpu.front();
                    open_gpu.pop();

                    GPUBatch& batch = gpu_batches[i];
                    // Wait for the batch to be Done (unless it's already done)
                    if (batch.state != GPUBatch::State::DONE) {
                        std::unique_lock<std::mutex> lock(gpu_batches[i].m);
                        batch.cv.wait(lock, [&] { return batch.state == GPUBatch::State::DONE ||!is_running; });
                    }

                    if (batch.state == GPUBatch::State::DONE) {
                        batch.results_start = high_resolution_clock::now();
                        // Read results, mark any found primes, and possible finalize m-interval
                        process_finished_batch(batch);

                        batch.results_end = high_resolution_clock::now();
                        //if (rand() % (1 * 1024) == 0) {
                        if (0) {
                            // TODO check if gpu times are the same.
                            // If so that means that they are running side by side which maybe isn't what we want.
                            printf("CPU: batch timing fill: %.4f, to gpu: %.4f, "
                                    "gpu: %.4f, to cpu: %.4f, process: %.4f\n",
                                   duration<double>(batch.fill_end - batch.fill_start).count(),
                                   duration<double>(batch.gpu_start - batch.fill_end).count(),
                                   duration<double>(batch.gpu_end - batch.gpu_start).count(),
                                   duration<double>(batch.results_start - batch.gpu_end).count(),
                                   duration<double>(batch.results_end - batch.results_start).count());
                        }


                        // Result batch to EMPTY
                        batch.state = GPUBatch::State::EMPTY;
                    }
                }
            }

            {
                bool any_pushed_to_overflow = false;
                // Process & remove finished intervals.
                // This is handled seperately for easier clean & single call to overflow.
                std::unique_lock<std::mutex> lock(interval_mtx);
                for (auto &row : processing) {
                    if (!row) {
                        continue;
                    }
                    DataM &interval = *row;

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
                            stats.s_gap_out_of_sieve_prev += 1;
                            interval.side = DataM::Side::PREV_P;
                            interval.state = DataM::State::OVERFLOWED;
                            overflowed.push_back(row.get());
                            //cout << "PREV_P for m=" << interval.m << endl;
                            any_pushed_to_overflow = true;
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

                    stats.process_results(og_config, interval.m, /* is_last */ false,
                        interval.unknowns[0].size(), interval.unknowns[1].size(),
                        prev_p, next_p,
                        interval.p_tests, interval.n_tests, merit);

                    distance_counts[next_p >> 1] += 1;

                    if (0 && (stats.s_tests % 1'000'000'000) == 100'000'000) {
                        size_t i = 0;  // max set index
                        for (size_t j = 0; j < 10000; j++) {
                            if (distance_counts[j] > 0) {
                                i = j;
                            }
                        }

                        cout << "COUNTS@" << stats.s_tests << ": " << distance_counts[0];
                        for (size_t j = 1; j <= i; j++) {
                            cout << "," << distance_counts[j];
                        }
                        cout << endl;
                    }

                    //mpz_clear(interval.center);
                    // Delete finished interval by reseting index.
                    row.reset();
                }
                if (any_pushed_to_overflow) {
                    // CPU sieving thread will start if unlocked and notified
                    lock.unlock();
                    overflow_cv.notify_one();
                }
            }

        }

        // ----- cleanup
        {
            mpz_clear(K);
            std::lock_guard lock(interval_mtx);
            for (auto& interval: processing) {
                if (interval) {
                    mpz_clear(interval->center);
                    interval.reset();
                }
            }
        }

        // Send notifies
        for (auto& gpu_batch : gpu_batches) {
            gpu_batch.cv.notify_all();
        }
        for (auto & gpu_thread : gpu_threads) {
            gpu_thread.join();
        }
    } catch (const std::exception &e) {
        cout << "ERROR in create_gpu_batches" << endl;
        cout << e.what() << endl;
        is_running = false;
    }
}

void coordinator_thread(struct Config global_config) {
    try {
        pthread_setname_np(pthread_self(), "MAIN_LOAD_AND_BATCH_THREAD");

        size_t total_ranges = 0;

        std::unique_lock<std::mutex> lock(sieve_mtx, std::defer_lock);
        while (is_running && (queue_new_work || sieveds.size())) {
            usleep(1'000'000); // 1,000ms

            std::unique_ptr<SieveResult> to_test = nullptr;
            lock.lock();
            if (queue_new_work) {  // Add new configs to sieveds queue.
                int needed = 2;
                int finished = 0;
                for (auto& sieved : sieveds) {
                    if (sieved->state == SieveResult::CONFIGURED) {
                        needed -= 1;
                    }
                    if (sieved->state == SieveResult::SIEVED) {
                        needed -= 1;
                        finished += 1;
                    }
                }
                if (finished == 0 && total_ranges > 2) {
                    printf("Currently no finished ranges!!!\n");
                }

                while (needed > 0) {
                    // Add new config to the queue
                    if (global_config.verbose >= 3) {
                        printf("\tQueued(%lu) %'lu to %'lu\n",
                                total_ranges,
                                global_config.mstart,
                                global_config.mstart + global_config.minc);
                    }

                    sieveds.emplace_back(std::make_unique<SieveResult>(global_config));
                    needed -= 1;
                    total_ranges += 1;

                    global_config.mstart += global_config.minc;
                }
            }

            // Remove any finished results.
            sieveds.erase(std::remove_if(std::begin(sieveds), std::end(sieveds),
                    [](const std::unique_ptr<SieveResult>& sieve) {
                        return sieve->state == SieveResult::DONE;
                    }), sieveds.end());
            if (sieveds.size() > (2+2)) {
                printf("%lu open ranges! `overflowed` might be falling behind!\n", sieveds.size());
            }


            lock.unlock();
        }
        cout << "\n\nCoordinator Done.\n" << endl;
    } catch (const std::exception &e) {
        cout << "ERROR in coordinator_thread" << endl;
        cout << e.what() << endl;
        is_running = false;
    }
}



void signal_callback_handler(int) {
    if (queue_new_work) {
       cout << endl;
       cout << "Caught CTRL+C stopping, winding down work." << endl;
       cout << endl;
       queue_new_work = false;
    } else {
       cout << endl;
       cout << "Caught 2nd CTRL+C, exit(2) now." << endl;
       cout << endl;
       is_running = false;
       exit(2);
    }
}

void prime_gap_test(struct Config config) {
    // Setup test runner
    printf("\n");
    printf("BITS=%d\tWINDOW_BITS=%d\n", BITS, WINDOW_BITS);
    printf("PRP/BATCH=%ld\n", GPU_BATCH_SIZE);
    printf("THREADS/PRP=%d\n", THREADS_PER_INSTANCE);

    assert( GPU_BATCH_SIZE == 1024 || GPU_BATCH_SIZE == 2048 || GPU_BATCH_SIZE == 4096 ||
            GPU_BATCH_SIZE == 8192 || GPU_BATCH_SIZE ==16384 || GPU_BATCH_SIZE ==32768 );

    mpz_t K;
    {
        init_K(config, K);
        // +4 is just is personal safety blanket buffer.
        size_t N_bits = mpz_sizeinbase(K, 2) + log2(config.mstart + 100ul * config.minc) + 4;

        // P# roughly 349, 709, 1063, 1447
        for (size_t bits : {512, 1024, 1536, 2048, 3036, 4096}) {
            if (N_bits <= bits) {
                if (bits < BITS) {
                    printf("\nFASTER WITH `make gap_test_gpu BITS=%ld` (may require `make clean`)\n\n", bits);
                    exit(1);
                }
                break;
            }
        }
        assert( (N_bits + 1) < BITS ); // See last debug line.
        assert( BITS <= (1 << (2 * WINDOW_BITS)) );
    }


    is_running     = true;
    queue_new_work = true;

    // Setup CTRL+C catcher
    signal(SIGINT, signal_callback_handler);

    /**
     * coordinator_thread: creates configs and adds them to queue
     * sieve_thread: reads configs from ^ queue and runs combined_sieve_small
     *
     * batch_thread: reads SieveResult's and creates GPU batches
     * gpu_thread: runs the GPU batches.
     *
     * overflow_sieve_thread: computes extra sieves for prev_p and overflowed next_p
     *
     */
    std::thread main_thread(coordinator_thread, config);
    std::thread sieve_thread(run_sieve_thread);

    std::thread batch_thread(create_gpu_batches, config);
    std::thread overflow_sieve_thread(run_overflow_thread, std::ref(K));

    main_thread.join();
    cout << "Joining threads" << endl;

    // Tell other threads to quit
    {
        is_running = false;
        sieve_thread.join();
        batch_thread.join();

        overflow_cv.notify_all();  // wake up all overflow thread
        overflow_sieve_thread.join();
    }

    mpz_clear(K);
}
