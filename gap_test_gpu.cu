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
#include <exception>
#include <fstream>
#include <iostream>
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
using std::vector;
using namespace std::chrono;

#ifdef GPU_BITS
#define BITS GPU_BITS
#else
#define BITS 1024
#endif

#define WINDOW_BITS ((BITS <= 1024) ? 5 : 6)

/**
 * GPU_BATCH_SIZE is 2^n | best is between 2K and 8K.
 * SEQUENTIAL_IN_BATCH = {1,2,4}
 *      1 => 0 overhead, 2 => 0.5 extra PRP/m, 4 => 1.5 extra PRP/M
 * BATCHED_M is number of M loaded at the same time
 */

const size_t GPU_BATCH_SIZE = 4 * 1024;
// 8 is best for BITS=1024, 4 is best for BITS=512
const int THREADS_PER_INSTANCE = 4;

// 701#
// 1080Ti - 16K,  8 -> 213K
// 1080Ti - 8K,   8 -> 220K -> 230K double batched
// 1080Ti - 4K,   8 -> 201K
// OLD 701# - A100  -> 327K!

// 347#
// 1080Ti - 4K, 8  -> 1034K
// 1080Ti - 2K, 8  -> 1054K!
// 1080Ti - 1K, 8  -> 780K
// 1080Ti - 8K, 4  -> 1300K
// 1080Ti - 4K, 4  -> 1448K!
// 1080Ti - 2K, 4  -> 1022K


// 1 is best, 2 if very large batch.
const size_t SEQUENTIAL_IN_BATCH = 1;
// Always use 1.
const int ROUNDS = 1;
// 20% extra for overflow is vary reasonable.
const size_t BATCHED_M = 160 / 100 * 2 * GPU_BATCH_SIZE / SEQUENTIAL_IN_BATCH;


// From 701# I believe.
// 1M -> 2000/second
// 200K -> 4000/second
const size_t CPU_SIEVE_LIMIT = 90'000;

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

    if (config.sieve_length < 6 * config.p || config.sieve_length > 22 * config.p) {
        int sl_min = ((config.p * 8 - 1) / 500 + 1) * 500;
        int sl_max = ((config.p * 20 - 1) / 500 + 1) * 500;
        printf("--sieve_length(%d) should be between [%d, %d]\n",
            config.sieve_length, sl_min, sl_max);
        return 1;
    }

    if (100'000 < config.max_prime && config.max_prime > 100'000'000) {
        printf("\tmax_prime(%'ld) should be between 100K and 100M\n", config.max_prime);
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


class GPUBatch {
    public:
        enum State : uint8_t { EMPTY, READY, DONE };
        std::atomic<State> state = EMPTY;

        // current index;
        size_t i;

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

            data_m.resize(n, -1);
            unknown_i.resize(n, -1);
        }

        ~GPUBatch() {
            cout << "~GPUBatch" << endl;
            for (size_t i = 0; i < elements; i++) {
                mpz_clear(z_array[i]);
            }
            free(z_array);
        }

        //GPUBatch(const GPUBatch&) = delete;
        GPUBatch& operator=(const GPUBatch&) = delete;

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


/** Shared state between threads */
std::atomic<bool> is_running;

/**
 * Note: Uses a double batched system
 * C++ Thread is preparing batch_a (even more m)
 * While GPU runs batch_b
 */
std::mutex batches_mtx;
std::array<GPUBatch, 2> gpu_batches = {GPUBatch(GPU_BATCH_SIZE), GPUBatch(GPU_BATCH_SIZE)};

// Don't read from sieveds without holding sieve_mtx
std::mutex sieve_mtx;
vector<std::unique_ptr<SieveResult>> sieveds;

// Don't mutate an interval (especially status) without holding interval_mtx
std::mutex interval_mtx;
std::condition_variable overflow_cv;
vector<DataM*> overflowed;

void run_gpu_thread(int verbose, int runner_num, GPUBatch& batch) {
    try {
        // TODO use runner in threadname.
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
            assert(std::count(batch.active.begin(), batch.active.end(), 1) == batch.i);
            lock.unlock();

            // Run batch on GPU and wait for results to be set
            if (1) {
                //printf("GPU(%d): Starting batch %lu\n", runner_num, processed_batches);
                runner.run_test(batch.z, batch.result);
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
        while (is_running) {
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
            // TODO set with const.
            to_sieve->config.threads = 8;
            to_sieve->config.verbose -= 3;
            auto result = prime_gap_parallel(to_sieve->config);
            to_sieve->config.verbose += 3;
            auto s_stop_t = high_resolution_clock::now();
            printf("\tCombined Sieve (%ldM to %ldM) took %.1f seconds\n",
                   to_sieve->config.mstart / 1'000'000, M_end / 1'000'000,
                   duration<double>(s_stop_t - s_start_t).count());

            lock.lock();
            to_sieve->state = SieveResult::SIEVED;
            to_sieve->last_m = result->m_start;
            to_sieve->result.swap(result);
            lock.unlock();
        }
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
                if (tested % 100'000 == 0 && overflowed.size() > 200) {
                    printf("CPU Sieve Queue: %lu open, %lu processed\n",
                            overflowed.size(), tested);
                }

                DataM& interval = *overflowed.back(); overflowed.pop_back();
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
                    // Consider not clearing and just adding to existing elements.
                    interval.unknowns[0].clear();
                    interval.p_index = 0;
                    assert(sieve_start + sl <= 0);
                    for (auto iter = unknowns.rbegin(); iter != unknowns.rend(); iter++) {
                        auto u = *iter;
                        assert(u >= 0);
                        assert( sieve_start + u < 0 );
                        interval.unknowns[0].push_back(sieve_start + u);
                    }
                } else {
                    // clear existing unknowns and replace
                    interval.unknowns[1].clear();
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
        std::unordered_map<int64_t, std::unique_ptr<DataM>> &processing,
        size_t n) {
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

    // We have a range of sieve results
    // Next item is sieve.current_index
    // have to first process m_inc to know m.
    size_t add = std::min(n, result.m_inc.size() -  sieve->current_index);

    //printf("Adding %lu from m_start: %lu @ index: %lu/%lu\n",
    //    add, result.m_start, sieve->current_index, result.m_inc.size());

    for (size_t i = 0; i < add; i++) {
        // Updated last_m (current_m during this loop)
        // Update current_index (next_index during the loop)

        const auto [m_add, found] = result.m_inc[sieve->current_index];
        const auto m_unknown_deltas = result.unknowns[sieve->current_index];
        assert( (size_t) found == m_unknown_deltas.size() );
        sieve->last_m += m_add;
        sieve->current_index++;

        uint64_t m = sieve->last_m;
        assert( gcd(m, D) == 1 );

        auto test = std::make_unique<DataM>(m, sl);
        test->unknowns[1].reserve(found);

        int32_t offset = 0;
        for (size_t j = 0; j < (unsigned) found; j++) {
            auto delta = m_unknown_deltas[j];
            offset += delta;
            test->unknowns[1].push_back(coprime_X[offset]);
        }

        assert( (size_t) offset < coprime_X.size() );
        assert( m_unknown_deltas.size() > 2 );
        assert( test->unknowns[1].size() == (size_t) found );
        assert( test->unknowns[1].back() <= sl );

        // No longer need a copy in SieveResult.
        sieve->result->unknowns[sieve->current_index - 1].clear();

        mpz_init(test->center);
        mpz_mul_ui(test->center, K, test->m);

        processing.emplace(test->m, std::move(test));
    }

    if (sieve->current_index == result.m_inc.size()) {
        std::lock_guard lock(sieve_mtx);
        sieve->state = SieveResult::DONE;
    }

    return add;
}


bool process_finished_batch(
        GPUBatch& batch,
        std::unordered_map<int64_t, std::unique_ptr<DataM>> &processing) {
    bool found_any_primes = false;
    std::lock_guard lock(interval_mtx);
    for (size_t i = 0; i < GPU_BATCH_SIZE; i++) {
        if (!batch.active[i]) {
            continue;
        }
        // Verify GPU really did write the result
        assert (batch.result[i] == 0 || batch.result[i] == 1);

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

            int offset_i = batch.unknown_i[i];
            int found_offset = interval.unknowns[side_i][offset_i];

            // Two primes happens fairly regularly with SEQUENTIAL_IN_BATCH > 1
            // Ignore the further away prime
            if (interval.side == DataM::Side::PREV_P) {
                assert( found_offset < 0 );
                if (interval.p_found) {
                    assert( SEQUENTIAL_IN_BATCH > 1 );
                    assert( abs(interval.prev_p) < abs(found_offset) );
                    continue;
                }
                assert(interval.prev_p == 0);
                interval.p_found = true;
                interval.prev_p = -found_offset;
            } else {
                if (interval.n_found) {
                    assert( SEQUENTIAL_IN_BATCH > 1 );
                    assert( interval.next_p < found_offset );
                    continue;
                }
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
        std::unordered_map<int64_t, std::unique_ptr<DataM>> &processing,
        StatsCounters& stats) {
    assert( batch.state == GPUBatch::State::EMPTY);

    // Not a full count but (only covers the portion till GPU_BATCH_SIZE is full).
    size_t currently_overflowed = 0;

    // Grap some entries from each item in M
    batch.i = 0;
    // Turn off all entries in batch
    std::fill_n(batch.active.begin(), GPU_BATCH_SIZE, false);
    // Mark all results as invalid
    std::fill_n(batch.result.begin(), GPU_BATCH_SIZE, -1);

    for (size_t i = 0; i < GPU_BATCH_SIZE; i++) {
        // TODO see if this can be moved into init or somewhere else.
        // This prevents the GPU from stalling in partial batches
        mpz_set_ui(*batch.z[i], 7);
    }

    bool any_pushed_to_overflow = false;
    {
        std::lock_guard lock(interval_mtx);
        for (auto& pair : processing) {
            auto& interval = *pair.second;
            currently_overflowed += (
                interval.state == DataM::State::OVERFLOWED ||
                interval.state == DataM::State::SIEVING
            );
            if (interval.state != DataM::State::READY) {
                continue;  // Not ready to be part of a batch.
            }

            int32_t last_offset = 0;

            auto side = interval.side;
            for (size_t j = 0; batch.i < GPU_BATCH_SIZE && j < SEQUENTIAL_IN_BATCH; j++) {
                int32_t index = 0;
                int32_t offset = 0;
                if (side == DataM::Side::NEXT_P) {
                    assert(!interval.n_found);
                    if (interval.n_index == interval.unknowns[1].size()) {
                        /* Don't push to overflow queue while GPU still has sieve to check.
                         * leads to race condition with overflow thread & uglier code.
                         */
                        if (j == 0) {
                            interval.state = DataM::State::OVERFLOWED;
                            stats.s_gap_out_of_sieve_next += 1;
                            overflowed.push_back(pair.second.get());
                            any_pushed_to_overflow = true;
                        }
                        break;
                    }
                    index = interval.n_index++;
                    offset = interval.unknowns[1][index];
                } else {
                    assert(interval.n_found);
                    assert(!interval.p_found);
                    if (interval.p_index == interval.unknowns[0].size()) {
                        if (j == 0) {
                            interval.state = DataM::State::OVERFLOWED;
                            stats.s_gap_out_of_sieve_prev += 1;
                            overflowed.push_back(pair.second.get());
                            any_pushed_to_overflow = true;
                        }
                        break;
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
                batch.data_m[gpu_i] = interval.m;  // [Data] index for GPU Batch
                batch.unknown_i[gpu_i] = index;
                batch.active[gpu_i] = true;

                batch.i++;
                last_offset = offset;
            }
            if (batch.i == GPU_BATCH_SIZE) break;
        }
    }
    if (any_pushed_to_overflow) {
        // CPU sieving thread will start if unlocked and notified
        overflow_cv.notify_one();
    }

    assert( batch.i <= GPU_BATCH_SIZE);

    // Batches should be full unless lots of overflowed results.
    if (batch.i > 0 && batch.i < GPU_BATCH_SIZE) {
        //printf("Partial load @ %lu -> %lu/%lu | %lu\n",
        //    batch.i > 0 ? batch.data_m[0] : 0,
        //    batch.i, GPU_BATCH_SIZE, currently_overflowed);
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

        std::unordered_map<int64_t, std::unique_ptr<DataM>> processing;

        const size_t QUEUE_SIZE = BATCHED_M;

        while (is_running) {
            usleep(50'000); // 50ms
            // Wait for first data before intitalizing StatsCounter to get more stable numbers.
            if (add_to_processing(K, D, processing, QUEUE_SIZE - processing.size())) {
                break;
            }
        }

        // Used for various stats
        StatsCounters stats(high_resolution_clock::now());

        std::thread gpu0(run_gpu_thread, og_config.verbose, 0, std::ref(gpu_batches[0]));
        std::thread gpu1(run_gpu_thread, og_config.verbose, 1, std::ref(gpu_batches[1]));

        // Silly but that's what life is.
        std::queue<int> open_gpu;

        //auto s_start_t = high_resolution_clock::now();
        //auto s_batch0_t = s_start_t;
        //auto s_batch1_t = s_start_t;

        // Main loop
        while (is_running) {
            /**
             * Try to fill all batches
             * send for async all ready batches
             * wait for result from the 1st batch sent
             * repeat?
             */

            // Add new DataM if free space in processing
            // TODO make use of return at some later point
            add_to_processing(K, D, processing, QUEUE_SIZE - processing.size());

            for (int i = 0; i < 2; i++) {
                GPUBatch& batch = gpu_batches[i];

                if (batch.state != GPUBatch::State::EMPTY)
                    continue;

                fill_batch(batch, processing, stats);

                if (batch.i > 0) {
                    batch.state = GPUBatch::State::READY;

                    /*
                    auto last = ((i == 0) ? s_batch0_t : s_batch1_t);
                    auto now = high_resolution_clock::now();
                    ((i == 0) ? s_batch0_t : s_batch1_t) = now;
                    printf("CPU: Batch(%d) ready to run @ %.5f | %.5f\n",
                        i, duration<double>(now - last).count(),
                        duration<double>(now - s_start_t).count());
                    // */

                    open_gpu.push(i);

                    // Start the batch
                    batch.cv.notify_one();
                }
            }

            // Wait for the next batch to be done.
            {
                if (open_gpu.empty()) {
                    usleep(5'000); // 5,000ms
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
                        /*
                        auto last = ((i == 0) ? s_batch0_t : s_batch1_t);
                        auto now = high_resolution_clock::now();
                        ((i == 0) ? s_batch0_t : s_batch1_t) = now;
                        printf("CPU: Processing batch from GPU(%d) took %.5f | %.5f\n",
                               i, duration<double>(now - last).count(),
                               duration<double>(now - s_start_t).count());
                        // */
                        // Read results, mark any found primes, and possible finalize m-interval
                        process_finished_batch(batch, processing);

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
                for (auto it = processing.begin(); it != processing.end(); /* increment in loop */) {
                    auto& interval = *it->second;

                    if (interval.state != DataM::State::PRIMED) {
                        it++;
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
                            overflowed.push_back(it->second.get());
                            //cout << "PREV_P for m=" << interval.m << endl;
                            any_pushed_to_overflow = true;
                            it++;
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

                    mpz_clear(interval.center);
                    it = processing.erase(it);  // Erase this element
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
            for (auto& [m, interval] : processing) {
                mpz_clear(interval->center);
            }
        }

        // Send notifies
        gpu_batches[0].cv.notify_all();
        gpu_batches[1].cv.notify_all();
        gpu0.join();
        gpu1.join();
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
        while (is_running) {
            usleep(1'000'000); // 1,000ms

            std::unique_ptr<SieveResult> to_test = nullptr;
            lock.lock();
            {  // Add new configs to sieveds queue.
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
            lock.unlock();
        }
        cout << "\n\nEND OF MAIN THREAD\n" << endl;
    } catch (const std::exception &e) {
        cout << "ERROR in coordinator_thread" << endl;
        cout << e.what() << endl;
        is_running = false;
    }
}



void signal_callback_handler(int) {
    if (is_running) {
       cout << "Caught CTRL+C stopping, exiting soon. " << endl;
       is_running = false;
    } else {
       cout << "Caught 2nd CTRL+C stopping now." << endl;
       exit(2);
    }
}

void prime_gap_test(struct Config config) {
    // Setup test runner
    printf("BITS=%d\tWINDOW_BITS=%d\n", BITS, WINDOW_BITS);
    printf("PRP/BATCH=%ld\tM/BATCH=%ld\n",
            GPU_BATCH_SIZE, GPU_BATCH_SIZE/SEQUENTIAL_IN_BATCH);
    printf("THREADS/PRP=%d\n", THREADS_PER_INSTANCE);

    assert( GPU_BATCH_SIZE == 1024 || GPU_BATCH_SIZE == 2048 || GPU_BATCH_SIZE == 4096 ||
            GPU_BATCH_SIZE == 8192 || GPU_BATCH_SIZE ==16384 || GPU_BATCH_SIZE ==32768 );
    assert( SEQUENTIAL_IN_BATCH == 1 || SEQUENTIAL_IN_BATCH == 2 || SEQUENTIAL_IN_BATCH == 4 );

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


    is_running = true;

    // Setup CTRL+C catcher
    // This isn't safe because some threads die faster and delete their references while other
    // threads might still make reads.
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
    // Tell other threads to quit
    is_running = false;
    cout << "Setting `is_running = false` and joining threads" << endl;

    sieve_thread.join();
    batch_thread.join();
    overflow_cv.notify_all();  // wake up all overflow thread
    overflow_sieve_thread.join();
    mpz_clear(K);
}
