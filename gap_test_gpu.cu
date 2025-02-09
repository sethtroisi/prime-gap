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

/**
 * Used for doing mpz_prevprime and mpz_nextprime on CPU
 * Given 40x slower, need a few to keep up with GPU
 */
const int CPU_THREADS = 6;

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

        // index into 'processing' (DataM)
        vector<int64_t> data_i;

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

        // If this entry needs to be handled on CPU with gmp prime finding.
        // One (and only one) of prev_p or next_p should be -1
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

std::mutex interval_mtx;
std::condition_variable overflow_cv;
vector<std::shared_ptr<DataM>> overflowed;

void run_gpu_thread(const struct Config config) {
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
                is_close_to_end |= batch.end_of_file;
                last_finished = high_resolution_clock::now();
            }
        }
        if (no_batch) {
            // Waiting doesn't count till 1st batch is ready
            auto now = high_resolution_clock::now();
            auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(
                    now - last_finished).count();
            if (config.verbose >= 0 && processed_batches > 0 && delta > 50) {
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

    if (config.verbose >= 1) {
        printf("Processed %'ld batches\n", processed_batches);
    }
}

void run_overflow_thread(const struct Config config) {
    pthread_setname_np(pthread_self(), "CPU_OVERFLOW_THREAD");
    mpz_t center, prime_test;
    mpz_init(center);
    mpz_init(prime_test);
    size_t tested = 0;

    std::unique_lock<std::mutex> lock(interval_mtx, std::defer_lock);

    while (is_running) {
        // Important so that overflow_cv / unlock waits correctly
        if (!lock.owns_lock())
            lock.lock();
        // Lock IS NOT held while waiting.
        overflow_cv.wait(lock, []{ return overflowed.size() || !is_running; });

        while (overflowed.size()) {
            DataM& interval = *overflowed.back(); overflowed.pop_back();
            assert (interval.overflow && interval.state == DataM::State::RUNNING);

            // NOTE: Overhead to doing this while GPU waits seems small (<1% of candidates)
            // But is actually A LOT because 40x slower. Becomes ~20-40% overhead quickly.
            assert(interval.overflow == 1);
            assert(interval.next_p == -1 || interval.prev_p == -1);
            mpz_set(center, interval.center);

            if (interval.next_p == -1) {
                assert(interval.n_tests > 0);

                lock.unlock();  // Allow main thread to add more things while we process
                mpz_add_ui(prime_test, center, config.sieve_length);
                mpz_nextprime(prime_test, prime_test);
                mpz_sub(prime_test, prime_test, center);
                lock.lock();

                // Atomic updated because we hold the lock
                interval.n_found = true;
                interval.next_p = mpz_get_ui(prime_test);
                //cout << "gap_out_of_sieve_next m=" << interval.m << " -> " << interval.next_p << endl;
            } else {
                assert(interval.prev_p == -1);
                assert(interval.p_tests == 0);

                lock.unlock();  // Allow main thread to add more things while we process
                mpz_prevprime(prime_test, center);
                mpz_sub(prime_test, center, prime_test);
                lock.lock();

                // Atomic updated because we hold the lock
                interval.p_found = true;
                interval.prev_p = mpz_get_ui(prime_test);
                //cout << "gap_out_of_sieve_prev m=" << interval.m << " -> " << interval.prev_p << endl;
            }

            interval.overflow = 0;
            interval.state = DataM::State::OVERFLOW_DONE;
            tested += 1;
        }
    }

    cout << "\tOverflowed " << tested << " intervals" << endl;
    mpz_clear(center);
    mpz_clear(prime_test);
}

void load_batch_thread(const struct Config config, const size_t QUEUE_SIZE) {
    pthread_setname_np(pthread_self(), "MAIN_LOAD_AND_BATCH_THREAD");

    mpz_t K;
    double K_log;
    std::ifstream unknown_file;

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
        K_log = prob_prime_and_stats(config, K);
        {
            float m_log = log(M_start);
            if (config.verbose >= 1) {
                printf("Min Gap ~= %d (for merit > %.1f)\n",
                    (int) (min_merit * (K_log + m_log)), min_merit);
                printf("Min Gap to continue ~= %d (for merit = %.1f)\n",
                    (int) (MIN_MERIT_TO_CONTINUE * (K_log + m_log)),
                    MIN_MERIT_TO_CONTINUE);
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

    // Used for various stats
    StatsCounters stats(high_resolution_clock::now());
    size_t pushed_batches = 0;

    // Main loop
    uint64_t mi = 0;
    while (mi < M_inc || !processing.empty()) {
        usleep(100); // 0.1ms
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
                //cout << "00 " << mi << " " << processing.size() << endl;

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
                            if (interval.state != DataM::State::READY || interval.overflow) {
                                currently_overflowed += interval.overflow;
                                // Already part of some other batch
                                continue;
                            }

                            for (size_t j = 0; j < SEQUENTIAL_IN_BATCH; j++) {
                                // interval shouldn't have any ends (as only 1 side is tested in this loop);
                                assert(!interval.n_found);
                                if (interval.n_tests == interval.unknowns[1].size()) {
                                    /* Don't push to overflow queue while GPU still has sieve to check.
                                     * leads to race condition with overflow thread & uglier code.
                                     */
                                    if (j == 0) {
                                        interval.next_p = -1;
                                        interval.overflow = 1; // Indicates next side has overflowed sieve
                                    }
                                    break;
                                }

                                interval.state = DataM::State::RUNNING;
                                int gpu_i = batch.i;  // [GPU] batch index
                                auto offset = interval.unknowns[1][interval.n_tests];
                                //printf("%d = %d,%d = %d*K+%d\n", gpu_i, i, j, interval.m, offset);
                                assert(offset > 0);
                                mpz_add_ui(*batch.z[gpu_i], interval.center, offset);
                                batch.data_i[gpu_i] = interval.m;  // [Data] index for GPU Batch
                                batch.unknown_i[gpu_i] = interval.n_tests++;
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

                // Mark batch as ready for GPU processing
                batch.end_of_file = (mi == M_inc);
                // Needs to be last so that GPU_thread doesn't read other parts early.
                batch.state = batch.i > 0 ? GPUBatch::State::READY : GPUBatch::State::EMPTY;
                //cout << "AA " << mi << " " << batch.i << endl;

                pushed_batches += batch.state == GPUBatch::State::READY;
                if (pushed_batches == 1) {
                    // Helpful for getting better tests/sec
                    stats.s_start_t = high_resolution_clock::now();
                }

            } else if (batch.state == GPUBatch::State::RESULT_WRITTEN) {
                //cout << "CC " << mi << " " << batch.i << endl;
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
                        DataM &interval = *processing.at(batch.data_i[i]);
                        // Mark interval as being ready again
                        interval.state = DataM::State::READY;

                        if (batch.result[i]) {
                            // Found prime in last partial batch of unknowns, no longer overflowed
                            assert(interval.overflow == 0 || SEQUENTIAL_IN_BATCH > 1);
                            interval.overflow = 0;

                            // next_prime found (and done)
                            int offset_i = batch.unknown_i[i];
                            assert(interval.n_tests > 0 );
                            int found_offset = interval.unknowns[1][offset_i];
                            if (interval.n_found && interval.next_p > 0 && interval.next_p < found_offset) {
                                // Two primes happens fairly regularly with SEQUENTIAL_IN_BATCH > 1
                                assert(SEQUENTIAL_IN_BATCH > 1);
                                //cerr << "\tFound two next primes for m=" << interval.m;
                                //cerr << " | " << interval.next_p;
                                //cerr << " vs " << found_offset << "(" << offset_i << ")" << endl;
                                continue;
                            }
                            assert(interval.next_p == 0);
                            interval.n_found = true;
                            interval.next_p = found_offset;
                        }
                    }
                }

                // Result batch to EMPTY
                batch.state = GPUBatch::State::EMPTY;
            }

            {
                bool any_pushed_to_overflow = false;
                // Push Out-Of-Sieve gaps to overflow queue and later notify that thread
                {
                    std::unique_lock<std::mutex> lock(interval_mtx);
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
                                overflowed.push_back(pair.second);
                                any_pushed_to_overflow = true;
                            }
                        }
                    }
                }

                if (any_pushed_to_overflow) {
                    // CPU overflow thread will start if unlocked and notified
                    overflow_cv.notify_one();
                    // TODO print warning if overflowed.size() is very large
                }
            }

            {
                // Process & remove finished intervals.
                std::unique_lock<std::mutex> lock(interval_mtx);
                for (auto it = processing.begin(); it != processing.end(); /* increment in loop */) {
                    auto& interval = *it->second;

                    int prev_p = interval.prev_p;
                    int next_p = interval.next_p;

                    // Potentially do Side-Skip if next_p is not very large (and just found)
                    if (interval.n_found && !interval.p_found && prev_p == 0) {
                        float next_merit = next_p / (K_log + log(interval.m));

                        if (next_merit < MIN_MERIT_TO_CONTINUE) {
                            stats.s_skips_after_one_side += 1;
                            goto next_process;
                        }

                        // Mark this for prev_p check (via overflow) and later alert
                        interval.state = DataM::State::READY;
                        interval.overflow = 1;
                        interval.prev_p = -1;
                        continue;
                    }

                    // Waiting on end primes.
                    if (!interval.p_found || !interval.n_found) {
                        ++it;
                        continue;
                    }

                    assert( prev_p > 0 && next_p > 0 );

next_process:
                    assert(interval.state == DataM::State::READY ||
                           interval.state == DataM::State::OVERFLOW_DONE);
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
        assert( (N_bits + 1) < BITS ); // See last debug line.
        assert( BITS <= (1 << (2 * WINDOW_BITS)) );
    }

    is_running = true;

    std::thread main_thread(load_batch_thread, config, BATCHED_M);
    std::thread gpu_thread(run_gpu_thread, config);
    vector<std::thread> overflow_threads;
    for (size_t t = 0; t < CPU_THREADS; t++) {
        overflow_threads.push_back(std::thread(run_overflow_thread, config));
    }

    main_thread.join();
    // Tell other threads to quit
    is_running = false;
    cout << "All done. Goodbye!" << endl;

    gpu_thread.join();
    overflow_cv.notify_all();  // wake up all overflow thread
    for (auto& overflow_thread : overflow_threads) {
        overflow_thread.join();
    }
}

