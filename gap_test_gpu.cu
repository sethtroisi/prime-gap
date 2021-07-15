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
#include <cmath>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
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


void test_interval_gpu(
        const uint64_t m, const mpz_t &K, const size_t SIEVE_LENGTH,
        vector<int32_t> (&unknowns)[2],
        size_t &p_tests, size_t &n_tests,
        int &prev_p, int &next_p) {
    /*
    size_t gap_out_of_sieve_count;
    test_interval_cpu(
            m, K, SIEVE_LENGTH,
            p_tests, gap_out_of_sieve_count, gap_out_of_sieve_count,
            unknowns,
            prev_p, next_p);
    // */
    // /*

    mpz_t center;
    mpz_init(center);
    mpz_mul_ui(center, K, m);

    //printf("%ld | %ld %ld\n", m, unknowns[0].size(), unknowns[1].size());

    /**
     * Originally 8 which has highest throuput but only if we have LOTS of instances
     * this helps reduce the number of parallel instances needed
     */
    const int THREADS_PER_INSTANCE = 32;
    const int ROUNDS = 1;

    typedef mr_params_t<THREADS_PER_INSTANCE, BITS, WINDOW_BITS> params;
    int prev_p_i = run_test<params>(center, -1, unknowns[0], ROUNDS, p_tests);
    int next_p_i = run_test<params>(center, +1, unknowns[1], ROUNDS, n_tests);
    assert(prev_p_i >= 0);
    assert(next_p_i >= 0);
    prev_p = unknowns[0][prev_p_i];
    next_p = unknowns[1][next_p_i];

//    printf("    | %d %d\t%d %d\n", prev_p_i, next_p_i, prev_p, next_p);
    // */
}

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

    for (uint32_t mi = 0; mi < M_inc; mi++) {
        long m = M_start + mi;
        if (gcd(m, D) > 1) {
            continue;
        }

        vector<int32_t> unknowns[2];

        load_and_verify_unknowns(
            compression, M_start + mi, SIEVE_LENGTH, unknown_file, unknowns);

        size_t unknown_l = unknowns[0].size();
        size_t unknown_u = unknowns[1].size();

        size_t p_tests = 0, n_tests = 0;
        int prev_p, next_p;

        test_interval_gpu(
            m, K, SIEVE_LENGTH,
            unknowns,
            p_tests, n_tests,
            prev_p, next_p);
        assert( prev_p > 0 && next_p > 0 );

        float merit = (next_p + prev_p) / (K_log + log(m));

        if (merit > min_merit)  {
            // TODO: write to file or database
            printf("%-5d %.4f  %ld * %ld#/%ld -%d to +%d\n",
                (next_p + prev_p), merit, m, P, D, prev_p, next_p);
            // TODO: Record finished mi in log file / db.
        }

        bool is_last = (mi == last_mi);
        stats.process_results(config, m, is_last, unknown_l, unknown_u, prev_p, next_p, p_tests, n_tests, merit);
    }

    // ----- cleanup

    mpz_clear(K);
}

