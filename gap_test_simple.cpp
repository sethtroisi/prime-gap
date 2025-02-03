// Copyright 2020 Seth Troisi
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
#include <sstream>
#include <string>
#include <vector>

#include <gmp.h>

#include "gap_common.h"
#include "gap_test_common.h"

using std::cout;
using std::endl;
using std::vector;
using namespace std::chrono;


void prime_gap_test(
        const struct Config config,
        std::ifstream &unknown_file);


int main(int argc, char* argv[]) {
    Config config = Args::argparse(argc, argv, Args::Pr::TEST_SIMPLE);
    if (config.valid == 0) {
        Args::show_usage(argv[0], Args::Pr::TEST_SIMPLE);
        return 1;
    }

    if (config.verbose >= 2) {
        printf("Compiled with GMP %d.%d.%d\n",
            __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);
    }

    if( !has_prev_prime_gmp() ) {
        cout << "Error: no mpz_prevprime!" << endl;
        cout << "See Notes in README.md for instructions on using dev GMPlib" << endl;
        return 1;
    }

    if (config.sieve_length == 0) {
        cout << "Must set sieve-length for " << argv[0] << endl;
        Args::show_usage(argv[0], Args::Pr::TEST_SIMPLE);
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

    // ----- Open unknown input file
    std::ifstream unknown_file;
    {
        std::string fn = Args::gen_unknown_fn(config, ".txt");
        if (config.verbose >= 1) {
            printf("\nReading unknowns from '%s'\n", fn.c_str());
        }
        unknown_file.open(fn, std::ios::in);
        assert( unknown_file.is_open() ); // Can't open save_unknowns file
        assert( unknown_file.good() );    // Can't open save_unknowns file

        // Ruins ability to make config const
        config.compression = Args::guess_compression(config, unknown_file);
    }

    prime_gap_test(config, unknown_file);
}


void prime_gap_test(
        const struct Config config,
        std::ifstream &unknown_file) {

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

    // For compressed lines
    BitArrayHelper helper(config, K);

    uint64_t valid_ms = 0;
    for (uint64_t mi = 0; mi < M_inc; mi++) {
        if (gcd(M_start + mi, D) == 1) {
            valid_ms++;
        }
    }
    assert(valid_ms > 0);

    uint64_t first_mi = 0;
    for (; first_mi > 0 && gcd(M_start + first_mi, D) > 1; first_mi++);
    assert(first_mi >= 0 && first_mi < M_inc);

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
        if (gcd(m, D) > 1) continue;

        std::string line;
        // Loop can be pragma omp parallel if this is placed in critical section
        std::getline(unknown_file, line);
        std::istringstream iss_line(line);

        if ((size_t) m < config.mskip) continue;

        vector<int32_t> unknowns[2];
        uint64_t m_parsed = parse_unknown_line(
            config, helper, m, iss_line, unknowns[0], unknowns[1]);
        assert(m_parsed == (uint64_t) m);

        size_t unknown_l = unknowns[0].size();
        size_t unknown_u = unknowns[1].size();

        size_t p_tests = 0, n_tests = 0;

        int prev_p = 0;
        int next_p = 0;
        size_t s_gap_out_of_sieve_prev, s_gap_out_of_sieve_next;

        test_interval_cpu(
            m, K, SIEVE_LENGTH,
            p_tests,
            s_gap_out_of_sieve_prev, s_gap_out_of_sieve_next,
            unknowns, prev_p, next_p);
        assert( prev_p > 0 && next_p > 0 );

        int gap  = next_p + prev_p;
        float merit = gap / (K_log + log(m));

        if (merit > min_merit)  {
            // TODO: write to file or database
            printf("%-5d %.4f  %ld * %ld#/%ld -%d to +%d\n",
                gap, merit, m, P, D, prev_p, next_p);
        }

        bool is_last = (mi == last_mi);
        stats.process_results(config, m, is_last, unknown_l, unknown_u, prev_p, next_p, p_tests, n_tests, merit);
    }

    // ----- cleanup

    mpz_clear(K);
}

