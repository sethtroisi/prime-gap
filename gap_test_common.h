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

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <chrono>

#include "gap_common.h"

using std::cout;
using std::endl;
using std::vector;

/** See gap_test_stats.py */
class StatsCounters {
    public:
        StatsCounters(std::chrono::high_resolution_clock::time_point now) : s_start_t(now) {}

        const std::chrono::high_resolution_clock::time_point s_start_t;

        uint32_t  s_tests     = 0;
        size_t    s_total_unknown = 0;
        size_t    s_t_unk_prev = 0;
        size_t    s_t_unk_next = 0;
        size_t    s_total_prp_tests = 0;
        size_t    s_gap_out_of_sieve_prev = 0;
        size_t    s_gap_out_of_sieve_next = 0;
        float     s_best_merit_interval = 0;
        size_t    s_best_merit_interval_m = 0;

        // This can change in const possible_print_stats
        mutable float     s_tests_per_second = 0;

        /** Return if stats were printed */
        void process_results(
            const Config &config,
            long m, bool is_last,
            size_t unknown_l, size_t unknown_u,
            int prev_p, int next_p,
            int p_tests, int n_tests,
            float merit);

        bool possible_print_stats(
            const Config &config,
            long m, bool is_last,
            size_t unknown_l, size_t unknown_u,
            int prev_p, int next_p) const;
};

void load_and_verify_unknowns(
        const int compression,
        const uint64_t m,
        const int SIEVE_LENGTH,
        std::ifstream &unknown_file,
        vector<int32_t> (&unknowns)[2]);

void test_interval_cpu(
        const uint64_t m, const mpz_t &K, const size_t SIEVE_LENGTH,
        size_t &s_total_prp_tests,
        size_t &s_gap_out_of_sieve_prev, size_t &s_gap_out_of_sieve_next,
        vector<int32_t> (&unknowns)[2],
        int &prev_p, int &next_p);
