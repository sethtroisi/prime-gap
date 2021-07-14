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

#include "gap_common.h"

using std::cout;
using std::endl;
using std::vector;


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
