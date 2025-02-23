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

#pragma once

#include <cstdio>
#include <memory>
#include <utility>
#include <vector>

#include "gap_common.h"

using std::vector;

class SieveOutput {
    public:
        // Prevent copying which would use lots of memory...
        SieveOutput(const SieveOutput&) = delete;
        void operator=(const SieveOutput&) = delete;

        SieveOutput(uint64_t m_start, int32_t sieve_length):
            m_start(m_start), sieve_length(sieve_length) {};

        /* Was used for debug.
        ~SieveOutput() {
            printf("~SieveOutput\n");
        }
        // */

        const uint64_t m_start;
        // Is there anything also needed from config?
        const int32_t sieve_length;

        // unknowns are delta run length encoded into this table.
        vector<uint16_t> coprime_X;

        // m_increment, unknown count
        vector<std::tuple<int8_t, uint8_t>> m_inc;

        vector<vector<uint16_t>> unknowns;
};

std::unique_ptr<SieveOutput> prime_gap_parallel(const struct Config& config);
