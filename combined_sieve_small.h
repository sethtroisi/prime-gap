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

#include <memory>

#include "gap_common.h"


class SieveOutput {
    public:
        // Prevent copying which would use lots of memory...
        SieveOutput(const SieveOutput&) = delete;
        void operator=(const SieveOutput&) = delete;

        SieveOutput(uint64_t m_start, int32_t sieve_length):
            m_start(m_start), sieve_length(sieve_length) {};

        const uint64_t m_start;
        // Is there anything also needed from config?
        const int32_t sieve_length;

        // m_increment, unknown count
        vector<std::tuple<int16_t, int16_t>> m_inc;

        // TODO change to uint8_t with 0 as a +255 jump without testing.
        vector<vector<uint16_t>> unknowns;
};

std::unique_ptr<SieveOutput> prime_gap_parallel(const struct Config& config);
