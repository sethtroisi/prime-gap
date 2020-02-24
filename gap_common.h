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

#pragma once

#include <vector>
#include <functional>
#include <map>

using std::vector;

struct Config {
    int valid   = 0;
    int mstart  = 0;
    int minc    = 0;
    int p       = 0;
    int d       = 0;
    float minmerit = 12;

    unsigned int sieve_length = 0;
    unsigned long sieve_range  = 0;

    bool run_prp = true;
    bool save_unknowns = false;
};

extern std::map<uint64_t,uint64_t> common_primepi;

void show_usage(char* name);
Config argparse(int argc, char* argv[]);
vector<uint32_t> get_sieve_primes(uint32_t n);
vector<uint64_t> get_sieve_primes_segmented(uint64_t n);
void             get_sieve_primes_segmented_lambda(uint64_t n, std::function<void (uint64_t)> lambda);
