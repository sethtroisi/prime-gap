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

const double GAMMA = 0.577215665;

struct Config {
    int valid   = 0;
    uint64_t mstart  = 0;
    uint64_t minc    = 0;
    uint32_t p       = 0;
    uint32_t d       = 0;
    float minmerit = 12;

    uint32_t sieve_length = 0;
    uint64_t sieve_range  = 0;

    bool run_prp = true;
    bool save_unknowns = false;

    bool method2 = false;
};

extern std::map<uint64_t,uint64_t> common_primepi;

std::string gen_unknown_fn(const struct Config& config, std::string suffix);

void show_usage(char* name);
Config argparse(int argc, char* argv[]);

uint32_t gcd(uint32_t a, uint32_t b);

double prop_gap_larger(
    const struct Config& config,
    double prob_prime,
    double *prob_prime_coprime,
    size_t *count_coprime);

vector<uint32_t> get_sieve_primes(uint32_t n);
vector<uint64_t> get_sieve_primes_segmented(uint64_t n);
void             get_sieve_primes_segmented_lambda(uint64_t n, std::function<void (uint64_t)> lambda);
