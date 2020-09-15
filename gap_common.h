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

#include <functional>
#include <map>
#include <string>
#include <vector>

#include <gmp.h>


using std::map;
using std::string;
using std::vector;

const double GAMMA = 0.577215665;

extern std::map<uint64_t,uint64_t> common_primepi;


/* Arg Parsing */

struct Config {
    int valid   = 0;
    uint64_t mstart  = 0;
    uint64_t minc    = 0;
    uint32_t p       = 0;
    uint32_t d       = 0;
    float min_merit = 12;

    uint32_t sieve_length = 0;
    uint64_t sieve_range  = 0;

    bool run_prp = false;
    bool save_unknowns = false;

    bool method1 = false;

    /**
     * -1: results only
     *  0: results & final stats
     *  1: stats
     *  2: stats, probs, debug
     *  3: ???
     */
    int verbose = 2;

    string unknown_filename;
};

void show_usage(char* name);
Config argparse(int argc, char* argv[]);

std::string gen_unknown_fn(const struct Config& config, std::string suffix);


/* Random Utils */

uint32_t gcd(uint32_t a, uint32_t b);

void K_stats(
        const struct Config& config,
        mpz_t &K, int *K_digits, double *K_log);

double prp_time_estimate(double K_log);

double prop_gap_larger(
    const struct Config& config,
    double prob_prime,
    double *prob_prime_coprime,
    size_t *count_coprime);



/* Prime Stuff */

vector<uint32_t> get_sieve_primes(uint32_t n);
vector<uint64_t> get_sieve_primes_segmented(uint64_t n);
void             get_sieve_primes_segmented_lambda(uint64_t n, std::function<bool (uint64_t)> lambda);
