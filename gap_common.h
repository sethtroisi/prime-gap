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

#include <cassert>
#include <cstdint>
#include <fstream>
#include <functional>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <gmp.h>
#include <sqlite3.h>


using std::map;
using std::string;
using std::vector;

const double GAMMA = 0.577215665;

/* Arg Parsing */

struct Config {
    int valid   = 0;
    uint32_t p       = 0;
    uint32_t d       = 0;

    uint64_t mstart  = 0;
    uint64_t minc    = 0;

    float min_merit = 12;

    uint32_t sieve_length = 0;
    uint64_t max_prime    = 0;

    bool save_unknowns = false;

    bool method1 = false;
#ifdef RLE
    bool rle = true;
#else
    bool rle = false;
#endif

    int threads = 1;

    /**
     * -1: results only
     *  0: results & final stats
     *  1: stats
     *  2: stats, probs, debug
     *  3: ???
     */
    int verbose = 2;

    string unknown_filename;

    string search_db  = "prime-gap-search.db";
    string records_db = "gaps.db";
};


class Args
{
    public:
        static void show_usage(char* name);
        static Config argparse(int argc, char* argv[]);
        static std::string gen_unknown_fn(const struct Config& config, std::string suffix);
        static bool is_rle_unknowns(std::ifstream& unknown_file);

    private:
        // Disallow creating instance
        Args() = default;
};


class DB
{
    public:
        explicit DB(const char* path);
        ~DB() { if (db != nullptr) sqlite3_close(db); };

        uint64_t    config_hash(const struct Config& config);
        sqlite3*    get_db() { assert(db != nullptr); return db; };

    private:
        sqlite3 *db = nullptr;
};


/* Random Utils */

bool has_prev_prime_gmp();

uint32_t gcd(uint32_t a, uint32_t b);


/* K stuff */
double calc_log_K(const struct Config& config);

void init_K(const struct Config& config, mpz_t &K);

void K_stats(
        const struct Config& config,
        mpz_t &K, int *K_digits, double *K_log);


/* Utils */
std::pair<uint32_t, uint32_t> calculate_thresholds_method2(
        struct Config config,
        size_t count_coprime_sieve,
        size_t valid_ms);


double prp_time_estimate_composite(double K_log, int verbose);

double combined_sieve_method2_time_estimate(
        const struct Config& config,
        const mpz_t &K,
        uint64_t valid_ms,
        double prp_time_est);



// Used to optimize d
std::tuple<double, uint32_t, double, double>
count_K_d(const struct Config& config);

double prob_prime_and_stats(const struct Config& config, mpz_t &K);


/* Prime Stuff */

bool is_prime_brute(uint32_t n);
vector<uint32_t> get_sieve_primes(uint32_t n);

size_t primepi_estimate(uint64_t max_prime);
