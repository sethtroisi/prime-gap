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
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "gap_common.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using namespace std::chrono;



std::map<uint64_t,uint64_t> common_primepi = {
    {     10'000'000,      664'579},
    {    100'000'000,    5'761'455},
    {    200'000'000,   11'078'937},
    {    400'000'000,   21'336'326},
    {    800'000'000,   41'146'179},
    {  1'000'000'000,   50'847'534},
    {  2'000'000'000,   98'222'287},
    {  3'000'000'000,   144'449'537},
    {  4'000'000'000,   189'961'812},
    {  5'000'000'000,   234'954'223},
    {  6'000'000'000,   279'545'368},
    { 10'000'000'000,   455'052'511},
    { 15'000'000'000,   670'180'516},
    { 20'000'000'000,   882'206'716},
    { 25'000'000'000, 1'091'987'405},
    { 30'000'000'000, 1'300'005'926},
    { 40'000'000'000, 1'711'955'433},
    { 50'000'000'000, 2'119'654'578},
    { 60'000'000'000, 2'524'038'155},
    {100'000'000'000, 4'118'054'813}
};


static void assert_file_exists(string path) {
    std::ifstream f(path);
    if (!f.good()) {
        printf("'%s' doesn't exist\n", path.c_str());
        exit(1);
    }
}

bool has_prev_prime_gmp() {
    return (
        (__GNU_MP_VERSION > 6) ||
        (__GNU_MP_VERSION == 6 && __GNU_MP_VERSION_MINOR > 3) ||
        (__GNU_MP_VERSION == 6 && __GNU_MP_VERSION_MINOR == 2 && __GNU_MP_VERSION_PATCHLEVEL == 99)
    );
}


uint32_t gcd(uint32_t a, uint32_t b) {
    if (b == 0) return a;
    return gcd(b, a % b);
}


void K_stats(
        const struct Config& config,
        mpz_t &K, int *K_digits, double *K_log) {
    mpz_init(K);
    mpz_primorial_ui(K, config.p);
    assert(0 == mpz_tdiv_q_ui(K, K, config.d));
    assert(mpz_cmp_ui(K, 1) > 0);  // K <= 1 ?!?

    *K_digits = mpz_sizeinbase(K, 10);

    long exp;
    double mantis = mpz_get_d_2exp(&exp, K);
    *K_log = log(mantis) + log(2) * exp;

    if (config.verbose >= 2) {
        int K_bits   = mpz_sizeinbase(K, 2);
        printf("K = %d bits, %d digits, log(K) = %.2f\n",
            K_bits, *K_digits, *K_log);
    }
}


/**
 * Return estimated time (in seconds) to PRP test a composite with no small factor
 */
double prp_time_estimate_composite(double K_log, int verbose) {
    // TODO: For large K_log, time smaller PRP then upscale with polynomial

    // Some rough estimates at
    // https://github.com/sethtroisi/misc-scripts/tree/master/prime-time

    float K_log_2 = K_log * K_log;
    float t_estimate_poly = -1.1971e-03
                        +  5.1072e-07 * K_log
                        +  9.4362e-10 * K_log_2
                        +  1.8757e-13 * K_log_2 * K_log
                        + -1.9582e-18 * K_log_2 * K_log_2;
    float t_estimate = std::max(1e-3f, t_estimate_poly);

    if (verbose >= 2) {
        if (t_estimate > 0.3) {
            printf("Estimated secs/PRP: %.1f\n", t_estimate);
        } else {
            // Benchmark in thread

            // Create some non-trivial semi-primes.
            mpz_t n, p, q;
            mpz_init(n);
            mpz_init(p);
            mpz_init(q);

            size_t bits = K_log * 1.442;
            assert( bits > 50 );

            // Larger static prime
            mpz_ui_pow_ui(p, 2, bits - 25);
            mpz_nextprime(p, p);

            // Smaller prime for fast nextprime.
            // Large enough to avoid being found with trial division.
            mpz_ui_pow_ui(q, 2, 25);
            mpz_nextprime(q, q);

            double t = 0;
            size_t count = 0;
            // A lot of overhead in nextprime so only time ~1 second of tests.
            for (; count < 10 || t < 1; count++) {
                mpz_mul(n, p, q);
                assert( mpz_sizeinbase(n, 2) >= bits );

                auto  s_start_t = high_resolution_clock::now();

                assert( mpz_probab_prime_p(n, 25) == 0 );

                t += duration<double>(high_resolution_clock::now() - s_start_t).count();
                mpz_nextprime(q, q);
            }

            printf("Estimating PRP/s: %ld / %.2f = %.1f/s vs polyfit estimate of %.1f/s\n",
                count, t, count / t, 1 / t_estimate);
            t_estimate = t / count;
        }
    }

    return t_estimate;
}


/**
 * Handles approx count of divisors by d
 * See "Optimizing Choice Of D" in THEORY.md for why this is required
 */
double prob_gap_larger(
        const struct Config& config, double prob_prime,
        double *prob_prime_coprime_p, size_t *count_coprime_p) {
    vector<uint32_t> P_primes = get_sieve_primes(config.p);
    assert( P_primes.back() == config.p );

    *prob_prime_coprime_p = 1;
    for (uint32_t prime : P_primes) {
        *prob_prime_coprime_p *= (1 - 1.0/prime);
    }

    // Count remaining coprime X for 0 <= X <= SL over SAMPLE_M m values
    const size_t SAMPLES_M = config.d == 1 ? 1 : 6;

    // Unoptimized but cleaner code is preferred here
    mpz_t P, N, temp;
    mpz_init(P);
    mpz_init(N);
    mpz_init(temp);

    mpz_primorial_ui(P, config.p);

    size_t m_counted = 0;
    for (uint64_t mi = 0; m_counted < SAMPLES_M; mi++) {
        uint64_t m = config.mstart + mi;
        if (gcd(m, config.d) > 1) continue;
        m_counted++;

        mpz_mul_ui(N, P, m);
        assert(0 == mpz_tdiv_q_ui(N, N, config.d));

        for (size_t X = 0; X <= config.sieve_length; X++) {
            mpz_gcd(temp, N, P);
            if (mpz_cmp_ui(temp, 1) == 0) {
                *count_coprime_p += 1;
            }
            mpz_add_ui(N, N, 1);
        }
    }
    mpz_clear(N);
    mpz_clear(P);
    mpz_clear(temp);

    *count_coprime_p /= SAMPLES_M;

    double chance_coprime_composite = 1 - prob_prime / *prob_prime_coprime_p;
    return pow(chance_coprime_composite, *count_coprime_p);
}


double prob_prime_and_stats(const struct Config& config, mpz_t &K) {

    int K_digits;
    double K_log;
    K_stats(config, K, &K_digits, &K_log);

    if (config.verbose >= 2) {
        // From Mertens' 3rd theorem
        double unknowns_after_sieve = 1 / (log(config.max_prime) * exp(GAMMA));
        const double N_log = K_log + log(config.mstart);
        const double prob_prime = 1 / N_log - 1 / (N_log * N_log);
        double prob_prime_after_sieve = prob_prime / unknowns_after_sieve;

        size_t count_coprime_p = 0;
        double prob_prime_coprime_p = 0;
        double prob_gap_hypothetical = prob_gap_larger(
            config, prob_prime, &prob_prime_coprime_p, &count_coprime_p);

        float expected = count_coprime_p * (unknowns_after_sieve / prob_prime_coprime_p);
        printf("\n");
        printf("\texpect %.0f left = 2 * %.0f (%.3f%%) of %u after %ldM\n",
                2 * expected, expected,  100.0 * expected / (config.sieve_length + 1),
                config.sieve_length, config.max_prime/1'000'000);
        printf("\t%.3f%% of %d digit numbers are prime\n",
                100 * prob_prime, K_digits);
        printf("\t%.3f%% of tests should be prime (%.1fx speedup)\n",
                100 * prob_prime_after_sieve, 1 / unknowns_after_sieve);
        printf("\t~2x%.1f = %.1f PRP tests per m\n",
                1 / prob_prime_after_sieve, 2 / prob_prime_after_sieve);
        printf("\tsieve_length=%d is insufficient ~~%.3f%% of time\n",
                config.sieve_length, 100 * prob_gap_hypothetical);
        printf("\n");
    }

    return K_log;
}


vector<uint32_t> get_sieve_primes(uint32_t n) {
    vector<uint32_t> primes = {2};
    uint32_t half_n = n >> 1;
    vector<bool> is_prime(half_n + 1, true);

    for (uint32_t p = 3; p <= n; p += 2) {
        if (is_prime[p >> 1]) {
            primes.push_back(p);
            uint64_t p2 = p * p;
            if (p2 > n) break;

            for (uint32_t m = p2 >> 1; m <= half_n; m += p)
                is_prime[m] = false;
        }
    }
    for (uint32_t p = primes.back() + 2; p <= n; p += 2) {
        if (is_prime[p >> 1])
            primes.push_back(p);
    }
    return primes;
}


// Faster because of better memory access patterns
vector<uint64_t> get_sieve_primes_segmented(uint64_t n) {
    assert( n > 10'000 );
    uint64_t sqrt_n = sqrt(n);
    while (sqrt_n * sqrt_n < n) sqrt_n++;

    const vector<uint32_t> small_primes = get_sieve_primes(sqrt_n);

    // First number in next block that primes[pi] divides.
    vector<int32_t> next_mod(small_primes.size(), 0);

    // Large enough to be fast and still fit in L1/L2 cache.
    uint32_t BLOCKSIZE = 1 << 16;
    uint32_t ODD_BLOCKSIZE = BLOCKSIZE >> 1;
    vector<char> is_prime(ODD_BLOCKSIZE, true);

    vector<uint64_t> primes = {2};

    uint32_t max_pi = 0;
    for (uint64_t B = 0; B < n; B += BLOCKSIZE) {
        uint64_t B_END = B + BLOCKSIZE - 1;
        if (B_END > n) {
            BLOCKSIZE = (n - B);
            ODD_BLOCKSIZE = (n - B + 1) >> 1;
            B_END = n;
        }

        while ((max_pi < small_primes.size()) &&
               small_primes[max_pi] * small_primes[max_pi] <= B_END) {
            uint64_t first = small_primes[max_pi] * small_primes[max_pi];
            next_mod[max_pi] = (first - B) >> 1;
            max_pi += 1;
        }

        // reset is_prime
        std::fill(is_prime.begin(), is_prime.end(), true);
        if (B == 0) is_prime[0] = 0; // Skip 1

        // Can skip some large pi up to certain B (would have to set next_mod correctly)
        for (uint32_t pi = 1; pi < max_pi; pi++) {
            const uint32_t prime = small_primes[pi];
            uint32_t first = next_mod[pi];
            for (; first < ODD_BLOCKSIZE; first += prime){
                is_prime[first] = false;
            }
            next_mod[pi] = first - ODD_BLOCKSIZE;
        }
        for (uint32_t prime = 0; prime < ODD_BLOCKSIZE; prime++) {
            if (is_prime[prime]) {
                primes.push_back(B + 2 * prime + 1);
            }
        }
    }
    return primes;
}


static bool isprime_brute(uint32_t n) {
    if ((n & 1) == 0)
        return false;
    for (uint32_t p = 3; p * p <= n; p += 2)
        if (n % p == 0)
            return false;
    return true;
}

void get_sieve_primes_segmented_lambda(uint64_t n, std::function<bool (uint64_t)> lambda) {
    // Large enough to be fast and still fit in L1/L2 cache.
    uint32_t BLOCKSIZE = 1 << 16;
    uint32_t ODD_BLOCKSIZE = BLOCKSIZE >> 1;
    vector<char> is_prime(ODD_BLOCKSIZE, true);

    lambda(2L);

    vector<int32_t> primes = {3};
    // First number in next block that primes[pi] divides.
    vector<int32_t> next_mod = {9 >> 1};

    uint32_t p_lim = 5;
    uint64_t p2_lim = p_lim * p_lim;

    for (uint64_t B = 0; B < n; B += BLOCKSIZE) {
        uint64_t B_END = B + BLOCKSIZE - 1;
        if (B_END > n) {
            BLOCKSIZE = (n - B);
            ODD_BLOCKSIZE = (n - B + 1) >> 1;
            B_END = n;
        }

        while (p2_lim <= B_END) {
            if (isprime_brute(p_lim)) {
                primes.push_back(p_lim);
                assert( p2_lim >= B );
                next_mod.push_back((p2_lim - B) >> 1);
            }
            p2_lim += 4 * p_lim + 4;
            p_lim += 2;
            //assert( p_lim * p_lim == p2_lim );
        }

        // reset is_prime
        std::fill(is_prime.begin(), is_prime.end(), true);
        if (B == 0) is_prime[0] = 0; // Skip 1

        // Can skip some large pi up to certain B (would have to set next_mod correctly)
        for (uint32_t pi = 0; pi < primes.size(); pi++) {
            const uint32_t prime = primes[pi];
            uint32_t first = next_mod[pi];
            for (; first < ODD_BLOCKSIZE; first += prime) {
                is_prime[first] = false;
            }
            next_mod[pi] = first - ODD_BLOCKSIZE;
        }
        for (uint32_t prime = 0; prime < ODD_BLOCKSIZE; prime++) {
            if (is_prime[prime]) {
                if (!lambda(B + 2 * prime + 1)) {
                    return;
                }
            }
        }
    }
}



void Args::show_usage(char* name) {
    cout << "Usage: " << name << endl;
    cout << "[REQUIRED]" << endl;
    cout << "  -p <p>" << endl;
    cout << "  -d <p>" << endl;
    cout << "  --mstart <start>" << endl;
    cout << "  --minc   <int>" << endl;
    cout << "OR" << endl;
    cout << "  --unknown-filename <filename>" << endl;
    cout << "    parse p, d, mstart, minc, sieve-length, max-prime from filename" << endl;
    cout  << endl;
    cout << "[OPTIONALLY]" << endl;
    cout << "  --min-merit <minmerit>" << endl;
    cout << "    only display prime gaps with merit >= minmerit" << endl;
    cout << "  --sieve-length" << endl;
    cout << "    how large the positive/negative sieve arrays should be" << endl;
    cout << "  --max-prime" << endl;
    cout << "    use primes <= max-prime (in millions) for checking composite" << endl;
    cout  << endl;
    //cout << "  --run-prp" << endl;
    //cout << "    run PRP tests" << endl;
    cout << "  --save-unknowns" << endl;
    cout << "    save unknowns to a temp file where they are processed in a 2nd pass." << endl;
    cout << endl;
    cout << "[OPTIONAL]" << endl;
    cout << "  --search-db" << endl;
    cout << "    Database for this project (Default: 'prime-gap-search.db')" << endl;
    cout << "  --prime-gaps-db" << endl;
    cout << "    Prime gap prime gap search db (Default: 'prime-gap-search.db')" << endl;
    cout << endl;
    cout << "  -q, --quiet" << endl;
    cout << "    suppress some status output (twice for more suppression)" << endl;
    cout << "  -h, --help" << endl;
    cout << "    print this help message" << endl;
    cout << endl;
    cout << "calculates prime_gaps for (mstart + mi) * p#/d, mi <= minc " << endl;
}


std::string Args::gen_unknown_fn(const struct Config& config, std::string suffix) {
    if (!config.unknown_filename.empty()) {
        // If trying to read the unknown_fn can cause issue (with having dropped basename)
        // handled specially in gap_stats / gap_test
        string t = config.unknown_filename;
        assert( 0 == t.compare(t.size() - 4, 4, ".txt") );
        return t.replace(t.size() - 4, 4, suffix);
    }

    return std::to_string(config.p) + "_" +
           std::to_string(config.d) + "_" +
           std::to_string(config.mstart) + "_" +
           std::to_string(config.minc) + "_s" +
           std::to_string(config.sieve_length) + "_l" +
           std::to_string(config.max_prime / 1'000'000) + "M" +
           (config.method1 ? ".m1" : "") +
           suffix;
}


Config Args::argparse(int argc, char* argv[]) {
    static struct option long_options[] = {
        {"p",                required_argument, 0,  'p' },
        {"d",                required_argument, 0,  'd' },

        {"mstart",           required_argument, 0,   1  },
        {"minc",             required_argument, 0,   2  },

        {"unknown-filename", required_argument, 0,  'u' },

        {"sieve-length",     required_argument, 0,   4  },
        {"max-prime",        required_argument, 0,   5  },

        {"min-merit",        required_argument, 0,   3  },
        {"run-prp",          no_argument,       0,   6  },
        {"save-unknowns",    no_argument,       0,   7  },

        {"search-db",        required_argument, 0,   9  },
        {"prime-gaps-db",    required_argument, 0,  10  },

        {"method1",          no_argument,       0,   8  },

        {"quiet",            no_argument,       0,  'q' },
        {"help",             no_argument,       0,  'h' },
        {0,                  0,                 0,   0  }
    };

    Config config;
    config.valid = 1;

    int option_index = 0;
    char c;
    while ((c = getopt_long(argc, argv, "qhp:d:", long_options, &option_index)) >= 0) {
        switch (c) {
            case 'h':
                show_usage(argv[0]);
                exit(0);

            case 'q':
                config.verbose--;
                break;

            case 'p':
                config.p = atoi(optarg);
                break;
            case 'd':
                config.d = atoi(optarg);
                break;
            case 1:
                config.mstart = atoi(optarg);
                break;
            case 2:
                config.minc = atoi(optarg);
                break;

            case 'u':
                {
                    char* t = optarg;
                    assert(*t != 0);

                    config.unknown_filename = t;

                    // Don't consider directory for validation.
                    t = basename(t);

                    assert( std::count(t, t + strlen(t), '_')  == 5);

                    config.p = atoi(t);
                    t = std::strchr(t, '_');
                    t++;

                    config.d = atoi(t);
                    t = std::strchr(t, '_');
                    t++;

                    config.mstart = atol(t);
                    t = std::strchr(t, '_');
                    t++;

                    config.minc = atol(t);
                    t = std::strchr(t, '_');
                    assert( t[0] == '_' && t[1] == 's' );
                    t += 2;

                    config.sieve_length = atoi(t);
                    t = std::strchr(t, '_');
                    assert( t[0] == '_' && t[1] == 'l' );
                    t += 2;

                    config.max_prime = atol(t) * 1'000'000;
                    t = std::strchr(t, 'M');

                    config.method1 = (t[3] == '1');

                    assert( std::strcmp(t, "M.txt") == 0 || std::strcmp(t, "M.m1.txt") == 0 );
                }
                break;

            case 3:
                config.min_merit = atof(optarg);
                break;
            case 4:
                config.sieve_length = atoi(optarg);
                break;
            case 5:
                config.max_prime = atol(optarg) * 1'000'000;
                break;

            case 6:
                config.run_prp = true;
                break;
            case 7:
                config.save_unknowns = true;
                break;
            case 8:
                config.method1 = true;
                break;

            case 9:
                config.search_db = optarg;
                assert_file_exists(optarg);
                break;
            case 10:
                config.records_db = optarg;
                assert_file_exists(optarg);
                break;

            case 0:
                printf("option %s arg %s\n", long_options[option_index].name, optarg);
                config.valid = 0;
                break;
            case '?':
                config.valid = 0;
                break;
            default:
                config.valid = 0;
                printf("getopt returned \"%d\"\n", c);
        }
    }

    if (optind < argc) {
        config.valid = 0;
        printf("unknown positional arguements: ");
        while (optind < argc) {
            printf("%s ", argv[optind++]);
        }
        printf("\n");
    }

    // ----- Validation

    if (config.mstart <= 0) {
        config.valid = 0;
        cout << "mstart must be greater than 0: " << config.mstart << endl;
    }

    uint64_t last_m = ((long) config.mstart) + config.minc;
    if (last_m >= 1'000'000'001 ) {
        config.valid = 0;
        cout << "mstart + minc must be < 1e9" << endl;
    }

    if (config.minc <= 0) {
        config.valid = 0;
        cout << "minc must be greater than 0: " << config.minc << endl;
    }

    if (config.minc >= 50'000'000) {
        config.valid = 0;
        cout << "minc > 50M will use to much memory" << endl;
    }

    if (config.max_prime > 1'000'000'000'000) {
        // This is kinda arbitrary.
        config.valid = 0;
        cout << "max_prime > 1000B not supported" << endl;
    }

    uint64_t max_m = (1UL << 62) / config.max_prime;
    if (max_m <= last_m) {
        config.valid = 0;
        cout << "max_prime * last_m(" << last_m << ") might overflow int64" << endl;
    }

    {
        // check if p is valid
        bool valid = config.p < 1'000'0000;
        for (size_t t = 2; valid && t*t <= config.p; t++) {
            valid = (config.p % t) > 0;
        }

        if (!valid) {
            config.valid = 0;
            cout << "p# not prime (p=" << config.p << ")" << endl;
            exit(1);
        }
    }

    {
        // Check that SL % 1000 = 0 or (SL, K) = 1
        size_t SL = config.sieve_length;
        if (SL % 1000 != 0) {
            for (size_t p = 2; p < config.p; p += 1 + (p > 2)) {
                if (isprime_brute(p)) {
                    if (SL % p == 0 && config.d % p != 0) {
                        config.valid = 0;
                        cout << "SL=" << SL << " not coprime (p=" << p << ")" << endl;
                    }
                }
            }
        }
    }

    if (config.d <= 0) {
        config.valid = 0;
        cout << "d must be greater than 0: " << config.d << endl;
    }

    if (config.valid == 0) {
        cout << endl;
    }

    return config;
}


DB::DB(const char* path) {
    assert_file_exists(path);

    if (sqlite3_open(path, &db) != SQLITE_OK) {
        printf("Can't open database(%s): %s\n", path, sqlite3_errmsg(db));
        exit(1);
    }
};


uint64_t DB::config_hash(const struct Config& config) {
    // Hash the config to a uint64
    uint64_t hash =    config.mstart;
    hash = hash * 31 + config.minc;
    hash = hash * 31 + config.p;
    hash = hash * 31 + config.d;
    hash = hash * 31 + config.sieve_length;
    hash = hash * 31 + config.max_prime;
    return hash;
}

