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
#include <cstdio>
#include <getopt.h>
#include <iostream>
#include <vector>

#include <gmp.h>

#include "gap_common.h"

using std::cout;
using std::endl;
using std::vector;


vector<uint32_t> get_sieve_primes(uint32_t n) {
    vector<uint32_t> primes = {2};
    vector<bool> is_prime((n+1) >> 1, true);
    uint32_t half_n = n >> 1;

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


bool isprime_brute(uint32_t n) {
    if ((n & 1) == 0)
        return false;
    for (uint32_t p = 3; p * p <= n; p += 2)
        if (n % p == 0)
            return false;
    return true;
}

void get_sieve_primes_segmented_lambda(uint64_t n, std::function<void (uint64_t)> lambda) {
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
                lambda(B + 2 * prime + 1);
            }
        }
    }
}


void show_usage(char* name) {
    cout << "Usage: " << name << endl;
    cout << "[REQUIRED]" << endl;
    cout << "  -p <p>" << endl;
    cout << "  -d <p>" << endl;
    cout << "  --mstart <start>" << endl;
    cout << "  --minc   <int>" << endl;
    cout << "[OPTIONALLY]" << endl;
    cout << "  --minmerit <minmerit>" << endl;
    cout << "    only display prime gaps with merit >= minmerit" << endl;
    cout << "  --sieve-length" << endl;
    cout << "    how large the positive/negative sieve arrays should be" << endl;
    cout << "  --sieve-range" << endl;
    cout << "    use primes <= sieve-range million for checking composite" << endl;
    cout  << endl;
    cout << "  --sieve-only" << endl;
    cout << "    only sieve ranges, don't run PRP. useful for benchmarking" << endl;
    cout << "  --save-unknowns" << endl;
    cout << "    save unknowns to be PRPed to a temp file (p_d_mstart_minc_sieve_range.txt)" << endl;
    cout << "    where they can be processed in a 2nd pass." << endl;
    cout << "  -h, --help" << endl;
    cout << "    print this help message" << endl;
    cout << endl;
    cout << "calculates prime_gaps for (mstart + mi) * p#/d, mi <= minc " << endl;
}


Config argparse(int argc, char* argv[]) {
    // TODO add print_interval option.

    static struct option long_options[] = {
        {"mstart",        required_argument, 0,   1  },
        {"minc",          required_argument, 0,   2  },
        {"p",             required_argument, 0,  'p' },
        {"d",             required_argument, 0,  'd' },

        {"sieve-length",  required_argument, 0,   4  },
        {"sieve-range",   required_argument, 0,   5  },

        {"min-merit",     required_argument, 0,   3  },
        {"sieve-only",    no_argument,       0,   6  },
        {"save-unknowns", no_argument,       0,   7  },

        {"help",          no_argument,       0,  'h' },
        {0,               0,                 0,   0  }
    };

    Config config;
    config.valid = 1;

    int option_index = 0;
    char c;
    while ((c = getopt_long(argc, argv, "hp:d:", long_options, &option_index)) >= 0) {
        switch (c) {
            case 'h':
                show_usage(argv[0]);
                exit(0);
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

            case 3:
                config.minmerit = atof(optarg);
                break;
            case 4:
                config.sieve_length = atoi(optarg);
                break;
            case 5:
                config.sieve_range = atol(optarg) * 1'000'000;
                break;

            case 6:
                config.run_prp = false;
                break;
            case 7:
                config.save_unknowns = true;
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
    if (last_m >= 0x7FFFFFFF ) {
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

    if (config.sieve_range > 100'000'000'000) {
        config.valid = 0;
        cout << "sieve_range > 100B not supported" << endl;
    }

    {
        mpz_t ptest;
        mpz_init_set_ui(ptest, config.p);
        bool valid = 0 != mpz_probab_prime_p(ptest, 25);
        mpz_clear(ptest);
        if (!valid) {
            config.valid = 0;
            cout << "p# not prime (p=" << config.p << ")" << endl;
        }
    }

    if (config.d <= 0) {
        config.valid = 0;
        cout << "d must be greater than 0: " << config.d << endl;
    }

    uint64_t max_m = (1UL << 63) / config.sieve_range;
    if (max_m <= last_m) {
        config.valid = 0;
        cout << "sieve_range * last_m(" << last_m << ") will overflow int64" << endl;
    }

    return config;
}
