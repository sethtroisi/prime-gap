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

#include <cassert>
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
    vector<bool> is_prime(n+1, true);
    for (uint32_t p = 3; p <= n; p += 2) {
        if (is_prime[p]) {
            primes.push_back(p);
            uint64_t p2 = p * p;
            if (p2 > n) break;

            for (uint32_t m = p2; m <= n; m += 2*p)
                is_prime[m] = false;
        }
    }
    for (uint32_t p = primes.back() + 2; p <= n; p += 2) {
        if (is_prime[p])
            primes.push_back(p);
    }
    return primes;
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
    cout << "    use primes <= sieve-range for checking composite" << endl;
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
                config.sieve_range = atol(optarg);
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

    if (((long) config.mstart + config.minc) >= MAX_INT) {
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

    if (config.sieve_range > 4'000'000'000) {
        config.valid = 0;
        cout << "sieve_range > 4B not supported" << endl;
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

    return config;
}
