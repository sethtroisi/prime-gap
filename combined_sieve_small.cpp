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

#include <clocale>
#include <cstdio>
#include <iostream>

#include <gmp.h>

#include "gap_common.h"
#include "sieve_small.h"

using std::cout;
using std::endl;


int main(int argc, char* argv[]) {
    // Display %'d with commas i.e. 12,345
    setlocale(LC_NUMERIC, "");

    Config config = Args::argparse(argc, argv, Args::Pr::SIEVE);

    if (config.verbose >= 2) {
        printf("\tCompiled with GMP %d.%d.%d\n\n",
            __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);
    }

    // More combined sieve specific validation
    {
        // Both shouldn't be true from gap_common.
        if (!config.testing) {
            cout << "Must set --testing" << endl;
            exit(1);
        }

        if (config.valid == 0) {
            Args::show_usage(argv[0], Args::Pr::SIEVE);
            exit(1);
        }

        if (config.max_prime > 100'000'000) {
            printf("\tmax_prime(%ldM) is probably too large\n",
                config.max_prime / 1'000'000);
        }
    }

    // Status lines
    if (config.verbose >= 0) {
        printf("\n");
        printf("Testing m * %u#/%u, m = %'ld + [0, %'ld)\n",
            config.p, config.d, config.mstart, config.minc);
    }

    if (config.verbose >= 2 && config.threads > 1) {
        printf("Running with %d threads\n", config.threads);
    }

    if (config.method1) {
        cout << "method1 not supported here" << endl;
    } else {
        prime_gap_parallel(config);
    }
}
