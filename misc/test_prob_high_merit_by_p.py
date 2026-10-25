#!/usr/bin/env python
#
# Copyright 2026 Seth Troisi
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Check probabilities of large record at various P#, d, merits.

Requires hand modifying sieve_limit at combined_sieve.cpp:97 limit
sqlite3 high-prob-stats.db < schema.sql
"""

import argparse
import os
import subprocess
import time

import gmpy2


def get_arg_parser():
    parser = argparse.ArgumentParser('Test prob of record with merit > X for various P#/d')

    parser.add_argument(
        '--start_p', type=int, default=300,
        help="Smallest P# to check.")

    parser.add_argument(
        '--stop_p', type=int, default=1000,
        help="Largest P# to check.")

    parser.add_argument(
        '-d', '--divisors', type=int, nargs='+', default=[30, 210],
        help='list of divisors to test fo each P#')

    parser.add_argument(
        '-m', '--min-merit', type=int, nargs='+', default=[30],
        help='Check for probability of record > this merit')

    return parser


def run_sieves(args):
    jobs = []
    for p in range(args.start_p, args.stop_p + 1):
        if not gmpy2.is_prime(p):
            continue

        for d in args.divisors:
            K = gmpy2.primorial(p) // d
            log_K = gmpy2.log(K)

            # Rounded up to next thousand
            sl = int(log_K * max(args.min_merit) - 1)
            sl += (-sl % 1000)

            # gap_stats looks at log(m) so use 9e9 to account for it's contribution.
            unknown_fn = f"{p}_{d}_900000000_50000_s{sl}_l10M.txt"
            jobs.append((p, d, unknown_fn))

            path = os.path.join("unknowns", unknown_fn);
            if os.path.exists(path):
                print(f"\t{unknown_fn} already sieved")
                continue

            t0 = time.time()
            subprocess.check_call(
                    ["./combined_sieve",
                     "--save", "-q", "-q",
                     "--search-db", "high-prob-stats.db",
                     "-u", unknown_fn])
            t1 = time.time()
            print(f"\t{unknown_fn} took {t1 - t0:.1f} seconds")

    return jobs


def gap_stats(args, jobs):
    probs = {}
    for p, d, unknown_fn in jobs:
        for m in args.min_merit:
            result = subprocess.check_output(
                    ["./gap_stats",
                     "--testing",
                     "--min-merit", str(m),
                     "-u", unknown_fn])
            print(f"Gap Stats for {unknown_fn!r}")
            prob = None
            for line in result.decode().split("\n"):
                if line.startswith("prob merit total:"):
                    prob = float(line.split(" ")[-1])
            assert prob is not None, result
            probs[(p, d, m)] = prob
            print(f"\t{p}#/{d} @ {m} = {prob:.2e}")

    print(probs)

    return probs


if __name__ == "__main__":
    parser = get_arg_parser()
    args = parser.parse_args()

    jobs = run_sieves(args)
    gap_stats(args, jobs)
