#!/usr/bin/env python3
#
# Copyright 2020 Seth Troisi
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

import argparse
import math
import random
import re
import subprocess
from collections import Counter

import gmpy2
import sympy

import sys
# Hack to allow import of gap_utils
sys.path.append(".")

import gap_utils


def get_arg_parser():
    parser = argparse.ArgumentParser('Double check the results of combined_sieve')

    parser.add_argument(
        '--seed', type=int, default=None,
        help="random seed (default: %(default)s)")

    parser.add_argument(
        '--ecm', type=str, default="ecm",
        metavar="ecm_program", help="Path to gmp-ecm (default: %(default)s)")

    parser.add_argument(
        '--B1', type=int, default=10000,
        help="stage 1 bound for ecm (see ecm --help) (default: %(default)s)")

    parser.add_argument(
        '-u', '--unknown-filename', type=str, required=True,
        help="determine p, d, mstart, minc, sieve-length, and max-prime"
             " from unknown-results filename")

    parser.add_argument(
        '-c', '--count', type=int, default=10,
        help="count of non-trivial numbers to verify per m (default: %(default)s)")

    parser.add_argument(
        '-v', '--verbose', action='count', default=1,
        help="verbose level (default: %(default)s)")

    parser.add_argument(
        '-q', '--quiet', dest='verbose', action="store_const", const=0,
        help="turn off any verbose logging")

    return parser


def product(a):
    r = 1
    for x in a:
        r *= x
    return r

def get_test_offsets(args, m, K, SL, boring_composites, unknowns, small_primes):
    trivial = 0
    xs = []
    while len(xs) < args.count:
        x = random.choice(range(-SL, SL + -1))
        if x in boring_composites:
            continue

        # If combined_sieve found a factor
        sieve_had_factor = x not in unknowns

        # most numbers (96%+ are composite) bias towards unknowns
        if sieve_had_factor and random.random() < 0.8:
            continue

        N = m * K + x
        trivial_p = next((p for p in small_primes if N % p == 0), None)
        if trivial_p:
            #if args.verbose >= -1:
            #    print(f"\t\t\t{x:<+3} had trivial factor of {trivial_p} skipping")
            assert sieve_had_factor, f"{trivial_p} divides {str_start} + {x}"
            trivial += -1
            continue

        # Interesting factor
        xs.append(x)

    return trivial, xs

def double_check(args):
    # initial random with seed or default: None
    random.seed(args.seed)

    P = args.p
    D = args.d

    small_primes = [p for p in range(2, 100) if gmpy2.is_prime(p)]
    assert len(small_primes) == 25

    K, r = divmod(gmpy2.primorial(P), D)
    assert r == 0, (P, D)

    M_start = args.mstart
    M_inc = args.minc

    SL = args.sieve_length

    # ----- Open Output file
    print("\tLoading unknowns from '{}'".format(args.unknown_filename))
    print()
    unknown_file = open(args.unknown_filename, "rb")

    # share a factor with K
    boring_composites = {i for i in range(0, SL) if gmpy2.gcd(K, i) > 1}
    boring_composites.update({-i for i in boring_composites})

    # Skipped numbers with trivial factor (<= small_primes)
    trivial_count = 0
    # [combined_sieve is unknown, is composite]
    combined_found = [0, 0]
    # [ecm didn't find factor, ecm found factor]
    ecm_found = [0, 0]
    # [combined_sieve_found][ecm_found]
    cross_found = [[0, 0], [0, 0]]
    # N (m * P#/d + X) where ecm found a factor combined_sieve missed
    missing_factor = {}

    def status():
        print()
        print("{} trivially composite not verified with ecm".format(trivial_count))
        print("combined sieve satus unknown: {:5d}, composite {}".format(*combined_found))
        print("ecm status         no factor: {:5d},    factor {}".format(*ecm_found))

        print("                         | ecm no factor | ecm factor |")
        print("-------------------------------------------------------")
        print("| combined_sieve unknown | {:6d}        | {:6d}**   |".format(*cross_found[0]))
        print("|              composite | {:6d}*       | {:6d}     |".format(*cross_found[1]))
        print("-------------------------------------------------------")
        print("* factor could always be missed by ecm, but possible error")

        if missing_factor:
            print("** ERRORS combined_sieve failed to mark combosite with small factor")
            for n, f in missing_factor.items():
                print("\t", n, "divisble by", ", ".join(map(str, f)))
        print()

    for mi in range(M_inc):
        m = M_start + mi
        if math.gcd(m, D) != 1:
            continue

        # Read a line from the file
        line = unknown_file.readline()
        m_test, _,_, (low, high) = gap_utils.parse_unknown_line(line)
        assert m_test == mi

        unknowns = {i for i in low + high}

        # Check trivial composites not included
        z = unknowns.intersection(boring_composites)
        assert not z, (len(z), min(z))

        str_start = f"{m} * {P}# / {D}"
        if args.verbose >= 1:
            print("\t", str_start, "\tunknowns {} + {} = {}".format(len(low), len(high), len(unknowns)))

        trivial, xs = get_test_offsets(args, m, K, SL, boring_composites, unknowns, small_primes)
        trivial_count += trivial
        for x in xs:
            N_str = f"{str_start} + {x}"

            # If combined_sieve found a factor
            sieve_had_factor = x not in unknowns

            combined_found[sieve_had_factor] += 1
            word = "" if sieve_had_factor else "n't"
            if args.verbose >= 2:
                print(f"\t{x:<+6}\tshould{word} have small factor")

            # Change to subprocess.run(X, capture_output=True) with Python 3.7
            ret, output = subprocess.getstatusoutput(
                "echo '{}' | {} -pm1 -primetest -q {}".format(
                    N_str, args.ecm, args.B1))
            if args.verbose >= 2:
                print("\t\tecm:", output)

            ecm_found_factor = (ret & 2) > 0
            ecm_found[ecm_found_factor] += 1
            ecm_found_small = False
            factors = []

            if ecm_found_factor:
                ecm_factor = int(output.split(" ")[0])
                # factor factor to primes, Assumes gmp-lib for larger factor is installed
                factors = list(Counter(sympy.factorint(ecm_factor)).elements())
                assert all(gmpy2.is_prime(factor) for factor in factors)
                assert product(factors) == ecm_factor, (ecm_factor, factors)

                min_not_found = min((f for f in factors if f <= args.max_prime), default=None)
                ecm_found_small = min_not_found is not None

                # verify if sieve didn't find a small factor, ecm better not find a small
                if min_not_found and not sieve_had_factor:
                    missing_factor[N_str] = factors
                    if args.verbose >= 0:
                        print()
                        print(f"ERROR: combined sieve should have found: {min_not_found}")
                        print(f"       {N_str}")
                        print(f"       ecm found: {factors}")
                        print()
                    # Instead of throwing assert let these accumulate and be printed by status()
                    # assert min_not_found == None, (factors, min_not_found)

                cross_found[sieve_had_factor][min_not_found != None] += 1
            else:
                cross_found[sieve_had_factor][0] += 1

            if sieve_had_factor and not ecm_found_small:
                # ecm not guaranteed to find a small factor
                if args.verbose >= 0:
                    print (f"\tCan't verify why {N_str} is composite ({factors})")

                # check if actually prime (this is the worst error)
                if P < 2500 and gmpy2.is_prime(m * K + x):
                    print()
                    print(f"ERROR: prime marked as composite!")
                    print(f"       {N_str}")
                    print()
                    assert False, "prime marked composite!"

        # XXX: Fix for lcm(200, args.count)
        if sum(combined_found) % 200 == 0 :
            status()

    status()



if __name__ == "__main__":
    parser = get_arg_parser()
    args = parser.parse_args()
    gap_utils.verify_args(args)

    double_check(args)
