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

import gmpy2

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
        '--unknown-filename', type=str, required=True,
        help="determine p, d, mstart, minc, sieve-length, and max-prime"
             " from unknown-results filename")

    parser.add_argument(
        '-c', '--count', type=int, default=10,
        help="count of non-trivial numbers to verify per m (default: %(default)s)")

    return parser


def product(a):
    r = 1
    for x in a:
        r *= x
    return r


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
    unknown_file = open(args.unknown_filename, "r")

    # share a factor with K
    boring_composites = {i for i in range(0, SL) if gmpy2.gcd(K, i) > 1}
    boring_composites.update({-i for i in boring_composites})

    tested = [0, 0, 0]
    ecm_found = [0, 0]
    for mi in range(M_inc):
        m = M_start + mi
        if math.gcd(m, D) != 1:
            continue

        # Read a line from the file
        line = unknown_file.readline()
        start, u_l, u_h = line.split("|")

        match = re.match(r"^([0-9]+) : -([0-9]+) \+([0-9]+)", start)
        assert match, (len(start), len(u_l), len(u_h))
        m_test, unknown_l, unknown_u = map(int, match.groups())
        assert m_test == mi

        low = list(map(int, u_l.strip().split(" ")))
        high = list(map(int, u_h.strip().split(" ")))

        unknowns = {i for i in low + high}

        # Check trivial composites not included
        z = unknowns.intersection(boring_composites)
        assert not z, (len(z), min(z))

        str_start = f"{m} * {P}# / {D}"
        print(str_start, "\tunknowns {} + {} = {}".format(len(low), len(high), len(unknowns)))

        # Choose some random numbers
        count = args.count
        while count > 0:
            t = random.choice(range(-SL, SL + 1))
            if t in boring_composites:
                tested[0] += 1
                continue

            # If sieve found a small factor
            sieve_had_factor = t not in unknowns

            N = (m * K + t)
            found = False
            for p in small_primes:
                if N % p == 0:
                    print(f"\t\t\t{t:<+6} had trivial factor of {p} skipping")
                    assert sieve_had_factor, f"{p} divides {str_start} + {t}"
                    found = True
                    break
            if found:
                continue

            tested[1 + sieve_had_factor] += 1
            word = "" if sieve_had_factor else "n't"
            print(f"\t{t:<+6}\tshould{word} have small factor")

            # Interesting factor
            count -= 1

            # Change to subprocess.run(X, capture_output=True) with Python 3.7
            ret, output = subprocess.getstatusoutput(
                "echo '{} + {}' | {} -one -primetest -q {}".format(
                    str_start, t, args.ecm, args.B1))

            found_factor = (ret & 2) > 0

            print("\t\tecm:", output)
            if found_factor:
                ecm_factor = int(output.split(" ")[0])

                # factor factor to primes, Assumes gmp-lib for larger factor is installed
                factor_output = subprocess.check_output(["factor", str(ecm_factor)])
                factors = [int(f) for f in factor_output.decode().split(":")[1].strip().split(" ")]
                if len(factors) > 1:
                    print("\t\tfactor", factor_output.decode().strip())
                assert all(gmpy2.is_prime(factor) for factor in factors)

                assert product(factors) == ecm_factor, (ecm_factor, factors)

                # verify if sieve didn't find a small ecm better not find a small
                if not sieve_had_factor:
                    assert all(f > args.max_prime for f in factors), (ecm_factor, factors)

                # Tempting to check if sieve_had_factor then any(f < args.max_prime)
                # but ecm not guaranteed to find that factor.

                # generally expected ecm to find the
                ecm_found[sieve_had_factor] += 1

    print("{} trivially composite, {} unknowns, {} known composites".format(*tested))
    print("ecm found {} composites: known {} composite, {} were unknown".format(
        sum(ecm_found), ecm_found[True], ecm_found[False]))


if __name__ == "__main__":
    parser = get_arg_parser()
    args = parser.parse_args()
    gap_utils.verify_args(args)

    double_check(args)
