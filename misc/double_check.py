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

import sys
# Hack to allow import of gap_utils
sys.path.append(".")

import gmpy2

import gap_utils



def get_arg_parser():
    parser = argparse.ArgumentParser('Double check the results of combined_sieve')

    parser.add_argument('--ecm', type=str, default="ecm",
        help="ecm tool")

    parser.add_argument('--unknown-filename', type=str, required=True,
        help="determine mstart, minc, p, d, sieve-length, and sieve-range"
             " from unknown-results filename")

    parser.add_argument('-c', '--count', type=int, default=20,
        help="count of non-trivial numbers to verify per m")

    return parser


def product(l):
    r = 1
    for x in l:
        r *= x
    return r

def double_check(args):
    P = args.p
    D = args.d

    small_primes = [p for p in range(2, 100) if gmpy2.is_prime(p)]
    assert len(small_primes) == 25

    K, r = divmod(gmpy2.primorial(P), D)
    assert r == 0, (P, D)

    M_start = args.mstart
    M_inc = args.minc

    SL = sieve_length = args.sieve_length
    sieve_range = args.sieve_range

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
        if math.gcd(m, D) != 1: continue

        # Read a line from the file
        line = unknown_file.readline()
        start, u_l, u_h = line.split("|")

        match = re.match(r"^([0-9]+) : -([0-9]+) \+([0-9]+)", start)
        assert match, (len(start), len(u_l), len(u_h))
        mtest, unknown_l, unknown_u = map(int, match.groups())
        assert mtest == mi

        low = list(map(int,u_l.strip().split(" ")))
        high = list(map(int,u_h.strip().split(" ")))

        unknowns = {i for i in low + high}

        # Check trivial composites not included
        z = unknowns.intersection(boring_composites)
        assert not z, (len(z), min(z))

        str_start = f"{m} * {P}# / {D}"
        print(str_start, "\tunknowns {} + {} = {}".format(len(low), len(high), len(unknowns)))

        # Choose some random numbers
        count = args.count
        while count > 0:
            t = random.choice(range(-SL+1, SL))
            if t in boring_composites:
                tested[0] += 1
                continue

            # If sieve found a small factor
            found_small = t not in unknowns

            N = (m * K + t)
            found = False
            for p in small_primes:
                if N % p == 0:
                    print(f"\t\t{t:<+6} had trivial factor of {p} skipping")
                    assert found_small, f"{p} divides {str_start} + {t}"
                    found = True
                    break
            if found:
                continue

            tested[1 + found_small] += 1
            word = "a" if found_small else "no"
            print (f"\t{t:<+6}\texpect {word} small composite factor")

            # Interesting factor
            count -= 1

            # Change to subprocess.run(X, capture_output=True) with Python 3.7
            ret, output = subprocess.getstatusoutput("echo '{} + {}' | {} -one -primetest -q 1e4".format(
                str_start, t,
                args.ecm))

            found_factor = (ret & 2) > 0

            print ("\t\t\tecm:", output)
            if found_factor:
                ecm_factor = int(output.split(" ")[0])

                # factor factor to primes, Assumes gmp-lib for larger factor is installed
                factor_output = subprocess.check_output(["factor", str(ecm_factor)])
                factors = [int(f) for f in factor_output.decode().split(":")[1].strip().split(" ")]
                if len(factors) > 1:
                    print("\t\t\tfactor", factor_output)
                assert all(gmpy2.is_prime(factor) for factor in factors)

                assert product(factors) == ecm_factor, (ecm_factor, factors)

                # verify if sieve didn't find a small ecm better not find a small
                if not found_small:
                    assert all(f > args.sieve_length for f in factors), (ecm_factor, factors)

                # Tempting to check if found_small then any(f < args.sieve_length)
                # but ecm not guarenteed to find that factor.

            # generally expected ecm to find the
            ecm_found[found_small] += 1

    print ("{} trivially composite, {} unknowns, {} known composites".format(*tested))
    print ("ecm found {} composites: not known {} composite, {} known".format(sum(ecm_found), *ecm_found))

if __name__ == "__main__":
    parser = get_arg_parser()
    args = parser.parse_args()
    gap_utils.verify_args(args, ".txt")

    double_check(args)

