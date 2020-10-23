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

import math
import time

import gmpy2
import sympy

def next_prime(start):
    return gmpy2.next_prime(start)


def sieve_percent(max_prime):
    GAMMA = 0.577215
    return 1 / (math.log(max_prime) * math.exp(GAMMA))


def prime_chance(ln_n):
    return 1 / ln_n - 1 / ln_n ** 2


def expected_PRP(max_prime, P):
    ln_n = float(gmpy2.log(gmpy2.primorial(P)))
    return ln_n, sieve_percent(max_prime) / prime_chance(ln_n)


def expected_PRP_gap(max_prime, gap):
    return gap * sieve_percent(max_prime)


def pgsurround_sieve_limit(ln):
    logn = ln
    log2n = ln / math.log(2)
    assert 900 <= log2n <= 203601, log2n
    return (0.05 + (log2n/8000.0)) * logn * logn * math.log(logn)

def ParameterSelection():
    for max_prime, primes in (
        (10**4, 1229), (10**5, 9592), (10**6, 78498),
        (10**8, 5761455), (4*10**9, 189961812), (10**10, 4118054813),
        (10**11, 4118054813), (10**12, 37607912018),
    ):

        power = math.log10(max_prime)
        primes = ("{}" if primes < 1e9 else "{:.1e}").format(primes)
        percent = sieve_percent(max_prime)
        print(f"$10^{{{power:<2}}}$  & {primes:8}\t& {percent:.6f} &\t\t\\\\")



def Appendix1():
    # 503, 1009, 1511, 5003, 10007
    for P, mp in [
            (1511, 250e9),
            (2111, 800e9),
            (4441, 400e9),
            (4441, 100e9),
            (5333, 5000e9),
            (8887, 1000e9),
    ]:
        ln, prps = expected_PRP(mp, P)
        pg_mp = pgsurround_sieve_limit(ln)
        speedup = sieve_percent(pg_mp) / sieve_percent(mp) - 1

        print (f"{P:5} & {ln:.0f}\t& {mp/1e9:4.0f}e9\t& {prps:.1f}\t& {pg_mp:.1e}\t& {speedup:.1%}\t\\\\")


def Trick2():
    # SLOW (~250s in python3, ~70s in pypy3)

    N = sympy.prod(sympy.primerange(2, 503+1)) // 210
    # N = gmpy2.primorial(503) // 210

    X = 100000

    first = N - X

    for low, high in [(200000, 10 ** 6), (10 ** 6, 1 * 10 ** 9)]:
        primes = 0
        count = 0
        t = time.time()
        for p in sympy.sieve.primerange(low, high):
            primes += 1
            count += ((first % p) + 2 * X) >= p
        print (f"{count}/{primes} = {count / primes:.2%}\t", time.time() - t)

#ParameterSelection()
Appendix1()
#Trick2()
