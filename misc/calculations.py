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

#print (expected_PRP_gap(40 * 10 ** 6, 155056))


def Appendix1():
    for P, mp in [
            (503, 10e9), (1009, 20e9), (1511, 80e9), (1999, 100e9),
            (5003, 200e9), (10007, 400e9), (20011, 500e9), (20011, 4e9)
    ]:
        print (f"{P:5} {mp:.1e}\t", expected_PRP(mp, P))


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

#Appendix1()
#Trick2()
