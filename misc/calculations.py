#!/usr/bin/env python3
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


def sieve_percent(max_prime):
    GAMMA = 0.577215
    return 1 / (math.log(max_prime) * math.exp(GAMMA))


def prime_chance(ln_n):
    return 1 / ln_n - 1 / ln_n ** 2


def expected_PRP(max_prime, P):
    import gmpy2
    ln_n = float(gmpy2.log(gmpy2.primorial(P)))
    return ln_n, sieve_percent(max_prime) / prime_chance(ln_n)


def expected_PRP_gap(max_prime, gap):
    return gap * sieve_percent(max_prime)


def pgsurround_sieve_limit(ln):
    log_n = ln
    log2_n = ln / math.log(2)
    assert 900 <= log2_n <= 203601, log2_n
    return (0.05 + (log2_n / 8000.0)) * log_n * log_n * math.log(log_n)


def Runtime():
    from math import log, log2
    import sympy
    import gmpy2

    M_m = 0.26149

    K = gmpy2.primorial(8887) // gmpy2.primorial(13)
    log2K = int(gmpy2.log2(K))
    M = 10000
    sl = 10 ** 10
    M_c = 1917
    c = 50
    S = 2 * 10 * (1000 * round(float(gmpy2.log(K)) / 1000)) + 1
    C_mod = 4 / 64

    C_eq3 = 30
    print(f"log2(K) = {log2K}, S = {S}")

    P = 455052511
    Psmall = int(sympy.primepi(S * c))
    Plarge = P - Psmall

    ll_sl = log(log(sl))

    A1 = P * log2K * C_mod
    B1 = S * (ll_sl + M_m)
    traditional = M_c * (A1 + B1)
    print(f"{M_c} \\times ({A1:.2e} + {B1:.2e})")
    print(f"{traditional=:.2e}")
    print()

    print(
        "{:,}*{:,}/{} + {:,}*{:,}({:.2f}+{:.2f}) + {:,}*{:,} + "
        "({:,}*{:,}*({:.2f}-{:.2f}) + {:,})*{}log2({}/{:,})".format(
            P, log2K, int(1 / C_mod),
            M_c, S, ll_sl, M_m,
            M_c, Psmall,
            M, S, ll_sl, log(log(c * S)), Plarge,
            C_eq3, sl, S))

    A = P * log2K * C_mod
    B = M_c * B1
    C = M_c * Psmall
    Da = M * S * (ll_sl - log(log(c * S))) + Plarge
    Db = 30 * log2(sl / S)

    total = A + B + C + Da * Db

    print(f"{A:.2e} + {B:.2e} + {C:.2e} + {Da * Db:.2e}".replace("e+", " \\times 10^"))
    print()
    print(f"{total=:.2e}".replace("e+", " \\times 10^"))

    print()
    print(f"{traditional / total:.0f}x speedup")


def MertensThird():
    for max_prime, primes in (
            (10 ** 4, 1229), (10 ** 5, 9592), (10 ** 6, 78498),
            (10 ** 8, 5761455), (4 * 10 ** 9, 189961812), (10 ** 10, 4118054813),
            (10 ** 11, 4118054813), (10 ** 12, 37607912018), (10 ** 13, 346965536839),
    ):
        power = math.log10(max_prime)
        primes = ("{}" if primes < 1e9 else "{:.1e}").format(primes)
        percent = sieve_percent(max_prime)
        print(f"$10^{{{power:<2}}}$  & {primes:8}\t& {percent:.6f} &\t\t\\\\")


def Speedup():
    import matplotlib.pyplot as plt
    import numpy as np

    # for mi in {100,300,1000,3000,10000,30000,100000,300000,1000000,3000000,10000000}; do \
    #   time ./combined_sieve --save -u 1511_2190_2000000_${mi}_s20000_l100000M.txt --search-db test.db; done
    # time ./large_sieve 73 1511 2190 -15000 30000 100''000''000''000 1>test.txt

    # Moved to colab


def Appendix1():
    # 503, 1009, 1511, 5003, 10007
    for P, mp in [
        (1511, 250e9),
        (1511, 15000e9),
        (2111, 1800e9),
        (4441, 2000e9),
        (5333, 2000e9),
        (8887, 1000e9),
    ]:
        ln, PRPs = expected_PRP(mp, P)
        pg_mp = pgsurround_sieve_limit(ln)
        speedup = sieve_percent(pg_mp) / sieve_percent(mp) - 1

        pg_mp_str = f"{pg_mp:.1e}".replace("e+0", " \\times 10^")

        print(f"{P:5} & {ln:.0f}\t& {mp / 1e9:4.0f} \\times 10^9\t& "
              f"{PRPs:.1f}\t& {pg_mp_str}\t& {speedup:.1%}\t\\\\")


def Trick2():
    # SLOW (~250s in python3, ~70s in pypy3)
    import sympy

    N = sympy.prod(sympy.primerange(2, 503 + 1)) // 210
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
        print(f"{count}/{primes} = {count / primes:.2%}\t", time.time() - t)


# Runtime()
# MertensThird()
# Speedup()
# Appendix1()
# Trick2()
