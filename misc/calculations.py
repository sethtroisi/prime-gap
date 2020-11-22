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
    logn = ln
    log2n = ln / math.log(2)
    assert 900 <= log2n <= 203601, log2n
    return (0.05 + (log2n/8000.0)) * logn * logn * math.log(logn)


def Runtime():
    from math import log, log2
    import sympy
    import gmpy2

    Mc = 0.26149

    logK = int(gmpy2.log(gmpy2.primorial(20011) // gmpy2.primorial(13)))
    M  = 10000
    sl = 10 ** 10
    Cf = 1/5
    c = 10 / Cf
    S = 2 * 20 * 20000 + 1
    C_mod = 4/64
    C_eq3 = 30

    P = 455052511
    Psmall = int(sympy.primepi(S * c))
    Plarge = P - Psmall

    llsl = log(log(sl))

    A1 = P * logK * C_mod
    B1 = S * (llsl + Mc)
    traditional = M * Cf * (A1 + B1)
    print (f"{M} \\times {Cf} \\times ({A1:.2e} + {B1:.2e})")
    print (f"{traditional:.2e}")
    print ()

    print ("{:,}*{:,}/{} + {:,}/{}*{:,}({:.2f}+{:.2f}) + {:,}/{}*{:,} + ({:,}*{:,}*({:.2f}-{:.2f}) + {:,})*{}log2({}/{:,})".format(
        P, logK, int(1/C_mod),
        M, 1/Cf, S, llsl, Mc,
        M, 1/Cf, Psmall,
        M, S, llsl, log(log(c * S)), Plarge,
        C_eq3, sl, S))

    A = P * logK * C_mod
    B = M * Cf * B1
    C = M * Cf * Psmall
    Da = M * S * (llsl - log(log(c * S))) + Plarge
    Db = 30 * log2(sl/S)

    total = A + B + C + Da * Db

    print (f"{A:.2e} + {B:.2e} + {C:.2e} + {Da*Db:.2e}")
    print (f"{total:.2e}")

    print ()
    print (f"{traditional/total:.0f}x speedup")


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


def Speedup():
    import matplotlib.pyplot as plt
    import numpy as np

    # for mi in {100,300,1000,3000,10000,30000,100000,300000,1000000}; do ./combined_sieve --save-unknowns --unknown-filename 1511_2190_2000000_${mi}_s20000_l100000M.txt; done
    # time ./large_sieve 73 1511 2190 -15000 30000 100''000''000''000 1>test.txt

    M_single = 224.6
    M_cost = {
        (1,1): M_single,
        (100, 26): 931,
        (300, 79): 945,
        (1000, 262): 934,
        (3000, 789): 944,
        (10000, 2629): 950,
        (30000, 7891): 1054.6,
        (100000, 26301): 1184.3,
        (300000, 78905): 1776.1,
        (1000000, 263013): 3966,
        (3000000, 789041): 10550,
        (10000000, 2630136): 34814,
    }

    line = "{:>10}\t& {:>7}\t& {:5.0f}\t& {:.4g} \t& {:.0f} \t&\\"

    first = line.format("M", "coprime", 0.0, 1.1, 2.2)
    first = (first.replace("0.0", "time(s)")
                  .replace("1.1", "time(s)/coprime")
                  .replace("2.2", "speedup (vs sequential sieve)")
                  .replace("\\", "\rule{0pt}{1em}\\"))

    for (m, coprime), time in M_cost.items():
      print(line.format(m, coprime, time, time / coprime, M_single / (time / coprime)))

    fig = plt.figure(figsize=(32/3,9/3), dpi=200)

    # M = coprime count.
    M, time = zip(*[(coprime,time) for (_, coprime), time in M_cost.items()])
    M = np.array(M)

    ax1 = fig.gca()
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    p1 = ax1.plot(M, time,                   marker="o", label="combined sieve time(s)")
    p2 = ax2.plot(M, M_single * M / time,    marker="o", label="speedup", color="C1")
    plt.grid(which="major", axis="y")
    plt.xscale("log")
    ax1.set_yscale("log")
    ax2.set_yscale("log")

    ax1.set_xlabel("Parallel M")

    ax1.set_ylabel("time(s)")
    ax2.set_ylabel("speedup")

    ax2.legend(p1 + p2, ["time(s)", "speedup"], loc="upper left", framealpha=0.9)
    plt.tight_layout()
    plt.show()


def Appendix1():
    # 503, 1009, 1511, 5003, 10007
    for P, mp in [
            (1511, 250e9),
            (2111, 800e9),
            (4441, 400e9),
            (4441, 100e9),
            (5333, 5000e9),
            (6663,  500e9),
            (8887, 1000e9),
    ]:
        ln, prps = expected_PRP(mp, P)
        pg_mp = pgsurround_sieve_limit(ln)
        speedup = sieve_percent(pg_mp) / sieve_percent(mp) - 1

        print (f"{P:5} & {ln:.0f}\t& {mp/1e9:4.0f}e9\t& {prps:.1f}\t& {pg_mp:.1e}\t& {speedup:.1%}\t\\\\")


def Trick2():
    # SLOW (~250s in python3, ~70s in pypy3)
    import sympy

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

#Runtime()
#ParameterSelection()
#Speedup()
#Appendix1()
#Trick2()
