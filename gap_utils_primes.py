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

import re
import subprocess
import time
import tempfile

import gmpy2

try:
    # primegapverify is only alpha so don't require it yet
    import primegapverify

    has_pgv = True
except ModuleNotFoundError:
    has_pgv = False


def is_prime(num, str_n, dist):
    # TODO print log of which library is being used.
    if gmpy2.num_digits(num, 2) > 8000:
        return openPFGW_is_prime(str_n + str(dist))

    return gmpy2.is_prime(num)


def openPFGW_is_prime(str_n):
    # XXX: getstatusoutput runs in a shell, does this double the overhead?
    # Overhead of subprocess calls seems to be ~0.03
    s = subprocess.getstatusoutput(f'./pfgw64 -f0 -q"{str_n}"')
    assert s[1].startswith('PFGW'), s
    return s[0] == 0


def openPFGW_ABC(m, K, offsets):
    """
    See pfgw's abcfileformats.txt for details

    Produces a temp file with this content
      ABC $a*1009#/3018+$b // {number_primes,$a,1}
      13 4
      13 12
      13 22
      13 26
      ...

    pfgw then tests till the first prime
    """

    assert len(offsets) < 50000, len(offsets)

    with tempfile.NamedTemporaryFile(mode="w") as abc:
        abc.write(f"ABC $a*{K}+$b // {{number_primes,$a,1}}\n")
        for offset in offsets:
            abc.write(f"{m} {offset}\n")
        abc.flush()

        s = subprocess.getstatusoutput(f"./pfgw64 -f0 {abc.name}")
        assert s[0] in (0, 1), s
        assert s[1].startswith('PFGW'), s

        # TODO how to validate doesn't skip because previous processed?
        # if "failed" in s[1] or "previous" in s[1]:
        #    print(s[1])
        #    assert False

        # Look for "1*1009#/3018+512 is 3-PRP! (0.0021s+0.0000s)"
        if 'is 3-PRP' in s[1]:
            offset = int(re.search(r'\+(-?[0-9]+) is 3-PRP', s[1]).group(1))
            index = offsets.index(offset)
            return offset, index + 1

    return None, len(offsets)


def determine_next_prime(m, str_n, K, unknowns, SL):
    center = m * K
    tests = 0

    for i in unknowns:
        assert i > 0
        tests += 1
        if is_prime(center + i, str_n, i):
            return tests, i

    # next_prime(...) outside of SL
    tests1, next_p = determine_next_prime_large(m, K, SL)
    return tests + tests1, next_p


def determine_prev_prime(m, str_n, K, unknowns, SL, primes, remainder):
    center = m * K
    tests = 0

    for i in unknowns:
        assert i < 0
        tests += 1
        if is_prime(center + i, str_n, i):
            return tests, -i

    # prev_prime(...) outside of SL need more code.
    tests1, prev_p = determine_prev_prime_large(m, str_n, K, SL, primes, remainder)
    return tests + tests1, prev_p


def determine_next_prime_large(m, K, SL):
    # XXX: PFGW fallback

    # XXX: parse to version and verify > 6.2.99
    assert gmpy2.mp_version() == 'GMP 6.2.99', gmpy2.mp_version()

    center = m * K
    # Double checks center + SL.
    next_p = int(gmpy2.next_prime(center + SL) - center)
    return 0, next_p


def determine_prev_prime_large(m, str_n, K, SL, primes, remainder):
    global has_pgv

    tests = 0
    center = m * K
    if has_pgv:
        t0 = time.time()
        # primegapverify sieve a big interval below
        interval = 4 * SL
        top = center - SL
        bottom = top - interval
        assert bottom > 0
        assert interval < 10 ** 7, SL
        # Cover [center - 5 * SL, center - SL]
        # technically last is included multiple times but safer is better.
        composites = primegapverify.sieve(bottom, 4 * SL)
        assert len(composites) == len(range(bottom, top + 1))

        for i, composite in enumerate(reversed(composites)):
            if not composite:
                tests += 1
                if is_prime(center - (SL + i), str_n, -(SL + i)):
                    t1 = time.time()
                    if (t1 - t0) > 60:
                        print("\tprimegapverify prev_prime({}{}) took {:.2f} second "
                              "({} tests, {:.3f}s/test".format(
                            str_n, -SL, t1 - t0, tests, (t1 - t0) / tests))
                    return tests, SL + i

        assert False, ("Huge prev_prime!", str_n, ">", 5 * SL)

    # Double checks center + SL.
    # Very ugly fallback.
    print("Falling back to slow prev_prime({}{})".format(str_n, -SL))
    t0 = time.time()
    tests0 = tests
    for i in range(SL, 5 * SL + 1):
        composite = False
        for prime, remain in zip(primes, remainder):
            modulo = (remain * m) % prime
            if i % prime == modulo:
                composite = True
                break
        if not composite:
            tests += 1
            if is_prime(center - i, str_n, -i):
                t1 = time.time()
                t = t1 - t0
                if t > 60:
                    num_tests = tests - tests0
                    print("\tfallback prev_prime({}{}) took {:.2f} second ({} tests, {:.3f}s/test"
                          .format(str_n, -SL, t, num_tests, t / num_tests))
                return tests, i

    assert False
