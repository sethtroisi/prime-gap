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

import math
import array


def calc_valid_mi(ms, mi, d):
    '''
    Equivilant to [mi for mi in range(M_inc) if math.gcd(ms + mi, D) == 1]
    '''
    valid_mi = array.array('q')  # signed long long AKA int64_t

    # Avoids allocating 100M+ byte array for valid
    interval_size = 1000000
    for start in range(0, mi, interval_size):
        end = min(mi, start + interval_size)
        size = end - start

        valid = array.array('b', (1 for _ in range(start, end)))
        factors_d = _factor_simple(d)
        for p in factors_d:
            first = -(ms + start) % p
            assert (ms + start + first) % p == 0
            for i in range(first, size, p):
                valid[i] = 0

        for mii, v in enumerate(valid, start):
            if v:
                valid_mi.append(mii)
    return valid_mi


def count_num_m(ms, mi, d):
    if d == 1:
        return mi

    if ms + mi < 10000:
        return sum(1 for m in range(ms, ms+mi) if math.gcd(m, d) == 1)

    factors_d = _factor_simple(d)
    return (_r_count_num_m(ms + mi - 1, factors_d, len(factors_d)-1) -
            _r_count_num_m(ms - 1,      factors_d, len(factors_d)-1))


def _factor_simple(d):
    # Most d are a primorial so this is quite quick
    factors = []
    for p in range(2, int(math.sqrt(d)+2)):
        if p*p > d:
            break
        while d % p == 0:
            factors.append(p)
            d //= p
    if d > 1:
        factors.append(d)
    return factors


def _r_count_num_m(n, factors_d, i):
    """Count of numbers coprime to d less than end; sum( gcd(m, d) == 1 for m in range(n, n+i) )

    Uses inclusion exclusion on prime factorization of d
    """

    if n == 0:
        return 0

    if i < 0:
        return n

    return _r_count_num_m(n, factors_d, i-1) - _r_count_num_m(n // factors_d[i], factors_d, i-1)
