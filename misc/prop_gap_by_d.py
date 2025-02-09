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

import array
import math

import gmpy2
import matplotlib.pyplot as plt
import numpy as np
import primesieve
import sympy


P = 347
PRIME_LIMIT = 10 ** 7
M_START = 10 ** 9
M_COUNT = 100

assert gmpy2.is_prime(P)
primorial = gmpy2.primorial(P)
ln_N = math.log(primorial) + math.log(M_START)

# Only interestd up to merit X
MAX_MERIT = 10
SL = int(MAX_MERIT * ln_N)
print()
print(f"ln({P}# * <M>) = {int(ln_N + .99)}")
print(f"\t{SL = } for merit {MAX_MERIT}")


other_primes = primesieve.primes(P+1, PRIME_LIMIT+1)
prob_prime = 1 / ln_N - 1 / (ln_N * ln_N)
prob_prime_after = prob_prime * math.log(PRIME_LIMIT) * math.exp(sympy.EulerGamma.evalf())
prob_prime_nth = [prob_prime_after * pow(1 - prob_prime_after, i) for i in range(SL+1)]
prob_greater_nth = 1 - np.cumsum(prob_prime_nth)
print(f"\ttakes ~{1/prob_prime_after:.1f} PRPs to find a prime")
print()

fig, (ax1, ax2) = plt.subplots(1, 2)
for D in [2, 3, 5, 11, 13]:
    assert gmpy2.is_prime(D)
    if D < 40:
        d = gmpy2.primorial(D)
    K = gmpy2.primorial(P) // d

    # Coprime X
    X = [i for i in range(1, SL+1) if gmpy2.gcd(i, K) == 1]

    print(d, len(X), X[:20])

    count0 = [gmpy2.gcd(x, d) == 1 for x in X]
    cum_count = np.cumsum(count0)

    # Probability of prime for each index

    factors_d = sympy.factorint(d)
    assert all(v == 1 for v in factors_d.values())
    factors_d = array.array('Q', factors_d.keys())

    all_primes = [p for p in factors_d + other_primes if p != 2]
    small_p_and_r = [(p, K % p) for p in all_primes if p < SL]
    large_p_and_r = [(p, K % p) for p in all_primes if p > SL]

    sieve_start = [0 for i in range(SL+1)]
    for x in X: sieve_start[x] = 1
    sieve_start = tuple(sieve_start)

    prob_gap = np.zeros(len(X))
    prob_greater = 0

    m = M_START - 1
    for i in range(M_COUNT):
        m += 1
        while math.gcd(m, d) > 1:
            m += 1

        # 1 unknown, 0 composite
        sieve = list(sieve_start)

        center = m * K
        assert center % 2 == 1

        # handle 2
        for z in range((-center) % 2, SL+1, 2):
            sieve[z] = 0

        for p, r in small_p_and_r:
            offset = (-r * m) % p
            first_odd = offset if (offset & 1 == 0) else offset + p
            for offset in range(first_odd, SL+1, 2*p):
                sieve[offset] = 0
        for p, r in large_p_and_r:
            offset = (-r * m) % p
            if offset <= SL:
                sieve[offset] = 0

        seen = 0
        # Find expected gap
        for i, x in enumerate(X):
            if sieve[x]:
                prob_gap[i] += prob_prime_nth[seen]
                seen += 1
        prob_greater = prob_greater_nth[seen]

        #print("\t", d, " ", m, seen, "remaining", prob_greater_nth[seen])

    cdf = np.cumsum(prob_gap)
    ax1.scatter(X, prob_gap / M_COUNT, label=str(d))
    ax2.plot(X, 1 - (cdf / M_COUNT), label=str(d))
    print(f"\tProb gap >= {SL} is {prob_greater / M_COUNT:.4f}")

    # Assume other side is roughly the same

plt.title(f'K={P}#/d')

ax1.set_ylabel("Probability of Gap")
ax1.set_xlabel("Gap Size")

ax2.set_ylabel("Prop of Gap > X")
ax2.set_xlabel("Gap Size")
ax2.set_yscale("log")

plt.legend(loc='upper right')
plt.show()
