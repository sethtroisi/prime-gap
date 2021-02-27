# Copyright 2021 Seth Troisi
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

def modulo_search(p, a, l, r):
    if l == 0:
        return 0

    if 2*a > p:
        return modulo_search(p, p - a, p - r, p - l)

    if a <= r:
        # mult = (l - 1) // a + 1;
        # test = mult * a
        # if test <= r:
        #   return mult

        mult = (l - 1) // a
        test = mult * a
        assert mult < l
        if a <= r - test:
            return mult + 1

    new_a = a - (p % a)
    assert 0 <= new_a < a, (a, new_a)
    k = modulo_search(a, new_a, l % a, r % a)

    tl = k * p + l
    mult = (tl - 1) // a + 1
    return mult

def unwind(k, stack):
    for p, a, l in reversed(stack):
        k = (k * p + l - 1) // a + 1
    return k

def modulo_search_stack(p, a, l, r, stack):
    if l == 0:
        return unwind(0, stack)

    if 2*a > p:
        a, l, r = p - a, p - r, p - l

    if a <= r:
        # mult = (l - 1) // a + 1;
        # test = mult * a
        # if test <= r:
        #   return mult

        mult = (l - 1) // a
        test = mult * a
        assert mult < l
        if a <= r - test:
            return unwind(mult + 1, stack)

    new_a = a - (p % a)
    assert 0 <= new_a < a, (a, new_a)
    # k = modulo_search(a, new_a, l % a, r % a, stack)
    # tl = k * p + l
    # return (tl - 1) // a + 1

    stack.append((p, a, l))
    return modulo_search_stack(a, new_a, l % a, r % a, stack)



import sympy
s = 0
for p in list(sympy.primerange(10 ** 4, 10 ** 5))[-100:]:
    s += modulo_search(p, 123, 1000, 1010)
print()
print(s, "vs 780012")


import time
primes = list(sympy.primerange(10 ** 4, 2 * 10 ** 6))
s = 0
t0 = time.time()
for p in primes:
    s += modulo_search(p, 123, 1000, 1010)

t1 = time.time()
for p in primes:
    s -= modulo_search_stack(p, 123, 1000, 1010, [])

t2 = time.time()

print (f"diff: {s} | {t1 - t0:.3f} vs {t2 - t1:.3f}")
