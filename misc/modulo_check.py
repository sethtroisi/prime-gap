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

    #assert l < r, (p, a, l, r)
    delta = abs(r - l)

    '''
    l_div, l_mod = divmod(l, a)
    r_div, r_mod = divmod(r, a)
    if l_div < r_div or l_mod == 0:
        assert l_mod == 0 or l_mod + delta >= a, (l_mod + delta, a)
        return l_div + (l_mod > 0)
    assert r_mod == l_mod + delta
    assert r_div == l_div
    '''
    l_div, l_mod = divmod(l - 1, a)
    l_mod += 1
    r_mod = l_mod + delta
    if r_mod >= a:
        return l_div + 1

    if 2*a > p:
        return modulo_search(p, p - a, p - r, p - l)

    new_a = a - (p % a)
    assert 0 <= new_a < a, (a, new_a)
    k = modulo_search(a, new_a, l_mod, r_mod)

    tl = k * p + l
    mult = (tl - 1) // a + 1
    return mult

def unwind(k, stack):
    for p, a, l in reversed(stack):
        #assert k < p
        #assert a < p
        #assert k < a
        #assert l > 0  # Helps guarentee r * k + l - 1 >= 0

        # XXX: these are probably already computed (for new_a)
        div, rem = divmod(p, a)
        k_1 = k * div + 1

        # No larger than a^2 + l
        k_2 = (rem * k + l - 1) // a
        k = k_1 + k_2

        #k = (k * p + l - 1) // a + 1
    return k

def modulo_search_stack(p, a, l, r, stack):
    if l == 0:
        return unwind(0, stack)

    if 2*a > p:
        a, l, r = p - a, p - r, p - l

    l_div, l_mod = divmod(l, a)
    r_div, r_mod = divmod(r, a)
    if l_div < r_div or l_mod == 0:
        return unwind(l_div + (l_mod > 0), stack)

    new_a = a - (p % a)
    assert 0 <= new_a < a, (a, new_a)
    # k = modulo_search(a, new_a, l % a, r % a, stack)
    # tl = k * p + l
    # return (tl - 1) // a + 1

    stack.append((p, a, l))
    return modulo_search_stack(a, new_a, l_mod, r_mod, stack)


def main():
    # Test case where m*K = L first
    z = modulo_search(10**9+7, 123, 5*123, 5*123+6)
    assert z == 5, z
    # Test case where m*K = R first
    z = modulo_search(10**9+7, 123, 5*123-6, 5*123)
    assert z == 5, z


    import sympy
    s = 0
    for p in list(sympy.primerange(10 ** 4, 10 ** 5))[-100:]:
        s1 = modulo_search(p, 123, 1000, 1010)
        s2 = modulo_search_stack(p, 123, 1000, 1010, [])
        if s1 != s2: print(f"\tError with {p} | {s1} != {s2}")
        s += s1
    print()
    print(s, "vs 780012")

    '''
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
    '''



if __name__ == "__main__":
    main()
