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

def modulo_search(p, a, l, r):
    if l == 0: return 0

    if 2*a > p:
        return modulo_search(p, p - a, p - r, p - l)

    if a <= r:
#        mult = (l - 1) // a + 1;
#        test = mult * a
#        if test <= r:
#            return mult

        mult = (l - 1) // a;
        test = mult * a
        assert mult < l
        if a <= r - test:
            return mult + 1

    print ("\t\tdelta:", r - l)

    new_a = a - (p % a)
    assert 0 <= new_a < a, (a, new_a)
    k = modulo_search(a, new_a, l % a, r % a)

    tl = k * p + l
    mult = (tl - 1) // a + 1
    return mult


import sympy
primes = list(sympy.primerange(2, 10 ** 5))

s = 0
for p in primes[-100:]:
    s += modulo_search(p, 123, 1000, 1010)
    print ("\tprime:", p, "search:", s)

print ()
print (s, "vs 780012")
