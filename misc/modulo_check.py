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
    assert 0 <= new_a < a
    k = modulo_search(a, new_a, l % a, r % a)

    tl = k * p + l
    mult = (tl - 1) // a + 1
    return mult


import MathLib
primes = MathLib.sieveOfErat(10 ** 5)

s = 0
for p in primes[-100:]:
    s += modulo_search(p, 123, 1000, 1010)
    print ("\tprime:", p, "search:", s)

print ()
print (s, "vs 780012")
