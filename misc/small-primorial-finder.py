z = """
2 3 7 23 89 139 199 113 1831 523 887 1129 1669 2477 2971 4297 5591 1327 9551
30593 19333 16141 15683 81463 28229 31907 19609 35617 82073 44293 43331
34061 89689 162143 134513 173359 31397 404597 212701 188029 542603 265621
461717 155921 544279 404851 927869 1100977 360653 604073 396733 1444309
1388483 1098847 2238823 1468277 370261 492113 5845193 1349533 1895359
3117299 6752623 1671781 3851459 5518687 1357201 6958667 6371401 3826019
7621259 10343761 11981443 6034247 2010733 13626257 8421251 4652353
17983717 49269581 33803689 39175217 20285099 83751121 37305713 27915737
38394127 52721113 38089277 39389989 17051707 36271601 79167733 147684137
134065829 142414669 123454691 166726367 70396393 46006769
"""

z = list(map(int, z.strip().split()))
print (z)
print()
import gmpy2
import sympy

for p in sympy.primerange(5, 30):
    prim = int(gmpy2.primorial(p))
    for prv in z:
        nxt = gmpy2.next_prime(prv)
        # check if a multiple of prim is in the middle of (k, n)
        mult = nxt // prim
        if mult > prv // prim and mult < 1000:
            mid = (nxt // prim) * prim
            print (f"{prv:<6}  {nxt-prv:2d}  |  {mid - prv} {nxt - mid}\t{nxt // prim} * {p}# = {mid}")
print("\n" * 3)



#mid = 16170
#X   = 30

#mid = int(8 * gmpy2.primorial(10))
#X   = 15

#mid = 210
#X   = 11

mid = 1680
X   = 15

start = mid - X
stop = mid + X
assert start % 2 == 1
assert stop % 2 == 1

search = 20

labels = {start, mid, stop}
factors = []

# Skip even numbers
for i in range(start, stop+1):
    if i % 2 == 0 and i != mid:
        continue

    for j in range(2, search):
        if i % j == 0:
            factors.append(j)
            break
    else:
        factors.append("?")
        labels.add(i)
        print ("\t", i, sympy.factorint(i))


names = []
for a, b in zip(sorted(labels), sorted(labels)[1:]):
    if a == mid:
        names.append("\\textbf{%d}" % mid)
    else:
        names.append(a)
    if a + 2 < b:
        names.append("\multicolumn{%d}{|c|}{...}" % (len(range(a+2, b, 2))))
names.append(b)

BEGIN_TAB = "\\begin{tabular}{%s}"

print (BEGIN_TAB % ("c" * (stop - start + 1)))
print ("\t", " & ".join(map(str, names)), "\\\\")
print ("\t", " & ".join(map(str, factors)), "\\\\")
print ("\\end{tabular}")
print ()
