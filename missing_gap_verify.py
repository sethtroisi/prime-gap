import json
import re
import time

import gmpy2

MISSING_GAP_BOTH_PRIME_FILE = "missing_gap_data.json"

with open(MISSING_GAP_BOTH_PRIME_FILE) as fn:
    data = json.load(fn)

success = data['success']
check   = data['check']
failed  = data['failed']

print("success: {}, checks: {}, failed: {}\n".format(
    len(success), len(check), len(failed)))

updates = []
for test in check:
    test = test.strip()
    if test == "":
        continue

    m, p, d, l, h = map(int, re.search(r"(\d+)\s?\*\s?(\d+)\#\s?\/\s*(\d+) (\d+) (\d+)", test).groups())
    N = m * gmpy2.primorial(p) // d
    low = N - l
    high = N + h

    print (f"Testing {m}*{p}#/{d} (-{l}, +{h})")

    t0 = time.time()
    assert gmpy2.is_prime(low)
    assert gmpy2.is_prime(high)
    t1 = time.time()

    print ("\tverified endpoints {:.2f} seconds".format(t1-t0))

    # TODO because many values on the low side are not-prime prev_prime is likely to be faster.
    # TODO if prev_prime was available, check next_prime(N) - N = h then check N - prev_prime(N) = l

    continue

    z = gmpy2.next_prime(low)
    t2 = time.time()

    print ("\t next_prime {}, {}   {:.1f} seconds".format(
        z == high, z - low, t2 - t1))

    update = f"\t{test}\t => gap = {z - low}"
    updates.append(update)
    print(update)
    if z == high:
        # Double print with lots of space for improved visibility
        print("\n"*3, update, "\n"*2)

if updates:
    print ("\n")
    for update in updates:
        print (update)
