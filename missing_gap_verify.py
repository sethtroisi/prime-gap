import json
import pprint
import re
import time

import gmpy2

"""
For missing gaps with prime end points (which would be a record)
Verify that no internal primes were missed (generally there will be missed records)

sqlite> select '"' || "BOTH SIDES PRIME: " || m || " * " || p || "# / " || d || " " || prev_p || " " || next_p || " " || (next_p + prev_p) || '",' from m_missing_stats where prev_p > 0 and next_p > 0 ORDER BY M;

"BOTH SIDES PRIME: 1001027 * 9511# / 310310 13520 116264 129784",
"BOTH SIDES PRIME: 1002873 * 9511# / 310310 34720 96682 131402",
...

Add these to missing_gap_data.json (under check)

"""


MISSING_GAP_BOTH_PRIME_FILE = "missing_gap_data.json"


def load_data():
    with open(MISSING_GAP_BOTH_PRIME_FILE) as f:
        data = json.loads(f.read().replace("'", '"'))

    return data['check'], data['valid'], data['fails']


def save_data(valid, fails):
    with open(MISSING_GAP_BOTH_PRIME_FILE, "w") as f:
        res = pprint.pformat({
            'check': [],
            'valid': normalize(valid),
            'fails': normalize(fails),
        })
        # HACK: fix quotes so json can load the file.
        f.write(res.replace("'", '"'))

def normalize(lines):
    def sort_key(line):
        m,p,d,*_ = parse_line(line)
        return p,d,m

    def change_to_new(line):
        # All success / fails have gap
        data = parse_line(line)
        assert data[-1] not in (None, ""), (line, data)

        new_format = "{:<7} * {:5}# / {:<6} {:<+6} {:+6}".format(*data[:-1])
        # manually gap align the final gap column
        return "{:50} gap {}".format(new_format, data[-1])

    # Sort by (p,d) then m
    ordered = sorted(lines, key=sort_key)
    return list(map(change_to_new, ordered))


def parse_line(line):
    LINE_RE = re.compile(r"(\d+)\s*\*\s*(\d+)\#\s*\/\s*(\d+)\s*(-\d+)\s*(\+\d+)"
                          "(?:\s*gap (\d+))?")

    match = LINE_RE.search(line)
    assert match, line
    groups = match.groups()

    gap = groups[-1] or "" # Replace None with ""
    groups = groups[:-1]

    m, p, d, l, h = map(int, groups)
    return (m, p, d, l, h, gap)


def test_records():
    check, valid, fails = load_data()

    print("checks: {}, valid: {}, fails: {}\n".format(
        len(check), len(valid), len(fails)))

    updates = []
    for test in check:
        test = test.strip()
        if test == "":
            continue

        m, p, d, l, h, gap = parse_line(test)
        assert gap == "", test

        N = m * gmpy2.primorial(p) // d
        low = N - l
        high = N + h

        print (f"Testing {m}*{p}#/{d} (-{l}, +{h})")

        t0 = time.time()
        assert gmpy2.is_prime(low)
        assert gmpy2.is_prime(high)
        t1 = time.time()

        print ("\tverified endpoints {:.2f} seconds".format(t1-t0))

        # XXX: because many low values checked prev_prime(high) is likely to be faster.
        # if prev_prime was available, check next_prime(N) - N = h
        # then check N - prev_prime(N) = l

        z = gmpy2.next_prime(low)
        t2 = time.time()

        print ("\t next_prime {}, {}   {:.1f} seconds".format(
            z == high, z - low, t2 - t1))

        update = f"\t{test}\t => gap = {z - low}"
        updates.append(update)
        print(update)
        if z == high:
            valid.append(update.replace("\t", "    "))
            # Double print with lots of space for improved visibility
            print("\n"*3, update, "\n"*2)
        else:
            fails.append(update.replace("\t", "    "))

    if updates:
        print ("\n")
        for update in updates:
            print (update)

    save_data(valid, fails)


if __name__ == "__main__":
    test_records()
