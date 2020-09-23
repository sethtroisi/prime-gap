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

import json
import pprint
import re
import time

import gmpy2

"""
For missing gaps with prime end points (which would be a record)
Verify that no internal primes were missed (generally there will be missed records)

1.
```
sqlite> select '"' || "BOTH SIDES PRIME: " || m || " * " || p || "# / " || d || " " || -prev_p || " +" || next_p || '",' from m_missing_stats where prev_p > 0 and next_p > 0 ORDER BY M;

"BOTH SIDES PRIME: 1001027 * 9511# / 310310 -13520 +116264",
"BOTH SIDES PRIME: 1002873 * 9511# / 310310 -34720 +96682",
...
```
2. Add these to missing_gap_data.json (under check)
3. time python missing_gap_verify.py

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
        f.write(res.replace("'", '"').strip())


def normalize_line_new_format(line):
    # All success / fails have gap
    data = parse_line(line)
    assert data[-1] not in (None, ""), (line, data)

    new_format = "{:<7} * {:5}# / {:<6} {:<+6} {:+6}".format(*data[:-1])
    # manually gap align the final gap column
    return "{:50} gap {}".format(new_format, data[-1])


def normalize(lines):
    def sort_key(line):
        m,p,d,*_ = parse_line(line)
        return p,d,m

    # Sort by (p,d) then m
    ordered = sorted(lines, key=sort_key)
    return list(map(normalize_line_new_format, ordered))


def parse_line(line):
    LINE_RE = re.compile(r"(\d+)\s*\*\s*(\d+)\#\s*\/\s*(\d+)\s*(-\d+)\s*(\+\d+)"
                          "(?:\s*gap (\d+))?")

    match = LINE_RE.search(line)
    assert match, line
    groups = match.groups()

    gap = groups[-1] or ""  # Replace None with ""
    groups = groups[:-1]

    m, p, d, l, h = map(int, groups)
    return (m, p, d, l, h, gap)


def check_check(test):
    m, p, d, l, h, gap = parse_line(test)
    assert gap == "", test

    N = m * gmpy2.primorial(p) // d
    assert l < 0 and h > 0, test
    low = N + l
    high = N + h

    print (f"Testing {m} * {p}# / {d} ({l}, {h})")

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

    success = z == high
    found_gap = z - low
    print ("\t next_prime {}, {}   {:.1f} seconds".format(
        success, found_gap, t2 - t1))

    return success, found_gap, test


def filter_checks(checks, valid, fails):
    processed = {(m,p,d) for m, p, d, *_ in map(parse_line, valid + fails)}
    already_processed = 0

    for test in checks:
        test = test.strip()
        if test == "":
            continue

        m, p, d, l, h, gap = parse_line(test)
        assert gap == "", test

        if (m,p,d) in processed:
            already_processed += 1
            continue
        processed.add((m,p,d))

        yield test

    if already_processed:
        print(f"\n\t{already_processed} checks were already processed")


def test_records():
    checks, valid, fails = load_data()
    print("checks: {}, valid: {}, fails: {}\n".format(
        len(checks), len(valid), len(fails)))

    updates = []
    for test_line in filter_checks(checks, valid, fails):
        success, found_gap, test = check_check(test_line)
        assert test == test_line

        update = normalize_line_new_format(test + " gap " + str(found_gap))
        updates.append(update)
        if success:
            valid.append(update)
            # Double print with lots of space for improved visibility
            print("\n"*3, update, "\n"*2)
        else:
            fails.append(update)

        print(update)

    if updates:
        print ("\n")
        for update in updates:
            print (update)

    save_data(valid, fails)


if __name__ == "__main__":
    test_records()
