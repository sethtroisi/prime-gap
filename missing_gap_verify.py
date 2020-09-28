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
import multiprocessing
import pprint
import re
import time

import gmpy2
import sqlite3

"""
For missing gaps with prime end points (which would be a record)
Verify that no internal (next) primes were missed (generally there will be missed records)
"""

SEARCH_GAPS_DB_FN           = "prime-gap-search.db"
MISSING_GAP_BOTH_PRIME_FILE = "missing_gap_data.json"


def load_data():
    with open(MISSING_GAP_BOTH_PRIME_FILE) as f:
        data = json.loads(f.read().replace("'", '"'))

    return data['valid'], data['fails']


def get_db():
    return sqlite3.connect(SEARCH_GAPS_DB_FN)


def load_sql():
    with get_db() as conn:
        conn.row_factory = sqlite3.Row

        rv = conn.execute(
            "SELECT m, P, D, prev_p, next_p, prp_prev, prp_next "
            "FROM m_stats "
            "WHERE prev_p > 0 AND next_p < 0")

        return tuple(dict(row) for row in rv)


def save_record(m, p, d, prev_p, next_p, test_time):
    with get_db() as conn:
        cur = conn.cursor()
        cur.execute(
            "UPDATE m_stats SET "
            "   next_p = ?, test_time = test_time + ?"
            "WHERE m=? AND P=? AND D=? AND prev_p=?",
            (next_p, test_time, m, p, d, prev_p))

        assert cur.rowcount == 1, (m, p, d, prev_p, next_p)
        conn.commit()


def save_data(valid, fails):
    with open(MISSING_GAP_BOTH_PRIME_FILE, "w") as f:
        res = pprint.pformat({
            'valid': normalize(valid),
            'fails': normalize(fails),
        })
        # HACK: fix quotes so json can load the file.
        output = res.replace("'", '"')
        # Add trailing newline so hand editing doesn't leave diff
        f.write(output.strip() + "\n")


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
    m, p, d = test['m'], test['P'], test['D']
    N = m * gmpy2.primorial(p) // d

    prev_p, next_p = test['prev_p'], test['next_p']

    assert prev_p > 0, test
    assert next_p < 0, test

    # set both positive to help out
    prev_p = abs(prev_p)
    next_p = abs(next_p)

    short_num = f"{m} * {p}# / {d} -{prev_p:<5}"
    print (f"Testing {short_num} to +{next_p}")

    low  = N - prev_p
    high = N + next_p

    t0 = time.time()
    assert gmpy2.is_prime(low)
    assert gmpy2.is_prime(high)
    t1 = time.time()

    print ("\tverified endpoints {:.2f} seconds".format(t1 - t0))

    found_prime = gmpy2.next_prime(N)
    gap = found_prime - low

    success = (gap == next_p + prev_p)
    if success:
        # Also verify prev_prime
        temp = gmpy2.next_prime(low)
        assert temp == found_prime, (prev_p, next_p, "|", gap, temp - low)

    t2 = time.time()

    print ("\nnext_prime({}) = {} {} {}\t{:.1f} seconds\n".format(
        short_num, gap, "=" if success else "!=", high - low, t2 - t1))

    save_record(m, p, d, prev_p, next_p, t2 - t0)

    update = normalize_line_new_format(f"{short_num} +{next_p} gap {gap}")
    return success, update, gap, t2 - t0


def verify_checks(checks, valid, fails):
    processed = {(m,p,d) for m, p, d, *_ in map(parse_line, valid + fails)}
    for row in checks:
        assert (row['m'], row['P'], row['D']) not in processed


def test_records():
    valid, fails = load_data()

    checks = load_sql()

    # Verify none have been previously processed
    verify_checks(checks, valid, fails)

    print("checks: {}, valid: {}, fails: {}\n".format(
        len(checks), len(valid), len(fails)))

    updates = []

    # XXX: make this configurable
    with multiprocessing.Pool(multiprocessing.cpu_count() // 4) as pool:
        # XXX: should write back to DB even in filtered case.
        for success, update, found_gap, time in pool.imap(check_check, checks):
            updates.append(update)
            if success:
                valid.append(update)
                # Double print with lots of space for improved visibility
                print("\n"*3, update, "\n"*2)
            else:
                fails.append(update)


    if updates:
        print ("\n")
        for update in updates:
            print (update)

    save_data(valid, fails)


if __name__ == "__main__":
    test_records()
