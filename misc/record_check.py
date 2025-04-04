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

"""Search logs directory for any record sized gaps

Replaces some awk hackery.
"""

import argparse
import glob
import os.path
import re
import sqlite3
import sys
import time
from collections import defaultdict

import gmpy2
import primegapverify


def get_arg_parser():
    parser = argparse.ArgumentParser('Search prime gaps logs')

    parser.add_argument(
        '--logs-directory', type=str,
        default=None,
        help="directory of logs")

    parser.add_argument(
        '--log-files', type=str, nargs='*',
        default=[],
        help="individual log files (e.g. logs/2221*)")

    parser.add_argument(
        '--search-db', type=str,
        default="prime-gap-search.db",
        help="Prime database from gap_test (or \"ignore\" to skip)")

    parser.add_argument(
        '--prime-gaps-db', type=str,
        default="gaps.db",
        help="Prime gap database see github.com/primegap-list-project/prime-gap-list")

    parser.add_argument(
        '--whoami', type=str,
        default="S.Troisi",
        help="Name in gap database to ignore for already submitted records")

    parser.add_argument(
        '--ignore-small', nargs='?',
        const=0, default=8, type=int,
        help="Ignore 'records' with small merit (default: 8 or passed arg)")

    parser.add_argument(
        '-s', '--sort', action="store_true",
        help="Sort new records by gapsize")

    parser.add_argument(
        '-v', '--verify', action="store_true",
        help="Verify gaps")

    return parser


def describe_found_gaps(gaps, verify):
    gaps_by_size = sorted(g[0] for g in gaps)
    gaps_by_merit = defaultdict(list)
    for gap in gaps:
        gaps_by_merit[int(gap[1])].append(gap)

    print("Found {} gaps  ({:<6} to {:>6})".format(
        len(gaps_by_size), gaps_by_size[0], gaps_by_size[-1]))
    print("      {} merit ({:.3f} to {:.3f})".format(
        " " * len(str(len(gaps_by_size))),
        min(g[1] for g in gaps), max(g[1] for g in gaps)))
    for int_merit, items in sorted(gaps_by_merit.items()):
        print("    Merit {:<2} x{:<4} | {}".format(
            int_merit,
            len(items),
            ", ".join([str(rec[0]) for rec in sorted(items)[:2]] +
                      ["..."] * (len(items) > 2))
        ))
    print()

    if verify:
        T0 = time.time()
        print(f"Verifying Gaps")
        for gap, merit, _, start, _ in gaps:
            n = primegapverify.parse(start)
            assert gmpy2.is_prime(n), (start)
            assert gmpy2.is_prime(n + gap), (start, gap)
            assert gmpy2.next_prime(n+1) == n + gap
            assert gmpy2.prev_prime(n+1) == n
        print(f"Verified Gaps ({time.time() - T0:.2f}s)")


def record_line_sort(line):
    return tuple(float(v) for v in line.split()[:3])


def print_record_gaps(args, gaps):
    with sqlite3.connect(args.prime_gaps_db) as conn:
        num_gaps = conn.execute('SELECT COUNT(*) FROM gaps').fetchone()[0]
        assert num_gaps > 50000, num_gaps
        min_merit, max_merit = conn.execute(
            'SELECT min(merit), max(merit) FROM gaps').fetchone()
        print("Records, {} gaps merit {:.2f} to {:.2f}".format(
            num_gaps, min_merit, max_merit))

        print("Checking", len(gaps), "gaps against records")
        small_merit = 0
        own_records = []
        record_lines = []
        improvements = defaultdict(float)
        for gap in gaps:
            # gapsize, merit, raw_data, startprime, "line"
            size = gap[0]
            new_merit = gap[1]
            raw_data = gap[2]
            startprime = gap[3]

            # See https://primegap-list-project.github.io/lists/prime-gaps-low-watermarks/
            if size <= 30000 and new_merit < 21.9:
                continue
            if size <= 50000 and new_merit < 19.8:
                continue
            if size <= 75000 and new_merit < 15:
                continue
            if size <= 100000 and new_merit < 10.7:
                continue

            existing = conn.execute(
                'SELECT merit,primedigits,startprime,discoverer FROM gaps WHERE'
                ' gapsize=?', (size,)).fetchone()

            if not existing:
                if new_merit < args.ignore_small and size < 1000000:
                    small_merit += 1
                    continue

                record_lines.append(raw_data)
                #print("\tRecord {:5} | {:70s} | Gap={:<6} (New!)".format(
                #    len(record_lines), raw_data.strip(), size))
                continue

            # Works most of the time, could have false positives
            is_same = existing[2].replace(" ", "") in startprime.replace(" ", "")
            is_own_record = existing[3] == args.whoami

            if is_same and not is_own_record:
                print("\tREDISCOVERED | {:70s} (old: {})".format(raw_data, existing))
                continue

            # If obvious not an improvement don't call parse(...)
            if existing[0] > new_merit + 0.04:
                improvement = new_merit - existing[0] + 6e-3
            else:
                old_n = primegapverify.parse(existing[2])
                if old_n:
                    # Strip to just the number section
                    new_n = primegapverify.parse(startprime)
                    assert new_n, new_n
                    improvement = float(size / gmpy2.log(new_n) - size / gmpy2.log(old_n))
                else:
                    improvement = new_merit - existing[0] + 6e-3

            if improvement >= 0:
                # if not is_same and improvement < 6e-3:
                #    print("Close:", existing[2], "vs newer", startprime)

                if is_same and is_own_record:
                    own_records.append(raw_data)
                    ith = len(own_records)
                    special = ith in (1, 2, 5, 10, 20, 50) or ith % 100 == 0
                    if not special:
                        continue
                else:
                    improvements[gap[0]] = max(improvements[gap[0]], improvement)
                    record_lines.append(raw_data)

                #print("\tRecord {:5} | {:70s} | Gap={:<6} (old: {:.2f}{} +{:.2f})".format(
                #    str(len(own_records)) + "*" if is_same else len(record_lines), gap[4], size,
                #    existing[0], " by you" * is_own_record, new_merit - existing[0]))

        if record_lines:
            if args.sort:
                record_lines.sort(key=record_line_sort)

            seen = set()
            for line in record_lines:
                gap = line.split()[0]
                if gap in seen:
                    continue
                seen.add(gap)
                if len(seen) % 50 == 1:
                    # Seperator before 0 and *after* 50 unique items, before the 51st
                    print()
                print(len(seen), "\t", line)

            print()
            print("Records {} unique {} {}".format(
                len(record_lines), len(seen),
                f"({len(own_records)} already submitted)" if own_records else ""))
            if improvements:
                print("Merit  : +{:.2f} total, +{:.2f} for {}".format(
                    sum(improvements.values()), *max((v, k) for k, v in improvements.items())))
            print("Smallst:", min(record_lines, key=record_line_sort))
            print("Largest:", max(record_lines, key=record_line_sort))
            if small_merit:
                print(f'\tHid {small_merit} new "records" with merit < {args.ignore_small}')
            print()


def search_logs(log_files):
    # 32280  10.9749  5641 * 3001#/2310 -18514 to +13766
    record_format = re.compile(
        r"(\d+)\s+(\d+\.\d+)\s+"
        r"(\d+)\s*(\*)\s*(\d+#)(/\d+#?)\s+"
        r"(-\d+)")

    assert any(map(os.path.isfile, log_files)), "Log files not present: " + ",".join(log_files)

    gaps = []
    seen = set()
    match_print = ""
    for log_fn in log_files:
        with open(log_fn, "r") as f:
            lines = f.readlines()

        print("Processing {:5d} lines from log: {}".format(
            len(lines), log_fn))

        file_match = 0
        for li, line in enumerate(lines):
            if ' ' not in line:
                continue

            if line in seen:
                # Ignore completely duplicate lines
                continue

            match = record_format.search(line)
            if match:
                file_match += 1
                partial_line = line.split(':')[-1].strip()
                match_print = f"    Match {file_match} at line {li}: {partial_line}"
                if file_match < 6:
                    print (match_print)

                seen.add(line)

                gaps.append([
                    # gap
                    int(match.group(1)),
                    # merit
                    float(match.group(2)),
                    # submit format
                    " ".join(match.groups()),
                    # number
                    " ".join(match.groups()[2:]),
                    # line
                    line.strip(),
                ])

        if file_match > 5:
            # Also print last match
            print (match_print)

    if gaps:
        describe_found_gaps(gaps, args.verify)
        print_record_gaps(args, gaps)
    else:
        print("Didn't find any gaps in logs files: %s".format(
            ", ".join(log_files)))


def search_db(args):
    assert os.path.exists(args.search_db)

    gaps = []
    with sqlite3.connect(args.search_db, timeout=30) as conn:
        conn.row_factory = sqlite3.Row

        num_gaps = conn.execute('SELECT COUNT(*) FROM result').fetchone()[0]
        assert num_gaps > 100, num_gaps
        print(f"{num_gaps} results in {args.search_db!r}")

        # Min gap for current record (filters 80% of results)
        existing = conn.execute(
            """SELECT p, d, m, next_p, prev_p, next_p + prev_p as gapsize, merit FROM result
            WHERE  (next_p > 0 AND prev_p > 0) AND
                    ((merit > 21.9) OR
                    (gapsize > 30000 AND merit > 18.3) OR
                    (gapsize > 50000 AND merit > 15.3) OR
                    (gapsize > 70000 AND merit > 10.4) OR
                    (gapsize > 100000))"""
        ).fetchall()
        for gap in existing:
            gapsize = gap['gapsize']
            merit = gap['merit']
            number = "{} * {}#/{} -{}".format(
                gap["m"], gap["P"], gap["D"], gap["prev_p"])

            submit = "{:6d}  {:.3f}  {} to +{}".format(
                gapsize, merit,
                number, gap["next_p"])

            gaps.append((
                gapsize, merit, submit, number,
                ", ".join(f"{k}={gap[k]}" for k in ('p', 'd', 'm', 'prev_p', 'next_p')),
            ))

    describe_found_gaps(gaps, args.verify)
    print_record_gaps(args, gaps)


if __name__ == "__main__":
    parser = get_arg_parser()
    args = parser.parse_args()

    assert os.path.exists(args.prime_gaps_db), (
        f"Prime gaps database ({args.prime_gaps_db!r}) doesn't exist")

    neither = True
    if args.logs_directory and os.path.isdir(args.logs_directory):
        neither = False
        args.log_files.extend(glob.glob(args.logs_directory + "/*.log*"))

    if any(map(os.path.isfile, args.log_files)):
        neither = False
        search_logs(args.log_files)

    if os.path.isfile(args.search_db):
        neither = False
        search_db(args)

    if neither:
        print("Must pass --log-files, --logs-directory, or --search-db 'gaps.db'")
        sys.exit(1)
