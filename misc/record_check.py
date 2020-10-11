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
import itertools
import os.path
import re
import sqlite3
import sys
from collections import Counter, defaultdict


def get_arg_parser():
    parser = argparse.ArgumentParser('Search prime gaps logs')

    parser.add_argument('--logs-directory', type=str,
        default="logs",
        help="directory of logs")

    parser.add_argument('--search-db', type=str,
        default="prime-gap-search.db",
        help="Prime database from gap_test")

    parser.add_argument('--prime-gaps-db', type=str,
        default="gaps.db",
        help="Prime gap database see github.com/primegap-list-project/prime-gap-list")

    parser.add_argument('--whoami', type=str,
        default="S.Troisi",
        help="Name in gap database to ignore for already submitted records")

    return parser


def describe_found_gaps(gaps):
    gaps_by_size  = sorted(g[0] for g in gaps)
    gaps_by_merit = defaultdict(list)
    for gap in gaps:
        gaps_by_merit[int(gap[1])].append(gap)

    print ("Found {} gaps  ({:<6} to {:>6})".format(
        len(gaps_by_size), gaps_by_size[0], gaps_by_size[-1]))
    print ("      {} merit ({:.3f} to {:.3f})".format(
        " "*len(str(len(gaps_by_size))),
        min(g[1] for g in gaps), max(g[1] for g in gaps)))
    for int_merit, items in sorted(gaps_by_merit.items()):
        print ("    Merit {:<2} x{:<4} | {}".format(
            int_merit,
            len(items),
            ", ".join([str(rec[0]) for rec in sorted(items)[:2]] +
                      ["..."] * (len(items) > 2))
        ))
    print ()


def print_record_gaps(args, gaps):
    with sqlite3.connect(args.prime_gaps_db) as conn:
        num_gaps = conn.execute('SELECT COUNT(*) FROM gaps').fetchone()[0]
        assert num_gaps > 50000, num_gaps
        min_merit, max_merit = conn.execute(
            'SELECT min(merit), max(merit) FROM gaps').fetchone()
        print ("Records, {} gaps merit {:.2f} to {:.2f}".format(
            num_gaps, min_merit, max_merit))

        print ("Checking", len(gaps), "gaps against records")
        own_records = []
        record_lines = []
        for gap in gaps:
                size = gap[0]
                new_merit = gap[1]
                # These filters generated with
                # sqlite3 gaps.db  "select min(merit) from gaps where gapsize BETWEEN 1500 AND 50000;"
                if size <= 30000 and new_merit < 21.9:
                    continue
                if size <= 50000 and new_merit < 18:
                    continue
                if size <= 10000 and new_merit < 10.4:
                    continue

                existing = conn.execute(
                    'SELECT merit,primedigits,startprime,discoverer FROM gaps WHERE'
                    ' gapsize=?', (size,)).fetchone()
                if (not existing) or (new_merit - existing[0] > 6e-3):
                    own_record = existing is not None and existing[3] == args.whoami
                    if own_record:
                        own_records.append(gap[2])
                    else:
                        record_lines.append(gap[2])
                    print ("Hi", existing, not existing, "{:.6e}".format(new_merit - existing[0]))

                    print ("\tRecord {:3d}{:1} | {}\tGap={} (old: {}{})".format(
                        len(record_lines), "*" * own_record, gap[3], size,
                        None if not existing else existing[0],
                        " by you" * own_record))

        if record_lines:
            print ()
            print ("Records({}) {}:".format(
                len(record_lines),
                f"({len(own_records)} already submitted)" if own_records else ""))
            for line in record_lines:
                print (line)
            print ()


def search_logs(args):
    # 32280  10.9749  5641 * 3001#/2310 -18514 to +13766
    record_format = re.compile(
        r"(\d+)\s+(\d+\.\d+)\s+"
        r"(\d+)\s*(\*)\s*(\d+\#)(/\d+\#?)\s+"
        r"(-\d+)")

    assert os.path.exists(args.logs_directory), (
        "Logs directory ({}) doesn't exist".format(args.logs_directory))

    gaps = []
    for log_fn in glob.glob(args.logs_directory + "/*.log"):
        with open(log_fn, "r") as f:
            lines = f.readlines()

        print ("Processing {:5d} lines from log: {}".format(
            len(lines), log_fn))

        file_match = 0
        for li, line in enumerate(lines):
            match = record_format.search(line)
            if match:
                file_match += 1
                if file_match < 6:
                    partial_line = line.split(':')[-1].strip()
                    print ("    Match {} at line {}: {}".format(
                        file_match, li, partial_line))
                gaps.append([
                    int(match.group(1)), # gap
                    float(match.group(2)), # merit
                    " ".join(match.groups()),
                    line,
                ])

    if gaps:
        describe_found_gaps(gaps)
        print_record_gaps(args, gaps)
    else:
        print ("Didn't find any gaps in logs diroctory({})".format(
            args.logs_directory))


def search_db(args):
    assert os.path.exists(args.search_db)

    gaps = []
    with sqlite3.connect(args.search_db, timeout=10) as conn:
        conn.row_factory = sqlite3.Row

        num_gaps = conn.execute('SELECT COUNT(*) FROM result').fetchone()[0]
        assert num_gaps > 1000, num_gaps
        print (f"{num_gaps} results in {args.search_db!r}")

        # Min gap for current record (filters 80% of results)
        existing = conn.execute(
            'SELECT p, d, m, next_p_i, prev_p_i, merit FROM result '
            'WHERE  merit > 18 or (next_p_i + prev_p_i) > 50000 '
#            'ORDER BY p, d, m'
        ).fetchall()
        for gap in existing:
            gapsize = gap['next_p_i'] + gap['prev_p_i']
            merit = gap['merit']
            line = "{:6d}  {:.3f}  {} * {}#/{} -{} to +{}".format(
                gapsize, merit,
                gap["m"], gap["P"], gap["D"],
                gap["prev_p_i"], gap["next_p_i"])

            gaps.append((
                gapsize, merit,
                line,
                ", ".join(f"{k}={gap[k]}" for k in gap.keys()),
            ))

    describe_found_gaps(gaps)
    print_record_gaps(args, gaps)



if __name__ == "__main__":
    parser = get_arg_parser()
    args = parser.parse_args()

    assert os.path.exists(args.prime_gaps_db), (
        f"Prime gaps database ({args.prime_gaps_db!r}) doesn't exist")

    neither = True
    if os.path.exists(args.logs_directory):
        neither = False
        search_logs(args)

    if os.path.exists(args.search_db):
        neither = False
        search_db(args)

    if neither:
        print ("Must pass --logs-directory or --search-db 'gaps.db'")
        sys.exit(1)

