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

import sys
# Hack to allow import of gap_utils
sys.path.append(".")

import argparse
import glob
import os.path
import math
import sqlite3
from collections import defaultdict

import gap_utils
import misc_utils


def get_arg_parser():
    parser = argparse.ArgumentParser('Search prime gaps logs')

    parser.add_argument('--search-db', type=str,
        default="prime-gap-search.db",
        help="Prime database from gap_test")

    return parser


def sql_and_file_ranges(sql_ranges, unknown_fns):
    # (P, D, ms, mi): {status, num_m, done}
    ranges = {}
    # key: [unknown_fn, sql_range]
    lookup = {}

    # ---- Add sql ranges
    for r in sql_ranges:
        key = tuple(r[key] for key in ('P', 'D', 'm_start', 'm_inc'))
        ranges[key] = [20, r['num_m'], 0]
        fn_fake = gap_utils.generate_unknown_filename(*key, "???", "???")
        lookup[key] = (fn_fake, r)
        assert r['num_m'] == misc_utils.count_num_m(r['m_start'], r['m_inc'], r['D'])

    # ---- Add file only ranges
    for unknown_fn in unknown_fns:
        parsed = gap_utils.parse_unknown_filename(unknown_fn)
        assert parsed, parsed

        p, d, ms, mi, sl, mp, m1 = parsed
        key = (p, d, ms, mi)
        if key in ranges:
            ranges[key][0] = 21
            lookup[key] = (unknown_fn, lookup[key][1])
        else:
            ranges[key] = [10, misc_utils.count_num_m(ms, mi, d), 0]
            lookup[key] = (unknown_fn, None)

    assert ranges.keys() == lookup.keys()
    return ranges, lookup


def build_and_count_pd_results(results, ranges, lookup):
    def add_new_range(p, d, start, end, count):
        count_m = misc_utils.count_num_m(start, end - start - 1, d)
        status = [41 if count_m == count else 30, count_m, count]
        k = (p, d, start, end)
        ranges[key] = status
        fn_fake = gap_utils.generate_unknown_filename(*key, "???", "???")
        lookup[key] = (fn_fake, None)
#        print(f"\t\tAdded range ({start}, {end}): {status}")

    pd_m = defaultdict(list)
    for P, D, m in results:
        pd_m[(P, D)].append(m)

    # Have a bunch of m, and ranges, find the uncovered ranges
    for pd, ms in pd_m.items():
        P, D = pd

        pd_ranges = {range(mstart, mstart + minc): 0
                for p, d, mstart, minc in ranges if P==p and D==d}

        ms.sort()
#        print (f"\tP: {P:5} D: {D:7}")
#        print (f"\t\tResults: {len(ms)} min: {min(ms)} max: {max(ms)}")
#        print ("\t\tKnown ranges:", ", ".join(map(str,pd_ranges.keys())))

        # In order find if uncovered
        new_range = [0, 0, 0]
        last_found = True
        for m in ms:
            found = False
            for r in pd_ranges:
                if m in r:
                    found = True
                    pd_ranges[r] += 1
            # Can't use for else, need to add +1 to possibly multiple ranges
            if found:
                if not last_found:
                    add_new_range(P, D, *new_range)
                    new_range = [0, 0, 0]
            if not found:
                # 0 if last_found, 1 if not last_found
                new_range[1 - last_found] = m
                new_range[2] += 1

            last_found = found

        if not found:
            add_new_range(P, D, *new_range)

        # Update ranges from pd_ranges
        for r, count in pd_ranges.items():
            key = (P, D, r.start, r.stop - r.start)
            ranges[key][2] = count
            status, count_m, _ = ranges[key]

            if count_m == count:
                ranges[key][0] = 41
            elif count > 0:
                ranges[key][0] = 30 + (status % 10)


def check_processed(args):
    """
    Check status of each UNKNOWN_FN, range, m_stats

    Find all pseudo ranges:
        range
        unknown_fn
        results (create a pseudo range)

    Classify as
        10: file

        20: range
        21: range + file

        30: partial_results
        31: partial_results + file,

        40: results + (range or file)
        41: results only
    """

    unknown_fns = glob.glob("*M.txt")

    with sqlite3.connect(args.search_db) as conn:
        conn.row_factory = sqlite3.Row
        sql_ranges = conn.execute('SELECT * FROM range').fetchall()

        # At some later point maybe don't load all (group by thousands or something)
        conn.row_factory = None
        results = conn.execute('SELECT P, D, m FROM result').fetchall()
        print (f"\tLoaded {len(results):,} results")

        # ---- Add file only ranges
        ranges, lookup = sql_and_file_ranges(sql_ranges, unknown_fns)

        # ---- Find results not belonging to any range
        build_and_count_pd_results(results, ranges, lookup)

        print_results(conn, ranges, lookup)

def check_mstatus(conn, p, d, ms, me):
    count = conn.execute(
        'SELECT COUNT(*) FROM m_stats WHERE p=? AND d=? AND m BETWEEN ? AND ?',
        (p, d, ms, me)).fetchone()
    return count[0]


def print_results(conn, ranges, lookup):
    # ---- Show Results
    def sort_by_status(range_kv):
        # Sort by status high to low
        # P, D, ms, mi
        return (-range_kv[1][0], range_kv[0])

    print (lookup)
    print()
    display_order = sorted(ranges.items(), key=sort_by_status)
    for key, value in display_order:
        p, d, ms, mi = key
        status, count_m, finished = value


        assert (status in (40, 41)) == (count_m == finished)
        print("P: {:5} D: {:7} M({:8}) {:9} to {:9} |".format(
            p, d, count_m, ms, ms + mi - 1), end=" ")
        print("Finished: {:8} {:6.1%} (filename: {})".format(
            finished, finished / count_m, lookup[key][0]))

    print("\n", "-" * 80, "\n", sep="")

    for key, value in display_order:
        p, d, ms, mi = key
        status, count_m, finished = value
        if status in (40, 41): continue

        fn = lookup[key][0]
        unk_fn_param = "--unknown-filename " + fn

        if status == 31:
            print("P: {:5} D: {:7} M({:8}) {:9} to {:9} |".format(
                p, d, count_m, ms, ms + mi - 1), end=" ")
            print("Finished: {:8}/{:<8} {:6.1%}".format(
                finished, count_m, finished/count_m))

            print ("Partially finished resume:")
            print (f"\t./gap_test.py {unk_fn_param}\n")
        elif status in (20, 30):
            print ("File missing {}(could be recreated with):".format(
                "(with partial results) " * (status == 30)))
            print (f"\t./combined_sieve --save-unknowns {unk_fn_param}\n")
            r = lookup[key][1]
            if r:
                print ("\tor deleted with")
                print (f"\tsqlite3 {args.search_db} 'PRAGMA foreign_keys=1; "\
                       f"DELETE FROM range WHERE rid={r['rid']}; SELECT TOTAL_CHANGES()'\n")
        elif status == 21:
            mstatus = check_mstatus(conn, p, d, ms, ms + mi - 1)
            if mstatus:
                print (f"Ready to start testing ({mstatus} m_stats):")
                print (f"\t./gap_test.py {unk_fn_param}\n")
            else:
                print ("Ready for gap_stats")
                print (f"\t./gap_stats --save-unknowns {unk_fn_param}\n")
        elif status == 10:
            print ("range missing:")
            print (f"\t./gap_stats --save-unknowns {unk_fn_param}\n")
        else:
            print (f"unknown status: {status}")


if __name__ == "__main__":
    parser = get_arg_parser()
    args = parser.parse_args()

    assert os.path.exists(args.search_db), (
        f"Search database ({args.search_db!r}) doesn't exist")

    check_processed(args)
