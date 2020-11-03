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
import csv
import glob
import os.path
import math
import sqlite3
import subprocess

import gap_utils
import misc_utils


LOGS_DIR = "logs"

def get_arg_parser():
    parser = argparse.ArgumentParser('Search prime gaps logs')

    parser.add_argument('--search-db', type=str,
        default="prime-gap-search.db",
        help="Prime database from gap_test")

    parser.add_argument('-p', type=int, required=True,
        help="Prime to finalize")
    parser.add_argument('-d', type=int, required=True,
        help="Divisor to finalize")

    parser.add_argument('--mstart', type=int, default=1,
        help="finalize m between mstart and mend")

    parser.add_argument('--mend', type=int, default=10 ** 20,
        help="finalize m between mstart and mend")

    parser.add_argument('--no-dump-csv', action="store_true",
        help="Don't save csv results to file in logs dir")

    #parser.add_argument('--delete-unknown-files', action='store_true',
    #    help="Should unknow files be deleted")

    return parser


def count_results(conn, p, d, ms, me):
    count_r = conn.execute(
        'SELECT COUNT(*) FROM result WHERE p=? AND d=? AND m BETWEEN ? AND ?',
        (p, d, ms, me)).fetchone()
    count_m_stats = conn.execute(
        'SELECT COUNT(*) FROM m_stats WHERE p=? AND d=? AND m BETWEEN ? AND ?',
        (p, d, ms, me)).fetchone()
    return count_r[0], count_m_stats[0]

def check_contained(args, ms, minc):
    # Don't allow range to partially intersect [mstart, mend]

    if ms < args.mstart and me >= args.mstart:
        print(f"[{ms},{me+1} = {ms} + {minc}) overlaps mstart({mstart})")
        exit(1)

    me = ms + minc - 1
    if ms <= args.mend and me > args.mend:
        print(f"[{ms},{me+1} = {ms} + {minc}) overlaps mend({mend})")
        exit(1)

    return args.mstart <= ms <= args.mend


def sql_ranges_and_files(args, sql_ranges, unknown_fns):
    ranges = []
    # ---- Add sql ranges
    for r in sql_ranges:
        if args.p == r['P'] and args.d == r['D']:
            # Don't allow mstart, mend to bridge result
            if check_contained(args, r['m_start'], r['m_inc']):
                ranges.append(r)

    # ---- Find matching files
    unknown_files = []
    for unknown_fn in unknown_fns:
        parsed = gap_utils.parse_unknown_filename(unknown_fn)
        assert parsed, (parsed, unknown_fn)

        p, d, ms, mi, sl, mp, m1 = parsed
        if args.p == p and args.d == d:
            unknown_files.append(unknown_fn)

    return ranges, unknown_files


def verify_all_results(conn, ranges):
    # Check that result exists for each valid m in the range
    # XXX: add some --force option later

    all_done = True

    for r in ranges:
        ms = r['m_start']
        mi = r['m_inc']
        num_m = misc_utils.count_num_m(ms, mi, args.d)

        me = ms + mi - 1
        count_r, count_m = count_results(conn, args.p, args.d, ms, me)

        done = (num_m == count_r) and (num_m == count_m)
        all_done &= done

        print("\t[{:<8d}, {:>8d}]  m({:6d})  results: {:6d}, m_stats: {:6d}  Good: {}".format(
            ms, me, num_m, count_r, count_m, done))

    return all_done


def dump_query_to_file(conn, path, query, options):
    with open(path, "w") as csv_file:
        csv_writer = csv.writer(csv_file)
        cursor = conn.cursor()
        results = cursor.execute(query, options).fetchall()

        header = [i[0] for i in cursor.description]
        csv_writer.writerow(header)
        csv_writer.writerows(results)
    return header, results


def size_or_zero(path):
    if not os.path.exists(path):
        return 0
    return os.path.getsize(path)


def merged_ms_mi(ranges, unknown_files):
    # elements are (ms, mi, unk)
    joined = []

    unknowns = set()
    for ufn in sorted(unknown_files):
        p, d, ms, mi, _, _, _ = gap_utils.parse_unknown_filename(ufn)
        assert args.p == p and args.d == d, ufn

        unknowns.add((ms, mi))
        joined.append((ms, mi, ufn))

    for row in ranges:
        assert args.p == row['p'] and args.d == row['D'], row

        if (row['m_start'], row['m_inc']) in unknowns:
            continue

        params = [str(row[key]) for key in (
            'P', 'D', 'm_start', 'm_inc', 'sieve_length', 'max_prime')]
        # divide max_prime by million
        params[-1] = str(row['max_prime'] // 1000000)
        fake_ufn = gap_utils.generate_unknown_filename(*params)

        joined.append((row['m_start'], row['m_inc'], fake_ufn))

    return joined


def dump_to_file(conn, args, ranges, unknown_files):
    finalize_range = range(args.mstart, args.mend + 1)

    # output csv for each unknown_filename
    # output csv for each range (without unknown_filename)

    # elements are (ms, mi, unk)
    joined = merged_ms_mi(args, ranges, unknown_files)
    for _, _, ufn in joined:
        if ufn not in unknown_files:
            print (f"\tExtra range: {ufn!r}")

    # TODO: consider dumping to seperate sqlite3 db
    p = args.p
    d = args.d
    for ms, mi, ufn in joined:
        me = ms + mi - 1
        assert ms in finalize_range and me in finalize_range, ufn

        base_fn = os.path.join(LOGS_DIR, "results_" + os.path.splitext(ufn)[0])
        results_csv_fn = base_fn + ".csv"
        stats_fn = base_fn + ".stats.txt"

        if not args.no_dump_csv:
            print(f"\tDumping m_stats to\t{results_csv_fn!r}")
            size = size_or_zero(results_csv_fn)
            if size > 10 * 1024:
                print (f"\tAlready Processed {ufn!r} ({size//1024:,} kb)")
            else:
                # Write m_stats to .csv file
                header, m_stats = dump_query_to_file(
                    conn, results_csv_fn,
                    'SELECT * FROM m_stats WHERE P=? AND D=? AND m BETWEEN ? AND ?',
                    (p, d, ms, me))
                # Verify results got written to m_stats
                next_p = header.index("next_p")
                prev_p = header.index("prev_p")
                for row in m_stats:
                    assert row[next_p] != 0 or row[prev_p] != 0, row

        size = size_or_zero(stats_fn)
        if size < 2 * 1024:
            # XXX: Ugly but easier than lots of python code
            print(f"\tDumping stats to\t{stats_fn!r}")
            subprocess.check_call([
                "sqlite3", args.search_db,
                "-readonly", "-header", "-csv", "-echo",
                "-cmd", ".output "+ repr(stats_fn),
                f"""
                SELECT * FROM range
                WHERE P={p} AND D={d} AND m_start = {ms} AND m_inc = {mi};

                SELECT avg(merit), avg(prob_record), avg(prob_missing), avg(prob_merit),
                       avg(prp_next), avg(prp_prev), avg(test_time) FROM m_stats
                WHERE P={p} AND D={d} AND m BETWEEN {ms} AND {me};

                SELECT max(merit), max(prob_record), max(prob_missing), max(prob_merit),
                       max(prp_next), max(prp_prev), max(test_time) FROM m_stats
                WHERE P={p} AND D={d} AND m BETWEEN {ms} AND {me};

                SELECT cast(merit as int), count(*) FROM m_stats
                WHERE P={p} AND D={d} AND m BETWEEN {ms} AND {me} GROUP BY 1;

                SELECT PRINTF("%.4f %d * %d#/%d -%d to +%d   ",
                              merit, m, P, D, prev_p, next_p), * FROM m_stats
                WHERE P={p} AND D={d} AND m BETWEEN {ms} AND {me} ORDER BY merit DESC LIMIT 50;
                """
                ])


def delete_range_and_low_merit(conn, args, ranges, unknown_files):
    print()
    cursor = conn.cursor()

    # Need to check for missing ranges?
    joined = merged_ms_mi(ranges, unknown_files)
    for ms, mi, ufn in joined:
        me = ms + mi - 1
        print(f"\tDeleting range / low_merit results from {ufn!r}")

        condition = """WHERE P=? AND D=? AND m BETWEEN ? AND ? AND
            ((next_p + prev_p < 30000 AND MERIT < 20) OR
             (next_p + prev_p < 50000 AND MERIT < 15) OR
             (next_p + prev_p < 100000 AND MERIT < 8))
            """

        cursor.execute(
            "DELETE FROM m_stats " + condition
            (args.p, args.d, ms, me))
        m_stat_deletes = cursor.rowcount

        cursor.execute(
            "DELETE FROM result " + condition
            (args.p, args.d, ms, me))
        m_stat_deletes = cursor.rowcount

        print (f"\tDeleted {m_stat_deletes} and {m_stat_deletes} rows")

    for r in ranges:
        cursor.execute("DELETE FROM range WHERE rid = ?", (r['rid'],))
        assert cursor.rowcount == 1, r

    print("rm", " ".join(unknown_files))


def check_processed(args):
    unknown_fns = glob.glob("*M.txt")

    with sqlite3.connect(args.search_db) as conn:
        conn.row_factory = sqlite3.Row
        sql_ranges = conn.execute('SELECT * FROM range ORDER BY m_start ASC').fetchall()

        print ("{} ranges, {} unknown files".format(
            len(sql_ranges), len(unknown_fns)))

        # ---- Find ranges & unknown_files involved
        ranges, unknown_files = sql_ranges_and_files(args, sql_ranges, unknown_fns)

        print ("{} matching ranges, {} matching files".format(
            len(ranges), len(unknown_files)))

        #assert verify_all_results(conn, ranges), "Not all results present"

        conn.row_factory = None
        # Always dump stats, sometimes dump files
        dump_to_file(conn, args, ranges, unknown_files)

        print()
        print("tar'ing results with:")
        name = f"results_{args.p}_{args.d}"
        print(f"\ttar cvzf logs/archive_{name}.tar.gz logs/{name}_*")

        # Delete low merit result & m_stat
        delete_range_and_low_merit(conn, args, ranges, unknown_files)



if __name__ == "__main__":
    parser = get_arg_parser()
    args = parser.parse_args()

    print(args)
    print()

    assert os.path.exists(args.search_db), (
        f"Search database ({args.search_db!r}) doesn't exist")

    check_processed(args)
