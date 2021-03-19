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

import argparse
import math
import multiprocessing
import os
import queue
import signal
import sqlite3
import time

import gmpy2

import gap_utils
import gap_utils_primes
import gap_test_plotting
import gap_test_stats
import misc.misc_utils as misc_utils


MIN_MERIT_DEFAULT = 18

def get_arg_parser():
    parser = argparse.ArgumentParser('Test prime gaps')

    parser.add_argument(
        '--search-db', type=str,
        default="prime-gap-search.db",
        help="Database for this project (default: %(default)s)")

    parser.add_argument(
        '--prime-gaps-db', type=str,
        default="gaps.db",
        help="Prime gap database (default: %(default)s) "
             "See github.com/primegap-list-project/prime-gap-list")

    parser.add_argument(
        '-u', '--unknown-filename', type=str, required=True,
        help="determine p, d, mstart, minc, sieve-length, and max-prime"
             " from unknown-results filename")

    parser.add_argument(
        '--min-merit', type=int, default=MIN_MERIT_DEFAULT,
        help="only display prime gaps with merit >= min merit")

    parser.add_argument(
        '--megagap', type=int, default=0,
        help="Search for megagaps (of this size or greater)")

    parser.add_argument(
        '-t', '--threads', type=int, default=1,
        help="Number of threads to use for searching (default: %(default)s)")

    parser.add_argument(
        '--taskset', action='store_true',
        help="taskset each thread, helps when using more than 50%% of CPU threads")

    parser.add_argument(
        '--num-plots', type=int, default=0,
        help="Show up to 3 plots about distributions")

    parser.add_argument(
        '--save-logs', action='store_true',
        help="Save logs and plots about distributions")

    parser.add_argument(
        '--prp-top-percent', type=int, default=None,
        help="Only test top X%% (sorted by prob(record))")

    parser.add_argument(
        '--no-one-side-skip', action='store_true',
        help="Don't skip when one prev_prime is very small")

    return parser


def K_and_stats(args):
    P = args.p
    D = args.d

    assert P <= 80000
    K = gmpy2.primorial(P)
    K, r = divmod(K, D)
    assert r == 0

    K_digits = gmpy2.num_digits(K, 10)
    K_bits = gmpy2.num_digits(K, 2)
    K_log = float(gmpy2.log(K))

    return K, K_digits, K_bits, K_log


# ---- gap_testing ---- #

# NOTE: Manager is prime_gap_test maintains workers who run process_line.
def process_line(
        early_stop_flag, work_q, results_q, prob_side_threshold, sides_tested,
        prob_prime, prob_prime_after_sieve, record_gaps, side_skip_enabled,
        thread_i, SL, K, P, D, megagap,
        primes, remainder):
    def cleanup(*_):
        work_q.close()
        work_q.join_thread()
        results_q.close()
        results_q.join_thread()
        print(f"\tThread {thread_i} stopped")
        exit(0)

    # Calculate extended gap (see gap_stats)
    prob_prime_coprime, coprime_extended = gap_test_stats.setup_extended_gap(SL, P, D, prob_prime)
    K_mod_d = K % D

    D_primes, is_offset_coprime, coprime_X = gap_test_stats.setup_bitarray(SL, P, D)

    # Ignore KeyboardInterrupt (Manager catches first and sets early_stop_flag)
    # Normal flow is receive sentinel (None) stop
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGTERM, cleanup)

    print(f"\tThread {thread_i} started")
    while True:
        # timeout should never happen, but might happen if iterating the file
        # without queuing anything on large file (or slow filesystem)
        work = work_q.get(timeout=600) # See note above

        # Always mark the item as done, progress is tracked with sc.will_test
        work_q.task_done()

        if early_stop_flag.is_set():
            assert work is None

        if work is None:
            cleanup()

        m, mi, prev_p, next_p, prob_record, log_n, line = work

        if b' || ' in line[:30]:
            m_test, unknown_l, unknown_u, unknowns = gap_utils.parse_compressed_line(
                SL, K, D, m,
                is_offset_coprime, coprime_X, D_primes,
                line)
        else:
            m_test, unknown_l, unknown_u, unknowns = gap_utils.parse_unknown_line(line)

        assert m == m_test

        # Used for openPFGW
        str_n = "{}*{}#/{}+".format(m, P, D)

        p_tests = 0
        n_tests = 0
        # next_p of 0 would mean next_prime was skipped
        next_p = 0
        # Prob(record) after finding prev_prime
        new_prob_record = 0
        # Will be updated based on side_skip_enabled / new_prob_record
        test_next = True

        if megagap:
            action = "Starting" if prev_p == next_p == 0 else "Resuming"
            print(f"\t{action} {str_n:25} on thread {thread_i}")

        t0 = time.time()

        # prev_p > 0 means we loaded a partial result
        if prev_p <= 0:
            p_tests, prev_p = gap_utils_primes.determine_prev_prime(
                    m, str_n, K, unknowns[0], SL, primes, remainder)

        prev_time = time.time() - t0

        if side_skip_enabled:
            # Check if slower to test a new record or next prime

            # On average sum(new_prob_record) should approx equal sum(prob_record)
            # benchmarked at 2.5% overheard of prev_prime for P=1511
            new_prob_record = gap_test_stats.prob_record_one_sided(
                    record_gaps, megagap, prev_p,
                    unknowns[1], prob_prime_after_sieve,
                    m, K_mod_d, D, coprime_extended, prob_prime_coprime)

            # This is probable the wrong logic for --megagap, because it requires gap_stats
            # to have been called with the correct --min-merit.

            with prob_side_threshold.get_lock():
                threshold = prob_side_threshold.value / sides_tested.value

                test_next = new_prob_record > threshold
                prob_side_threshold.value += prob_record - (0 if test_next else new_prob_record)
                sides_tested.value += 1 + test_next

        should_save_partial = megagap or (p_tests and prev_time > 100)
        if should_save_partial:
            print("\t{:5} | m: {:8} | count: {}/{} {} | {:.3g} => {:.3g} | {:.0f} secs".format(
                    prev_p, m, p_tests, unknown_l, unknown_u,
                    prob_record, new_prob_record, prev_time))

            if test_next:
                # Save partial before starting next_prime
                # XXX: consider updating prob_record to new_prob_record
                results_q.put((
                    m, log_n,
                    unknown_l, unknown_u,
                    0, -1,  # next_p = -1, indicates partial result
                    p_tests, prev_p,
                    0, 0, prev_time))

                # Can take A LONG TIME, so allow early quitting if partial is saved
                if early_stop_flag.is_set():
                    # will cleanup at top
                    continue

        if test_next:
            n_tests, next_p = gap_utils_primes.determine_next_prime(m, str_n, K, unknowns[1], SL)

        test_time = time.time() - t0

        results_q.put((
            m, log_n,
            unknown_l, unknown_u,
            n_tests, next_p,
            p_tests, prev_p,
            prob_record, new_prob_record,
            test_time,
        ))


def run_in_parallel(
        args, conn, unknown_file, record_gaps,
        prob_prime, prob_prime_after_sieve,
        existing,
        K, K_log,
        data
):
    # XXX: Cleanup after gmpy2.prev_prime.
    # Remainders of (p#/d) mod prime
    primes = tuple([p for p in range(3, 80000+1) if gmpy2.is_prime(p)])
    remainder = tuple([K % prime for prime in primes])

    # Based on args
    prob_threshold, m_probs = gap_test_stats.determine_test_threshold(args, data)
    if prob_threshold >= 0:
        print("Testing {} m where prob(record) >= {:.3g}".format(
            len(m_probs), prob_threshold))

    # Are there any non-processed (or partial) m_probs.
    if all(m in existing and existing[m][1] >= 0 for m in m_probs):
        print(f"All prp-top-percent({len(m_probs)}) already processed!")
        print()
        return

    # Worker setup
    assert args.threads in range(1, 65), args.threads

    # Try to keep at least this many in the queue
    min_work_queued = 3 * (args.threads + 1)

    early_stop_flag = multiprocessing.Event()
    work_q = multiprocessing.JoinableQueue()
    results_q = multiprocessing.Queue()

    # Used to dynamically set prob_threshold, See THEORY.md#one-sided-tests
    prob_side_threshold = multiprocessing.Value('d', 5 * prob_threshold)
    sides_tested = multiprocessing.Value('l', 10)

    shared_args = (
        early_stop_flag,
        work_q,
        results_q,
        prob_side_threshold,
        sides_tested,
    )

    one_sided_args = (
        prob_prime,
        prob_prime_after_sieve,
        record_gaps,
        not args.no_one_side_skip,
    )

    print()

    processes = []
    for i in range(args.threads):
        process = multiprocessing.Process(
            target=process_line,
            args=(
                *shared_args,
                *one_sided_args,
                i, args.sieve_length, K, args.p, args.d,
                args.megagap,
                primes, remainder
            )
        )
        process.start()
        if args.taskset:
            os.system("taskset -p -c {} {}".format(i, process.pid))
        time.sleep(0.01)
        processes.append(process)

    time.sleep(1)
    print()

    # This sets the start time which can be hundreds of seconds off if most of
    # valid_mi has already been processed. In that case some stats take a while
    # to stablize.
    sc = gap_test_stats.StatCounters(time.time(), time.time())

    try:
        for mi in data.valid_mi:
            m = args.mstart + mi
            log_n = (K_log + math.log(m))

            # Read a line from the file (can be the slowest part of restarting)
            line = unknown_file.readline()
            assert line.startswith(str(m).encode()), (
                "file format has changed, lines should start with <mstart + mi> not <mi>"
                "you can add <m> to each line, recreate, or git checkout <old commit>"
                "sorry")

            exist = existing.get(m)
            is_partial = exist and exist[1] < 0
            if is_partial:
                # Should always run partial results even if don't match prob_threshold.
                if m not in m_probs or m_probs[m][1] < prob_threshold:
                    m_probs[m] = (1.01 * prob_threshold for _ in range(2))

            # Check if prob_record < threshold
            if m not in m_probs or m_probs[m][1] < prob_threshold:
                continue

            prev_p, next_p = 0, 0
            if exist:
                prev_p, next_p = exist
                # next_p < -1, related to missing_gap
                # next_p = -1, is partial result (should be continued)
                if next_p >= 0:
                    # next_p = 0, is side skip
                    # next_p > 0, result
                    continue

            sc.will_test += 1

            # Should never wait because of lower min_work_queued blocking below
            work_q.put_nowait((m, mi, prev_p, next_p, m_probs[m][1], log_n, line))

            # Process any finished results
            while (sc.will_test - sc.tested) >= min_work_queued:
                result = results_q.get()
                # Process result contains save to DB call see gap_test_stats.py
                gap_test_stats.process_result(conn, args, record_gaps, m_probs, data, sc, result)

    except (KeyboardInterrupt, queue.Empty) as e:
        print("Received first  Ctrl+C | Waiting for current work to finish")
        time.sleep(0.1)

        early_stop_flag.set()

        # Flush queue and wait on current results
        try:
            while True:
                work_q.get_nowait()
                work_q.task_done()
                sc.will_test -= 1
        except queue.Empty:
            pass

        time.sleep(0.01)

    print("No more work pushing NONE")
    for _ in processes:
        work_q.put(None)

    while sc.tested < sc.will_test:
        # Partial results cause multiple prints of "waiting on X..."
        print(f"Waiting on {sc.will_test-sc.tested} of {sc.will_test} results")
        try:
            result = results_q.get()

            # partial (-1) doesn't increment tested in process_result
            if early_stop_flag.is_set() and result[5] == -1:
                sc.tested += 1

            gap_test_stats.process_result(conn, args, record_gaps, m_probs, data, sc, result)
        except KeyboardInterrupt:
            # commit any potentially uncommitted save() results.
            conn.commit()
            print("Received second Ctrl+C | Terminating now")
            for i, p in enumerate(processes):
                if p.is_alive():
                    print("\tTerminating thread", i)
                    p.terminate()
            exit(1)

    # Make sure everything was committed
    conn.commit()

    print("Joining work_q (should be instant)")
    work_q.join()
    work_q.close()
    work_q.join_thread()

    print(f"Joining {len(processes)} processes")
    for i, process in enumerate(processes):
        process.join(timeout=0.1)

    print("Done!")


def prime_gap_test(args):
    P = args.p
    D = args.d

    M = args.mstart
    M_inc = args.minc

    SL = args.sieve_length
    max_prime = args.max_prime

    K, K_digits, K_bits, K_log = K_and_stats(args)
    M_log = K_log + math.log(M)
    print("K = {} bits, {} digits, log(K) = {:.2f}".format(K_bits, K_digits, K_log))

    # ----- Check if pfgw64 is present for large P
    if K_bits > 7800:
        try:
            assert gap_utils_primes.openPFGW_is_prime("2 ^ 8000 + 3081")
        except AssertionError:
            print()
            print("pfgw64 suggested for P# > 7800#")
            print("please download and place executable in prime-gap/ directory")
            print()
            exit(1)

    # ----- Open prime-gap-search.db
    # Longer timeout so that record_checking doesn't break when saving
    conn = sqlite3.connect(args.search_db, timeout=60)

    # Try to load range data for args
    range_data_db = gap_test_stats.load_range(conn, args)
    if (range_data_db and
            'min_merit' in range_data_db and
            args.min_merit == MIN_MERIT_DEFAULT):
        args.min_merit = range_data_db['min_merit']


    # ----- Load Record sized gaps
    with sqlite3.connect(args.prime_gaps_db) as prime_db_conn:
        record_gaps = gap_test_stats.load_records(prime_db_conn, M_log)

        # Calculate min_merit record, min_merit for 50% record, min_merit for 90% record
        min_record = record_gaps[0] / M_log
        # Look 250 records later and check if contains more than half of numbers
        avg_record = next(record_gaps[r] for r in range(len(record_gaps) - 250)
            if record_gaps[r + 250] - record_gaps[r] < 1000)

        min_merit_gap = int(args.min_merit * M_log)
        print("Min Gap ~= {} (for merit > {:.1f})\n".format(
            min_merit_gap, args.min_merit))

        print("\tLoaded {} records ({} to {}) from {!r}".format(
            len(record_gaps), record_gaps[0], record_gaps[-1], args.prime_gaps_db))
        print("\tMin merit for record {:.2f}, gapsize for likely record {:.2f} {}".format(
            min_record, avg_record / M_log, avg_record))
        if args.megagap:
            print("\tMin merit for megagap({}) {:.2f}\n".format(
                args.megagap, args.megagap / M_log))

        record_gaps = set(record_gaps)


    # ----- Open Output file
    print("\tLoading unknowns from {!r}".format(args.unknown_filename))
    print()
    fn_with_unk = gap_utils.transform_unknown_filename(
            args.unknown_filename, "unknowns", "txt")
    if os.path.exists(fn_with_unk):
        args.unknown_filename = fn_with_unk
    unknown_file = open(args.unknown_filename, "rb")

    count_m = misc_utils.count_num_m(M, M_inc, D)

    # ----- Load data from prime-gap-search.db
    t0 = time.time()
    existing = gap_test_stats.load_existing(conn, args)
    t1 = time.time()
    n_exist = len(existing)
    print(f"Found {n_exist:,} ({n_exist/count_m:.1%}) results ({t1-t0:.1f} sec)")

    if len(existing) == count_m:
        print(f"All processed!")
    elif args.prp_top_percent == 0:
        print("--prp-top-percent=0, skipping testing")

    # ----- Allocate memory for a handful of utility functions.

    # ----- Sieve stats
    prob_prime = 1 / M_log - 1 / (M_log * M_log)
    prob_prime_after_sieve = gap_test_stats.prob_prime_sieve_length(
        M, K, D, P, prob_prime, K_digits, SL, max_prime)

    # ----- Main sieve loop.

    valid_mi = misc_utils.calc_valid_mi(M, M_inc, D)

    data = gap_test_stats.GapData()

    data.first_m = M + valid_mi[0]
    data.last_m = M + valid_mi[-1]
    data.valid_mi = valid_mi

    print(f"\nStarting m({len(valid_mi):,}) {data.first_m:,} to {data.last_m:,}")
    print()

    # Load stats for prob_record
    temp = gap_test_stats.load_probs_only(conn, args)
    assert len(temp.prob_merit_gap) == len(valid_mi), \
            ("run ./gap_stats first", len(valid_mi), len(temp.prob_merit_gap))
    assert len(temp.prob_record_gap) == len(valid_mi), \
            ("run ./gap_stats first", len(valid_mi), len(temp.prob_record_gap))
    data.prob_merit_gap = temp.prob_merit_gap
    data.prob_record_gap = temp.prob_record_gap
    del temp

    run_in_parallel(
        args, conn, unknown_file, record_gaps,
        prob_prime, prob_prime_after_sieve,
        existing,
        K, K_log,
        data
    )

    # ----- Plots
    if args.num_plots:
        # Load stats from gap_stats
        data_db, misc_db = gap_test_stats.load_stats(conn, args)
        data_db.valid_mi = valid_mi

        assert len(data_db.experimental_gap) >= len(data.experimental_gap)
        assert misc_db.prob_gap_comb, len(misc_db.prob_gap_comb)

        del data

        if False:
            # VERY slowly validates gap_stats results.
            gap_test_stats.validate_prob_record_merit(
                    args, data_db, K_log, prob_prime_after_sieve)

        # XXX: move inside plot_stuff
        with open(args.unknown_filename, "rb") as unknown_file_repeat:
            if False:
                # First three lines is easy.
                unk_mi_of_interest = valid_mi[:3]
            else:
                # Requires reading more of the file,
                # generates more interesting result (best, worst, average)
                pi = sorted(zip(data_db.prob_record_gap, valid_mi), reverse=True)
                unk_mi_of_interest = [pi[0][1], pi[-1][1], pi[len(pi)//2][1]]
                print("Best, Worst, Avg m:", [M + mi for mi in unk_mi_of_interest])

            for mi in valid_mi:
                line = unknown_file_repeat.readline()
                assert line, mi
                if mi in unk_mi_of_interest:
                    m, _, _, unknowns = gap_utils.parse_unknown_line(line)
                    assert m == config.mstart + mi, (m, config.mstart, mi)

                    # TODO append valid_mi[k] and prob
                    misc_db.test_unknowns[mi] = unknowns

                    if mi == max(unk_mi_of_interest):
                        break

        gap_test_plotting.plot_stuff(
            args, data_db, misc_db,
            min_merit_gap, record_gaps, prob_prime_after_sieve)


if __name__ == "__main__":
    parser = get_arg_parser()
    args = parser.parse_args()
    gap_utils.verify_args(args)

    # TeeLogger context if args.save_logs
    with gap_utils.logger_context(args):
        prime_gap_test(args)

