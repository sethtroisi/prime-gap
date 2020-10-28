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
import itertools
import math
import multiprocessing
import os
import signal
import sqlite3
import time
from collections import defaultdict
from dataclasses import dataclass

import gmpy2

import gap_utils
import gap_test_plotting


def get_arg_parser():
    parser = argparse.ArgumentParser('Test prime gaps')

    parser.add_argument('--search-db', type=str,
        default="prime-gap-search.db",
        help="Database for this project (default: %(default)s)")

    parser.add_argument('--prime-gaps-db', type=str,
        default="gaps.db",
        help="Prime gap database (default: %(default)s) see github.com/primegap-list-project/prime-gap-list")

    parser.add_argument('--unknown-filename', type=str,
        help="determine p, d, mstart, minc, sieve-length, and max-prime"
             " from unknown-results filename")

    parser.add_argument('--min-merit', type=int, default=12,
        help="only display prime gaps with merit >= minmerit")

    parser.add_argument('--threads', type=int, default=1,
        help="Number of threads to use for searching (default: %(default)s)")

    parser.add_argument('--num-plots', type=int, default=0,
        help="Show plots about distributions")

    parser.add_argument('--stats', action='store_true',
        help="SLOWLY Calculate stats (as doublecheck of gap_stats)")

    parser.add_argument('--save-logs', action='store_true',
        help="Save logs and plots about distributions")

    parser.add_argument('--prp-top-percent', type=int, default=None,
        help="Only test top X%% (sorted by prob(record))")

    parser.add_argument('--no-one-side-skip', action='store_true',
        help="Don't skip when one prev_prime is very small")

    return parser


#---- gap_testing ----#


def config_hash(config):
    h =          config.mstart
    h = h * 31 + config.minc;
    h = h * 31 + config.p;
    h = h * 31 + config.d;
    h = h * 31 + config.sieve_length;
    h = h * 31 + config.max_prime;
    return h % (2 ** 64)


def load_records(conn, log_N):
    # Find records merit <= gap * log_N
    # Stop at merit = 35

    # NOTE: to make sure 'count(records)' still counts after a
    # found record is submitted subtract a small amound from merit
    rv = conn.execute(
        "SELECT gapsize FROM gaps "
        "WHERE gapsize > merit * ? - 0.001"
        "  AND gapsize < 35 * ?",
        (log_N, log_N))
    records = [g[0] for g in rv]

    # Handle missing gaps
    SMALLEST_MISSING = 113326
    largest_likely = int(log_N * 30)
    if SMALLEST_MISSING < largest_likely:
        rv = conn.execute(
            "SELECT gapsize FROM gaps "
            "WHERE gapsize BETWEEN ? AND ?",
            (SMALLEST_MISSING, largest_likely))
        # Any number not present in this list is 'missing'
        all_gaps = range(SMALLEST_MISSING, largest_likely+1, 2)
        missing = list(set(all_gaps) - set(g[0] for g in rv))
        records.extend(missing)

    return tuple(sorted(records))


def load_existing(conn, args):
    rv = conn.execute(
        "SELECT m, next_p_i, prev_p_i FROM result"
        " WHERE P = ? AND D = ? AND m BETWEEN ? AND ?",
        (args.p, args.d, args.mstart, args.mstart + args.minc - 1))
    return {row['m']: [row['prev_p_i'], row['next_p_i']] for row in rv}


@dataclass
class Datas:
    first_m = -1
    last_m = -1

    valid_m = []

    # Values for each m in valid_m
    expected_prev = []
    expected_next = []
    expected_gap  = []
    experimental_side = []
    experimental_gap = []
    prob_merit_gap = []
    prob_record_gap = []

@dataclass
class Misc:
    prob_gap_side  = defaultdict(float)
    prob_gap_comb  = defaultdict(float)

    test_unknowns = []


def load_stats(conn, args):
    data = Datas()
    misc = Misc()

    rv = conn.execute(
        "SELECT e_gap_prev, e_gap_next, e_gap_prev + e_gap_next,"
        "       prev_p, next_p,"
        "       prob_merit, prob_record"
        " FROM m_stats WHERE P = ? AND D = ? AND m BETWEEN ? AND ?",
        (args.p, args.d, args.mstart, args.mstart + args.minc - 1))

    # Will fail if non-present
    (data.expected_prev, data.expected_next, data.expected_gap,
     exp_prev, exp_next,
     data.prob_merit_gap, data.prob_record_gap)= zip(*[tuple(row) for row in rv])

    # Interweave these
    data.experimental_side = [s for s in exp_prev + exp_next if s > 0]
    data.experimental_gap = [(p + n) for p, n in zip(exp_prev, exp_next) if p > 0 and n > 0]

    # Need dictionary
    m_values = len(data.prob_merit_gap)

    rv = conn.execute(
        "SELECT gap, prob_combined, prob_low_side, prob_high_side "
        "FROM range_stats where rid = ?",
        (config_hash(args),))
    for row in rv:
        gap = row['gap']
        # TODO: test if this works with out the normalization
        # Values were normalized by / m_values in gap_stats
        misc.prob_gap_comb[gap] += row['prob_combined'] * m_values
        misc.prob_gap_side[gap] += row['prob_low_side'] / 2 * m_values
        misc.prob_gap_side[gap] += row['prob_high_side'] / 2 * m_values

    return data, misc


def save(conn, p, d, m, next_p_i, prev_p_i, merit,
         n_tests, p_tests, test_time):
    assert p in range(201, 80000, 2), (p, d, m)
    conn.execute(
        "INSERT INTO result(P,D,m,next_p_i,prev_p_i,merit)"
        "VALUES(?,?,?,  ?,?,  ?)",
        (p, d, m, next_p_i, prev_p_i,  round(merit,4)))

    conn.execute(
        "UPDATE m_stats "
        "SET next_p=?, prev_p=?, merit=?,"
        "    prp_next=?, prp_prev=?, test_time=?"
        "WHERE p=? AND d=? AND m=?",
        (next_p_i, prev_p_i, round(merit,4),
         p_tests, n_tests, test_time, p, d, m))
    conn.commit()


def prob_prime_sieve_length(M, K, D, prob_prime, K_digits, P_primes, SL, max_prime):
    assert max_prime >= 10 ** 6, max_prime

    # From Mertens' 3rd theorem
    gamma = 0.577215665
    unknowns_after_sieve = 1.0 / (math.log(max_prime) * math.exp(gamma))

    prob_prime_after_sieve = prob_prime / unknowns_after_sieve

    prob_prime_coprime_p = 1
    for prime in P_primes:
        prob_prime_coprime_p *= (1 - 1/prime)

    # TODO can improve by doing gmpy2.gcd(K, i) in setup.

    # See "Optimizing Choice Of D" in THEORY.md for why this is required
    count_coprime_p = 0
    m_tests = [M+i for i in range(100) if math.gcd(M+i, D) == 1][:6]
    for m in m_tests:
        N0 = m * K
        for i in range(1, SL+1):
            N = N0 + i
            if math.gcd(N, D) == 1 and gmpy2.gcd(N, K) == 1:
                count_coprime_p += 1

    count_coprime_p //= len(m_tests)

    chance_coprime_composite = 1 - prob_prime / prob_prime_coprime_p
    prob_gap_shorter_hypothetical = chance_coprime_composite ** count_coprime_p

    expected = count_coprime_p * (unknowns_after_sieve / prob_prime_coprime_p)
    print("\texpect {:.0f} left, {:.3%} of SL={} after {}M".format(
        expected, expected / (SL + 1), SL, max_prime//1e6))
    print("\t{:.3%} of {} digit numbers are prime".format(
        prob_prime, K_digits))
    print("\t{:.3%} of tests should be prime ({:.1f}x speedup)".format(
        prob_prime_after_sieve, 1 / unknowns_after_sieve))
    print("\t~2x{:.1f} = {:.1f} PRP tests per m".format(
        1 / prob_prime_after_sieve, 2 / prob_prime_after_sieve))
    print("\tsieve_length={} is insufficient ~{:.2%} of time".format(
        SL, prob_gap_shorter_hypothetical))
    print()

    return prob_prime_after_sieve


def calculate_expected_gaps(
        SL, min_merit_gap, prob_nth, prob_longer, log_n, unknowns,
        data, misc):

    for i, side in enumerate(unknowns):
        expected_length = 0
        for v, prob in zip(side, prob_nth):
            expected_length += abs(v) * prob
            # Normalized to 1 by matplotlib.hist later.
            misc.prob_gap_side[abs(v)] += prob

        # expected to encounter a prime at distance ~= ln(n)
        prob_gap_longer = prob_nth[len(side)]
        assert prob_gap_longer < 0.01, (prob_gap_longer, len(side))
        expected_length += (SL + log_n) * prob_gap_longer

        if i == 0:
            data.expected_prev.append(expected_length)
        else:
            data.expected_next.append(expected_length)

    data.expected_gap.append(data.expected_prev[-1] + data.expected_next[-1])

    p_merit = 0
    for prob_i, lower in zip(prob_nth, unknowns[0]):
        for prob_j, upper in zip(prob_nth, unknowns[1]):
            prob_joint = prob_i * prob_j
            gap = -lower + upper
            if gap >= min_merit_gap:
                p_merit += prob_joint

            misc.prob_gap_comb[gap] += prob_joint
            if prob_joint < 1e-9:
                break
        else:
            # gap = (K + SL+1) - (K - lower) = SL+1 - lower
            if -lower + SL+1 > min_merit_gap:
                p_merit += prob_i * prob_longer[len(unknowns[1])]

    for prob_j, upper in zip(prob_nth, unknowns[1]):
        if SL+1 + upper > min_merit_gap:
            p_merit += prob_j * prob_longer[len(unknowns[0])]

    if 2*SL+1 > min_merit_gap:
        p_merit += prob_longer[len(unknowns[0])] * prob_longer[len(unknowns[1])]

    assert 0 <= p_merit <= 1.00, (p_merit, unknowns)
    data.prob_merit_gap.append(p_merit)


def setup_extended_gap(SL, K, P, D, prob_prime):
    K_primes = [p for p in range(2, P+1) if gmpy2.is_prime(p)]
    prob_prime_coprime = prob_prime
    for p in K_primes:
        # multiples of d are handled in prob_record_one_side
        prob_prime_coprime /= (1 - 1 / p)

    is_coprime = [True for i in range(2 * SL)]
    for p in K_primes:
        if D % p == 0:
            continue

        for i in range(0, len(is_coprime), p):
            is_coprime[i] = False

    coprime_x = [i for i, v in enumerate(is_coprime) if v and i > SL]
    return prob_prime_coprime, coprime_x


def prob_record_one_side(
        record_gaps, other_side,
        unknowns_high, prob_prime_after_sieve,
        m, K_mod_d, d, coprime_extended, prob_prime_coprime):
    assert other_side > 0
    prob_record = 0

    prob_nth = 1
    for unknown in unknowns_high:
        if unknown + other_side in record_gaps:
            prob_record += prob_nth * prob_prime_after_sieve
        prob_nth *= 1 - prob_prime_after_sieve

    # Use -m if wanted to consider unknowns_low
    for x in coprime_extended:
        if math.gcd(m * K_mod_d + x, d) > 1:
            continue

        if x + other_side in record_gaps:
            prob_record += prob_nth * prob_prime_coprime

        prob_nth *= 1 - prob_prime_coprime

    # Anything larger than 2*SL is likely record
    return prob_record + prob_nth


def determine_test_threshold(args, valid_mi, data):
    percent = args.prp_top_percent
    if not percent or percent == 100:
        return 0, {
            mi: (p_merit, p_record) for mi, p_merit, p_record in zip(
                valid_mi, data.prob_merit_gap, data.prob_record_gap)
        }

    assert 1 <= percent <= 99
    # Could be several million datapoints.
    # XXX: debug/print probs to see they match gap_stats
    best_probs = sorted(data.prob_record_gap, reverse=True)
    prob_threshold = best_probs[round(len(best_probs) * percent / 100)]
    return prob_threshold, {
        mi: (p_merit, p_record) for mi, p_merit, p_record in zip(
                valid_mi, data.prob_merit_gap, data.prob_record_gap)
        if p_record >= prob_threshold
    }


def should_print_stats(
        args, sc, data, mi, m,
        unknown_l, unknown_u,
        prev_p_i, next_p_i):

    stop_t = time.time()
    print_secs = stop_t - sc.last_print_t

    # Print a little bit if we resume but mostly as we test.
    if sc.tested in (1,10,30,100,300,1000,3000) or (sc.tested and sc.tested % 5000 == 0) \
            or m == data.last_m or print_secs > 1200:
        secs = stop_t - sc.start_t

        print("\t{:3d} {:4d} <- unknowns -> {:-4d}\t{:4d} <- gap -> {:-4d}".format(
            m,
            unknown_l, unknown_u,
            prev_p_i, next_p_i))
        if mi <= 10 and secs < 6:
            return False

        def roundSig(n, sig):
            return '{:g}'.format(float('{:.{p}g}'.format(n, p=sig)))

        def rate(count, time):
            # Want 3 sig figs which is hard in python
            if count and count < time:
                return "{} secs/test".format(roundSig(time / count, 3))
            else:
                return "{}/sec".format(roundSig(count / max(1e-5, time), 3))

        timing = rate(sc.tested, secs)
        timing_threads = rate(sc.tested, sc.test_time)

        # Stats!
        if sc.tested:
            print("\t    tests      {:<9d} ({}, {})  {:.0f}, {:.0f} secs".format(
                sc.tested, timing, timing_threads, secs, sc.test_time))

            print("\t    sum(prob_minmerit):  {:6.2g}, {:.3g}/day\tfound: {}".format(
                sc.prob_minmerit, 86400 / secs * sc.prob_minmerit, sc.count_minmerit))
            print("\t    sum(prob_record):    {:6.2g}, {:.3g}/day\tfound: {}".format(
                sc.prob_record, 86400  / secs * sc.prob_record, sc.count_record))

        total_unknown = sc.t_unk_low + sc.t_unk_hgh
        if total_unknown:
            print("\t    unknowns   {:<9d} (avg: {:.2f}), {:.2f}% composite  {:.2f}% <- % -> {:.2f}%".format(
                total_unknown, total_unknown / sc.tested,
                100 * (1 - total_unknown / ((2 * args.sieve_length + 1) * sc.tested)),
                100 * sc.t_unk_low / total_unknown,
                100 * sc.t_unk_hgh / total_unknown))
        prp_tests = sc.total_prp_tests
        if sc.tested and prp_tests:
            print("\t    prp tests  {:<9d} (avg: {:.2f}) ({:.3f} tests/sec)".format(
                prp_tests, prp_tests / sc.tested, prp_tests / secs))
            if True or sc.one_side_skips:
                print("\t    one side skips: {} ({:.1%}) (delta: {:.1} avg: {:.3}".format(
                    sc.one_side_skips, sc.one_side_skips / sc.tested,
                    sc.prob_one_side_delta, sc.prob_one_side_delta/sc.tested))
            if sc.gap_out_of_sieve_prev or sc.gap_out_of_sieve_next:
                print("\t    fallback prev_gap {} ({:.1%}), next_gap {} ({:.1%})".format(
                    sc.gap_out_of_sieve_prev, sc.gap_out_of_sieve_prev / sc.tested,
                    sc.gap_out_of_sieve_next, sc.gap_out_of_sieve_next / sc.tested))
            print("\t    merit {:.3f} (at m={})".format(
                sc.best_merit_interval, sc.best_merit_interval_m))

        return True
    return False


def determine_next_prime_i(m, strn, K, unknowns, SL):
    center = m * K
    tests = 0

    for i in unknowns:
        tests += 1;
        if gap_utils.is_prime(center + i, strn, i):
            return tests, i

    # XXX: parse to version and verify > 6.2.99
    assert gmpy2.mp_version() == 'GMP 6.2.99'
    # Double checks center + SL.
    next_p_i = int(gmpy2.next_prime(center + SL) - center)
    return tests, next_p_i


def determine_prev_prime_i(m, strn, K, unknowns, SL, primes, remainder):
    center = m * K
    tests = 0

    for i in unknowns:
        assert i < 0
        tests += 1;
        if gap_utils.is_prime(center + i, strn, i):
            return tests, -i

    # Double checks center + SL.
    # Medium ugly fallback.
    for i in range(SL, 5*SL+1):
        composite = False
        for prime, remain in zip(primes, remainder):
            modulo = (remain * m) % prime
            if i % prime == modulo:
                composite = True
                break
        if not composite:
            if gap_utils.is_prime(center - i, strn, -i):
                return tests, i

    assert False


def handle_result(
        args, record_gaps, data, sc,
        mi, m, log_n, prev_p_i, next_p_i, unknown_l, unknown_u):
    '''Called for existing and new records'''
    assert prev_p_i > 0, (mi, m, next_p_i, prev_p_i)

    gap = next_p_i + prev_p_i
    data.experimental_gap.append(gap)
    data.valid_m.append(m)
    data.experimental_side.append(next_p_i)
    data.experimental_side.append(prev_p_i)

    merit = gap / log_n
    if gap in record_gaps or merit > args.min_merit:
        print("{}  {:.4f}  {} * {}#/{} -{} to +{}".format(
            gap, merit, m, args.p, args.d, prev_p_i, next_p_i))

    if merit > sc.best_merit_interval:
        sc.best_merit_interval = merit
        sc.best_merit_interval_m = m

    if should_print_stats(
            args, sc, data,
            mi, m,
            unknown_l, unknown_u,
            prev_p_i, next_p_i,
            ):
        sc.best_merit_interval = 0
        sc.best_merit_interval_m = -1
        sc.last_print_t = time.time()


# NOTE: Manager is prime_gap_test maintains workers wh run process_line.
def process_line(
        done_flag, work_q, results_q, thread_i,
        SL, K, P, D,
        prob_prime, prob_prime_after_sieve, record_gaps, one_side_skip_threshold,
        primes, remainder):
    def cleanup(*_):
        print (f"Thread {thread_i} stopping")
        work_q.close()
        work_q.join_thread()
        results_q.close()
        results_q.join_thread()
        exit(0)

    # Calculate extended gap (see gap_stats)
    prob_prime_coprime, coprime_extended = setup_extended_gap(SL, K, P, D, prob_prime)
    K_mod_d = K % D

    # Ignore KeyboardInterrupt and let Manager terminate us.
    signal.signal(signal.SIGINT, cleanup)

    print (f"Thread {thread_i} started")
    while True:
        work = work_q.get()
        if work == "STOP":
            assert done_flag.is_set()
            cleanup()
            return

        m, mi, prob_record, log_n, line = work

        mtest, unknown_l, unknown_u, unknowns = gap_utils.parse_unknown_line(line)
        assert mi == mtest

        # Used for openPFGW
        strn = "{}*{}#/{}+".format(m, P, D)

        t0 = time.time()

        p_tests, prev_p_i = determine_prev_prime_i(m, strn, K, unknowns[0], SL, primes, remainder)

        test_next = True
        if one_side_skip_threshold:
            # Check if slower to test a new record or this one

            new_prob_record = prob_record_one_side(
                    record_gaps, prev_p_i,
                    unknowns[1], prob_prime_after_sieve,
                    m, K_mod_d, D, coprime_extended, prob_prime_coprime)

            #print ("{:5} | m: {:8} | count: {} {} | {:.3g} => {:.3g} | {:.3g}".format(
            #        prev_p_i, m, unknown_l, unknown_u,
            #        prob_record, new_prob_record,
            #        one_side_skip_threshold))

            # XXX: this increases effective prob/test (recursievly)
            #      which is hard to factor here.
            # Takes 2x (see XXX note above) time to test a new interval.
            #       (2 - skip_percent)/(tested_prob) x time to test a new interval
            if 3 * new_prob_record < one_side_skip_threshold:
                test_next = False

            #if new_prob_record > 10 * one_side_skip_threshold:
            #    print ("{:5} | m: {:8} | count: {} {} | {:.3g} => {:.3g}".format(
            #        prev_p_i, m, unknown_l, unknown_u,
            #        prob_record, new_prob_record))

            # On average sum(new_prob_record) ~= sum(prob_record)
        else:
            new_prob_record = 0

        n_tests, next_p_i = 0, 0
        if test_next:
            n_tests, next_p_i = determine_next_prime_i(m, strn, K, unknowns[1], SL)

        test_time = time.time() - t0

        results_q.put((
            m, mi, log_n,
            unknown_l, unknown_u,
            n_tests, next_p_i,
            p_tests, prev_p_i,
            prob_record, new_prob_record,
            test_time,
        ))


def process_result(conn, args, record_gaps, mi_probs, data, sc, result):
    ''' Handles new results '''

    (m, mi, r_log_n, unknown_l, unknown_u,
     n_tests, next_p_i,
     p_tests, prev_p_i,
     prob_record, one_sided_prob_record,
     test_time) = result

    sc.tested += 1
    sc.test_time += test_time

    sc.t_unk_low += unknown_l
    sc.t_unk_hgh += unknown_u

    gap = next_p_i + prev_p_i
    merit = gap / r_log_n

    sc.prob_minmerit           += mi_probs[mi][0]
    if merit > args.min_merit:
        sc.count_minmerit += 1

    if gap in record_gaps:
        sc.count_record += 1

    if next_p_i:
        assert mi_probs[mi][1] == prob_record
        sc.prob_record += prob_record
    else:
        sc.prob_record += prob_record - one_sided_prob_record
        sc.one_side_skips += 1

    sc.prob_one_side_delta += prob_record - one_sided_prob_record

    sc.total_prp_tests += p_tests
    sc.gap_out_of_sieve_prev += prev_p_i > args.sieve_length

    sc.total_prp_tests += n_tests
    sc.gap_out_of_sieve_next += next_p_i > args.sieve_length

    save(conn,
        args.p, args.d, m, next_p_i, prev_p_i, merit,
        n_tests, p_tests, test_time)

    handle_result(
        args, record_gaps, data, sc,
        mi, m, r_log_n, prev_p_i, next_p_i, unknown_l, unknown_u)


def run_in_parallel(
        args, conn, unknown_file, record_gaps,
        prob_nth, prob_prime, prob_prime_after_sieve,
        existing, valid_mi,
        K, K_log,
        data, sc, misc
):
    # XXX: Cleanup after gmpy2.prev_prime.
    # Remainders of (p#/d) mod prime
    primes = tuple([p for p in range(3, 80000+1) if gmpy2.is_prime(p)])
    remainder = tuple([K % prime for prime in primes])

    # Based on args
    prob_threshold, mi_probs = determine_test_threshold(args, valid_mi, data)
    if prob_threshold >= 0:
        print ("Testing {} m where prob(record) >= {:.3g}".format(
            len(mi_probs), prob_threshold))
    # XXX: how to set one_side_skip when prp_top_percent == 100

    # Any non-processed mi_probs?
    if all(args.mstart + mi in existing for mi in mi_probs):
        print(f"All prp-top-percent({len(mi_probs)}) already processed!")
        print()
        return

    # Worker setup
    done_flag = multiprocessing.Event()
    work_q    = multiprocessing.Queue(40)
    results_q = multiprocessing.Queue()

    one_sided_args = (
        prob_prime,
        prob_prime_after_sieve,
        record_gaps,
        0 if args.no_one_side_skip else prob_threshold,
    )

    assert args.threads in range(1, 65), args.threads
    processes = []
    for i in range(args.threads):
        process = multiprocessing.Process(
            target=process_line,
            args=(
                done_flag, work_q, results_q,
                i, args.sieve_length, K, args.p, args.d,
                *one_sided_args,
                primes, remainder
            )
        )
        process.start()
        time.sleep(0.01)
        processes.append(process)
    print()
    time.sleep(0.5)

    for mi in valid_mi:
        m = args.mstart + mi
        log_n = (K_log + math.log(m))

        # Read a line from the file
        line = unknown_file.readline()

        if len(data.valid_m) <= 3:
            _, _, _, unknowns = gap_utils.parse_unknown_line(line)
            misc.test_unknowns.append(unknowns)


        if args.stats and (args.save_logs or args.num_plots):
            # NOTES: calculate_expected_gaps is really slow,
            # only real use is to doublecheck gap_stats.cpp
            _, _, _, unknowns = gap_utils.parse_unknown_line(line)
            calculate_expected_gaps(
                SL, min_merit_gap, prob_nth, prob_longer,
                log_n, unknowns,
                # Results saved to data / misc
                data, misc)
            mi_probs[mi] = (data.prob_merit_gap[-1], data.prob_record_gap[-1])

        # Check if this prob_record high enough to run
        if mi not in mi_probs or mi_probs[mi][1] < prob_threshold:
            continue

        if m in existing:
            # XXX: missing_gap_test only sets one side.
            prev_p_i, next_p_i = existing[m]
            handle_result(
                args, record_gaps, data, sc,
                mi, m, log_n, prev_p_i, next_p_i, 0, 0)

        else:
            sc.will_test += 1
            try:
               work_q.put((m, mi, mi_probs[mi][1], log_n, line), True)
            except KeyboardInterrupt:
                print(" Breaking from Ctrl+C ^C")
                for p in processes:
                    p.terminate()
                return

        # Process any finished results
        while not results_q.empty():
            result = results_q.get(block=False)
            process_result(conn, args, record_gaps, mi_probs, data, sc, result)


    print("Everything Queued, done.set() & pushing STOP")
    done_flag.set()

    # Push "STOP" for every worker
    for i in range(args.threads):
        work_q.put("STOP")
    work_q.close()

    while sc.tested < sc.will_test:
        print(f"Waiting on {sc.will_test - sc.tested} of {sc.will_test} results")
        result = results_q.get(block=True)
        process_result(conn, args, record_gaps, mi_probs, data, sc, result)

    print("Joining work_q")
    work_q.join_thread()
    time.sleep(0.2)

    for i, process in enumerate(processes):
        print(f"Joining process({i})")
        process.join(timeout=0.1)
    print ("Done!")


def prime_gap_test(args):
    P = args.p
    D = args.d

    M = args.mstart
    M_inc = args.minc

    SL = args.sieve_length
    max_prime = args.max_prime

    K, K_digits, K_bits, K_log = gap_utils.K_and_stats(args)
    M_log    = K_log + math.log(M)
    min_merit_gap = int(args.min_merit * M_log)
    print("K = {} bits, {} digits, log(K) = {:.2f}".format(
        K_bits, K_digits, K_log))
    print("Min Gap ~= {} (for merit > {:.1f})\n".format(
        min_merit_gap, args.min_merit))

    # ----- Load Record sized gaps

    with sqlite3.connect(args.prime_gaps_db) as conn:
        record_gaps = load_records(conn, M_log)
        print ("\tLoaded {} records ({} to {}) from {!r}".format(
            len(record_gaps), record_gaps[0], record_gaps[-1], args.prime_gaps_db))

        # Calculate min_merit record, min_merit for 50% record, min_merit for 90% record
        min_record = record_gaps[0] / M_log
        # Look 250 records later and check if contains more than half of numbers
        avg_record = next(record_gaps[r] for r in range(len(record_gaps) - 250)
                if record_gaps[r + 250] - record_gaps[r] < 1000)
        print ("\tMin merit for record {:.2f}, gapsize for likely record {:.2f} {}".format(
            min_record, avg_record / M_log, avg_record))

        record_gaps = set(record_gaps)

    # ----- Open Output file
    print("\tLoading unknowns from {!r}".format(args.unknown_filename))
    print()
    folder_unk = gap_utils.transform_unknown_filename(
            args.unknown_filename, "unknowns", "txt")
    if os.path.exists(folder_unk):
        unknown_file = open(folder_unk, "r")
    else:
        unknown_file = open(args.unknown_filename, "r")

    # ----- Open Prime-Gap-Search DB
    # Longer timeout so that record_checking doesn't break saving
    conn = sqlite3.connect(args.search_db, timeout=15)
    conn.row_factory = sqlite3.Row
    existing = load_existing(conn, args)
    print (f"Found {len(existing)} existing results")

    # used in next_prime
    assert P <= 80000
    # Very slow but works.
    P_primes = [p for p in range(2, P+1) if gmpy2.is_prime(p)]

    # ----- Allocate memory for a handful of utility functions.

    # ----- Sieve stats
    prob_prime = 1 / M_log - 1 / (M_log * M_log)
    prob_prime_after_sieve = prob_prime_sieve_length(
        M, K, D, prob_prime, K_digits, P_primes, SL, max_prime)

    # Geometric distribution
    prob_nth = []
    prob_longer = []
    prob_gap_longer = 1
    while prob_gap_longer > 1e-13:
        prob_nth.append(prob_gap_longer * prob_prime_after_sieve)
        prob_longer.append(prob_gap_longer)
        prob_gap_longer *= (1 - prob_prime_after_sieve)
    assert min(prob_nth) > 0
    assert min(prob_longer) > 0
    #print (f"|prob_nth| = {len(prob_nth)}")


    # ----- Main sieve loop.

    @dataclass
    class StatCounters:
        start_t: float
        last_print_t: float
        will_test = 0
        tested = 0

        test_time = 0
        t_unk_low = 0
        t_unk_hgh = 0
        total_prp_tests = 0
        one_side_skips = 0
        gap_out_of_sieve_prev = 0
        gap_out_of_sieve_next = 0

        best_merit_interval = 0
        best_merit_interval_m = -1

        prob_record = 0.0
        prob_minmerit = 0.0

        # To validate one sided prob is working
        prob_one_side_delta = 0.0

        count_record = 0
        count_minmerit = 0


    sc = StatCounters(time.time(), time.time())
    data = Datas()
    misc = Misc()

    valid_mi = [mi for mi in range(M_inc) if math.gcd(M + mi, D) == 1]
    data.first_m = M + valid_mi[0]
    data.last_m = M + valid_mi[-1]

    if len(existing) == len(valid_mi):
        print(f"All processed!")
    else:
        print(f"\nStarting m({len(valid_mi)}) {data.first_m} to {data.last_m}")
        print()

        # Load stats for prob_record
        if not args.stats and (args.num_plots or args.prp_top_percent):
            data_db, misc_db = load_stats(conn, args)
            assert len(data_db.prob_merit_gap)  == len(valid_mi), "run ./gap_stats first"
            assert len(data_db.prob_record_gap) == len(valid_mi), "run ./gap_stats first"
            data.prob_merit_gap = data_db.prob_merit_gap
            data.prob_record_gap = data_db.prob_record_gap

        run_in_parallel(
            args, conn, unknown_file, record_gaps,
            prob_nth, prob_prime, prob_prime_after_sieve,
            existing, valid_mi,
            K, K_log,
            data, sc, misc
        )

    # ----- Plots
    if args.num_plots:
        gap_test_plotting.plot_stuff(
            args, conn, data, sc, misc,
            min_merit_gap, record_gaps, prob_nth)



if __name__ == "__main__":
    parser = get_arg_parser()
    args = parser.parse_args()
    gap_utils.verify_args(args)

    # TeeLogger context if args.save_logs
    with gap_utils.logger_context(args):
        prime_gap_test(args)

