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
import queue
from collections import defaultdict
from dataclasses import dataclass

import gmpy2

import gap_utils
import gap_utils_primes
import gap_test_plotting
import misc.misc_utils as misc_utils


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

    prob_record_all = 0.0
    prob_minmerit = 0.0

    # handles one sided testes
    prob_record_processed = 0.0
    sides_processed = 0.0
    # To validate one sided prob is working
    prob_one_side_delta = 0.0

    count_record = 0
    count_minmerit = 0


@dataclass
class GapData:
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


def config_hash(config):
    """Shared with gap_stats to load/save range rid"""
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
        "SELECT m, next_p, prev_p FROM result"
        " WHERE P = ? AND D = ? AND m BETWEEN ? AND ?",
        (args.p, args.d, args.mstart, args.mstart + args.minc - 1))
    return {row['m']: [row['prev_p'], row['next_p']] for row in rv}


def load_stats(conn, args):
    data = GapData()
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


def save(conn, p, d, m, next_p, prev_p, merit,
         n_tests, p_tests, test_time):
    assert p in range(201, 80000, 2), (p, d, m)
    conn.execute(
        "REPLACE INTO result(P,D,m,next_p,prev_p,merit)"
        "VALUES(?,?,?,  ?,?,  ?)",
        (p, d, m, next_p, prev_p,  round(merit,4)))

    conn.execute(
        """UPDATE m_stats
        SET next_p=?, prev_p=?, merit=?,
            prp_next=prp_next+?, prp_prev=prp_prev+?,
            test_time=test_time+?
        WHERE p=? AND d=? AND m=?""",
        (next_p, prev_p, round(merit,4),
         n_tests, p_tests, test_time, p, d, m))
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

    K_coprime = [math.gcd(K, i) == 1 for i in range(SL+1)]

    # See "Optimizing Choice Of D" in THEORY.md for why this is required
    count_coprime_p = 0
    m_tests = [M+i for i in range(100) if math.gcd(M+i, D) == 1][:6]
    for m in m_tests:
        N_mod_D = m * K % D
        for i in range(1, SL+1):
            if K_coprime[i] and math.gcd(N_mod_D + i, D) == 1:
                # assert gmpy2.gcd(N, K*D) == 1,
                count_coprime_p += 1

    count_coprime_p //= len(m_tests)

    chance_coprime_composite = 1 - prob_prime / prob_prime_coprime_p
    prob_gap_shorter_hypothetical = chance_coprime_composite ** count_coprime_p

    expected = count_coprime_p * (unknowns_after_sieve / prob_prime_coprime_p)
    print("\texpect {:.0f} left, {:.3%} of SL={} after {}M".format(
        expected, expected / (SL + 1), SL, max_prime//10 ** 6))
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
    """
    Extends from SL to 2*SL based on coprime K

    Used in prob_record_one_sided
    """
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


def prob_record_one_sided(
        record_gaps, megagap, other_side,
        unknowns_high, prob_prime_after_sieve,
        m, K_mod_d, d, coprime_extended, prob_prime_coprime):
    "Probability of record given one side"
    assert other_side > 0
    prob_record = 0

    prob_nth = 1
    for unknown in unknowns_high:
        test_gap = unknown + other_side
        if (megagap and test_gap > megagap) or (not megagap and test_gap in record_gaps):
            prob_record += prob_nth * prob_prime_after_sieve
        prob_nth *= 1 - prob_prime_after_sieve

    # Use -m if wanted to consider unknowns_low
    for x in coprime_extended:
        test_gap = x + other_side
        if math.gcd(m * K_mod_d + x, d) > 1:
            if (megagap and test_gap > megagap) or (not megagap and test_gap in record_gaps):
                prob_record += prob_nth * prob_prime_coprime

            prob_nth *= 1 - prob_prime_coprime

    # Anything larger than 2*SL is likely record
    return prob_record + prob_nth


def determine_test_threshold(args, valid_mi, data):
    percent = args.prp_top_percent
    valid_m = (args.mstart + mi for mi in valid_mi)
    if not percent or percent == 100:
        return 0, {
            m: (p_merit, p_record) for m, p_merit, p_record in zip(
                valid_m, data.prob_merit_gap, data.prob_record_gap)
        }

    if args.megagap:
        data.prob_record_gap = data.prob_merit_gap

    assert 1 <= percent <= 99
    # Could be several million datapoints.
    best_probs = sorted(data.prob_record_gap, reverse=True)
    prob_threshold = best_probs[round(len(best_probs) * percent / 100)]
    return prob_threshold, {
        m: (p_merit, p_record) for m, p_merit, p_record in zip(
                valid_m, data.prob_merit_gap, data.prob_record_gap)
        if p_record >= prob_threshold
    }


def should_print_stats(
        args, sc, data, mi, m,
        unknown_l, unknown_u,
        prev_p, next_p):

    stop_t = time.time()
    print_secs = stop_t - sc.last_print_t

    # Print a little bit if we resume but mostly as we test.
    if (args.megagap and sc.tested % 10 == 0) or sc.tested in (1,10,30,100,300,1000,3000) \
            or (sc.tested and sc.tested % 5000 == 0) \
            or m == data.last_m or print_secs > 1200:
        secs = stop_t - sc.start_t

        print("\t{:3d} {:4d} <- unknowns -> {:-4d}\t{:4d} <- gap -> {:-4d}".format(
            m, unknown_l, unknown_u, prev_p, next_p))
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

        total_unknown = sc.t_unk_low + sc.t_unk_hgh
        if total_unknown:
            print("\t    unknowns   {:<9d} (avg: {:.2f}), {:.2f}% composite  {:.2f}% <- % -> {:.2f}%".format(
                total_unknown, total_unknown / sc.tested,
                100 * (1 - total_unknown / ((2 * args.sieve_length + 1) * sc.tested)),
                100 * sc.t_unk_low / total_unknown,
                100 * sc.t_unk_hgh / total_unknown))
        prp_tests = sc.total_prp_tests
        if sc.tested and prp_tests:
            print("\t    prp tests  {:<9d} (avg: {:.2f}/side, {:.2f}/m) ({:.3f} tests/sec)".format(
                prp_tests, prp_tests / sc.sides_processed, prp_tests / sc.tested, prp_tests / secs))
            if sc.one_side_skips:
                avg_prob_side = sc.prob_record_processed / sc.sides_processed
                avg_prob_no_skip = sc.prob_record_all / (2 * sc.tested)
                print("\t    side skips {:<9d} ({:.1%}) (avg prob_record {:.3g}/side, {:.3g}/m/2, {:+.1%})".format(
                    sc.one_side_skips, sc.one_side_skips / sc.tested,
                    avg_prob_side, avg_prob_no_skip, avg_prob_side / avg_prob_no_skip - 1))

            print("\t    sum(prob_minmerit):  {:7.5g}, {:.3g}/day\tfound: {}".format(
                sc.prob_minmerit, 86400 / secs * sc.prob_minmerit, sc.count_minmerit))
            print("\t    sum(prob_record):    {:7.5g}, {:.3g}/day\tfound: {}".format(
                sc.prob_record_processed, sc.prob_record_processed / secs * 86400, sc.count_record))

            if sc.gap_out_of_sieve_prev or sc.gap_out_of_sieve_next:
                print("\t    fallback prev_gap {} ({:.1%}), next_gap {} ({:.1%})".format(
                    sc.gap_out_of_sieve_prev, sc.gap_out_of_sieve_prev / sc.tested,
                    sc.gap_out_of_sieve_next, sc.gap_out_of_sieve_next / sc.tested))
            print("\t    merit      {:<6.3f}    (at m={})".format(
                sc.best_merit_interval, sc.best_merit_interval_m))

        return True
    return False


def handle_result_for_plots(args, data, mi, m, prev_p, next_p):
    '''Called for existing and new records'''
    assert prev_p > 0, (mi, m, next_p, prev_p)

    if args.num_plots:
        gap = next_p + prev_p
        data.experimental_gap.append(gap)
        data.valid_m.append(m)
        data.experimental_side.append(next_p)
        data.experimental_side.append(prev_p)


def process_result(conn, args, record_gaps, m_probs, data, sc, result):
    ''' Handles new results '''
    (m, mi, r_log_n, unknown_l, unknown_u,
     n_tests, next_p,
     p_tests, prev_p,
     prob_record, new_prob_record,
     test_time) = result

    gap = next_p + prev_p
    merit = gap / r_log_n

    save(conn,
        args.p, args.d, m, next_p, prev_p, merit,
        n_tests, p_tests, test_time)

    if next_p < 0 or prev_p < 0:
        # partial result don't print anything
        return

    sc.tested += 1
    sc.test_time += test_time

    sc.t_unk_low += unknown_l
    sc.t_unk_hgh += unknown_u

    sc.prob_minmerit += m_probs[m][0]
    if merit > args.min_merit:
        sc.count_minmerit += 1

    if gap in record_gaps:
        sc.count_record += 1

    assert m_probs[m][1] == prob_record
    sc.prob_record_all += prob_record

    # Should mirror logic in process_line
    if next_p:
        sc.prob_record_processed += prob_record
        sc.sides_processed += 2
    else:
        sc.one_side_skips += 1
        sc.prob_record_processed += prob_record - new_prob_record
        sc.sides_processed += 1

    sc.prob_one_side_delta += prob_record - new_prob_record

    sc.total_prp_tests += p_tests
    sc.gap_out_of_sieve_prev += prev_p > args.sieve_length

    sc.total_prp_tests += n_tests
    sc.gap_out_of_sieve_next += next_p > args.sieve_length

    handle_result_for_plots(args, data, mi, m, prev_p, next_p)

    is_record = gap in record_gaps
    if is_record or merit > args.min_merit:
        print("{}  {:.4f}  {} * {}#/{} -{} to +{}{}".format(
            gap, merit, m, args.p, args.d, prev_p, next_p,
            "\tRECORD!" if is_record else ""))

    if merit > sc.best_merit_interval:
        sc.best_merit_interval = merit
        sc.best_merit_interval_m = m

    if should_print_stats(
            args, sc, data,
            mi, m,
            unknown_l, unknown_u,
            prev_p, next_p,
            ):
        sc.best_merit_interval = 0
        sc.best_merit_interval_m = -1
        sc.last_print_t = time.time()
