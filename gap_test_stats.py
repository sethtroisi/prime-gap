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

import array
import itertools
import math
import random
import time
from collections import defaultdict
from dataclasses import dataclass

import gmpy2
import sqlite3

import gap_utils


@dataclass
class StatCounters:
    start_t: float
    last_print_t: float
    will_test = 0
    tested = 0

    test_time = 0
    # Total unknowns on previous / next side
    t_unk_prev = 0
    t_unk_next = 0
    # Total PRP tests performed
    total_prp_tests = 0
    # Number of times one side was skipped
    one_side_skips = 0
    # Number of prev_prime / next_prime that exceeded SIEVE_LENGTH
    gap_out_of_sieve_prev = 0
    gap_out_of_sieve_next = 0

    # Highest merit seen recently, m of that merit
    best_merit_interval = 0
    best_merit_interval_m = -1

    prob_record_all = 0.0
    prob_min_merit = 0.0

    # handles one sided testes
    prob_record_processed = 0.0
    sides_processed = 0.0

    # To validate one sided prob is working
    # TODO: validate in some way (value trends to zero?)
    prob_one_side_delta = 0.0

    count_record = 0
    count_min_merit = 0


@dataclass
class GapData:
    first_m = -1
    last_m = -1

    valid_mi = []

    # Values for each m in valid_mi
    expected_prev = []
    expected_next = []
    expected_gap = []
    prob_merit_gap = []
    prob_record_gap = []

    # Because of --prp-top-percent may be less than |valid_mi|
    experimental_side = []
    experimental_gap = []


@dataclass
class Misc:
    prob_gap_side = defaultdict(float)
    prob_gap_comb = defaultdict(float)

    test_unknowns = {}


def config_hash(config):
    """Shared with gap_stats to load/save range rid"""
    h = config.mstart
    h = h * 31 + config.minc
    h = h * 31 + config.p
    h = h * 31 + config.d
    h = h * 31 + config.sieve_length
    h = h * 31 + config.max_prime
    return h % (2 ** 64)


def load_range(conn, args):
    """Load range data for args"""
    cursor = conn.execute(
            "SELECT * FROM range WHERE rid = ?",
            (config_hash(args),))
    # Set row factor for just this cursor
    cursor.row_factory = sqlite3.Row
    return cursor.fetchone()

def load_records(conn, log_N):
    # Find records merit <= gap * log_N
    # Stop at merit = 35

    # NOTE: to make sure 'count(records)' still counts after a
    # found record is submitted subtract a small amount from merit
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
        all_gaps = range(SMALLEST_MISSING, largest_likely + 1, 2)
        missing = list(set(all_gaps) - set(g[0] for g in rv))
        records.extend(missing)

    return tuple(sorted(records))


def load_existing(conn, args):
    rv = conn.execute(
        "SELECT m, prev_p, next_p FROM result"
        " WHERE P = ? AND D = ? AND m BETWEEN ? AND ?",
        (args.p, args.d, args.mstart, args.mstart + args.minc - 1))
    return {row[0]: (row[1], row[2]) for row in rv}


def _zip_to_array(cursor):
    """convert rows to columns"""
    init = True
    for row in cursor:
        if init:
            init = False
            arrays = [array.array('f') for _ in row]

        for i, v in enumerate(row):
            arrays[i].append(v)

    if init:
        exit("Empty m_stats for this range, have you run gap_stats?")

    return arrays


def load_probs_only(conn, args):
    data = GapData()

    rv = conn.execute(
        "SELECT prob_merit, prob_record"
        " FROM m_stats WHERE P = ? AND D = ? AND m BETWEEN ? AND ?",
        (args.p, args.d, args.mstart, args.mstart + args.minc - 1))

    data.prob_merit_gap, data.prob_record_gap = _zip_to_array(rv)

    return data


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
     gap_prev, gap_next,
     data.prob_merit_gap, data.prob_record_gap) = _zip_to_array(rv)

    interleaved = itertools.chain(*itertools.zip_longest(gap_prev, gap_next))
    data.experimental_side = array.array('i', (int(s) for s in interleaved if s and s > 0))
    zipped_sum_gen = (int(p + n) for p, n in zip(gap_prev, gap_next) if p > 0 and n > 0)
    data.experimental_gap = array.array('i', zipped_sum_gen)

    rv = conn.execute(
        "SELECT gap, prob_combined, prob_low_side, prob_high_side "
        "FROM range_stats where rid = ?",
        (config_hash(args),))
    for row in rv:
        gap = row[0]
        misc.prob_gap_comb[gap] += row[1]
        misc.prob_gap_side[gap] += row[2] / 2
        misc.prob_gap_side[gap] += row[3] / 2

    return data, misc


def save(conn, p, d, m, next_p, prev_p, merit,
         n_tests, p_tests, test_time):
    assert p in range(201, 80000, 2), (p, d, m)
    conn.execute(
        "REPLACE INTO result(P,D,m,next_p,prev_p,merit)"
        "VALUES(?,?,?,  ?,?,  ?)",
        (p, d, m, next_p, prev_p, round(merit, 4)))

    conn.execute(
        """UPDATE m_stats
        SET next_p=?, prev_p=?, merit=?,
            prp_next=prp_next+?, prp_prev=prp_prev+?,
            test_time=test_time+?
        WHERE p=? AND d=? AND m=?""",
        (next_p, prev_p, round(merit, 4),
         n_tests, p_tests, test_time, p, d, m))

    # Commit on average every 200 seconds (always commit if test_time > 200).
    # Keeps number of raw commits down when testing small p#
    # conn.commit() in gap_test.py should prevent losing results during cleanup
    if 200 * random.random() < test_time:
        conn.commit()


def prob_prime_sieve_length(M, K, D, P, prob_prime, K_digits, SL, max_prime):
    assert max_prime >= 10 ** 6, max_prime

    # From Mertens' 3rd theorem
    gamma = 0.577215665
    unknowns_after_sieve = 1.0 / (math.log(max_prime) * math.exp(gamma))

    prob_prime_after_sieve = prob_prime / unknowns_after_sieve

    P_primes = [p for p in range(2, P+1) if gmpy2.is_prime(p)]

    prob_prime_coprime_p = 1
    for prime in P_primes:
        prob_prime_coprime_p *= (1 - 1 / prime)

    K_coprime = [math.gcd(K, i) == 1 for i in range(SL + 1)]

    # See "Optimizing Choice Of D" in THEORY.md for why this is required
    count_coprime_p = 0
    m_tests = [M + i for i in range(100) if math.gcd(M + i, D) == 1][:6]
    for m in m_tests:
        N_mod_D = m * K % D
        for i in range(1, SL + 1):
            if K_coprime[i] and math.gcd(N_mod_D + i, D) == 1:
                # assert gmpy2.gcd(N, K*D) == 1,
                count_coprime_p += 1

    count_coprime_p //= len(m_tests)

    chance_coprime_composite = 1 - prob_prime / prob_prime_coprime_p
    prob_gap_shorter_hypothetical = chance_coprime_composite ** count_coprime_p

    # Printed by C code, less interesting here
    #print("\t{:.3%} of {} digit numbers are prime".format(
    #    prob_prime, K_digits))
    #print("\t{:.3%} of 0tests should be prime ({:.1f}x speedup)".format(
    #    prob_prime_after_sieve, 1 / unknowns_after_sieve))
    print("\t~2x{:.1f} = {:.1f} PRP tests per m".format(
        1 / prob_prime_after_sieve, 2 / prob_prime_after_sieve))
    print("\tsieve_length={} is insufficient ~{:.2%} of time".format(
        SL, prob_gap_shorter_hypothetical))
    print()

    return prob_prime_after_sieve


def calculate_expected_gaps(
        SL, min_merit_gap, prob_nth, prob_longer, log_n, unknowns, data):
    for i, side in enumerate(unknowns):
        expected_length = 0
        for v, prob in zip(side, prob_nth):
            expected_length += abs(v) * prob

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

            if prob_joint < 1e-9:
                break
        else:
            # gap = (K + SL+1) - (K - lower) = SL+1 - lower
            if -lower + SL + 1 > min_merit_gap:
                p_merit += prob_i * prob_longer[len(unknowns[1])]

    for prob_j, upper in zip(prob_nth, unknowns[1]):
        if SL + 1 + upper > min_merit_gap:
            p_merit += prob_j * prob_longer[len(unknowns[0])]

    if 2 * SL + 1 > min_merit_gap:
        p_merit += prob_longer[len(unknowns[0])] * prob_longer[len(unknowns[1])]

    assert 0 <= p_merit <= 1.00, (p_merit, unknowns)
    data.prob_merit_gap.append(p_merit)


def validate_prob_record_merit(
        args, data_db, K_log, prob_prime_after_sieve):
    """Validate gap_stats with slow python calculation"""
    print("\tSLOWLY calculating prob_merit_gap for validation purpose")

    import scipy.stats

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

    data = GapData()
    with open(args.unknown_filename, "rb") as unknown_file_repeat:
        for line in unknown_file_repeat.readlines():
            m, _, _, unknowns = gap_utils._parse_unknown_line(line)
            log_n = K_log + math.log(m)
            min_merit_gap = int(args.min_merit * log_n)

            calculate_expected_gaps(
                args.sieve_length, min_merit_gap, prob_nth, prob_longer,
                log_n, unknowns,
                # Results saved to data / misc
                data)

    def lin_reg(name, a, b):
        # p_val is "two-sided p-val for null hypothesis that the slope is zero."
        slope, intercept, r_val, p_val, std_err = scipy.stats.linregress(a, b)
        print("\t{:15} | scale: {:.3f}, r^2: {:.2f}".format(name, slope, r_val))
        print("\t\t", p_val, std_err)

    # Validate several distributions have high correlation & r^2
    lin_reg("expected gap ", data_db.expected_gap, data.expected_gap)
    lin_reg("expected prev", data_db.expected_prev, data.expected_prev)
    lin_reg("expected next", data_db.expected_next, data.expected_next)

    lin_reg("expected prev", data_db.prob_merit_gap, data.prob_merit_gap)


def _sieve_side(SL, primes):
    is_coprime = [True for _ in range(SL+1)]
    for p in primes:
        for i in range(0, len(is_coprime), p):
            is_coprime[i] = False

    return is_coprime


def setup_extended_gap(SL, P, D, prob_prime):
    """
    Extends from SL to 2*SL based on coprime K

    Used in prob_record_one_sided
    """
    P_primes = [p for p in range(2, P + 1) if gmpy2.is_prime(p)]
    prob_prime_coprime = prob_prime
    for p in P_primes:
        # multiples of d are handled in prob_record_one_side
        prob_prime_coprime /= (1 - 1 / p)

    is_coprime = _sieve_side(2 * SL, (p for p in P_primes if D % p != 0))
    for p in P_primes:
        if D % p == 0:
            continue

        for i in range(0, len(is_coprime), p):
            is_coprime[i] = False

    coprime_X = [i for i, v in enumerate(is_coprime) if v and i > SL]
    return prob_prime_coprime, coprime_X


def setup_bitarray(SL, P, D):
    """
    Extends from [-SL, SL] based on coprime K

    Used in parse_compressed_line for --bitcompress'ed lines
    """
    D_primes = [p for p in range(2, P + 1) if D % p == 0 and gmpy2.is_prime(p)]
    K_primes = [p for p in range(2, P + 1) if gmpy2.is_prime(p) and D % p != 0]

    # [0, SL]
    is_coprime = _sieve_side(SL, K_primes)
    # [-Sl, SL]
    is_offset_coprime= is_coprime[::-1] + is_coprime[1:]
    coprime_X = [x for x, v in zip(range(-SL, SL+1), is_offset_coprime) if v]

    return D_primes, is_offset_coprime, coprime_X


def prob_record_one_sided(
        record_gaps, megagap, other_side,
        unknowns_next, prob_prime_after_sieve,
        m, K_mod_d, d, coprime_extended, prob_prime_coprime):
    """Probability of record given one side"""
    assert other_side > 0
    prob_record = 0

    prob_nth = 1
    for unknown in unknowns_next:
        test_gap = unknown + other_side
        if (megagap and test_gap > megagap) or (not megagap and test_gap in record_gaps):
            prob_record += prob_nth * prob_prime_after_sieve
        prob_nth *= 1 - prob_prime_after_sieve

    for x in coprime_extended:
        test_gap = x + other_side
        # Would use -m if wanted to consider unknowns_prev
        if math.gcd(m * K_mod_d + x, d) > 1:
            if (megagap and test_gap > megagap) or (not megagap and test_gap in record_gaps):
                prob_record += prob_nth * prob_prime_coprime

            prob_nth *= 1 - prob_prime_coprime

    # Anything larger than 2*SL is likely record
    return prob_record + prob_nth


def determine_test_threshold(args, data):
    percent = args.prp_top_percent
    valid_m = (args.mstart + mi for mi in data.valid_mi)
    if not percent or percent == 100:
        return 0, {
            m: (p_merit, p_record) for m, p_merit, p_record in zip(
                valid_m, data.prob_merit_gap, data.prob_record_gap)
        }

    if args.megagap:
        data.prob_record_gap = data.prob_merit_gap

    assert 1 <= percent <= 99
    # Could be several million data points.
    best_probs = sorted(data.prob_record_gap, reverse=True)
    prob_threshold = best_probs[round(len(best_probs) * percent / 100)]
    del best_probs

    # XXX: Could be array of (m, prob, prob) which would be MUCH smaller
    # Would make run_in_parallel loop slightly harder
    m_probs = {
        m: (p_merit, p_record) for m, p_merit, p_record in zip(
            valid_m, data.prob_merit_gap, data.prob_record_gap)
        if p_record >= prob_threshold
    }

    # No longer needed
    del data.prob_merit_gap[:]
    del data.prob_record_gap[:]

    return prob_threshold, m_probs


def should_print_stats(
        args, sc, data, m,
        unknown_l, unknown_u,
        prev_p, next_p):
    stop_t = time.time()
    print_secs = stop_t - sc.last_print_t

    print_small = (1, 10, 30, 100, 300, 1000, 3000, 5000, 10000, 30000, 100000, 300000)
    # <1000 => 1e6, <2000 => 1e5, <
    print_every = 10 ** max(0, 16 - int(math.log2(args.p)))
    assert 1 < print_every <= 10 ** 7, args.p

    # Print a little bit if we resume but mostly as we test.
    if ((args.megagap and sc.tested % 10 == 0)
            or sc.tested in print_small
            or (sc.tested and sc.tested % print_every == 0)
            or m == data.last_m or print_secs > 3600): # 24/day
        secs = stop_t - sc.start_t

        print("\t{:3d} {:4d} <- unknowns -> {:-4d}\t{:4d} <- gap -> {:-4d}".format(
            m, unknown_l, unknown_u, prev_p, next_p))

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

        total_unknown = sc.t_unk_prev + sc.t_unk_next
        if total_unknown:
            print("\t    unknowns   {:<9d} (avg: {:.2f}), {:.2%} composite  "
                  "{:.2%} <- % -> {:.2%}".format(
                total_unknown, total_unknown / sc.tested,
                (1 - total_unknown / ((2 * args.sieve_length + 1) * sc.tested)),
                sc.t_unk_prev / total_unknown,
                sc.t_unk_next / total_unknown))
        prp_tests = sc.total_prp_tests
        if sc.tested and prp_tests:
            print("\t    prp tests  {:<9d} (avg: {:.2f}/side, {:.2f}/m) ({:.3f} tests/sec)".format(
                prp_tests, prp_tests / sc.sides_processed, prp_tests / sc.tested, prp_tests / secs))
            if sc.one_side_skips:
                avg_prob_side = sc.prob_record_processed / sc.sides_processed
                avg_prob_no_skip = sc.prob_record_all / (2 * sc.tested)
                print("\t    side skips {:<9d} ({:.1%}) "
                      "(prob_record {:.3g}/side, {:.3g}/m/2, {:+.1%})".format(
                    sc.one_side_skips, sc.one_side_skips / sc.tested,
                    avg_prob_side, avg_prob_no_skip,
                    avg_prob_side / avg_prob_no_skip - 1))

            print("\t    sum(prob_min_merit):  {:7.5g}, {:.3g}/day\tfound: {}".format(
                sc.prob_min_merit, 86400 / secs * sc.prob_min_merit, sc.count_min_merit))
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


def process_result(conn, args, record_gaps, m_probs, data, sc, result):
    """Handles new results"""
    (m, r_log_n, unknown_l, unknown_u,
     n_tests, next_p,
     p_tests, prev_p,
     prob_record, new_prob_record,
     test_time) = result

    gap = next_p + prev_p
    merit = gap / r_log_n

    save(conn, args.p, args.d, m, next_p, prev_p, merit, n_tests, p_tests, test_time)

    if next_p < 0 or prev_p < 0:
        # partial result don't print anything
        return

    sc.tested += 1
    sc.test_time += test_time

    sc.t_unk_prev += unknown_l
    sc.t_unk_next += unknown_u

    sc.prob_min_merit += m_probs[m][0]
    if merit > args.min_merit:
        sc.count_min_merit += 1

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
            m,
            unknown_l, unknown_u,
            prev_p, next_p,
    ):
        sc.best_merit_interval = 0
        sc.best_merit_interval_m = -1
        sc.last_print_t = time.time()
