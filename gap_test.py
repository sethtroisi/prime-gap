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
import re
import sqlite3
import time
from collections import defaultdict

import gmpy2

import gap_utils



def get_arg_parser():
    parser = argparse.ArgumentParser('Test prime gaps')

    parser.add_argument('--search-db', type=str,
        default="prime-gap-search.db",
        help="Prime database from gap_test")

    parser.add_argument('--unknown-filename', type=str,
        help="determine mstart, minc, p, d, sieve-length, and sieve-range"
             " from unknown-results filename")

    parser.add_argument('--min-merit',    type=int, default=10,
        help="only display prime gaps with merit >= minmerit")
    parser.add_argument('--run-prp',  action='store_true',
        help="Run PRP, leave off to benchmarking and print some stats")

    parser.add_argument('--plots', action='store_true',
        help="Show plots about distributions")

    parser.add_argument('--save-logs', action='store_true',
        help="Save logs and plots about distributions")

    return parser


#---- Stats to plot ----#
def stats_plots(
        args,
        min_merit_gap,

        valid_m,
        s_expected_gap,
        s_expected_prev, s_expected_next,
        s_experimental_gap, s_experimental_side,
        p_gap_side, p_gap_comb, p_gap_merit
    ):


    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats

    if args.run_prp:
        slope, _, R, _, _ = stats.linregress(s_expected_gap, s_experimental_gap)
        print ("R^2 for expected gap: {:.3f}, corr: {:.3f}".format(R**2, slope))

    # Set up subplots.
    fig3 = plt.figure(constrained_layout=True, figsize=(8, 8))
    gs = fig3.add_gridspec(6, 2)

    for plot_i, (d, label, color) in enumerate([
            (s_expected_prev, 'prev', 'lightskyblue'),
            (s_expected_next, 'next', 'darksalmon'),
            (s_expected_gap,  'expected', 'seagreen'),
            (s_experimental_side, 'next/prev gap', 'sandybrown'),
            (s_experimental_gap, 'gap', 'peru'),
            (p_gap_merit, f'P(gap > min_merit({args.min_merit}))', 'dodgerblue'),
    ]):
        if not d: continue

        if plot_i == 0: # Not called again on plot_i == 1
            fig3.add_subplot(gs[0:2, 0])
        elif plot_i == 2:
            fig3.add_subplot(gs[2:4, 0])
        elif plot_i == 3:
            fig3.add_subplot(gs[0:2, 1])
        elif plot_i == 4:
            fig3.add_subplot(gs[2:4, 1])
        elif plot_i == 5:
            fig3.add_subplot(gs[4:6, 0])

        hist_data = np.histogram(d, bins=100, density=True)
        plt.scatter(hist_data[1][:-1], hist_data[0], color=color,
                    marker='o' if plot_i in (3,4) else 'x', s=8)
        max_y = hist_data[0].max()

        if plot_i == 3:
            trend, _ = np.polyfit([mi for mi in valid_m for d in ['l','up']], d, 1)
        else:
            trend, _ = np.polyfit(valid_m, d, 1)

        if trend > 2e-3:
            print("\n")
            print("NON-ZERO TREND: ", trend) # Verify that expected value doesn't vary with M.
            print("\n")

        E = np.mean(d)

        if plot_i != 5:
            plt.axvline(x=E, ymax=1.0/1.2, color=color, label=f"E({label}) = {E:.0f}")
            plt.legend(loc='upper right')
        else:
            plt.legend([label])

        gap_span = np.linspace(0.95 * np.percentile(d, 1), 1.05 * np.percentile(d, 99), 400)
        if plot_i == 5:
            gap_span = np.linspace(0, max(d) * 1.1, 400)

        if 'gap' not in label:
            mu, std = stats.norm.fit(d)
            p = stats.norm.pdf(gap_span, mu, std)
            plt.plot(gap_span, p, color=color)

        plt.xlim(np.percentile(d, 0.01), np.percentile(d, 99.9))
        plt.ylim(2e-5, 1.2 * max_y)

    for d, color in [
            (p_gap_side, 'blueviolet'),
            (p_gap_comb, 'seagreen'),
    ]:
        if color == 'blueviolet':
            fig3.add_subplot(gs[0:2, 1])
        else:
            fig3.add_subplot(gs[2:4, 1])

        d_x, d_w = zip(*sorted(d.items()))

        plt.hist(d_x, weights=d_w, bins=100, density=True,
                 label='Theoretical P(gap)', color=color, alpha=0.4)
        plt.legend(loc='upper right')

    assert len(p_gap_merit) == len(s_experimental_gap)
    zipped = list(itertools.zip_longest(p_gap_merit, s_experimental_gap))

    # P(gap > min_merit_gap) & Count(gap > min_merit_gap)
    # sorted and unsorted order
    fig3.add_subplot(gs[4:6, 1])

    p_gap_merit_sorted, _ = zip(*sorted(zipped, reverse=True))
    p_gap_merit_ord, gap_real_ord = zip(*zipped)

    tests = list(range(1, len(p_gap_merit_ord)+1))

    # Theoretical
    cumsum_p = np.cumsum(p_gap_merit_ord)
    cumsum_p_sorted = np.cumsum(p_gap_merit_sorted)

    # Theoretical with restart
    def cumsum_restarts(restarts):
        part_size = len(tests) // (restarts + 1)
        t = []
        for i in range(restarts):
            t.extend(p_gap_merit_sorted[:part_size])
        t.extend(p_gap_merit_sorted[:len(tests) - len(t)])
        return np.cumsum(t)

    cumsum_p_restart = cumsum_restarts(1)
    cumsum_p_restart_freq = cumsum_restarts(9)

    # Experimental
    cumcount_large = np.cumsum(np.array(gap_real_ord) > min_merit_gap)

    # Want this one below next graph
    z3 = plt.plot(tests, cumsum_p_restart_freq, label='10 runs of top 10%')
    z2 = plt.plot(tests, cumsum_p_restart, label='P(gap > min_merit) (top 50% of two sieves)')
    z  = plt.plot(tests, cumsum_p_sorted, label='P(gap > min_merit) (sorting by best first)')
    plt.plot(tests, cumsum_p, label='P(gap > min_merit)')
    plt.plot(tests, cumcount_large, label='Count gap > min_merit')
    plt.xlabel(" # of m's tests")
    plt.ylabel(f'Sum(P(gap > min_merit)')

    # Plot speedup at 50th percentile
    mid_t = len(tests) // 2

    y = [cumsum_p[mid_t], cumsum_p_sorted[mid_t]]
    plt.plot([tests[mid_t], tests[mid_t]], y, c=z[0].get_color(),
                label="+{:.1%} by sorting".format(y[1] / y[0] - 1))

    y = [cumsum_p[-1], cumsum_p_restart[-1]]
    plt.plot([tests[-1], tests[-1]], y, c=z2[0].get_color(),
                label="+{:.1%} using top 50% & sorting".format(y[1] / y[0] - 1))

    plt.legend(loc='upper left')


    # TODO plot number of PRP tests per M
    # TODO(P(record))
    # TODO(P(gap > X))

    # Another plot needed is sieve_depth
    # warm audiance up to time on Y axis with pause will sieving

    # Another plot needed is
    #   P(record gap) sorted vs unsorted
    #   P(rg) sorted 50% then new gap



    if args.save_logs:
        plt.savefig(args.unknown_filename + ".png", dpi=1080//8)

    if args.plots:
        plt.show()

    plt.close()


#---- gap_testing ----#


def load_existing(conn, args):
    # TODO to_process_range
    rv = conn.execute(
        "SELECT m, next_p_i, prev_p_i FROM result"
        " WHERE P = ? and D = ? and m BETWEEN ? AND ?",
        (args.p, args.d, args.mstart, args.mstart + args.minc - 1))
    return {row['m']: [row['prev_p_i'], row['next_p_i']] for row in rv}


def save(conn, m, p, d, n_p_i, p_p_i, merit):
    conn.execute(
        "INSERT INTO result VALUES(null, ?, ?, ?, ?, ?, ?)",
        (m, p, d, n_p_i, p_p_i, round(merit,3)))
    conn.commit()


def prob_prime_sieve_length(M, D, prob_prime, K_digits, K_primes, SL, sieve_range):
    assert sieve_range >= 10 ** 6, sieve_range

    # From Mertens' 3rd theorem
    gamma = 0.577215665
    unknowns_after_sieve = 1.0 / (math.log(sieve_range) * math.exp(gamma))

    prob_prime_coprime = 1
    prob_prime_after_sieve = prob_prime / unknowns_after_sieve

    for prime in K_primes:
        if D % prime != 0:
            prob_prime_coprime *= (1 - 1/prime)

    count_coprime = SL
    for i in range(1, SL+1):
        for prime in K_primes:
            if (i % prime) == 0 and (D % prime) != 0:
                count_coprime -= 1
                break

    chance_coprime_composite = 1 - prob_prime / prob_prime_coprime
    prob_gap_shorter_hypothetical = chance_coprime_composite ** count_coprime

    # count_coprime already includes some parts of unknown_after_sieve
    print("\t{:.3f}% of SL should be unknown ({}M) ~= {:.0f}".format(
        100 * unknowns_after_sieve,
        sieve_range//1e6,
        count_coprime * (unknowns_after_sieve / prob_prime_coprime)))
    print("\t{:.3f}% of {} digit numbers are prime".format(
        100 * prob_prime, K_digits))
    print("\t{:.3f}% of tests should be prime ({:.1f}x speedup)".format(
        100 * prob_prime_after_sieve, 1 / unknowns_after_sieve))
    print("\t~2x{:.1f} = {:.1f} PRP tests per m".format(
        1 / prob_prime_after_sieve, 2 / prob_prime_after_sieve))
    print("\tsieve_length={} is insufficient ~{:.2f}% of time".format(
        SL, 100 * prob_gap_shorter_hypothetical))
    print()

    return prob_prime_after_sieve


def calculate_expected_gaps(composites, SL, prob_prime_after_sieve, log_m,
                            p_gap_side, p_gap_comb, min_merit_gap):
    expected_side = []

    # Geometric distribution (could be cached)
    probs = []
    prob_longer = []
    prob_gap_longer = 1
    for i in range(max(len(composites[0]), len(composites[1]))+1):
        probs.append(prob_gap_longer * prob_prime_after_sieve)
        prob_longer.append(prob_gap_longer)

        prob_gap_longer *= (1 - prob_prime_after_sieve)

    assert min(probs) > 0

    for side in composites:
        expected_length = 0
        for v, prob in zip(side, probs):
            expected_length += abs(v) * prob
            p_gap_side[abs(v)] += prob

        # expected to encounter a prime at distance ~= ln(start)
        prob_gap_longer = probs[len(side)]
        assert prob_gap_longer < 0.01, (prob_gap_longer, len(side))
        expected_length += (SL + log_m) * prob_gap_longer
        expected_side.append(expected_length)

    # TODO really slow, gap_test.cpp much faster.
    # TODO lookup best merit of every gap.
    p_merit = 0
    for prob_i, lower in zip(probs, composites[0]):
        for prob_j, upper in zip(probs, composites[1]):
            prob_joint = prob_i * prob_j
            gap = -lower + upper
            if gap >= min_merit_gap:
                p_merit += prob_joint

            p_gap_comb[gap] += prob_joint
            if prob_joint < 1e-8:
                break
        else:
            # gap = (K + SL+1) - (K - lower) = SL+1 - lower
            if -lower + SL+1 > min_merit_gap:
                p_merit += prob_i * prob_longer[len(composites[1])]

    for prob_j, upper in zip(probs, composites[1]):
        if SL+1 + upper > min_merit_gap:
            p_merit += prob_j * prob_longer[len(composites[0])]

    if 2*SL+1 > min_merit_gap:
        p_merit += prob_longer[len(composites[0])] * prob_longer[len(composites[1])]

    assert 0 <= p_merit <= 1.00, (p_merit, composites)
    return expected_side + [p_merit]


def determine_next_prime_i(m, strn, K, composites, SL):
    center = m * K
    tests = 0

    for i in composites:
        tests += 1;
        if gap_utils.is_prime(center + i, strn, i):
            next_p_i = i
            break
    else:
        # XXX: parse to version and verify > 6.2.99
        assert gmpy2.mp_version() == 'GMP 6.2.99'
        # Double checks center + SL.
        next_p_i = int(gmpy2.next_prime(center + SL) - center)

    return tests, next_p_i


def determine_prev_prime_i(m, strn, K, composites, SL, primes, remainder):
    center = m * K
    tests = 0

    for i in composites:
        assert i < 0
        tests += 1;
        if gap_utils.is_prime(center + i, strn, i):
            prev_p_i = -i
            break
    else:
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
                    prev_p_i = i
                    break

    return tests, prev_p_i


def prime_gap_test(args):
    P = args.p
    D = args.d

    M = args.mstart
    M_inc = args.minc

    SL = sieve_length = args.sieve_length
    sieve_range = args.sieve_range

    min_merit = args.min_merit

    K, K_digits, K_bits, K_log = gap_utils.K_and_stats(args)
    M_log    = K_log + math.log(M)
    min_merit_gap = int(min_merit * M_log)
    print("K = {} bits, {} digits, log(K) = {:.2f}".format(
        K_bits, K_digits, K_log))
    print("Min Gap ~= {} (for merit > {:.1f})\n".format(
        min_merit_gap, min_merit))

    # ----- Open Output file
    print("\tLoading unknowns from '{}'".format(args.unknown_filename))
    print()
    unknown_file = open(args.unknown_filename, "r")

    # ----- Open Prime-Gap-Search DB
    conn = sqlite3.connect(args.search_db)
    conn.row_factory = sqlite3.Row
    existing = {} #load_existing(conn, args)
    print (f"Found {len(existing)} existing results")

    # used in next_prime
    assert P <= 80000
    # Very slow but works.
    primes = [2] + [p for p in range(3, 80000+1, 2) if gmpy2.is_prime(p)]
    K_primes = [p for p in primes if p <= P]

    # ----- Allocate memory for a handful of utility functions.

    # Remainders of (p#/d) mod prime
    remainder   = [K % prime for prime in primes]

    # ----- Sieve stats
    prob_prime = 1 / M_log - 1 / (M_log * M_log)
    prob_prime_after_sieve = prob_prime_sieve_length(
        M, D, prob_prime, K_digits, K_primes, SL, sieve_range)

    # ----- Main sieve loop.
    print("\nStarting m={}".format(M))
    print()

    # Used for various stats
    s_start_t = time.time()
    s_last_print_t = time.time()
    s_total_unknown = 0
    s_t_unk_low = 0
    s_t_unk_hgh = 0
    s_total_prp_tests = 0
    s_gap_out_of_sieve_prev = 0
    s_gap_out_of_sieve_next = 0
    s_best_merit_interval = 0
    s_best_merit_interval_m = 0

    valid_m = []
    s_expected_prev = []
    s_expected_next = []
    s_expected_gap  = []
    s_experimental_side = []
    s_experimental_gap = []
    p_gap_side  = defaultdict(float)
    p_gap_comb  = defaultdict(float)
    p_gap_merit = []

    last_mi = M_inc - 1
    while math.gcd(M + last_mi, D) != 1:
        last_mi -= 1

    tested = 0
    for mi in range(M_inc):
        m = M + mi
        if math.gcd(m, D) != 1: continue

        valid_m.append(m)

        log_m = (K_log + math.log(m))

        unknowns = [[], []]

        prev_p_i = 0
        next_p_i = 0

        # Read a line from the file
        line = unknown_file.readline()
        start, c_l, c_h = line.split(" | ")

        match = re.match(r"^([0-9]+) : -([0-9]+) \+([0-9]+)", start)
        assert match, start
        mtest, unknown_l, unknown_u = map(int, match.groups())
        assert mtest == mi

        unknowns[0] = list(map(int,c_l.split(" ")))
        unknowns[1] = list(map(int,c_h.split(" ")))

        unknown_l_test = len(unknowns[0])
        unknown_u_test = len(unknowns[1])
        assert unknown_l == unknown_l_test, (unknown_l, unknown_l_test, "\t", start)
        assert unknown_u == unknown_u_test, (unknown_u, unknown_u_test, "\t", start)

        s_total_unknown += unknown_l + unknown_u
        s_t_unk_low += unknown_l
        s_t_unk_hgh += unknown_u

        if args.save_logs or args.plots:
            e_prev, e_next, p_merit = calculate_expected_gaps(
                unknowns, SL, prob_prime_after_sieve, log_m,
                p_gap_side, p_gap_comb, min_merit_gap)
            s_expected_prev.append(e_prev)
            s_expected_next.append(e_next)
            s_expected_gap.append(e_prev + e_next)
            p_gap_merit.append(p_merit)

        if args.run_prp:

            if m in existing:
                prev_p_i, next_p_i = existing[m]
            else:
                tested += 1
                # Used for openPFGW
                strn = "{}*{}#/{}+".format(m, P, D)

                p_tests, prev_p_i = determine_prev_prime_i(m, strn, K, unknowns[0],
                                                         SL, primes, remainder)
                s_total_prp_tests += p_tests
                s_gap_out_of_sieve_prev += prev_p_i > SL

                n_tests, next_p_i = determine_next_prime_i(m, strn, K, unknowns[1], SL)
                s_total_prp_tests += n_tests
                s_gap_out_of_sieve_next += next_p_i > SL

            assert prev_p_i > 0 and next_p_i > 0
            gap = int(next_p_i + prev_p_i)
            s_experimental_gap.append(gap)
            s_experimental_side.append(next_p_i)
            s_experimental_side.append(prev_p_i)
            assert next_p_i > 0 and prev_p_i > 0, (m, next_pi, prev_p_i)

            merit = gap / log_m
            if m not in existing:
                #save(conn, m, P, D, next_p_i, prev_p_i, merit)
                if merit > min_merit:
                    print("{}  {:.4f}  {} * {}#/{} -{} to +{}".format(
                        gap, merit, m, P, D, prev_p_i, next_p_i))

            if merit > s_best_merit_interval:
                s_best_merit_interval = merit
                s_best_merit_interval_m = m

        s_stop_t = time.time()
        print_secs = s_stop_t - s_last_print_t
        if len(valid_m) in (1,10,30,100,300,1000) or len(valid_m) % 5000 == 0 \
                or mi == last_mi or print_secs > 1200:
            secs = s_stop_t - s_start_t

            print("\t{:3d} {:4d} <- unknowns -> {:-4d}\t{:4d} <- gap -> {:-4d}".format(
                m,
                unknown_l, unknown_u,
                prev_p_i, next_p_i))
            if mi <= 10 and secs < 6: continue
            s_last_print_t = s_stop_t

            def roundSig(n, sig):
                return '{:g}'.format(float('{:.{p}g}'.format(n, p=sig)))

            # Want 3 sig figs which is hard in python
            if tested and tested < secs:
                timing = "{} secs/test".format(roundSig(secs / tested, 3))
            else:
                timing = "{}/sec".format(roundSig(tested / secs, 3))

            # Stats!
            print("\t    tests     {:<10d} ({})  {:.0f} seconds elapsed".format(
                tested, timing, secs))
            print("\t    unknowns  {:<10d} (avg: {:.2f}), {:.2f}% composite  {:.2f}% <- % -> {:.2f}%".format(
                s_total_unknown, s_total_unknown / len(valid_m),
                100 * (1 - s_total_unknown / ((2 * sieve_length + 1) * len(valid_m))),
                100 * s_t_unk_low / s_total_unknown,
                100 * s_t_unk_hgh / s_total_unknown))
            if tested and args.run_prp:
                print("\t    prp tests {:<10d} (avg: {:.2f}) ({:.3f} tests/sec)".format(
                    s_total_prp_tests, s_total_prp_tests / tested, s_total_prp_tests / secs))
                print("\t    fallback prev_gap {} ({:.1f}%), next_gap {} ({:.1f}%)".format(
                    s_gap_out_of_sieve_prev, 100 * s_gap_out_of_sieve_prev / tested,
                    s_gap_out_of_sieve_next, 100 * s_gap_out_of_sieve_next / tested))
                print("\t    best merit this interval: {:.3f} (at m={})".format(
                    s_best_merit_interval, s_best_merit_interval_m))

            s_best_merit_interval = 0
            s_best_merit_interval_m = -1

    if args.plots or args.save_logs:
        stats_plots(
            args,
            min_merit_gap,

            valid_m,
            s_expected_gap,
            s_expected_prev, s_expected_next,
            s_experimental_gap, s_experimental_side,
            p_gap_side, p_gap_comb, p_gap_merit
        )


if __name__ == "__main__":
    parser = get_arg_parser()
    args = parser.parse_args()
    gap_utils.verify_args(args, ".txt")

    # TeeLogger context if args.save_logs
    with gap_utils.logger_context(args):
        prime_gap_test(args)

