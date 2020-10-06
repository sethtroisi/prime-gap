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
import re
import signal
import sqlite3
import sys
import time
from collections import defaultdict
from dataclasses import dataclass

import gmpy2

import gap_utils



def get_arg_parser():
    parser = argparse.ArgumentParser('Test prime gaps')

    parser.add_argument('--search-db', type=str,
        default="prime-gap-search.db",
        help="Prime database from gap_test")

    parser.add_argument('--prime-gaps-db', type=str,
        default="gaps.db",
        help="Prime gap database see github.com/primegap-list-project/prime-gap-list")

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

    return parser


#---- Stats to plot ----#
def stats_plots(
        args,
        min_merit_gap, record_gaps, prob_nth,
        valid_m,
        data, misc):

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats

    # XXX: How to handle when only partial results

    def verify_no_trend(valid_m, data):
        """There should be no trend based on m value"""
        # Have to adjust for Expected gap which has two data points for each m
        m_values = valid_m
        if len(d) == 2 * len(valid_m):
            m_values = [mi for mi in valid_m for side in ['l', 'r']]

        trend, _ = np.polyfit(m_values, d, 1)
        if trend > 2e-3:
            print("\n")
            print("NON-ZERO TREND: ", trend) # Verify that expected value doesn't vary with M.
            print("\n")

    def plot_hist(axis, data, color, marker):
        hist_data = np.histogram(data, bins=100, density=True)
        axis.scatter(hist_data[1][:-1], hist_data[0], color=color, marker=marker, s=8, label='Observed')
        max_y = hist_data[0].max()

        axis.set_xlim(np.percentile(data, 0.01), np.percentile(data, 99.9))
        axis.set_ylim(top=1.2 * max_y)

    def plot_cdf(axis, data, color, label):
        n = len(data)
        d_sorted = np.sort(data)
        dist_label = f"Empirical CDF({label})" if 'P(' not in label else "CDF"
        axis.plot(
                d_sorted, np.arange(1, n+1) / n,
                color=color, label=dist_label)
        axis.set_xlim(np.percentile(data, 0.01), np.percentile(data, 99.9))

        # Draw some lines for 50th, 90th, 95th percentile
        for percent in (50, 90, 95):
            percentile = np.percentile(d_sorted, percent)
            if percentile > 100:
                dist_label=f"{percent}th percentile = {percentile:.0f}"
            elif percentile > .01:
                dist_label=f"{percent}th percentile = {percentile:.4f}"
            else:
                dist_label=f"{percent}th percentile = {percentile:.2e}"
            axis.plot(
                    [0, percentile, percentile],
                    [percent/100, percent/100, 0],
                    color="sandybrown",
                    label=dist_label)
        axis.legend(loc='upper left')

    def plot_prob_hist(axis, probs, color):
        x, w = zip(*sorted(((g, v) for g,v in probs.items() if v > 0)))
        axis.hist(x, weights=w, bins=100, density=True,
                  label='Theoretical P(gap)', color=color, alpha=0.4)
        print (f"|P(gap)| = {len(x)}, Sum(P(gap)) = {sum(w):.1f}")

    def prob_histogram_all(axis, probs, experimental, label):
        plot_prob_hist(axis, probs,  'blueviolet')
        # Expected value
        add_expected_value(axis, experimental, 'peru', label)
        # Experimental values
        plot_hist(axis, experimental, 'peru', 'x')

        # XXX: can I query this from experimental? axis.hist?
        min_y = 0.8 * min(v for v in probs.values() if v > 0) / sum(probs.values())
        max_y = 1.2 * max(probs.values()) / sum(probs.values())
        axis.set_yscale('log')
        axis.legend(loc='upper right')

        return min_y, max_y

    def fit_normal_dist(axis, data):
        x_start = max(0, 0.95 * np.percentile(data, 1))
        x_end   = 1.05 * np.percentile(data, 99)
        gap_span = np.linspace(x_start, x_end, 400)

        mu, std = stats.norm.fit(data)
        p = stats.norm.pdf(gap_span, mu, std)
        axis.plot(gap_span, p, color=color)

    def add_expected_value(axis, data, color, label):
        E = np.mean(data)
        label_e = f"E({label}) = {E:.0f}"
        axis.axvline(x=E, ymax=1.0/1.2, color=color, label=label_e)

    egap_n = len(data.experimental_gap)
    if len(data.expected_gap) != egap_n:
        print("experimental_gap size mismatch", len(data.expected_gap), egap_n)
    slope, _, R, _, _ = stats.linregress(data.expected_gap[:egap_n], data.experimental_gap)
    print ()
    print ("R^2 for expected gap: {:.3f}, corr: {:.3f}".format(R**2, slope))
    print ()

    if args.num_plots > 0:
        # Plot 1: Gap(side):
        #   [ prev mi,    next mi]
        #   [ BLANK,   prob all m]
        #   [ expected,   cdf    ]

        # Set up subplots.
        fig = plt.figure(
            "Per Side Statistics",
            constrained_layout=True,
            figsize=(8, 12))
        gs = fig.add_gridspec(3, 2)

        if misc.test_unknowns:
            axis_prev = fig.add_subplot(gs[0, 0])
            axis_next = fig.add_subplot(gs[0, 1])
        axis_prob_gap     = fig.add_subplot(gs[1, 1])
        axis_expected_gap = fig.add_subplot(gs[2, 0])
        axis_cdf_gap      = fig.add_subplot(gs[2, 1])

        def plot_prob_nth(axis, unknowns, color, label):
            n = min(len(prob_nth), len(unknowns))
            axis.scatter(
                unknowns[:n], prob_nth[:n],
                marker='.', s=12, color=color, label=label)
            # calculate expected value = sum(i * prob(i))
            E = sum(u * p for u, p in zip(unknowns, prob_nth))
            axis.axvline(
                x=E, ymax=1.0/1.2,
                color=color, label=f"E({label}) = {E:.0f}")

        if misc.test_unknowns:
            # prob_prev, prev_next for individual m
            # See Prob_nth in gap_stats
            colors = ['lightskyblue', 'tomato', 'seagreen']
            for m, c, (u_p, u_n) in zip(valid_m, colors, misc.test_unknowns):
                label = f"m={m}"
                plot_prob_nth(axis_prev, u_p, c, label)
                plot_prob_nth(axis_next, u_n, c, label)

            axis_prev.legend(loc='upper left')
            axis_next.legend(loc='upper right')
            axis_prev.set_yscale('log')
            axis_next.set_yscale('log')
            axis_prev.set_xlim(-args.sieve_length, 0)
            axis_next.set_xlim(0, args.sieve_length)

        min_y, max_y = prob_histogram_all(
                axis_prob_gap, misc.prob_gap_side, data.experimental_side, 'next')
        axis_prob_gap.set_xlim(0, args.sieve_length)
        axis_prob_gap.set_ylim(bottom=min_y / 10)
        print(f"Min Prob(gap side): {min_y:.2e}")

        for e_data, color, label in (
                (data.expected_prev, 'lightskyblue', 'prev'),
                (data.expected_next, 'tomato', 'next'),
        ):
            plot_hist(axis_expected_gap,          e_data, color, 'x')
            add_expected_value(axis_expected_gap, e_data, color, label)
            fit_normal_dist(axis_expected_gap, e_data)
            axis_expected_gap.legend(loc='upper left')

            # CDF of gap <= x
            plot_cdf(axis_cdf_gap, e_data, color, label)

    if args.num_plots > 1:
        # Plot 2: Gap(combined):
        #   [ prob all m,               expected,   cdf ]
        #   [ sum(prob) & count,  dist,       cdf ] (for >minmerit)
        #   [ sum(prob) & count,  dist,       cdf ] (for record)

        # Set up subplots.
        fig = plt.figure(
            "Combined Gap Statistics",
            constrained_layout=True,
            figsize=(12, 12))
        gs = fig.add_gridspec(3, 3)

        axis_prob_comb     = fig.add_subplot(gs[0, 0])
        axis_expected_comb = fig.add_subplot(gs[0, 1])
        axis_cdf_comb      = fig.add_subplot(gs[0, 2])

        # Combining all probs for a pseudo distribution of P(gap combined)
        min_y, max_y = prob_histogram_all(
                axis_prob_comb, misc.prob_gap_comb, data.experimental_gap, 'gap')
        print(f"Min Prob(gap comb): {min_y:.2e}")
        min_y = max(1e-7, min_y)
        axis_prob_comb.set_ylim(bottom=min_y, top=max_y)

        color = 'seagreen'
        axis = axis_expected_comb
        e_gap = data.expected_gap
        plot_hist(axis,          e_gap, color, 'x')
        add_expected_value(axis, e_gap, color, 'combined gap')
        fit_normal_dist(axis,    e_gap)
        axis.legend(loc='upper left')

#        # CDF of gap <= x
        plot_cdf(axis_cdf_comb,  e_gap, color, 'combined gap')

        # Sum(prob) vs m's tested
        # Hist of Prob(record/merit)
        # CDF of Prob(record/merit)
        for row, prob_data, label, color in (
            (1, data.prob_merit_gap, f'gap > min_merit({args.min_merit})', 'dodgerblue'),
            (2, data.prob_record_gap, 'gap = record', 'lightskyblue'),
        ):
            axis = fig.add_subplot(gs[row, 0])

            # P(gap > min_merit_gap) & Count(gap > min_merit_gap)
            # sorted and unsorted order
            axis.set_xlabel(" # of m's tests")
            axis.set_ylabel(f'Sum(P(gap {label})')

            #assert len(prob_data) == len(data.experimental_gap)
            zipped = list(zip(prob_data, data.experimental_gap))

            p_gap_merit_sorted, _ = zip(*sorted(zipped, reverse=True))
            p_gap_merit_ord, gap_real_ord = zip(*zipped)

            print(f"{label:20} | sum(P) = {sum(prob_data):.3f}")

            #Experimental
            if row == 1:
                cumcount = np.cumsum(np.array(gap_real_ord) > min_merit_gap)
            else:
                cumcount = np.cumsum(np.array([g in record_gaps for g in gap_real_ord]))

            tests = list(range(1, len(p_gap_merit_ord)+1))
            axis.plot(tests, cumcount, label='Count ' + label)

            # Theoretical
            cumsum_p = np.cumsum(p_gap_merit_ord)
            cumsum_p_sorted = np.cumsum(p_gap_merit_sorted)

            axis.plot(tests, cumsum_p, label=f'Sum(P({label}))')
            z  = axis.plot(tests, cumsum_p_sorted, label=f'Sum(P({label})) (best first)')
            axis.legend(loc='upper left')

            # Hist
            axis = fig.add_subplot(gs[row, 1])
            plot_hist(axis, prob_data, color, 'x')

            # CDF
            axis = fig.add_subplot(gs[row, 2])
            plot_cdf(axis,  prob_data, color, label)


    if args.num_plots > 2:
        # Older 1-page 3x3 layout

        fig = plt.figure(
            "Prime Gap Statistics",
            constrained_layout=True,
            figsize=(12, 12))
        gs = fig.add_gridspec(3, 3)
        axis_one_gap      = fig.add_subplot(gs[0, 2])
        axis_combined_gap = fig.add_subplot(gs[1, 2])

        for plot_i, (d, label, color) in enumerate((
                (data.expected_prev, 'prev', 'lightskyblue'),
                (data.expected_next, 'next', 'tomato'),
                (data.expected_gap,  'expected', 'seagreen'),
                (data.experimental_side, 'next/prev gap', 'sandybrown'),
                (data.experimental_gap, 'gap', 'peru'),
                (data.prob_merit_gap, f'P(gap > min_merit({args.min_merit}))', 'dodgerblue'),
        )):
            if not d: continue

            if plot_i == 0: # Not called again on plot_i == 1
                axis = fig.add_subplot(gs[0, 0])
                dist_axis = fig.add_subplot(gs[0, 1])
            elif plot_i == 2:
                axis = fig.add_subplot(gs[1, 0])
                dist_axis = fig.add_subplot(gs[1, 1])
            elif plot_i == 3:
                axis = axis_one_gap
                dist_axis = None
            elif plot_i == 4:
                axis = axis_combined_gap
                dist_axis = None
            elif plot_i == 5:
                axis = fig.add_subplot(gs[2, 0])
                dist_axis = fig.add_subplot(gs[2, 1])

            verify_no_trend(valid_m, d)

            marker='o' if plot_i in (3,4) else 'x'
            plot_hist(axis, d, color, marker)

            if 'gap' not in label:
                # Fit normal distribution to data
                fit_normal_dist(axis, d)

            if plot_i != 5:
                # Plot a line for expected value
                add_expected_value(axis, d, color, label)
                axis.legend()
            else:
                axis.legend([label])

            if dist_axis:
                # Cumulative sum of probability by gap
                plot_cdf(dist_axis, d, color, label)

        for d, color in [
                (misc.prob_gap_side, 'blueviolet'),
                (misc.prob_gap_comb, 'seagreen'),
        ]:
            if color == 'blueviolet':
                axis = axis_one_gap
            else:
                axis = axis_combined_gap

            plot_prob_hist(axis, d,  color)
            axis.legend(loc='upper right')

        assert len(data.prob_merit_gap) == len(data.experimental_gap)
        zipped = list(itertools.zip_longest(data.prob_merit_gap, data.experimental_gap))

        # P(gap > min_merit_gap) & Count(gap > min_merit_gap)
        # sorted and unsorted order
        fig.add_subplot(gs[2, 2])
        plt.xlabel(" # of m's tests")
        plt.ylabel(f'Sum(P(gap > min_merit)')

        p_gap_merit_sorted, _ = zip(*sorted(zipped, reverse=True))
        p_gap_merit_ord, gap_real_ord = zip(*zipped)

        tests = list(range(1, len(p_gap_merit_ord)+1))

        #Experimental
        cumcount_large = np.cumsum(np.array(gap_real_ord) > min_merit_gap)
        plt.plot(tests, cumcount_large, label='Count gap > min_merit')

        # Theoretical
        cumsum_p = np.cumsum(p_gap_merit_ord)
        cumsum_p_sorted = np.cumsum(p_gap_merit_sorted)

        plt.plot(tests, cumsum_p, label='Sum(P(gap > min_merit))')
        z  = plt.plot(tests, cumsum_p_sorted, label='Sum(P(gap > min_merit)) (best first)')

        '''
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

        # Want this one below next graph
        z3 = plt.plot(tests, cumsum_p_restart_freq, label='(top 10% of 10x larger run)')
        z2 = plt.plot(tests, cumsum_p_restart, label='P(gap > min_merit) (top 50% of two sieves)')

        # Plot speedup at 50th percentile
        mid_t = len(tests) // 2

        y = [cumsum_p[mid_t], cumsum_p_sorted[mid_t]]
        plt.plot([tests[mid_t], tests[mid_t]], y, c=z[0].get_color(),
                    label="+{:.1%} by sorting at midpoint".format(y[1] / y[0] - 1))

        y = [cumsum_p[-1], cumsum_p_restart[-1]]
        plt.plot([tests[-1], tests[-1]], y, c=z2[0].get_color(),
                    label="+{:.1%} using top 50% & sorting".format(y[1] / y[0] - 1))
        '''

        #plt.legend(loc='upper left')

    if args.save_logs:
        plt.savefig(args.unknown_filename + ".png", dpi=1080//8)

    if args.num_plots:
        plt.show()

    plt.close()


def plot_stuff(
        args, conn, data, sc, misc,
        min_merit_gap, record_gaps, prob_nth):
    if args.stats:
        # Not calculated
        data.prob_record_gap = data.prob_merit_gap

        stats_plots(
            args, min_merit_gap, record_gaps, prob_nth,
            data.valid_m, data, misc)

    # Load stats from gap_stats
    data_db, misc_db = load_stats(conn, args)
    assert data_db.expected_prev
    assert misc_db.prob_gap_comb, len(misc.prob_gap_comb)

    # test_unknowns come from unknown-file not DB.
    misc_db.test_unknowns = misc.test_unknowns

    stats_plots(
        args, min_merit_gap, record_gaps, prob_nth,
        data.valid_m, data_db, misc_db)


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

    # NOTE: to make sure 'count(records) still count after a gap
    # record is submitted subtract a small amound from merit
    rv = conn.execute(
        "SELECT gapsize FROM gaps "
        "WHERE gapsize > merit * ? - 0.001"
        "  AND gapsize < 35 * ?",
        (log_N, log_N))
    return tuple(g[0] for g in sorted(rv))


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
    print("\texpect {:.0f} left, {:.3f}% of SL={} after {}M".format(
        expected, 100.0 * expected / (SL + 1), SL, max_prime//1e6))
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


def should_print_stats(
        args, sc, data, mi, m,
        unknown_l, unknown_u,
        prev_p_i, next_p_i):

    stop_t = time.time()
    print_secs = stop_t - sc.last_print_t

    N = len(data.valid_m)
    if N in (1,10,30,100,300,1000) or N % 5000 == 0 \
            or sc.tested in (1,10,30,100,300,1000) or (sc.tested and sc.tested % 5000 == 0) \
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
        print("\t    tests      {:<9d} ({}, {})  {:.0f}, {:.0f} secs".format(
            sc.tested, timing, timing_threads, secs, sc.test_time))

        print("\t    sum(prob_minmerit):  {:.2g}, {:.3g}/day\tfound: {}".format(
            sc.prob_minmerit, 86400 / secs * sc.prob_minmerit, sc.count_minmerit))
        print("\t    sum(prob_record):    {:.2g}, {:.3g}/day".format(
            sc.prob_record, 86400  / secs * sc.prob_record))

        total_unknown = sc.t_unk_low + sc.t_unk_hgh
        if total_unknown:
            print("\t    unknowns   {:<9d} (avg: {:.2f}), {:.2f}% composite  {:.2f}% <- % -> {:.2f}%".format(
                total_unknown, total_unknown / N,
                100 * (1 - total_unknown / ((2 * args.sieve_length + 1) * N)),
                100 * sc.t_unk_low / total_unknown,
                100 * sc.t_unk_hgh / total_unknown))
        prp_tests = sc.total_prp_tests
        if sc.tested and prp_tests:
            print("\t    prp tests  {:<9d} (avg: {:.2f}) ({:.3f} tests/sec)".format(
                prp_tests, prp_tests / sc.tested, prp_tests / secs))
            if sc.gap_out_of_sieve_prev or sc.gap_out_of_sieve_next:
                print("\t    fallback prev_gap {} ({:.1f}%), next_gap {} ({:.1f}%)".format(
                    sc.gap_out_of_sieve_prev, 100 * sc.gap_out_of_sieve_prev / sc.tested,
                    sc.gap_out_of_sieve_next, 100 * sc.gap_out_of_sieve_next / sc.tested))
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
    assert next_p_i > 0 and prev_p_i > 0, (mi, m, next_pi, prev_p_i)

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


# NOTE: Manager is prime_gap_test, Worker is maintains process_line workers
def process_line(done_flag, work_q, results_q, thread_i, SL, K, P, D, primes, remainder):
    def cleanup(*_):
        print (f"Thread {thread_i} stopping")
        work_q.close()
        work_q.join_thread()
        results_q.close()
        results_q.join_thread()
        exit(0)

    # Ignore KeyboardInterrupt and let Manager terminate us.
    signal.signal(signal.SIGINT, cleanup)

    print (f"Thread {thread_i} started")
    while True:
        work = work_q.get()
        if work == "STOP":
            assert done_flag.is_set()
            cleanup()
            return

        #print (f"Work({thread_i}) done({done_flag.is_set()}): {work if work[0] == 'S' else work[:3]}")

        m, mi, log_n, line = work

        mtest, unknown_l, unknown_u, unknowns = gap_utils.parse_unknown_line(line)
        assert mi == mtest

        # Used for openPFGW
        strn = "{}*{}#/{}+".format(m, P, D)

        t0 = time.time()

        p_tests, prev_p_i = determine_prev_prime_i(m, strn, K, unknowns[0], SL, primes, remainder)
        n_tests, next_p_i = determine_next_prime_i(m, strn, K, unknowns[1], SL)

        test_time = time.time() - t0

        results_q.put((
            m, mi, log_n,
            unknown_l, unknown_u,
            n_tests, next_p_i,
            p_tests, prev_p_i,
            test_time,
        ))


def process_result(conn, args, record_gaps, mi_probs, data, sc, result):
    ''' Handles new results '''

    (m, mi, r_log_n, unknown_l, unknown_u,
     n_tests, next_p_i,
     p_tests, prev_p_i, test_time) = result

    sc.tested += 1
    sc.test_time += test_time

    sc.t_unk_low += unknown_l
    sc.t_unk_hgh += unknown_u

    gap = next_p_i + prev_p_i
    merit = gap / r_log_n

    sc.prob_minmerit           += mi_probs[mi][0]
    if merit > args.min_merit:
        sc.count_minmerit += 1
    sc.prob_record             += mi_probs[mi][1]

    save(conn,
        args.p, args.d, m, next_p_i, prev_p_i, merit,
        n_tests, p_tests, test_time)

    sc.total_prp_tests += p_tests
    sc.gap_out_of_sieve_prev += prev_p_i > args.sieve_length

    sc.total_prp_tests += n_tests
    sc.gap_out_of_sieve_next += next_p_i > args.sieve_length

    handle_result(
        args, record_gaps, data, sc,
        mi, m, r_log_n, prev_p_i, next_p_i, unknown_l, unknown_u)

def run_in_parallel(
        args, conn, unknown_file, record_gaps,
        prob_nth,
        existing, valid_mi,
        K, K_log,
        data, sc, misc
):

    # XXX: Cleanup after gmpy2.prev_prime.
    # Remainders of (p#/d) mod prime
    primes = tuple([p for p in range(3, 80000+1) if gmpy2.is_prime(p)])
    remainder = tuple([K % prime for prime in primes])

    done_flag = multiprocessing.Event()
    work_q    = multiprocessing.Queue(20)
    results_q = multiprocessing.Queue()
    #import queue
    #work_q    = queue.Queue(20)
    #results_q = queue.Queue()

    assert args.threads in range(1, 65), args.threads
    processes = []
    for i in range(args.threads):
        process = multiprocessing.Process(
            target=process_line,
            args=(
                done_flag, work_q, results_q,
                i, args.sieve_length, K, args.p, args.d,
                primes, remainder))
        process.start()
        processes.append(process)
    print()
    time.sleep(0.2)

    mi_probs = {
        mi: (p_merit, p_record) for mi, p_merit, p_record in zip(
            valid_mi, data.prob_merit_gap, data.prob_record_gap)
    }

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


        if m in existing:
            prev_p_i, next_p_i = existing[m]
            handle_result(
                args, record_gaps, data, sc,
                mi, m, log_n, prev_p_i, next_p_i, 0, 0)

        else:
            #print (f"Adding {m} to Queue")
            sc.will_test += 1
            try:
               work_q.put((m, mi, log_n, line), True)
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

    # ----- Open Output file
    print("\tLoading unknowns from {!r}".format(args.unknown_filename))
    print()
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
        gap_out_of_sieve_prev = 0
        gap_out_of_sieve_next = 0

        best_merit_interval = 0
        best_merit_interval_m = -1

        prob_record = 0.0
        prob_minmerit = 0.0
        count_minmerit = 0


    sc = StatCounters(time.time(), time.time())
    data = Datas()
    misc = Misc()

    valid_mi = [mi for mi in range(M_inc) if math.gcd(M + mi, D) == 1]
    data.first_m = M + valid_mi[0]
    data.last_m = M + valid_mi[-1]
    print(f"\nStarting m({len(valid_mi)}) {data.first_m} to {data.last_m}")
    print()

    # XXX: consider using all of data_db for resuming partially complete runs

    # Load stats for prob_record
    if not args.stats:
        data_db, misc_db = load_stats(conn, args)
        assert len(data_db.prob_merit_gap)  == len(valid_mi), "run ./gap_stats first"
        assert len(data_db.prob_record_gap) == len(valid_mi), "run ./gap_stats first"
        data.prob_merit_gap = data_db.prob_merit_gap
        data.prob_record_gap = data_db.prob_record_gap

    run_in_parallel(
        args, conn, unknown_file, record_gaps,
        prob_nth,
        existing, valid_mi,
        K, K_log,
        data, sc, misc
    )

    # ----- Stats and Plots
    if args.num_plots or args.save_logs:
        plot_stuff(
            args, conn, data, sc, misc,
            min_merit_gap, record_gaps, prob_nth)



if __name__ == "__main__":
    parser = get_arg_parser()
    args = parser.parse_args()
    gap_utils.verify_args(args)

    # TeeLogger context if args.save_logs
    with gap_utils.logger_context(args):
        prime_gap_test(args)

