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

    parser.add_argument('--prime-gaps-db', type=str,
        default="gaps.db",
        help="Prime gap database see github.com/primegap-list-project/prime-gap-list")

    parser.add_argument('--unknown-filename', type=str,
        help="determine mstart, minc, p, d, sieve-length, and max-prime"
             " from unknown-results filename")

    parser.add_argument('--min-merit', type=int, default=12,
        help="only display prime gaps with merit >= minmerit")
    parser.add_argument('--run-prp',  action='store_true',
        help="Run PRP, leave off to benchmarking and print some stats")

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
        min_merit_gap, record_gaps,

        valid_m,
        s_test_unknowns, prob_nth,
        s_expected_gap,
        s_expected_prev, s_expected_next,
        s_experimental_gap, s_experimental_side,
        p_gap_side, p_gap_comb,
        p_merit_gap, p_record_gap
    ):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats

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
        dist_label = f"Empirical CDF({label} gap <= X)" if 'P(' not in label else "CDF"
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

    slope, _, R, _, _ = stats.linregress(s_expected_gap, s_experimental_gap)
    print ()
    print ("R^2 for expected gap: {:.3f}, corr: {:.3f}".format(R**2, slope))
    print ()

    if args.num_plots > 1:
        # Set up subplots.
        fig = plt.figure(
            "Per Side Statistics",
            constrained_layout=True,
            figsize=(8, 12))
        gs = fig.add_gridspec(3, 2)

        # Plot 1: Gap(side):
        #   [ prev mi,    next mi]
        #   [ BLANK,   prob all m]
        #   [ expected,   cdf    ]

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

        # prob_prev, prev_next for individual m
        # See Prob_nth in gap_stats
        colors = ['lightskyblue', 'tomato', 'seagreen']
        for m, c, (u_p, u_n) in zip(valid_m, colors, s_test_unknowns):
            label = f"m={m}"
            plot_prob_nth(axis_prev, u_p, c, label)
            plot_prob_nth(axis_next, u_n, c, label)

        axis_prev.legend(loc='upper left')
        axis_next.legend(loc='upper right')
        axis_prev.set_yscale('log')
        axis_next.set_yscale('log')
        axis_prev.set_xlim(-args.sieve_length, 0)
        axis_next.set_xlim(0, args.sieve_length)

        # Combining all probs for a pseudo distribution of P(next_gap) / P(prev_gap)
        plot_prob_hist(axis_prob_gap, p_gap_side,  'blueviolet')
        # Expected value
        add_expected_value(axis_prob_gap, s_experimental_side, 'peru', 'next')
        # Experimental values
        plot_hist(axis_prob_gap, s_experimental_side, 'peru', 'x')
        min_y = 0.8 * min(v for v in p_gap_side.values() if v > 0) / sum(p_gap_side.values())
        print(f"Min Prob(gap side): {min_y:.2e}")
        axis_prob_gap.set_xlim(0, args.sieve_length)
        axis_prob_gap.set_yscale('log')
        axis_prob_gap.set_ylim(bottom=min_y / 10)
        axis_prob_gap.legend(loc='upper right')

        for data, color, label in (
                (s_expected_prev, 'lightskyblue', 'prev'),
                (s_expected_next, 'tomato', 'next'),
        ):
            plot_hist(axis_expected_gap,          data, color, 'x')
            add_expected_value(axis_expected_gap, data, color, label)
            fit_normal_dist(axis_expected_gap, data)
            axis_expected_gap.legend(loc='upper left')

            # CDF of gap <= x
            plot_cdf(axis_cdf_gap, data, color, label)

    if args.num_plots > 0:
        # Set up subplots.
        fig = plt.figure(
            "Combined Gap Statistics",
            constrained_layout=True,
            figsize=(12, 12))
        gs = fig.add_gridspec(3, 3)

        # Plot 2: Gap(combined):
        #   [ prob all m,               expected,   cdf ]
        #   [ sum(prob) & count,  dist,       cdf ] (for highmerit)
        #   [ sum(prob) & count,  dist,       cdf ] (for record)

        axis_prob_comb     = fig.add_subplot(gs[0, 0])
        axis_expected_comb = fig.add_subplot(gs[0, 1])
        axis_cdf_comb      = fig.add_subplot(gs[0, 2])

        # Combining all probs for a pseudo distribution of P(gap combined)
        plot_prob_hist(axis_prob_comb, p_gap_comb,  'blueviolet')
        # Expected value
        add_expected_value(axis_prob_comb, s_experimental_gap, 'peru', 'next')
        # Experimental values
        plot_hist(axis_prob_comb, s_experimental_gap, 'peru', 'x')
        min_y = 0.8 * min(v for v in p_gap_comb.values() if v > 0) / sum(p_gap_comb.values())
        min_y = max(1e-7, min_y)
        max_y = 1.2 * max(p_gap_comb.values()) / sum(p_gap_comb.values())
        print(f"Min Prob(gap comb): {min_y:.2e}")
        axis_prob_comb.set_yscale('log')
        axis_prob_comb.set_ylim(bottom=min_y, top=max_y)
        axis_prob_comb.legend(loc='upper right')

        color = 'seagreen'
        axis = axis_expected_comb
        data = s_expected_gap
        plot_hist(axis,          data, color, 'x')
        add_expected_value(axis, data, color, 'combined gap')
        fit_normal_dist(axis,    data)
        axis.legend(loc='upper left')

#        # CDF of gap <= x
        plot_cdf(axis_cdf_comb,  data, color, 'combined gap')

        # Sum(prob) vs m's tested
        # Hist of Prob(record/merit)
        # CDF of Prob(record/merit)
        for row, data, label, color in (
            (1, p_merit_gap, f'gap > min_merit({args.min_merit})', 'dodgerblue'),
            (2, p_record_gap, 'gap = record', 'lightskyblue'),
        ):
            axis = fig.add_subplot(gs[row, 0])

            # P(gap > min_merit_gap) & Count(gap > min_merit_gap)
            # sorted and unsorted order
            axis.set_xlabel(" # of m's tests")
            axis.set_ylabel(f'Sum(P(gap {label})')

            assert len(data) == len(s_experimental_gap)
            zipped = list(itertools.zip_longest(data, s_experimental_gap))

            p_gap_merit_sorted, _ = zip(*sorted(zipped, reverse=True))
            p_gap_merit_ord, gap_real_ord = zip(*zipped)

            print(f"{label:20} | sum(P) = {sum(data):.3f}")

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
            plot_hist(axis, data, color, 'x')

            # CDF
            axis = fig.add_subplot(gs[row, 2])
            plot_cdf(axis,  data, color, label)


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
                (s_expected_prev, 'prev', 'lightskyblue'),
                (s_expected_next, 'next', 'tomato'),
                (s_expected_gap,  'expected', 'seagreen'),
                (s_experimental_side, 'next/prev gap', 'sandybrown'),
                (s_experimental_gap, 'gap', 'peru'),
                (p_merit_gap, f'P(gap > min_merit({args.min_merit}))', 'dodgerblue'),
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
                (p_gap_side, 'blueviolet'),
                (p_gap_comb, 'seagreen'),
        ]:
            if color == 'blueviolet':
                axis = axis_one_gap
            else:
                axis = axis_combined_gap

            plot_prob_hist(axis, d,  color)
            axis.legend(loc='upper right')

        assert len(p_merit_gap) == len(s_experimental_gap)
        zipped = list(itertools.zip_longest(p_merit_gap, s_experimental_gap))

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
    rv = conn.execute(
        "SELECT gapsize FROM gaps "
        "WHERE gapsize > merit * ?"
        "  AND gapsize < 35 * ?",
        (log_N, log_N))
    return tuple(g[0] for g in sorted(rv))

def load_existing(conn, args):
    rv = conn.execute(
        "SELECT m, next_p_i, prev_p_i FROM result"
        " WHERE P = ? AND D = ? AND m BETWEEN ? AND ?",
        (args.p, args.d, args.mstart, args.mstart + args.minc - 1))
    return {row['m']: [row['prev_p_i'], row['next_p_i']] for row in rv}


def load_stats(conn, args):
    rv = conn.execute(
        "SELECT e_gap_prev, e_gap_next, e_gap_prev + e_gap_next,"
        "       prob_merit, prob_record"
        " FROM m_stats WHERE P = ? AND D = ? AND m BETWEEN ? AND ?",
        (args.p, args.d, args.mstart, args.mstart + args.minc - 1))

    # Will fail if non-present
    (s_expected_prev, s_expected_next, s_expected_gap,
     p_merit_gap, p_record_gap)= zip(*[tuple(row) for row in rv])

    # Need dictionary
    p_gap_comb  = defaultdict(float)
    p_gap_side  = defaultdict(float)

    m_values = len(p_merit_gap)

    rv = conn.execute(
        "SELECT gap, prob_combined, prob_low_side, prob_high_side "
        "FROM range_stats where rid = ?",
        (config_hash(args),))
    for row in rv:
        gap = row['gap']
        # Values were normalized by / m_values in gap_stats
        p_gap_comb[gap] += row['prob_combined'] * m_values
        p_gap_side[gap] += row['prob_low_side'] / 2 * m_values
        p_gap_side[gap] += row['prob_high_side'] / 2 * m_values

    return (
        s_expected_prev, s_expected_next, s_expected_gap,
        p_merit_gap, p_record_gap,
        p_gap_comb, p_gap_side
    )

def save(conn, m, p, d, n_p_i, p_p_i, merit):
    conn.execute(
        "INSERT INTO result(m,P,D,next_p_i,prev_p_i,merit)"
        "VALUES(?,?,?,  ?,?,  ?)",
        (m, p, d,  n_p_i, p_p_i,  round(merit,4)))
    conn.commit()


def prob_prime_sieve_length(M, K, D, prob_prime, K_digits, K_primes, SL, max_prime):
    assert max_prime >= 10 ** 6, max_prime

    # From Mertens' 3rd theorem
    gamma = 0.577215665
    unknowns_after_sieve = 1.0 / (math.log(max_prime) * math.exp(gamma))

    prob_prime_coprime = 1
    prob_prime_after_sieve = prob_prime / unknowns_after_sieve

    for prime in K_primes:
        if D % prime != 0:
            prob_prime_coprime *= (1 - 1/prime)

    count_coprime = 0
    for i in range(1, SL+1):
        count_coprime += math.gcd(K, i) == 1
    print (f"{count_coprime=}")

    chance_coprime_composite = 1 - prob_prime / prob_prime_coprime
    prob_gap_shorter_hypothetical = chance_coprime_composite ** count_coprime

    # count_coprime already includes some parts of unknown_after_sieve
    print("\t{:.3f}% of SL should be unknown ({}M) ~= {:.0f}".format(
        100 * unknowns_after_sieve,
        max_prime//1e6,
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


def calculate_expected_gaps(
        SL, min_merit_gap, prob_nth, prob_longer,
        log_n, unknowns, p_gap_side, p_gap_comb):
    expected_side = []

    for side in composites:
        expected_length = 0
        for v, prob in zip(side, prob_nth):
            expected_length += abs(v) * prob
            # Normalized to 1 by matplotlib.hist later.
            p_gap_side[abs(v)] += prob

        # expected to encounter a prime at distance ~= ln(n)
        prob_gap_longer = prob_nth[len(side)]
        assert prob_gap_longer < 0.01, (prob_gap_longer, len(side))
        expected_length += (SL + log_n) * prob_gap_longer
        expected_side.append(expected_length)

    p_merit = 0
    for prob_i, lower in zip(prob_nth, composites[0]):
        for prob_j, upper in zip(prob_nth, composites[1]):
            prob_joint = prob_i * prob_j
            gap = -lower + upper
            if gap >= min_merit_gap:
                p_merit += prob_joint

            p_gap_comb[gap] += prob_joint
            if prob_joint < 1e-9:
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


def should_print_stats(
        s_start_t, s_last_print_t,
        valid_m, mi, m, last_mi,
        unknown_l, unknown_u,
        prev_p_i, next_p_i,
        tested,
        s_total_unknown,
        s_t_unk_low, s_t_unk_hgh,
        s_total_prp_tests,
        s_gap_out_of_sieve_prev, s_gap_out_of_sieve_next,
        s_best_merit_interval, s_best_merit_interval_m):
    s_stop_t = time.time()
    print_secs = s_stop_t - s_last_print_t
    if len(valid_m) in (1,10,30,100,300,1000) or len(valid_m) % 5000 == 0 \
            or mi == last_mi or print_secs > 1200:
        secs = s_stop_t - s_start_t

        print("\t{:3d} {:4d} <- unknowns -> {:-4d}\t{:4d} <- gap -> {:-4d}".format(
            m,
            unknown_l, unknown_u,
            prev_p_i, next_p_i))
        if mi <= 10 and secs < 6:
            return False

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
            100 * (1 - s_total_unknown / ((2 * args.sieve_length + 1) * len(valid_m))),
            100 * s_t_unk_low / s_total_unknown,
            100 * s_t_unk_hgh / s_total_unknown))
        if tested and args.run_prp:
            print("\t    prp tests {:<10d} (avg: {:.2f}) ({:.3f} tests/sec)".format(
                s_total_prp_tests, s_total_prp_tests / tested, s_total_prp_tests / secs))
            if s_gap_out_of_sieve_prev or s_gap_out_of_sieve_next:
                print("\t    fallback prev_gap {} ({:.1f}%), next_gap {} ({:.1f}%)".format(
                    s_gap_out_of_sieve_prev, 100 * s_gap_out_of_sieve_prev / tested,
                    s_gap_out_of_sieve_next, 100 * s_gap_out_of_sieve_next / tested))
            print("\t    best merit this interval: {:.3f} (at m={})".format(
                s_best_merit_interval, s_best_merit_interval_m))

        return True
    return False

def prime_gap_test(args):
    P = args.p
    D = args.d

    M = args.mstart
    M_inc = args.minc

    SL = sieve_length = args.sieve_length
    max_prime = args.max_prime

    min_merit = args.min_merit

    K, K_digits, K_bits, K_log = gap_utils.K_and_stats(args)
    M_log    = K_log + math.log(M)
    min_merit_gap = int(min_merit * M_log)
    print("K = {} bits, {} digits, log(K) = {:.2f}".format(
        K_bits, K_digits, K_log))
    print("Min Gap ~= {} (for merit > {:.1f})\n".format(
        min_merit_gap, min_merit))

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
    conn = sqlite3.connect(args.search_db)
    conn.row_factory = sqlite3.Row
    existing = load_existing(conn, args)
    print (f"Found {len(existing)} existing results")

    # used in next_prime
    assert P <= 80000
    # Very slow but works.
    primes = [2] + [p for p in range(3, 80000+1, 2) if gmpy2.is_prime(p)]
    K_primes = [p for p in primes if p <= P]

    # ----- Allocate memory for a handful of utility functions.

    # XXX: Cleanup after gmpy2.prev_prime.
    # Remainders of (p#/d) mod prime
    remainder   = [K % prime for prime in primes]

    # ----- Sieve stats
    prob_prime = 1 / M_log - 1 / (M_log * M_log)
    prob_prime_after_sieve = prob_prime_sieve_length(
        M, K, D, prob_prime, K_digits, K_primes, SL, max_prime)

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
    s_test_unknowns = []
    s_expected_prev = []
    s_expected_next = []
    s_expected_gap  = []
    s_experimental_side = []
    s_experimental_gap = []
    p_gap_side  = defaultdict(float)
    p_gap_comb  = defaultdict(float)
    p_merit_gap = []

    last_mi = M_inc - 1
    while math.gcd(M + last_mi, D) != 1:
        last_mi -= 1

    tested = 0
    for mi in range(M_inc):
        m = M + mi
        if math.gcd(m, D) != 1: continue

        valid_m.append(m)

        log_n = (K_log + math.log(m))

        # Read a line from the file
        line = unknown_file.readline()

        mtest, unknown_l, unknown_u, unknowns = gap_utils.parse_unknown_line(line)
        assert mtest == mi
        if len(valid_m) <= 3:
            s_test_unknowns.append(unknowns)

        prev_p_i = 0
        next_p_i = 0

        s_total_unknown += unknown_l + unknown_u
        s_t_unk_low += unknown_l
        s_t_unk_hgh += unknown_u

        if args.stats and (args.save_logs or args.num_plots):
            # NOTES: calculate_expected_gaps is really slow, only used to
            # doublecheck gap_stats (with --stats).
            e_prev, e_next, p_merit = calculate_expected_gaps(
                SL, min_merit_gap, prob_nth, prob_longer,
                log_n, unknowns, p_gap_side, p_gap_comb)
            s_expected_prev.append(e_prev)
            s_expected_next.append(e_next)
            s_expected_gap.append(e_prev + e_next)
            p_merit_gap.append(p_merit)

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

            merit = gap / log_n
            if m not in existing:
                save(conn, m, P, D, next_p_i, prev_p_i, merit)

            if gap in record_gaps or merit > min_merit:
                print("{}  {:.4f}  {} * {}#/{} -{} to +{}".format(
                    gap, merit, m, P, D, prev_p_i, next_p_i))

            if merit > s_best_merit_interval:
                s_best_merit_interval = merit
                s_best_merit_interval_m = m

            if should_print_stats(
                    s_start_t, s_last_print_t,
                    valid_m, mi, m, last_mi,
                    unknown_l, unknown_u,
                    prev_p_i, next_p_i,
                    tested,
                    s_total_unknown,
                    s_t_unk_low, s_t_unk_hgh,
                    s_total_prp_tests,
                    s_gap_out_of_sieve_prev, s_gap_out_of_sieve_next,
                    s_best_merit_interval, s_best_merit_interval_m):
                s_best_merit_interval = 0
                s_best_merit_interval_m = -1
                s_last_print_t = time.time()

    if args.num_plots or args.save_logs:
        if args.stats:
            # Not calculated
            p_record_gap = p_merit_gap

            stats_plots(
                args,
                min_merit_gap, record_gaps,

                valid_m,
                s_test_unknowns, prob_nth,
                s_expected_gap,
                s_expected_prev, s_expected_next,
                s_experimental_gap, s_experimental_side,
                p_gap_side, p_gap_comb,
                p_merit_gap, p_record_gap
            )

        # Load stats from gap_stats (fails if empty)
        try:
            (s_expected_prev_db, s_expected_next_db, s_expected_gap_db,
             p_merit_gap_db, p_record_gap_db,
             p_gap_comb_db, p_gap_side_db) = load_stats(conn, args)
            assert s_expected_prev_db and p_gap_comb_db
        except:
            # Failed to load
            if not args.stats and args.num_plots:
                print("Failed to load from DB so no plots.")
                exit(1)

        stats_plots(
            args,
            min_merit_gap, record_gaps,

            valid_m,
            s_test_unknowns, prob_nth,
            s_expected_gap_db,
            s_expected_prev_db, s_expected_next_db,
            s_experimental_gap, s_experimental_side,
            p_gap_side_db, p_gap_comb_db,
            p_merit_gap_db, p_record_gap_db
        )



if __name__ == "__main__":
    parser = get_arg_parser()
    args = parser.parse_args()
    gap_utils.verify_args(args, ".txt")

    # TeeLogger context if args.save_logs
    with gap_utils.logger_context(args):
        prime_gap_test(args)

