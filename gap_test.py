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
import contextlib
import logging
import math
import os.path
import re
import sqlite3
import sys
import time

import gmpy2


class TeeLogger:
    def __init__(self, fn, sysout):
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)-15s %(levelname)s: %(message)s",
            handlers=[
                logging.StreamHandler(sys.stdout),
                logging.FileHandler(fn, mode="w"),
            ]
        )
        self.logger = logging.getLogger()

    def write(self, msg):
        if msg and not msg.isspace():
            self.logger.info(msg)

    def flush(self):
        self.logger.flush()


def get_arg_parser():
    parser = argparse.ArgumentParser('Test prime gaps generated by gap_search')
    parser.add_argument('--mstart',  type=int)
    parser.add_argument('--minc',    type=int)
    parser.add_argument('-p',       type=int)
    parser.add_argument('-d',       type=int)

    parser.add_argument('--sieve-length', type=int,
        help="how large the search window was (must match gap_search param)")
    parser.add_argument('--sieve-range',  type=int,
        help="Use primes <= sieve-range for checking composite (must match gap_search param)")

    parser.add_argument('--unknown-filename', type=str,
        help="determine mstart, minc, p, d, sieve-length, and sieve-range"
             " from unknown-results filename")

    parser.add_argument('--min-merit',    type=int, default=10,
        help="only display prime gaps with merit >= minmerit")
    parser.add_argument('--sieve-only',  action='store_true',
        help="only sieve ranges, don't run PRP. useful for benchmarking")

    parser.add_argument('--plots', action='store_true',
        help="Show plots about distributions")

    parser.add_argument('--save-logs', action='store_true',
        help="Save logs and plots about distributions")

    return parser

def verify_args(args):
    if args.unknown_filename:
        fn = args.unknown_filename
        if not os.path.exists(fn):
            print ("\"{}\" doesn't exist".format(fn))
            sys.exit(1)
        match = re.match(
            "^(\d+)_(\d+)_(\d+)_(\d+)_s(\d+)_l(\d+)M.txt",
            os.path.basename(fn))
        if not match:
            print ("\"{}\" doesn't match unknown file format".format(fn))
            sys.exit(1)

        ms, p, d, mi, sl, sr = map(int, match.groups())
        args.mstart = ms
        args.minc = mi
        args.p = p
        args.d = d
        args.sieve_length = sl
        args.sieve_range = sr

    if args.sieve_range <= 4000:
        args.sieve_range *= 10 ** 6

    for arg in ('mstart', 'minc', 'p', 'd', 'sieve_length', 'sieve_range'):
        if arg not in args or args.__dict__[arg] in (None, 0):
            print ("Missing required argument", arg)
            sys.exit(1)

    fn = "{}_{}_{}_{}_s{}_l{}M.txt".format(
        args.mstart, args.p, args.d, args.minc,
        args.sieve_length, args.sieve_range // 10 ** 6)

    if args.unknown_filename:
        assert fn == os.path.basename(args.unknown_filename), (fn, args.unknown_filename)
    else:
        args.unknown_filename = fn


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

    count_coprime = SL-1
    for i in range(1, SL):
        for prime in K_primes:
            if (i % prime) == 0 and (D % prime) != 0:
                count_coprime -= 1
                break

    chance_coprime_composite = 1 - prob_prime / prob_prime_coprime
    prob_gap_shorter_hypothetical = chance_coprime_composite ** count_coprime

    # count_coprime already includes some parts of unknown_after_sieve
    print("\t{:.3f}% of sieve should be unknown ({}M) ~= {:.0f}".format(
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

def calculate_expected_gap(composites, SL, prob_prime_after_sieve, log_m):
    expected_length = 0
    prob_gap_longer = 1
    for v in composites:
        expected_length += abs(v) * prob_gap_longer * prob_prime_after_sieve
        prob_gap_longer *= 1 - prob_prime_after_sieve

    # expected to encounter a prime at distance ~= ln(start)
    expected_length += (SL + log_m) * prob_gap_longer
    return expected_length


def determine_next_prime_i(m, K, composites, SL):
    center = m * K
    tests = 0

    for i in composites:
        tests += 1;
        if gmpy2.is_prime(center + i):
            next_p_i = i
            break
    else:
        # Using fallback to slower gmp routine
        next_p_i = int(gmpy2.next_prime(center + SL - 1) - center)

    return tests, next_p_i


def determine_prev_prime_i(m, K, composites, SL, primes, remainder):
    center = m * K
    tests = 0

    for i in composites:
        assert i < 0
        tests += 1;
        if gmpy2.is_prime(center + i):
            prev_p_i = -i
            break
    else:
        # Medium ugly fallback.
        for i in range(SL, 5*SL+1):
            composite = False
            for prime, remain in zip(primes, remainder):
                modulo = (remain * m) % prime
                if i % prime == modulo:
                    composite = True
                    break
            if not composite:
                if gmpy2.is_prime(center - i):
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

    run_prp = not args.sieve_only

    min_merit = args.min_merit

    K = gmpy2.primorial(P)
    assert K % D == 0
    K //= D

    K_digits = gmpy2.num_digits(K, 10)
    K_bits   = gmpy2.num_digits(K, 2)
    K_log    = float(gmpy2.log(K))
    M_log    = K_log + math.log(M)
    print("K = {} bits, {} digits, log(K) = {:.2f}".format(
        K_bits, K_digits, K_log))
    print("Min Gap ~= {} (for merit > {:.1f})\n".format(
        int(min_merit * M_log), min_merit))

    # ----- Open Output file
    print("\tLoading unknowns from '{}'".format(args.unknown_filename))
    print()
    unknown_file = open(args.unknown_filename, "r")

    # used in next_prime
    assert P <= 80000
    # Very slow but works.
    primes = [2] + [p for p in range(3, 80000+1, 2) if gmpy2.is_prime(p)]
    K_primes = [p for p in primes if p <= P]

    # ----- Allocate memory for a handful of utility functions.

    # Remainders of (p#/d) mod prime
    remainder   = [K % prime for prime in primes]

    # ----- Sieve stats
    prob_prime_after_sieve = prob_prime_sieve_length(
        M, D, 1 / M_log, K_digits, K_primes, SL, sieve_range)

    # ----- Main sieve loop.
    print("\nStarting m={}".format(M))
    print()

    # Used for various stats
    s_start_t = time.time()
    s_total_unknown = 0
    s_t_unk_low = 0
    s_t_unk_hgh = 0
    s_total_prp_tests = 0
    s_gap_out_of_sieve_prev = 0
    s_gap_out_of_sieve_next = 0
    s_best_merit_interval = 0
    s_best_merit_interval_m = 0

    s_expected_prev = []
    s_expected_next = []
    s_expected_gap  = []
    s_experimental_side = []
    s_experimental_gap = []

    for mi in range(M_inc):
        m = M + mi
        # TODO if gcd(m, d) != 1 continue?

        log_m = (K_log + math.log(m))

        composite = [[], []]

        unknown_l = -1
        unknown_u = -1
        prev_p_i = 0
        next_p_i = 0

        # Read a line from the file
        line = unknown_file.readline()
        start, c_l, c_h = line.split("|")

        match = re.match(r"^([0-9]+) : -([0-9]+) \+([0-9]+)", start)
        assert match, start
        mtest, unknown_l, unknown_u = map(int, match.groups())

        composite[0] = list(map(int,c_l.strip().split(" ")))
        composite[1] = list(map(int,c_h.strip().split(" ")))

        unknown_l_test = len(composite[0])
        unknown_u_test = len(composite[1])
        assert unknown_l == unknown_l_test, (unknown_l, unknown_l_test, "\t", start)
        assert unknown_u == unknown_u_test, (unknown_u, unknown_u_test, "\t", start)

        s_total_unknown += unknown_l + unknown_u
        s_t_unk_low += unknown_l
        s_t_unk_hgh += unknown_u

        if args.save_logs or args.plots:
            e_prev = calculate_expected_gap(composite[0], SL, prob_prime_after_sieve, log_m)
            e_next = calculate_expected_gap(composite[1], SL, prob_prime_after_sieve, log_m)
            s_expected_prev.append(e_prev)
            s_expected_next.append(e_next)
            s_expected_gap.append(e_prev + e_next)

        if run_prp:
            tests, prev_p_i = determine_prev_prime_i(m, K, composite[0], SL, primes, remainder)
            s_total_prp_tests += tests
            s_gap_out_of_sieve_prev += prev_p_i >= SL

            tests, next_p_i = determine_next_prime_i(m, K, composite[1], SL)
            s_total_prp_tests += tests
            s_gap_out_of_sieve_next += next_p_i >= SL

            assert prev_p_i > 0 and next_p_i > 0
            gap = int(next_p_i + prev_p_i)
            s_experimental_gap.append(gap)
            s_experimental_side.append(next_p_i)
            s_experimental_side.append(prev_p_i)

            merit = gap / log_m
            if merit > min_merit:
                # TODO write to file.
                print("{}  {:.4f}  {} * {}#/{} -{} to +{}".format(
                    gap, merit, m, P, D, prev_p_i, next_p_i))

            if merit > s_best_merit_interval:
                s_best_merit_interval = merit
                s_best_merit_interval_m = m

        if mi in (1,10,100,500,1000, M_inc-1) or m % 5000 == 0:
            s_stop_t = time.time()
            secs = s_stop_t - s_start_t

            print("\t{:3d} {:4d} <- unknowns -> {:-4d}\t{:4d} <- gap -> {:-4d}".format(
                m,
                unknown_l, unknown_u,
                prev_p_i, next_p_i))
            if mi <= 10: continue

            # Stats!
            tests = mi + 1
            print("\t    tests     {:<10d} ({:.2f}/sec)  {:.0f} seconds elapsed".format(
                tests, tests / secs, secs))
            print("\t    unknowns  {:<10d} (avg: {:.2f}), {:.2f}% composite  {:.2f}% <- % -> {:.2f}%".format(
                s_total_unknown, s_total_unknown / tests,
                100 * (1 - s_total_unknown / (2 * (sieve_length - 1) * tests)),
                100 * s_t_unk_low / s_total_unknown,
                100 * s_t_unk_hgh / s_total_unknown))
            if run_prp:
                print("\t    prp tests {:<10d} (avg: {:.2f}) ({:.1f} tests/sec)".format(
                    s_total_prp_tests, s_total_prp_tests / tests, s_total_prp_tests / secs))
                print("\t    fallback prev_gap {} ({:.1f}%), next_gap {} ({:.1f}%)".format(
                    s_gap_out_of_sieve_prev, 100 * s_gap_out_of_sieve_prev / tests,
                    s_gap_out_of_sieve_next, 100 * s_gap_out_of_sieve_next / tests))
                print("\t    best merit this interval: {:.3f} (at m={})".format(
                    s_best_merit_interval, s_best_merit_interval_m))

            s_best_merit_interval = 0
            s_best_merit_interval_m = -1

    if args.plots or args.save_logs:
        import numpy as np
        import matplotlib.pyplot as plt
        from scipy import stats

        if run_prp:
            corr, _ = stats.pearsonr(s_expected_gap, s_experimental_gap)
            print ("Pearson's correlation for expected gap: {:.3f}".format( corr))

        x = list(range(M, M + M_inc))
        rows = 5 if run_prp else 4

        def plot_gaps_with_trendline(data, ylabel):
            trendline = np.polyfit(x, data, 1)
            trend = np.poly1d(trendline)

            plt.plot(x, data)
            trend_equation = str(trend)
            plt.plot(x, trend(x), label=trend_equation,
                color='lightblue', linestyle='dashed')

            # Histogram of gap sizes
            #hist_data = np.histogram(data, bins=50)
            #plt.plot(4 * hist_data[0], hist_data[1][:-1],
            #    color='lightblue', marker='x', linestyle='');

            plt.ylabel(ylabel)
            plt.legend(loc='upper right')

        # Set up subplots.
        fig3 = plt.figure(constrained_layout=True)
        gs = fig3.add_gridspec(rows, 2)

        ax1 = fig3.add_subplot(gs[0, 0])
        plot_gaps_with_trendline(s_expected_prev, "E(previous gap)")

        ax2 = fig3.add_subplot(gs[0, 1])
        ax1.get_shared_y_axes().join(ax1, ax2)
        plot_gaps_with_trendline(s_expected_next, "E(next gap)")

        fig3.add_subplot(gs[1, :])
        plot_gaps_with_trendline(s_expected_gap, "E(combined gap)")

        if run_prp:
            fig3.add_subplot(gs[2, :])
            plot_gaps_with_trendline(s_experimental_gap, "gap")

        fig3.add_subplot(gs[2 + run_prp:, :])
        for d, label, color in [
                (s_expected_prev, 'prev', 'blueviolet'),
                (s_expected_next, 'next', 'peru'),
                (s_expected_gap,  'expected', 'dodgerblue'),
                (s_experimental_side, 'one side gap', 'darkorange'),
                (s_experimental_gap, 'gap', 'forestgreen')]:

            if not d:
                continue
            max_gap = max(max(s_expected_gap), max(s_experimental_gap, default=0))

            hist_data = np.histogram(d, bins=80, density=True)
            plt.plot(hist_data[1][:-1], hist_data[0], color=color, marker='x')

            gap_span = np.linspace(0, max_gap, 400)
            if 'gap' in label:
                mu, std = stats.norm.fit(d)
                p = stats.norm.pdf(gap_span, mu, std)
            else:
                a = stats.gamma.fit(d)
                p = stats.gamma.pdf(gap_span, a[0], loc=a[1], scale=a[2])

            plt.plot(gap_span, p, label=label, color=color)

        # Don't show as many zeros
        plt.ylim(bottom=1e-5)
        plt.legend(loc='upper left')

        if args.save_logs:
            plt.savefig(args.unknown_filename + ".png", dpi=200)

        if args.plots:
            plt.show()

        plt.close()


if __name__ == "__main__":
    parser = get_arg_parser()
    args = parser.parse_args()
    verify_args(args)

    context = contextlib.suppress()
    if args.save_logs:
        assert args.unknown_filename
        log_fn = args.unknown_filename + '.log'
        context = contextlib.redirect_stdout(TeeLogger(log_fn, sys.stdout))

    with context:
        prime_gap_test(args)

