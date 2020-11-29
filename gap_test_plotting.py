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

import itertools

import gap_utils



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
    slope, intercept, R, _, _ = stats.linregress(data.expected_gap[:egap_n], data.experimental_gap)
    print ()
    print ("R^2 for expected gap: {:.3f}, gap = {:.1f} + {:.3f} * expected".format(
        R**2, intercept, slope))
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
                marker='.', s=12, color=color)
            # calculate expected value = sum(i * prob(i))
            E = sum(u * p for u, p in zip(unknowns, prob_nth))
            axis.axvline(
                x=E, ymax=1.0/1.2,
                color=color, label=f"E({label}) = {E:.0f}")

        if misc.test_unknowns:
            # prob_prev, prev_next for individual m
            # See Prob_nth in gap_stats
            colors = plt.cm.tab10
            for i, (u_p, u_n) in enumerate(misc.test_unknowns):
                label = f"m={valid_m[i]}"
                plot_prob_nth(axis_prev, u_p, colors(i), label)
                plot_prob_nth(axis_next, u_n, colors(i), label)

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
        #print(f"Min Prob(gap side): {min_y:.2e}")

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
            figsize=(12, 8), dpi=150)
        gs = fig.add_gridspec(3, 3)

        axis_prob_comb     = fig.add_subplot(gs[0, 0])
        axis_expected_comb = fig.add_subplot(gs[0, 1])
        axis_cdf_comb      = fig.add_subplot(gs[0, 2])

        # Combining all probs for a pseudo distribution of P(gap combined)
        min_y, max_y = prob_histogram_all(
                axis_prob_comb, misc.prob_gap_comb, data.experimental_gap, 'gap')
        #print(f"Min Prob(gap comb): {min_y:.2e}")
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
            axis.set_xlabel(" # of m's tested")
            axis.set_ylabel(f'Sum(P(gap {label})')

            #assert len(prob_data) == len(data.experimental_gap)
            zipped = list(zip(prob_data, data.experimental_gap))

            p_gap_merit_sorted, _ = zip(*sorted(zipped, reverse=True))
            p_gap_merit_ord, gap_real_ord = zip(*zipped)

            print(f"{label:20} | sum(P) = {sum(prob_data):.4f}")

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
        plt.xlabel(" # of m's tested")
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
        png_path = gap_utils.transform_unknown_filename(
            args.unknown_filename, "logs/", ".png")
        if not png_path.startswith("logs/"):
            print ("No 'logs/' directory to save figure into")

        plt.savefig(png_path, dpi=1080//8)

    if args.num_plots:
        plt.show()

    plt.close()


def plot_stuff(
        args, conn, sc, data, misc,
        data_db, misc_db,
        min_merit_gap, record_gaps, prob_nth):

    assert data_db.expected_prev
    assert misc_db.prob_gap_comb, len(misc.prob_gap_comb)

    if args.stats:
        # Not calculated
        data.prob_record_gap = data.prob_merit_gap

        stats_plots(
            args, min_merit_gap, record_gaps, prob_nth,
            data.valid_m, data, misc)

    # test_unknowns come from unknown-file not DB.
    misc_db.test_unknowns = misc.test_unknowns

    stats_plots(
        args, min_merit_gap, record_gaps, prob_nth,
        data.valid_m, data_db, misc_db)
