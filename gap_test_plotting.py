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

import gap_utils



#---- Stats to plot ----#
def stats_plots(
        args,
        min_merit_gap, record_gaps, prob_nth,
        data, misc):

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats

    # XXX: How to handle when only partial results

    def verify_no_trend(x, y):
        """There should be no trend based on m value"""
        # Have to adjust for Expected gap which has two data points for each m
        if len(y) == 2 * len(x):
            x = [xi for xi in x for side in ['l', 'r']]

        trend, _ = np.polyfit(x, y, 1)
        if trend > 2e-3:
            print("\n")
            print("NON-ZERO TREND: ", trend) # Verify that expected value doesn't vary with M.
            print("\n")

    def plot_hist(axis, d, color, marker='x', label='Observed', mult=1):
        hist_data = np.histogram(d, bins=100, density=True)
        axis.scatter(hist_data[1][:-1], hist_data[0]*mult, color=color, marker=marker, s=8, label=label)
        max_y = hist_data[0].max()

        axis.set_xlim(np.percentile(d, 0.01), np.percentile(d, 99.9))
        axis.set_ylim(top=1.2 * max_y)

    def plot_cdf(axis, d, color, label):
        n = len(d)
        d_sorted = np.sort(d)
        dist_label = f"Empirical CDF({label})" if 'P(' not in label else "CDF"
        axis.plot(
                d_sorted, np.arange(1, n+1) / n,
                color=color, label=dist_label)
        axis.set_xlim(np.percentile(d, 0.01), np.percentile(d, 99.9))

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

    def plot_prob_hist(axis, label, probs, max_x, color):
        x, w = zip(*sorted(((g, v) for g,v in probs.items() if v > 0 and g <= max_x)))
        n, _, _ = axis.hist(x, weights=w, bins=100, density=True,
                  label='Theoretical P(gap)', color=color, alpha=0.4)
        print (f"|P({label})| = {len(x)}, Sum(P({label})) = {sum(w):.1f}")
        return n

    def prob_histogram_all(axis, probs, experimental, label, c1='blueviolet', c2='peru'):
        max_x = 0.95 * max(probs.keys())

        if experimental:
            # If this is P(combined gap) and used --one-side-skipped
            # It's possible this should normalize by len(experimental) but I'm not sure
            no_skips = 2 * len(experimental) >= len(data.experimental_side)
            m_tested = len(data.experimental_side) - len(experimental)
            mult = 1 if no_skips else len(experimental) / m_tested
            label_hist = 'Observed' if no_skips else 'Observed (--one-side-skipped used)'

            # Experimental values
            plot_hist(axis, experimental, 'peru', label=label_hist, mult=mult)
            if not no_skips:
                print("\nNormalizing {} gaps by {} m's tested (--one-side-skipped)\n".format(
                    len(experimental), m_tested))

            # Expected value
            add_expected_value(axis, experimental, c2, label)

            # Need to balance x which is generally longer for probs then 99.9 percentile of experimental
            max_x = min(1.1 * np.percentile(experimental, 99.9), max_x)

        # Theoretical probabilities
        tn = plot_prob_hist(axis, label, probs, max_x, c1)

        axis.set_xlim(left=0, right=max_x)
        axis.set_yscale('log')

        axis.legend(loc='upper right')

        return 0.9 * np.min(tn), 1.2 * np.max(tn)

    def fit_normal_dist(axis, d):
        x_start = max(0, 0.95 * np.percentile(d, 1))
        x_end   = 1.05 * np.percentile(d, 99)
        gap_span = np.linspace(x_start, x_end, 400)

        mu, std = stats.norm.fit(d)
        p = stats.norm.pdf(gap_span, mu, std)
        axis.plot(gap_span, p, color=color)

    def add_expected_value(axis, d, color, label):
        E = np.mean(d)
        label_e = f"E({label}) = {E:.0f}"
        axis.axvline(x=E, ymax=1.0/1.2, color=color, label=label_e)

    e_gap_n = len(data.experimental_gap)
    if len(data.expected_gap) != e_gap_n:
        print("experimental_gap size mismatch", len(data.expected_gap), e_gap_n)
    if e_gap_n:
        slope, intercept, R, _, _ = stats.linregress(
            data.expected_gap[:e_gap_n], data.experimental_gap)
        print()
        print("R^2 for expected gap: {:.3f}, gap = {:.1f} + {:.3f} * expected".format(
            R**2, intercept, slope))
        print()

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
            axis_prev = fig.add_subplot(gs[0, 0])
            axis_next = fig.add_subplot(gs[0, 1])
            axis_xnth = fig.add_subplot(gs[1, 0])

            # prob_prev, prev_next for individual m
            # See Prob_nth in gap_stats
            colors = plt.cm.tab10
            for i, (mi, (u_p, u_n)) in enumerate(misc.test_unknowns.items()):
                label = f"m={args.mstart + mi}"
                color = colors(i)
                plot_prob_nth(axis_prev, u_p, color, label)
                plot_prob_nth(axis_next, u_n, color, label)

                x_p = range(1, len(u_p) + 1)
                x_n = range(1, len(u_n) + 1)
                axis_xnth.plot(x_p, u_p, marker='.', color=color, label=label)
                axis_xnth.plot(x_n, u_n, marker='.', color=color, label=label)

            axis_xnth.legend()
            axis_prev.legend(loc='upper left')
            axis_next.legend(loc='upper right')
            axis_prev.set_yscale('log')
            axis_next.set_yscale('log')
            axis_prev.set_xlim(-args.sieve_length, 0)
            axis_next.set_xlim(0, args.sieve_length)

        min_y, max_y = prob_histogram_all(
                axis_prob_gap, misc.prob_gap_side, data.experimental_side, 'next')
        axis_prob_gap.set_xlim(0, args.sieve_length)
        axis_prob_gap.set_ylim(bottom=max(10 ** -8, min_y / 10))
        #print(f"Min Prob(gap side): {min_y:.2e}")

        for e_data, color, label in (
                (data.expected_prev, 'lightskyblue', 'prev'),
                (data.expected_next, 'tomato', 'next'),
        ):
            plot_hist(axis_expected_gap,          e_data, color)
            add_expected_value(axis_expected_gap, e_data, color, label)
            fit_normal_dist(axis_expected_gap, e_data)
            axis_expected_gap.legend(loc='upper left')

            # CDF of gap <= x
            plot_cdf(axis_cdf_gap, e_data, color, label)

    if args.num_plots > 1:
        # Plot 2: Gap(combined):
        #   [ prob all m,               expected,   cdf ]
        #   [ sum(prob) & count,  dist,       cdf ] (for >min merit)
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
        min_y = max(10 ** -8, min_y)
        axis_prob_comb.set_ylim(bottom=min_y, top=max_y)

        color = 'seagreen'
        axis = axis_expected_comb
        e_gap = data.expected_gap
        plot_hist(axis,          e_gap, color)
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

            # This assumes that experimental_gap is indexed the same as prob_data
            # This is not true unless --prp-top-percent is 100
            if len(prob_data) == len(data.experimental_gap):
                zipped = list(zip(prob_data, data.experimental_gap))
            else:
                print("Not all gaps are present (--prp-top-percent < 100) Not able to show 'Sum(P(...))'")
                zipped = list(zip(prob_data, [0 for i in range(len(prob_data))]))

            p_gap_merit_sorted, _ = zip(*sorted(zipped, reverse=True))
            p_gap_merit_ord, gap_real_ord = zip(*zipped)

            #Experimental
            if row == 1:
                cum_count = np.cumsum(np.array(gap_real_ord) > min_merit_gap)
            else:
                cum_count = np.cumsum(np.array([g in record_gaps for g in gap_real_ord]))

            print(f"{label:20} | sum(P) = {sum(prob_data):.4f}, count(experimental) = {cum_count[-1]}")

            tests = list(range(1, len(p_gap_merit_ord)+1))
            if cum_count[-1] > 0:
                axis.plot(tests, cum_count, label='Count ' + label)

            # Theoretical
            cum_sum_p = np.cumsum(p_gap_merit_ord)
            cum_sum_p_sorted = np.cumsum(p_gap_merit_sorted)

            axis.plot(tests, cum_sum_p, label=f'Sum(P({label}))')
            axis.plot(tests, cum_sum_p_sorted, label=f'Sum(P({label})) (best first)')
            axis.legend(loc='upper left')

            # Hist
            axis = fig.add_subplot(gs[row, 1])
            plot_hist(axis, prob_data, color)

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
                (data.prob_merit_gap, f'P(gap > min_merit({args.min_merit}))', 'dodgerblue'),
        )):
            if not d: continue

            if plot_i == 0:
                axis = fig.add_subplot(gs[0, 0])
                dist_axis = fig.add_subplot(gs[0, 1])
            elif plot_i == 1:
                # Reuse from plot_i == 0
                pass
            elif plot_i == 2:
                axis = fig.add_subplot(gs[1, 0])
                dist_axis = fig.add_subplot(gs[1, 1])
            elif plot_i == 3:
                axis = fig.add_subplot(gs[2, 0])
                dist_axis = fig.add_subplot(gs[2, 1])
            else:
                assert False

            assert not label.endswith("gap"), label
            verify_no_trend(data.valid_mi, d)

            plot_hist(axis, d, color)

            if plot_i <= 2:
                # Fit normal distribution to data
                fit_normal_dist(axis, d)
                # Plot a line for expected value
                add_expected_value(axis, d, color, label)
                axis.legend()
            else:
                axis.legend([label])

            # Cumulative sum of probability by gap
            plot_cdf(dist_axis, d, color, label)

        # Probability of Gap=X | Theoretical + Experimental
        for axis, label, theory, c1, experimental, c2 in [
                (axis_one_gap, 'next/prev gap',
                    misc.prob_gap_side, 'blueviolet',
                    data.experimental_side, 'sandybrown'),
                (axis_combined_gap, 'gap',
                    misc.prob_gap_comb, 'seagreen',
                    data.experimental_gap, 'peru'),
        ]:
            min_y, max_y = prob_histogram_all(axis, theory, experimental, label, c1, c2)
            axis.set_ylim(bottom=max(10 ** -8, min_y / 10))
            axis.legend(loc='upper right')

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
        args, conn, sc, data_db, misc_db,
        min_merit_gap, record_gaps, prob_prime_after_sieve):

    assert data_db.expected_prev

    # Geometric distribution
    prob_nth = []
    prob_gap_longer = 1
    while prob_gap_longer > 1e-13:
        prob_nth.append(prob_gap_longer * prob_prime_after_sieve)
        prob_gap_longer *= (1 - prob_prime_after_sieve)
    assert min(prob_nth) > 0

    stats_plots(args, min_merit_gap, record_gaps, prob_nth, data_db, misc_db)
