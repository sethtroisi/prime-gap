import matplotlib.pyplot as plt
import os

def plot_record_vs_plimit():
    """
    Plots probability of record as P_limit increases

    Takes [((prime1, prob_record1), (prime2, prob_record2), ...), ...]
    """
    """
    Commands run were:

    for m in 293609 811207 183047; do
        time primegapverify/large_sieve $m 1511 312270 -50000 100000 100''000''000''000 > \
            ../prime-gap/unknowns/1511_312270_${m}_1_s50000_l100000M.txt;
    done

    make gap_stats
    for m in 293609 811207 183047; do
        time ./gap_stats --unknown-filename unknowns/1511_312270_${m}_1_*.txt | tee ${m}_test.txt;
    done

    python -c 'import misc.paper_plots; misc.paper_plots.plot_record_vs_plimit()'
    """

    plt.rcParams.update({'mathtext.default':  'regular' })

    # Set up subplots.
    fig = plt.figure(
        "Probability of record gap vs P_limit",
        constrained_layout=True,
        figsize=(8, 5))

    def plot_record_probs(plimit, prob, color, label, marker):
        plt.scatter(plimit, prob, marker=marker, s=18, color=color, label=label)
        plt.plot(plimit, prob, color=color)

    # prob_prev, prev_next for individual m
    # See Prob_nth in gap_stats
    colors = plt.cm.tab10
    #for i, (u_p, u_n) in enumerate(misc.test_unknowns):
    #    label = f"m={valid_mi[i]}"
    #    plot_prob_nth(axis_prev, u_p, colors(i), label)
    #    plot_prob_nth(axis_next, u_n, colors(i), label)

    plt.xscale('log')

    plt.xlabel("$P_{limit}$")
    plt.ylabel("P(record gap|{partial sieve, $P_{limit}$})")

    for i, m in enumerate([293609, 811207, 183047]):
        for j, sl in enumerate([15, 20, 25]):
            fn = f"data/{m}_{sl}_test.txt"
            if not os.path.exists(fn):
                continue
            with open(fn) as data_file:
                lines = [line for line in data_file.readlines() if line and line[0].isdigit()]
                if not lines:
                    print(fn, "is empty")
                    continue

                p_limit, probs = zip(*[list(map(float, line.split(", "))) for line in lines])
                plot_record_probs(p_limit, probs, colors(i), f"{m=} ({sl=}", "Dhp"[j])

    plt.legend(loc='upper left')

    #if args.save_logs:
    #    png_path = gap_utils.transform_unknown_filename(
    #        args.unknown_filename, "logs/", ".png")
    #    if not png_path.startswith("logs/"):
    #        print ("No 'logs/' directory to save figure into")
    #    plt.savefig(png_path, dpi=1080//8)

    plt.show()
    plt.close()


def plot_record_vs_sl():
    """
    Plots probability of record at several different sieve-lengths

    make clean combined_sieve gap_stats
    for sl in 15000 30000 50000; do
    for sl in {1500..15000..1500}; do
        FN=1511_312270_1_1000000_s${sl}_l10000M.txt
        time ./combined_sieve --save-unknowns --unknown-filename $FN
        DB=data/test_1511_${sl}.db;
        rm -i "$DB";
        sqlite3 "$DB" < schema.sql;
        time ./gap_stats --save-unknowns --search-db "$DB" --unknown-filename $FN
    done

    python -c 'import misc.paper_plots; misc.paper_plots.plot_record_vs_sl()'
    """
    from collections import defaultdict

    m = None
    data = []

    for sl in (15000, 30000, 50000):
        with sqlite3.connect(f"DB=data/test_1511_{sl}.db") as conn:
            rv = conn.execute("SELECT m,prob_record FROM m_stats")
            ms, probs = zip(*[row for row in rv])
            if m is None:
                m = ms
            else:
                assert m == ms, "m mismatch between {sl} and other"
            data.pushback(probs)

