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
        DB=data/test_1511_${sl}.db;
        time ./combined_sieve --save-unknowns --unknown-filename $FN
        rm -f "$DB";
        sqlite3 "$DB" < schema.sql;
        time ./gap_stats --save-unknowns --search-db "$DB" --unknown-filename $FN
        #time ./gap_stats --unknown-filename $FN -q -q

        sqlite3 "$DB" <<EOL
        SELECT sieve_length,COUNT(*),SUM(prob_record),mean,
               AVG((prob_record-mean)*(prob_record-mean)) as variance
        FROM m_stats,
                (SELECT AVG(prob_record) as mean FROM m_stats),
                (SELECT sieve_length FROM range)
EOL
    done

    python -c 'import misc.paper_plots; misc.paper_plots.plot_record_vs_sl()'
    """

    import sqlite3

    PERCENT = 20
    P = 1511
    #P = 1667

    m_values = None
    data = []

    for fn in glob.glob(f"data/test_{P}_[0-9]*.db"):
        sl = int(re.match(f"data/test_{P}_([0-9]*).db", fn).group(1))
        if sl % 1500 != 0:
            continue
        print ("\tSL:", sl)
        with sqlite3.connect(fn) as conn:
            rv = conn.execute("SELECT m,prob_record FROM m_stats")
            ms, probs = zip(*[row for row in rv if len(row) == 2])
            if m_values is None:
                m_values = ms
            else:
                assert m_values == ms, "m mismatch between {sl} and other"
            data.append((sl, probs))

    # Could consider counting inversions but O(n^2) and this metric also works

    def top_n_percent(probs, percent):
        best_probs = sorted(probs, reverse=True)
        prob_threshold = best_probs[round(len(best_probs) * percent / 100)]
        del best_probs
        return [m for m, prob in zip(m_values, probs) if prob >= prob_threshold]

    index = [i for i in range(len(data)) if data[i][0] == 45000][0]
    assert data[index][0] == 45000
    best_guess_probs = {m: prob for m, prob in zip(m_values, data[index][1])}

    sum_top_percent = []
    for sl, probs in data:
        # Determine what m would be used and what
        top = top_n_percent(probs, PERCENT)
        sum_prob = sum(best_guess_probs[m] for m in top) / len(top)
        sum_top_percent.append((sl, sum_prob))
    sum_top_percent.sort()

    print (sum_top_percent)

    for row in sum_top_percent:
        print (row[0], ",",  row[1])
