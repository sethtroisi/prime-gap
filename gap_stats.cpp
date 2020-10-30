// Copyright 2020 Seth Troisi
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>

#include <gmp.h>
#include <sqlite3.h>

#include "gap_common.h"


using std::cout;
using std::endl;
using std::pair;
using std::tuple;
using std::vector;
using namespace std::chrono;


typedef const double cdouble;


// Limits the size of record list
const uint32_t MAX_GAP = 1'000'000;
const float    GAP_INF = std::numeric_limits<float>::max();   // log(starting_prime)

// Generated from https://primegap-list-project.github.io/lists/missing-gaps/
// Range of missing gap to search, values are loaded from records_db.
const uint32_t MISSING_GAPS_LOW  = 113326;
const uint32_t MISSING_GAPS_HIGH = 132928;


class ProbNth {
    public:
        // See `prob_nth_prime`
        // Probability that i'th unknown is prime (inside sieve)
        vector<float> prime_nth_sieve;
        // Probability that prime is >= i'th (inside sieve)
        vector<float> great_nth_sieve;

        /**
         * Probability that prev_prime & next_prime have X unknown composites in middle
         * prob_combined_sieve[i+j] = prime * (1 - prime)^i * (1 - prime)^j * prime
         *                          = prime ^ 2 * (1 - prime)^(i+j)
         */
        vector<float> combined_sieve;

        /**
         * Probability of gap[i] on prev side, next gap > SL, and is record.
         * Sum(prob_combined_sieve[i-1 + unknowns[side] + j-1,
         *      where unknowns[i] + exentended[j] is record)
         *
         * Because this uses m, need to handle prev and next side differently
         */
        int wheel_d;
        map<int, vector<float>> extended_record_high;
        float average_coprime;
        float prob_extended;
};


void prime_gap_stats(struct Config config);
bool is_range_already_processed(const struct Config& config);

static double average_v(vector<float> &a) {
    return std::accumulate(a.begin(), a.end(), 0.0) / a.size();
}

static void prob_stats(char const *name, vector<float> &probs) {
    vector<float> sorted = probs;
    std::sort(sorted.begin(), sorted.end(), std::greater<>());

    for (auto percent : {1, 5, 10, 20, 50, 100}) {
        size_t count = probs.size() * percent / 100;
        if (count == 0)
            continue;

        double sum_prob = std::accumulate(sorted.begin(), sorted.begin() + count, 0.0);
        printf("\t%-7s: top %3d%% (%5ld) => sum(prob) = %.2e (avg: %.2e)\n",
            name, percent, count, sum_prob, sum_prob / count);

        if (sorted[count-1] == 0)
            break;
    }
}


int main(int argc, char* argv[]) {

    Config config = Args::argparse(argc, argv);

    if (config.verbose >= 3) {
        printf("\tCompiled with GMP %d.%d.%d\n\n",
            __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);
    }

    if (config.valid == 0) {
        Args::show_usage(argv[0]);
        return 1;
    }

    if (config.save_unknowns == 0) {
        printf("Not saving unknowns (--save-unknowns=0)\n");
    } else if (is_range_already_processed(config)) {
        cout << "Range already processed!" << endl;
        return 1;
    }

    prime_gap_stats(config);
}


//---------------------------------------------------------------------------//


vector<float> get_record_gaps(const struct Config& config) {
    vector<float> records(MAX_GAP, GAP_INF);

    DB db(config.records_db.c_str());

    /* Create SQL statement */
    char sql[] = "SELECT gapsize, merit FROM gaps";
    char *zErrMsg = 0;

    /* Execute SQL statement */
    int rc = sqlite3_exec(db.get_db(), sql, [](void* recs, int argc, char **argv, char **azColName)->int {
        uint64_t gap = atol(argv[0]);
        if (gap < MAX_GAP) {
            // Recover log(startprime)
            (*static_cast<vector<float>*>(recs))[gap] = gap / atof(argv[1]);
        }
        return 0;
    }, (void*)&records, &zErrMsg);

    if( rc != SQLITE_OK ) {
        printf("SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
    }

    return records;
}


void load_possible_records(
        const double N_log,
        const vector<float> &records,
        vector<uint32_t> &poss_record_gaps) {
    for (size_t g = 2; g < MAX_GAP; g += 2) {
        // Ignore the infintesimal odds of finding >merit 35 gap.
        if (g / N_log > 35) {
            break;
        }

        if (records[g] > N_log) {
            poss_record_gaps.push_back(g);
        }
    }

    assert(is_sorted(poss_record_gaps.begin(), poss_record_gaps.end()));
    assert(poss_record_gaps.front() <= MISSING_GAPS_LOW);
}


bool is_range_already_processed(const struct Config& config) {
    DB db_helper(config.search_db.c_str());
    sqlite3 *db = db_helper.get_db();

    uint64_t hash = db_helper.config_hash(config);
    char sql[200];
    sprintf(sql, "SELECT count(*) FROM range WHERE rid = %ld and time_stats > 0", hash);
    char *zErrMsg = 0;

    int count = 0;
    int rc = sqlite3_exec(db, sql, [](void* data, int argc, char **argv, char **azColName)->int {
        assert( argc == 1 );
        *static_cast<int*>(data) = atoi(argv[0]);
        return 0;
    }, &count, &zErrMsg);

    if (rc != SQLITE_OK) {
        printf("\nrange SELECT failed '%s' | %d: '%s'\n",
            zErrMsg, rc, sqlite3_errmsg(db));
        exit(1);
    }
    return count > 0;
}


double get_range_time(const struct Config& config) {
    DB db_helper(config.search_db.c_str());
    sqlite3 *db = db_helper.get_db();

    uint64_t hash = db_helper.config_hash(config);
    char sql[200];
    sprintf(sql, "SELECT time_sieve + time_stats FROM range WHERE rid = %ld and time_sieve > 0", hash);
    char *zErrMsg = 0;

    double time = 0;
    int rc = sqlite3_exec(db, sql, [](void* data, int argc, char **argv, char **azColName)->int {
        assert( argc == 1 );
        *static_cast<double*>(data) = atof(argv[0]);
        return 0;
    }, &time, &zErrMsg);

    if (rc != SQLITE_OK) {
        printf("\nrange SELECT failed '%s' | %d: '%s'\n",
            zErrMsg, rc, sqlite3_errmsg(db));
    }

    return time;
}


void store_stats(
        const struct Config& config,
        float K_log,
        double time_stats,
        /* Over all M values */
        vector<float>& prob_gap_norm,
        vector<float>& prob_gap_low,
        vector<float>& prob_gap_high,
        /* Per m value */
        vector<uint64_t>& M_vals,
        vector<float>& expected_prev,
        vector<float>& expected_next,
        vector<float>& probs_seen,
        vector<float>& probs_record,
        vector<float>& probs_missing,
        vector<float>& probs_highmerit) {

    assert( M_vals.size() == expected_prev.size() );
    assert( M_vals.size() == expected_next.size() );
    assert( M_vals.size() == probs_seen.size() );
    assert( M_vals.size() == probs_record.size() );
    assert( M_vals.size() == probs_missing.size() );
    assert( M_vals.size() == probs_highmerit.size() );

    assert( !is_range_already_processed(config) );

    DB db_helper(config.search_db.c_str());
    sqlite3 *db = db_helper.get_db();

    char *zErrMsg = 0;
    if (sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &zErrMsg) != SQLITE_OK) {
        printf("BEGIN TRANSACTION failed: %s\n", zErrMsg);
        exit(1);
    }

    const uint64_t rid = db_helper.config_hash(config);
    const size_t num_rows = M_vals.size();
    char sSQL[500];
    sprintf(sSQL,
        "INSERT INTO range(rid, P, D, m_start, m_inc,"
                          "sieve_length, max_prime,"
                          "min_merit,"
                          "num_m, num_remaining,"
                          "time_stats)"
         "VALUES(%ld,  %d,%d, %ld,%ld,"
                "%d,%ld,  %.3f,"
                "%ld,%ld,  %.2f)"
        "ON CONFLICT(rid) DO UPDATE SET time_stats=%.2f",
            rid,  config.p, config.d, config.mstart, config.minc,
            config.sieve_length, config.max_prime,
            config.min_merit,
            num_rows, num_rows,
            time_stats, time_stats);

    int rc = sqlite3_exec(db, sSQL, NULL, NULL, &zErrMsg);
    if (rc != SQLITE_OK) {
        printf("\nrange INSERT/UPDATE failed %d: %s\n",
            rc, sqlite3_errmsg(db));
        exit(1);
    }

#define BIND_OR_ERROR(func, stmt, index, value)                             \
    if (func(stmt, index, value) != SQLITE_OK) {                            \
        printf("Failed to bind param %d: %s\n", index, sqlite3_errmsg(db)); \
        break;                                                              \
    }

    /* Create SQL statement to INSERT into range_stats. */
    char insert_range_stats[] = (
        "INSERT INTO range_stats(rid, gap, prob_combined, prob_low_side, prob_high_side)"
        " VALUES(?,?, ?,?,?)"
    );
    sqlite3_stmt *insert_range_stmt;
    /* Prepare SQL statement */
    rc = sqlite3_prepare_v2(db, insert_range_stats, -1, &insert_range_stmt, 0);
    if (rc != SQLITE_OK) {
        printf("Could not prepare statement: '%s'\n", insert_range_stats);
        exit(1);
    }

    assert( prob_gap_norm.size() == prob_gap_low.size() );
    assert( prob_gap_norm.size() == prob_gap_high.size() );
    size_t skipped_gap_stats = 0;
    for (size_t g = 1; g < prob_gap_norm.size(); g ++) {
        if (prob_gap_norm[g] < 1e-10 &&
            prob_gap_low[g]  < 1e-10 &&
            prob_gap_high[g] < 1e-10) {
            // XXX: Consider summing the missing prob at g=0.
            skipped_gap_stats += 1;
            continue;
        }

        BIND_OR_ERROR(sqlite3_bind_int64, insert_range_stmt, 1, rid);

        BIND_OR_ERROR(sqlite3_bind_int,    insert_range_stmt, 2, g);
        BIND_OR_ERROR(sqlite3_bind_double, insert_range_stmt, 3, prob_gap_norm[g]);
        BIND_OR_ERROR(sqlite3_bind_double, insert_range_stmt, 4, prob_gap_low[g]);
        BIND_OR_ERROR(sqlite3_bind_double, insert_range_stmt, 5, prob_gap_high[g]);

        int rc = sqlite3_step(insert_range_stmt);
        if (rc != SQLITE_DONE) {
            printf("\nrange_stats insert failed (%ld): %d: %s\n", g, rc, sqlite3_errmsg(db));
            break;
        }

        if (sqlite3_reset(insert_range_stmt) != SQLITE_OK) {
            printf("Failed to reset statement\n");
        }

        if (sqlite3_clear_bindings(insert_range_stmt) != SQLITE_OK) {
            printf("Failed to clear bindings\n");
        }
    }
    if (config.verbose >= 0) {
        printf("Saved %ld rows to 'gap_stats' table\n", prob_gap_norm.size() - skipped_gap_stats);
    }

    /* Create SQL statement to INSERT into m_stats. */
    // NOTE: IGNORE so that can rerun with different max-prime/sieve-length
    char insert_m_stats[] = (
            "INSERT OR IGNORE INTO m_stats"
            "(rid, P, D, m, "
            " prob_record, prob_missing, prob_merit,"
            " e_gap_next, e_gap_prev)"
            "VALUES"
            "(?, ?, ?, ?,"
            " ?, ?, ?,"
            " ?, ?)"
            );

    sqlite3_stmt *stmt;
    /* Prepare SQL statement */
    rc = sqlite3_prepare_v2(db, insert_m_stats, -1, &stmt, 0);
    if (rc != SQLITE_OK) {
        printf("Could not prepare statement: '%s'\n", insert_m_stats);
        exit(1);
    }

    if (config.verbose >= 2) {
        printf("\n");
    }

    for (size_t i = 0; i < num_rows; i++) {
        uint64_t m = M_vals[i];

        float e_next = expected_next[i];
        float e_prev = expected_prev[i];

        size_t r = i + 1;
        if (config.verbose >= 2 && (
                    (r <= 5) || (r <= 500 && r % 100 == 0) ||
                    (r <= 10000 && r % 1000 == 0) ||
                    (r <= 10000 && r % 10000 == 0) ||
                    (r % 100000 == 0) ||
                    (r == num_rows))) {
            printf("Saving Row: %6ld/%ld %6ld: %.1f, %.1f | R: %.1e M: %.1e HM(%.1f): %.1e\n",
                    r, num_rows, m,
                    e_next, e_prev,
                    probs_record[i], probs_missing[i],
                    config.min_merit, probs_highmerit[i]);
        }

        BIND_OR_ERROR(sqlite3_bind_int64, stmt, 1, rid);

        // P, D, m
        BIND_OR_ERROR(sqlite3_bind_int, stmt, 2, config.p);
        BIND_OR_ERROR(sqlite3_bind_int, stmt, 3, config.d);
        BIND_OR_ERROR(sqlite3_bind_int, stmt, 4, m);

        // prob_record, prob_missing, prob_merit
        BIND_OR_ERROR(sqlite3_bind_double, stmt, 5, probs_record[i]);
        BIND_OR_ERROR(sqlite3_bind_double, stmt, 6, probs_missing[i]);
        BIND_OR_ERROR(sqlite3_bind_double, stmt, 7, probs_highmerit[i]);

        // e_next, e_prev
        BIND_OR_ERROR(sqlite3_bind_double, stmt, 8, e_next);
        BIND_OR_ERROR(sqlite3_bind_double, stmt, 9, e_prev);

        int rc = sqlite3_step(stmt);
        if (rc != SQLITE_DONE) {
            printf("\nm_stats insert failed %ld: (%d, %d, %ld): %d: %s\n",
                i, config.p, config.d, m, rc, sqlite3_errmsg(db));
            break;
        }

        if (sqlite3_reset(stmt) != SQLITE_OK) {
            printf("Failed to reset statement\n");
        }

        if (sqlite3_clear_bindings(stmt) != SQLITE_OK) {
            printf("Failed to clear bindings\n");
        }
    }

    if (sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &zErrMsg) != SQLITE_OK) {
        printf("END TRANSACTION failed: %s\n", zErrMsg);
        exit(1);
    }

    if (config.verbose >= 0) {
        printf("Saved %ld rows to 'm_stats' table\n",
                num_rows);
    }
}


/**
 * Precalculate and cache two calculations
 *
 * Change that ith unknown number (in sieve) is the first prime
 * Prob_prime_nth[i] = (1 - prob_prime)^(i-1) * prob_prime
 *
 * Change that the first prime is later than i'th unknown
 * prod_great_nth[i] = (1 - prob_prime)^i
 */
cdouble DOUBLE_NTH_PRIME_CUTOFF = 1e-13;

void prob_nth_prime(
        double prob_prime,
        vector<float>& prob_prime_nth,
        vector<float>& prob_great_nth) {
    double prob_still_prime = 1.0;
    for (; prob_still_prime > DOUBLE_NTH_PRIME_CUTOFF;) {
        prob_prime_nth.push_back(prob_still_prime * prob_prime);
        prob_great_nth.push_back(prob_still_prime);
        prob_still_prime *= 1 - prob_prime;
    }
}


void prob_combined_gap(
        double prob_prime,
        vector<float>& prob_combined) {
    double prob = prob_prime * prob_prime;
    // Want error < 1e-9 | unknown_i * unkown_j * 1e-15 ~= 2000 * 2000 * 2.5e-6 = 1e-9
    for (; prob > 2.5e-16;) {
        prob_combined.push_back(prob);
        prob *= 1 - prob_prime;
    }
}


void prob_extended_gap(
        const struct Config& config,
         double PROB_PRIME,
        const vector<uint32_t>& poss_record_gaps,
        ProbNth &gap_probs) {

    const unsigned int SL = config.sieve_length;

    // ----- Generate primes for P
    vector<uint32_t> K_primes = get_sieve_primes(config.p);
    assert( K_primes.back() == config.p);

    // Need to correct for gcd_ui(K, i) below
    double prob_prime_coprime = PROB_PRIME;
    for (auto prime : K_primes) {
        if (config.d % prime != 0) {
            prob_prime_coprime /= (1 - 1.0 / prime);
        }
    }

    const size_t EXT_SIZE = 2 * SL;
    // Sieve a lot more (it's fast)
    vector<char> is_coprime(EXT_SIZE, true);
    {
        for (auto prime : K_primes) {
            if (config.d % prime == 0)
                continue;

            for (size_t i = 0; i < EXT_SIZE; i += prime)
                is_coprime[i] = false;
        }
    }

    const vector<uint32_t> wheel_primes = {2, 3, 5, 7};
	map<uint32_t, uint32_t> k_mod_p;

    int wheel = 1;
    for (auto p : wheel_primes) {
        if (config.d % p == 0) {
            wheel *= p;
            prob_prime_coprime /= (1 - 1.0 / p);

			uint32_t k_mod = 1;
			for (auto k : K_primes) {
            	if (config.d % k != 0)
					k_mod = (k_mod * k) % p;
			}
			k_mod_p[p] = k_mod;
			//printf("\tK %% %d = %d\n", p, k_mod);
        }
    }
    gap_probs.wheel_d = wheel;

    // Same as prob_prime_nth_sieve but not considering sieving (because outside of interval).
    vector<float> prob_prime_nth;
    vector<float> prob_great_nth;
    prob_nth_prime(prob_prime_coprime, prob_prime_nth, prob_great_nth);

    // For each wheel mark off a divisors of small primes in d.
	map<int, vector<char>> coprime_ms;
    float average_inner_coprime = 0;
    float average_extended_coprime = 0;

    for (int m = 0; m < wheel; m++) {
        // Don't need any where m is coprime with d
        if (gcd(m, wheel) > 1) continue;

        // Copy is_coprime
        vector<char> is_coprime_m(is_coprime);

        // Mark off multiples of d primes
        for (auto p : wheel_primes) {
            if (config.d % p != 0) continue;

			// (m * K) % p;
			uint32_t first = (m * k_mod_p[p]) % p;

            // first multiple on the positive side: -m % p
            for (size_t i = p - first; i < EXT_SIZE; i += p) {
                is_coprime_m[i] = false;
            }
        }
		coprime_ms[m] = is_coprime_m;

        size_t inner_coprime    = std::count(is_coprime_m.begin(), is_coprime_m.begin() + SL, true);
        size_t extended_coprime = std::count(is_coprime_m.begin() + SL, is_coprime_m.end(), true);
        average_inner_coprime += inner_coprime;
        average_extended_coprime += extended_coprime;

        //printf("\tWheel: %-3d %ld/%d inner, %ld/%d extended coprime)\n",
		//	m, inner_coprime, SL, extended_coprime, SL);
	}
    average_inner_coprime /= coprime_ms.size();
    average_extended_coprime /= coprime_ms.size();
    gap_probs.average_coprime = average_extended_coprime;
    gap_probs.prob_extended = prob_great_nth[(int) average_extended_coprime];

    if (config.verbose >= 2) {
        printf("Using Wheel: %d for extended probs\n", wheel);
        printf("\tAverage %5.0f inner    coprimes => %.3g%% prob_greater\n",
            average_inner_coprime, 100 * prob_great_nth[(int) average_inner_coprime]);
        printf("\tAverage %5.0f extended coprimes => %.3g%% prob_greater\n",
            gap_probs.average_coprime, 100 * gap_probs.prob_extended);
    }

    for (int m = 0; m < wheel; m++) {
        if (gcd(m, wheel) > 1) continue;

		vector<char> &is_coprime_m = coprime_ms.at(m);
		vector<char> &is_coprime_m_prev = coprime_ms.at(wheel - m);

        vector<uint32_t> count_coprime_m(EXT_SIZE, 0);
        {
            size_t count = 0;
            // partial_sum starting at SL+1
            for (size_t i = SL+1; i < EXT_SIZE; i++) {
                count += is_coprime_m[i];
                count_coprime_m[i] = count;
				//if (count < 10 && is_coprime[i])
				//	cout << "\t" << i << endl;
            }
        }

        vector<float> prob_extended_record(SL+1, 0.0);
        for (size_t gap_prev = 1; gap_prev <= SL; gap_prev++) {
            // only needed for values that can be coprime with K
            if (!is_coprime_m_prev[gap_prev]) {
                prob_extended_record[gap_prev] = std::nan("");
                continue;
            }

            if (gap_prev + EXT_SIZE < poss_record_gaps.front()) {
                continue;
            }

            double prob_record = 0;
            for (uint32_t record_gap : poss_record_gaps ) {
                uint32_t dist = record_gap - gap_prev;
                if (dist <= SL) continue;

                if (dist >= is_coprime_m.size()) break;

                // dist can never be prime.
                if (!is_coprime_m[dist]) continue;

                // This is the nth possible prime after SL
                uint32_t num_coprime = count_coprime_m[dist];
                if (num_coprime >= prob_prime_nth.size()) break;

                // chance of dist_after being first prime.
                prob_record += prob_prime_nth[num_coprime];
            }

            // Prob record gap, with 1 <= gap_prev <= SL, SL <= gap_next
			assert(prob_record >= 0 && prob_record < 1);
            prob_extended_record[gap_prev] = prob_record;
        }
        gap_probs.extended_record_high[m] = prob_extended_record;
    }
}


void setup_probnth(
        const struct Config &config,
        double N_log,
        const vector<uint32_t> &poss_record_gaps,
        ProbNth &gap_probs) {

    // ----- Sieve stats
    cdouble PROB_PRIME = 1 / N_log - 1 / (N_log * N_log);
    cdouble UNKNOWNS_AFTER_SIEVE = 1 / (log(config.max_prime) * exp(GAMMA));
    cdouble UNKNOWNS_AFTER_COPRIME = 1 / (log(config.p) * exp(GAMMA));
    cdouble PROB_PRIME_AFTER_SIEVE = PROB_PRIME / UNKNOWNS_AFTER_SIEVE;
    if (config.verbose >= 2) {
        printf("prob prime             : %.7f\n", PROB_PRIME);
        printf("prob prime coprime     : %.7f\n", PROB_PRIME / UNKNOWNS_AFTER_COPRIME);
        printf("prob prime after sieve : %.5f\n\n", PROB_PRIME_AFTER_SIEVE);
    }

    prob_nth_prime(
        PROB_PRIME_AFTER_SIEVE,
        gap_probs.prime_nth_sieve, gap_probs.great_nth_sieve);

    prob_combined_gap(
        PROB_PRIME_AFTER_SIEVE,
        gap_probs.combined_sieve);

    // Prob record with gap[i] and other gap > SL
    {
        auto s_start_t = high_resolution_clock::now();

        prob_extended_gap(
            config,
            PROB_PRIME,
            poss_record_gaps,
            gap_probs
        );

        if (config.verbose >= 1) {
            auto s_stop_t = high_resolution_clock::now();
            double   secs = duration<double>(s_stop_t - s_start_t).count();
            printf("Extended prob records considered (%.2f seconds)\n\n", secs);
        }
    }
}

/** Parse line (potentially with rle) to two positive lists */
void read_unknown_line(
        const struct Config& config,
        uint64_t mi,
        std::ifstream& unknown_file,
        vector<uint32_t>& unknown_low,
        vector<uint32_t>& unknown_high) {

    int unknown_l = 0;
    int unknown_u = 0;

    // Read a line from the file
    {
        int mtest;
        unknown_file >> mtest;
        assert( mtest >= 0 );
        assert( (size_t) mtest == mi );

        std::string delim;
        char delim_char;
        unknown_file >> delim;
        assert( delim == ":" );

        unknown_file >> unknown_l;
        unknown_l *= -1;
        unknown_file >> unknown_u;

        unknown_file >> delim;
        assert( delim == "|" );
        delim_char = unknown_file.get(); // get space character
        assert( delim_char == ' ');

        unsigned char a, b;
        int c = 0;
        for (int k = 0; k < unknown_l; k++) {
            if (config.rle) {
                // Read bits in pairs (see save_unknowns_method2)
                a = unknown_file.get();
                b = unknown_file.get();
                c += (a - 48) * 75 + (b - 48);
            } else {
                unknown_file >> c;
                c *= -1;
            }
            unknown_low.push_back(c);
        }

        unknown_file >> delim;
        assert( delim == "|" );
        delim_char = unknown_file.get(); // get space character
        assert( delim_char == ' ');

        c = 0;
        for (int k = 0; k < unknown_u; k++) {
            if (config.rle) {
                a = unknown_file.get();
                b = unknown_file.get();
                c += (a - 48) * 75 + (b - 48);
            } else {
                unknown_file >> c;
            }
            unknown_high.push_back(c);
        }
    }
}

void run_gap_file(
        /* input */
        const struct Config& config,
        const float K_log,
        const vector<float>& records,
        const uint32_t min_record_gap,
        const uint32_t min_gap_min_merit,
        const ProbNth &gap_probs,
        vector<uint32_t>& valid_m,
        std::ifstream& unknown_file,
        /* output */
        vector<float>& prob_gap_norm,
        vector<float>& prob_gap_low,
        vector<float>& prob_gap_high,
        vector<uint64_t>& M_vals,
        vector<float>& expected_prev,
        vector<float>& expected_next,
        vector<float>& probs_seen,
        vector<float>& probs_record,
        vector<float>& probs_missing,
        vector<float>& probs_highmerit) {

    auto  s_start_t = high_resolution_clock::now();

    prob_gap_norm.clear();
    prob_gap_low.clear();
    prob_gap_high.clear();

    // NOTE: prob_gap_low only use values <=  SL but helps with store_stats
    prob_gap_norm.resize(2*config.sieve_length+1, 0);
    prob_gap_low .resize(2*config.sieve_length+1, 0);
    prob_gap_high.resize(2*config.sieve_length+1, 0);

    // sum prob_record_inside sieve
    // sum prob_record_outer (extended)
    float sum_prob_inner = 0.0;
    float sum_prob_outer = 0.0;

    // max_prob_record and max_prob_missing_record
    float max_p_record = 1e-10;
    float max_m_record = 1e-10;

    for (uint32_t mi : valid_m) {
        uint64_t m = config.mstart + mi;

        vector<uint32_t> unknown_low, unknown_high;
        read_unknown_line(config, mi, unknown_file, unknown_low, unknown_high);

        if (config.verbose >= 2 && probs_seen.empty()) {
            for (size_t i = 0; i < unknown_low.size(); i += 1) {
                if (2 * unknown_low[i] > MISSING_GAPS_LOW) {
                    printf("MISSING_GAP prob pair: ~%.2e (%ld in middle) "
                            "if both prime: ~%.2e (~%.0f tests per record)\n",
                            gap_probs.combined_sieve[2*i],
                            2 * i - 2,
                            gap_probs.great_nth_sieve[2*i-2],
                            1 / gap_probs.great_nth_sieve[2*i-2]);
                    break;
                }
            }
        }

        // Note slightly different from N_log
        float log_start_prime = K_log + log(m);

        cdouble prob_prev_end = gap_probs.great_nth_sieve[unknown_low.size()];
        cdouble prob_next_end = gap_probs.great_nth_sieve[unknown_high.size()];
        cdouble prob_extended = gap_probs.prob_extended;

        // Directly examined (1 - prob_prev_end) * (1 - prob_next_end)
        // +
        // Extended examined (1 - prob_prev_end) * (prob_next_end * (1 - prob_extended))
        // =                 (1 - prob_prev_end) * (1 - prob_next_end * prob_extended)
        cdouble prob_seen =
            ((1 - prob_prev_end) * (1 - prob_next_end) +
             (1 - prob_prev_end) * (prob_next_end * (1 - prob_extended)) +
             (1 - prob_next_end) * (prob_prev_end * (1 - prob_extended)));

        double prob_record = 0;
        double prob_is_missing_gap = 0;
        double prob_highmerit = 0;

        size_t min_j = unknown_high.size() - 1;
        uint32_t min_interesting_gap = std::min(min_gap_min_merit, min_record_gap);
        for (size_t i = 0; i < unknown_low.size(); i++) {
            uint32_t gap_low = unknown_low[i];
            while ((min_j > 0) && (gap_low + unknown_high[min_j] > min_interesting_gap)) {
                min_j -= 1;
            }

            size_t max_j = std::min(unknown_high.size(), gap_probs.combined_sieve.size() - i);

            for (size_t j = 0; j < max_j; j++) {
                uint32_t gap_high = unknown_high[j];
                uint32_t gap = gap_low + gap_high;

                // Same as prob_prime_nth[i] * prob_prime_nth[j];
                float prob_this_gap = gap_probs.combined_sieve[i + j];

                // XXX: Costs some performance to calculate all of these
                prob_gap_norm[gap] += prob_this_gap;

                if (gap >= min_gap_min_merit) {
                    prob_highmerit += prob_this_gap;
                }

                if (gap >= min_record_gap && records[gap] > log_start_prime) {
                    prob_record += prob_this_gap;

                    if (MISSING_GAPS_LOW <= gap && gap <= MISSING_GAPS_HIGH &&
                            records[gap] == GAP_INF) {
                        prob_is_missing_gap += prob_this_gap;
                    }
                }
            }
        }

        // expected_gap_low | expected_gap_high
        // prob_gap_low     | prob_gap_high

        double e_prev = 0, e_next = 0;
        double prob_record_outer = 0;

        const float PROB_HIGH_GREATER = gap_probs.great_nth_sieve[unknown_high.size()];
        const float PROB_LOW_GREATER  = gap_probs.great_nth_sieve[unknown_low.size()];

        // See `prob_extended_gap`
        int m_high = m % gap_probs.wheel_d;
        const vector<float> &extended_record_high =
            gap_probs.extended_record_high.at(m_high);
        // want -m % wheel_d => wheel_d - m
        const vector<float> &extended_record_low =
            gap_probs.extended_record_high.at(gap_probs.wheel_d - m_high);

		//printf("m: %ld | wheeled: %d %d |\n", m, m_high, gap_probs.wheel_d - m_high);

        for (size_t i = 0; i < std::max(unknown_low.size(), unknown_high.size()); i++) {
            float prob_i = gap_probs.prime_nth_sieve[i];

            // unknown[i'th] is prime, on the otherside have prime be outside of sieve.
            if (i < unknown_low.size()) {
                float conditional_prob = extended_record_high[unknown_low[i]];
                assert(conditional_prob >= 0);
				assert(conditional_prob < 1);

                prob_record_outer += prob_i * PROB_HIGH_GREATER * conditional_prob;
                e_prev += unknown_low[i] * prob_i;

                prob_gap_low[unknown_low[i]] += prob_i;
            }
            if (i < unknown_high.size()) {
                float conditional_prob = extended_record_low[unknown_high[i]];
                assert(conditional_prob >= 0);
				assert(conditional_prob < 1);

                prob_record_outer += prob_i * PROB_LOW_GREATER * conditional_prob;
                e_next += unknown_high[i] * prob_i;

                prob_gap_high[unknown_high[i]] += prob_i;
            }
        }

        // Combination of observed (0 < i, j <= SL) + extended (i or j > SL)
        double prob_record_combined = prob_record + prob_record_outer;

        sum_prob_inner += prob_record;
        sum_prob_outer += prob_record_outer;

        M_vals.push_back(m);
        expected_prev.push_back(e_prev);
        expected_next.push_back(e_next);
        probs_seen.push_back(prob_seen);
        probs_record.push_back(prob_record_combined);
        probs_missing.push_back(prob_is_missing_gap);
        probs_highmerit.push_back(prob_highmerit);

        if (prob_record_combined > max_p_record) {
            max_p_record = prob_record_combined;
            if (config.verbose >= 1) {
                printf("RECORD :%-6ld line %-6ld  unknowns: %3ld, %3ld "
                        "| e: %.1f, %.1f\t| "
                        "prob record: %.2e (%.2e + %.2e)\t| %.7f\n",
                        m, M_vals.size(),
                        unknown_low.size(), unknown_high.size(),
                        e_prev, e_next,
                        prob_record_combined, prob_record, prob_record_outer,
                        prob_seen);
            }
        }

        if (prob_is_missing_gap > max_m_record) {
            max_m_record = prob_is_missing_gap;
            if (config.verbose >= 2) {
                printf("MISSING:%-6ld line %-6ld  unknowns: %3ld, %3ld "
                        "|\t\t\t| prob record: %.2e  missing: %.4e\t| %.7f\n",
                        m, M_vals.size(),
                        unknown_low.size(), unknown_high.size(),
                        prob_record_combined, prob_is_missing_gap, prob_seen);
            }
        }
    }

    // Normalize the probability of gap (across all m) to per m
    for (size_t i = 0; i < prob_gap_norm.size(); i++) {
        prob_gap_norm[i] /= valid_m.size();
        prob_gap_low[i]  /= valid_m.size();
        prob_gap_high[i] /= valid_m.size();
    }

    if (config.verbose >= 0) {
        long  s_tests = probs_seen.size();
        auto s_stop_t = high_resolution_clock::now();
        double   secs = duration<double>(s_stop_t - s_start_t).count();
        printf("%ld m's processed in %.2f seconds (%.2f/sec)\n",
                s_tests, secs, s_tests / secs);

        if (config.verbose >= 1)
            cout << endl;
    }
    if (config.verbose >= 2) {
        printf("prob record inside sieve: %.5f   prob outside: %.5f\n\n",
                sum_prob_inner, sum_prob_outer);
        printf("\tsum(prob(gap[X])): %.5f\n", average_v(prob_gap_norm) * prob_gap_norm.size());
        printf("\tavg seen prob    : %.7f\n", average_v(probs_seen));
    }
}


void calculate_prp_top_percent(
        struct Config& config,
        uint64_t valid_ms,
        double N_log,
        vector<float> &probs_record) {

    printf("\n");

    // Determine PRP time, time per m
    cdouble prp_time_est = prp_time_estimate_composite(N_log, 2 /* verbose */);
    cdouble prob_prime = 1 / N_log - 1 / (N_log * N_log);
    cdouble estimated_prp_per_m = 1 / (prob_prime * log(config.max_prime) * exp(GAMMA));
    cdouble time_per_side = prp_time_est * estimated_prp_per_m;

    // Try to load combined_time from db, fallback to estimate.
    double combined_time = get_range_time(config);
    bool exact = combined_time > 0;
    if (combined_time <= 0) {
        // Calculate combined_sieve time
        mpz_t K;
        init_K(config, K);
        config.verbose = 0;
        // Inflate slightly to account for gap_stat, starting up...
        combined_time = 1.05 * combined_sieve_method2_time_estimate(
            config, K, valid_ms, 0.0 /* prp_time_est */);
        mpz_clear(K);
    }

    printf("\n");
    printf("%sieve time: %.0f seconds (%.2f hours)\n",
        exact ? "S" : "Estimated s", combined_time, combined_time / 3600);
    printf("Estimated time/m: 2 * (%.1f PRP/m / %.1f PRP/s) = %.2f seconds\n",
        estimated_prp_per_m, 1 / prp_time_est, 2 * time_per_side);
    printf("\n");

    // Sort probs, greater first
    vector<float> sorted = probs_record;
    std::sort(sorted.begin(), sorted.end(), std::greater<>());

    vector<size_t> print_points;
    for (auto percent : {1, 5, 10, 20, 30, 50, 100}) {
        size_t count = sorted.size() * percent / 100;
        if (count == 0)
            continue;
        print_points.push_back(count);
    }

    printf("Sum(prob(record)) at different --prp-top-percent.\n");
    printf("\tUses estimate for combined_sieve timing.\n");
    printf("\tEstimate of optimal printed with *\n");
    printf("\n");

    // Both sides & One sided at 20% (assume 80% of prob also)
    for (size_t side_percent : {100, 20}) {

        double sum_prob = 0.0;
        double time = combined_time;

        bool max_happened = false;
        for (size_t i = 1; i <= sorted.size(); i++) {
            // Print at print_points (1, 5, 10, ... 100% and also at "optimal" percent)
            double sides_tested = 1 + side_percent / 100.0;
            double add_t = sides_tested * time_per_side;
            double add_p = sorted[i-1] * (side_percent == 100 ? 1.0 : 0.8);

            time += add_t;
            sum_prob += add_p;

            double avg = sum_prob / time;
            bool is_below_avg = (add_p / add_t) < avg;
            bool first_below = !max_happened && is_below_avg;
            max_happened |= first_below;

            if (first_below || std::count(print_points.begin(), print_points.end(), i)) {
                cdouble percent = i * 100.0 / sorted.size();

                // testing one side and other side smaller percent
                printf("\t%7ld %c(%5.1f%%) | sum(prob) = %9.5f / (%.0f + %6ld * %3g * %.2f) => %6.3f/%.1f hr = %.6f prob/hour\n",
                    i, " *"[first_below], percent,
                    sum_prob, combined_time, i, sides_tested, time_per_side,
                    sum_prob, time / 3600, avg);
            }
        }

        printf("\n");
    }
}


void prime_gap_stats(struct Config config) {
    const unsigned int SIEVE_LENGTH = config.sieve_length;
    const unsigned int SL = SIEVE_LENGTH;
    assert( SL > 1000 );

    auto  s_start_t = high_resolution_clock::now();

    // ----- Read from unknown file
    std::ifstream unknown_file;
    {
        std::string fn = Args::gen_unknown_fn(config, ".txt");
        if (config.verbose >= 0) {
            printf("\nReading from '%s'\n\n", fn.c_str());
        }
        unknown_file.open(fn, std::ios::in);
        assert( unknown_file.is_open() ); // Can't open save_unknowns file
        assert( unknown_file.good() );    // Can't open save_unknowns file
    }

    // ----- Merit Stuff
    mpz_t K;

    int K_digits;
    double K_log;
    K_stats(config, K, &K_digits, &K_log);
    double N_log = K_log + log(config.mstart);
    mpz_clear(K);

    uint32_t min_gap_min_merit = std::ceil(config.min_merit * N_log);
    if (config.verbose >= 2) {
        printf("Min Gap ~= %d (for merit > %.1f)\n\n", min_gap_min_merit, config.min_merit);
    }

    // ----- Get Record Prime Gaps
    vector<float> records = get_record_gaps(config);

    // gap that would be a record with m*P#/d
    vector<uint32_t> poss_record_gaps;
    load_possible_records(N_log, records, poss_record_gaps);
    assert( poss_record_gaps.size() >= 2);
    {
        if (config.verbose >= 1) {
            printf("Found %ld possible record gaps (%d to %d)\n",
                poss_record_gaps.size(), poss_record_gaps.front(), poss_record_gaps.back());

            #if SEARCH_MISSING_GAPS
            {
                size_t missing_gaps = 0;
                for (size_t gap = MISSING_GAPS_LOW; gap <= MISSING_GAPS_HIGH; gap += 2) {
                    missing_gaps += (records[gap] == GAP_INF);
                }
                printf("Found %ld missing gaps between %d and %d\n",
                    missing_gaps, MISSING_GAPS_LOW, MISSING_GAPS_HIGH);
            }
            #endif  // SEARCH_MISSING_GAPS
        }
        if (config.verbose >= 2) {
            for (int gi = 0; gi <= 2; gi++) {
                size_t g = poss_record_gaps[gi];
                printf("\tIf found Gap: %ld (current: %.2f) would improve to %.3f\n",
                    g, g / records[g], g / N_log);
            }
        }

        if (poss_record_gaps.front() > 2 * SIEVE_LENGTH) {
            printf("\n\n\n");
            printf("\tHard to determine record prob, 2 * sieve_length < min_record_gap");
            printf("\n\n\n");
        }
        if (config.verbose >= 1) {
            cout << endl;
        }
    }

    ProbNth gap_probs;
    setup_probnth(
        config, N_log,
        poss_record_gaps,
        gap_probs);

    vector<uint32_t> valid_m;
    for (uint64_t mi = 0; mi < config.minc; mi++) {
        if (gcd(config.mstart + mi, config.d) == 1) {
            valid_m.push_back(mi);
        }
    }

    /* Over all m values */
    vector<float> prob_gap_norm;
    vector<float> prob_gap_low;
    vector<float> prob_gap_high;

    /* Per m stats */
    vector<uint64_t> M_vals;
    vector<float> expected_prev;
    vector<float> expected_next;
    /* Per m probabilities */
    vector<float> probs_seen;
    vector<float> probs_record;
    vector<float> probs_missing;
    vector<float> probs_highmerit;

    // ----- Main calculation
    run_gap_file(
        /* Input */
        config, K_log,
        records, poss_record_gaps.front(), min_gap_min_merit,
        gap_probs,
        valid_m,
        /* sieve input */
        unknown_file,
        /* output */
        prob_gap_norm, prob_gap_low, prob_gap_high,
        M_vals,
        expected_prev, expected_next,
        probs_seen,
        probs_record, probs_missing, probs_highmerit
    );

    // Compute sum_missing_prob, sum_record_prob @1,5,10,20,50,100%
    if (config.verbose >= 1) {
        vector<float> expected_gap;
        for (size_t i = 0; i < expected_prev.size(); i++) {
            expected_gap.push_back(expected_prev[i] + expected_next[i]);
        }

        prob_stats("EXPECTED GAP", expected_gap);
        printf("\n");

        prob_stats("RECORD", probs_record);

        double avg_missing = average_v(probs_missing);
        double avg_record  = average_v(probs_record);
        // missing mostly incluses > 3 * SL, which is likely to be a record.
        double uncertainty = avg_missing / (avg_missing + avg_record);

        if (uncertainty > 1e-5) {
            printf("\tRECORD : avg: %.2e | missing: %.2e | uncertainty: %.4f%% \n",
                avg_record, avg_missing, 100 * uncertainty);
        }

        if (config.verbose >= 2) {
            if (avg_missing > 0) {
                printf("\n");
                prob_stats("MISSING", probs_missing);
            }
        }
        printf("\n");
    }

    // Make sure probs_record didn't sort our copy.
    assert(!std::is_sorted(probs_record.begin(), probs_record.end(), std::greater<>()));


    if (config.save_unknowns) {
        auto s_stop_t = high_resolution_clock::now();
        double   secs = duration<double>(s_stop_t - s_start_t).count();

        store_stats(
            config, K_log,
            secs,
            prob_gap_norm, prob_gap_low, prob_gap_high,
            M_vals,
            expected_prev, expected_next,
            probs_seen,
            probs_record, probs_missing, probs_highmerit
        );
    }

    if (config.verbose >= 1) {
        calculate_prp_top_percent(config, valid_m.size(), N_log, probs_record);
    }
}
