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


const char records_db[] = "gaps.db";
const char gaps_db[]    = "prime-gap-search.db";

// Limits the size of record list
const uint32_t MAX_GAP = 1'000'000;
const float    GAP_INF = std::numeric_limits<float>::max();   // log(starting_prime)

// Generated from https://primegap-list-project.github.io/lists/missing-gaps/
// Range of missing gap to search, values are loaded from records_db.
const uint32_t MISSING_GAPS_LOW  = 113326;
const uint32_t MISSING_GAPS_HIGH = 132928;

// TODO make this an option (enabled or percent).
const float MISSING_GAP_SAVE_PERCENT = 0.05;


void prime_gap_stats(const struct Config config);
bool is_range_already_processed(const struct Config& config);

static double average_v(vector<float> &a) {
    return std::accumulate(a.begin(), a.end(), 0.0) / a.size();
}

static void prob_stats(char const *name, vector<float> probs) {
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

    Config config = argparse(argc, argv);

    if (config.verbose >= 3) {
        printf("\tCompiled with GMP %d.%d.%d\n\n",
            __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);
    }

    if (config.valid == 0) {
        show_usage(argv[0]);
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


sqlite3* get_db(const char* path) {
    {
        std::ifstream f(path);
        if (!f.good()) {
            printf("database(%s) doesn't exist\n", path);
            exit(1);
        }
    }

    sqlite3 *db;
    if (sqlite3_open(path, &db) != SQLITE_OK) {
        printf("Can't open database(%s): %s\n", path, sqlite3_errmsg(db));
        exit(1);
    }
    return db;
}


vector<float> get_record_gaps() {
    vector<float> records(MAX_GAP, GAP_INF);

    // TODO accept db file as param.
    sqlite3 *db = get_db(records_db);

    /* Create SQL statement */
    char sql[] = "SELECT gapsize, merit FROM gaps";
    char *zErrMsg = 0;

    /* Execute SQL statement */
    int rc = sqlite3_exec(db, sql, [](void* recs, int argc, char **argv, char **azColName)->int {
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
    sqlite3_close(db);

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


uint64_t config_hash(const struct Config& config) {
    // Hash the config to a uint64
    uint64_t hash =    config.mstart;
    hash = hash * 31 + config.minc;
    hash = hash * 31 + config.p;
    hash = hash * 31 + config.d;
    hash = hash * 31 + config.sieve_length;
    hash = hash * 31 + config.max_prime;
    return hash;
}


bool is_range_already_processed(const struct Config& config) {
    uint64_t hash = config_hash(config);
    char sql[200];
    sprintf(sql, "SELECT count(*) FROM range WHERE rid = %ld", hash);
    char *zErrMsg = 0;

    sqlite3 *db = get_db(gaps_db);

    int count = 0;
    int rc = sqlite3_exec(db, sql, [](void* data, int argc, char **argv, char **azColName)->int {
        assert( argc == 1 );
        *static_cast<int*>(data) = atoi(argv[0]);
        return 0;
    }, &count, &zErrMsg);
    sqlite3_close(db);

    if (rc != SQLITE_OK) {
        printf("\nrange SELECT failed '%s' | %d: '%s'\n",
            zErrMsg, rc, sqlite3_errmsg(db));
        exit(1);
    }
    return count > 0;
}


void store_stats(
        const struct Config& config,
        float K_log,
        double time_stats,
        vector<float>& prob_gap_norm,
        vector<uint64_t>& M_vals,
        vector<float>& expected_prev,
        vector<float>& expected_next,
        vector<float>& probs_seen,
        vector<float>& probs_record) {

    assert( M_vals.size() == expected_prev.size() );
    assert( M_vals.size() == expected_next.size() );
    assert( M_vals.size() == probs_seen.size() );
    assert( M_vals.size() == probs_record.size() );

    sqlite3 *db = get_db(gaps_db);

    assert( !is_range_already_processed(config) );

    char *zErrMsg = 0;
    if (sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &zErrMsg) != SQLITE_OK) {
        printf("BEGIN TRANSACTION failed: %s\n", zErrMsg);
        exit(1);
    }

    const uint64_t rid = config_hash(config);
    const size_t num_rows = M_vals.size();
    char sSQL[300];
    sprintf(sSQL,
        "INSERT INTO range(rid, m_start, m_inc, P, D,"
                          "sieve_length, max_prime,"
                          "min_merit,"
                          "num_m,num_processed, num_remaining,"
                          "time_sieve,time_stats,time_tests)"
         "VALUES(%ld,  %ld,%ld, %d,%d,"
                "%d,%ld,  %.3f,"
                "%ld,%d,%ld,  %.2f,%.2f,%.2f)",
            rid,  config.mstart, config.minc,  config.p, config.d,
            config.sieve_length, config.max_prime,
            config.min_merit,
            num_rows, 0, num_rows,
            0.0, time_stats, 0.0);

    int rc = sqlite3_exec(db, sSQL, NULL, NULL, &zErrMsg);
    if (rc != SQLITE_OK) {
        printf("\nrange INSERT failed %d: %s\n",
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
        "INSERT INTO range_stats(rid, gap, prob) VALUES(?,?,?)"
    );
    sqlite3_stmt *insert_range_stmt;
    /* Prepare SQL statement */
    rc = sqlite3_prepare_v2(db, insert_range_stats, -1, &insert_range_stmt, 0);
    if (rc != SQLITE_OK) {
        printf("Could not prepare statement: '%s'\n", insert_range_stats);
        exit(1);
    }

    size_t skipped_gap_stats = 0;
    for (size_t g = 2; g < prob_gap_norm.size(); g += 2) {
        if (prob_gap_norm[g] < 1e-10) {
            // XXX: Consider summing the misc prob at g=0.
            skipped_gap_stats += 1;
            continue;
        }

        BIND_OR_ERROR(sqlite3_bind_int64, insert_range_stmt, 1, rid);

        BIND_OR_ERROR(sqlite3_bind_int,    insert_range_stmt, 2, g);
        BIND_OR_ERROR(sqlite3_bind_double, insert_range_stmt, 3, prob_gap_norm[g]);

        int rc = sqlite3_step(insert_range_stmt);
        if (rc != SQLITE_DONE) {
            // Allow BUSY with retry.
            printf("\nCould not execute stmt (%ld): %d: %s\n", g, rc, sqlite3_errmsg(db));
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
    char insert_m_stats[] = (
            "INSERT INTO m_stats"
            "(rid, m, P, D, next_p, prev_p, merit,"
            " prob_record, prob_missing, prob_merit,"
            " e_gap_next, e_gap_prev,"
            " prp_next, prp_prev, test_time)"
            "VALUES"
            "(?, ?, ?, ?, 0, 0, 0,"
            " ?, ?, ?,"
            " ?, ?,"
            " 0, 0, 0)"
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
                    (r <= 5) || (r < 500 && r % 100 == 0) ||
                    (r % 2000 == 0) || r == num_rows)) {
            printf("Saving Row: %6ld/%ld %ld: %.1f, %.1f | %.1e\n",
                    r, num_rows, m,
                    e_next, e_prev, probs_record[i]);
        }

        BIND_OR_ERROR(sqlite3_bind_int64, stmt, 1, rid);

        // m, P, D
        BIND_OR_ERROR(sqlite3_bind_int, stmt, 2, m);
        BIND_OR_ERROR(sqlite3_bind_int, stmt, 3, config.p);
        BIND_OR_ERROR(sqlite3_bind_int, stmt, 4, config.d);

        // prob_record, prob_missing, prob_merit
        BIND_OR_ERROR(sqlite3_bind_double, stmt, 5, probs_record[i]);
        BIND_OR_ERROR(sqlite3_bind_double, stmt, 6, probs_record[i]);
        BIND_OR_ERROR(sqlite3_bind_double, stmt, 7, probs_record[i]);

        // e_next, e_prev
        BIND_OR_ERROR(sqlite3_bind_double, stmt, 8, e_next);
        BIND_OR_ERROR(sqlite3_bind_double, stmt, 9, e_prev);

        int rc = sqlite3_step(stmt);
        if (rc != SQLITE_DONE) {
            // Allow BUSY with retry.
            printf("\nCould not execute stmt (%ld): %d: %s\n", i, rc, sqlite3_errmsg(db));
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

    sqlite3_close(db);
}


/**
 * Precalculate and cache two calculations
 *
 * Change that ith unknown number (in sieve) is the first prime
 * Prob_prime_nth[i] = (1 - prob_prime)^(i-1) * prob_prime
 *
 * Change that the first prime is past ith unknown
 * prod_great_nth[i] = (1 - prob_prime)^i
 */
const double DOUBLE_NTH_PRIME_CUTOFF = 1e-13;

void prob_nth_prime(
        const double prob_prime,
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
        const double prob_prime,
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
        const mpz_t &K,
        const double PROB_PRIME,
        const double PROB_GREATER_CUTOFF,
        const double CUTOFF_CONTRIB,
        vector<uint32_t>& poss_record_gaps,
        vector<float>& prob_record_extended_gap) {

    const unsigned int SL = config.sieve_length;

    // ----- Generate primes for P
    vector<uint32_t> K_primes = get_sieve_primes(config.p);
    assert( K_primes.back() == config.p);

    // Need to correct for gcd_ui(K, i) below
    double prob_prime_coprime = PROB_PRIME;
    for (auto prime : K_primes) {
        if (prime % config.d != 0) {
            prob_prime_coprime /= (1 - 1.0 / prime);
        }
    }

    // Same as prob_prime_nth_sieve but not considering sieving (because outside of interval.
    vector<float> prob_prime_nth;
    vector<float> prob_great_nth;
    prob_nth_prime(prob_prime_coprime, prob_prime_nth, prob_great_nth);

    // Count of numbers [SL, SL+j] < i coprime to D
    // only (K, SL+i) == 1 can be prime
    vector<bool>     is_coprime;
    vector<uint32_t> count_coprime;
    {
        uint32_t count = 0;
        // TODO: verify it doesn't make sense to check larger gaps.
        for (size_t i = SL; i <= poss_record_gaps.back() ; i++) {
            bool coprime = mpz_gcd_ui(NULL, K, i) == 1;
            count += coprime;
            is_coprime.push_back(coprime);
            count_coprime.push_back(count);
            if (prob_great_nth[count] < PROB_GREATER_CUTOFF) {
                break;
            }
        }

        if (config.verbose >= 2) {
            printf("Considering SL + %ld (%d coprime) after record extended gap\n",
                    count_coprime.size(), count);
            printf("\tUsing PROB_GREATER_CUTOFF %.2e * %.2e\n",
                    prob_great_nth[count], CUTOFF_CONTRIB);
        }
        assert( count <= prob_prime_nth.size() );
    }

    for (size_t gap_one = 1; gap_one <= SL; gap_one++) {
        // only needed for values that can be coprime with K
        if (mpz_gcd_ui(NULL, K, gap_one) > 1) {
            prob_record_extended_gap[gap_one] = std::nan("");
            continue;
        }

        double prob_record = 0;
        for (uint32_t record_gap : poss_record_gaps ) {
            uint32_t dist = record_gap - gap_one;
            if (dist <= SL) continue;

            uint32_t dist_after = dist - SL;
            if (dist_after >= is_coprime.size()) break;

            // dist can never be prime.
            if (!is_coprime[dist_after]) continue;

            // This is the nth possible prime after SL
            uint32_t num_coprime = count_coprime[dist_after];

            // chance of dist_after being first prime.
            prob_record += prob_prime_nth[num_coprime];
        }

        // Prob record gap, with 1 <= gap_one <= SL, SL <= X
        prob_record_extended_gap[gap_one] = prob_record;
    }
}


void read_unknown_line(
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
        unknown_file >> delim;
        assert( delim == ":" );

        unknown_file >> unknown_l;
        unknown_l *= -1;
        unknown_file >> unknown_u;

        unknown_file >> delim;
        assert( delim == "|" );

        int c;
        for (int k = 0; k < unknown_l; k++) {
            unknown_file >> c;
            unknown_low.push_back(-c);
        }

        unknown_file >> delim;
        assert( delim == "|" );

        for (int k = 0; k < unknown_u; k++) {
            unknown_file >> c;
            unknown_high.push_back(c);
        }
    }
}

void run_gap_file(
        /* input */
        const struct Config config,
        const float K_log,
        const vector<float>& records,
        const uint32_t min_record_gap,
        const vector<float>& prob_prime_nth,
        const vector<float>& prob_great_nth,
        const vector<float>& prob_combined,
        const vector<float>& prob_record_extended_gap,
        std::ifstream& unknown_file,
        /* output */
        vector<float>& prob_gap_norm,
        vector<uint64_t>& M_vals,
        vector<float>& expected_prev,
        vector<float>& expected_next,
        vector<float>& probs_seen,
        vector<float>& probs_record,
        vector<float>& probs_missing,
        vector<float>& probs_highmerit) {

    auto  s_start_t = high_resolution_clock::now();

    vector<uint32_t> valid_m;
    for (uint64_t mi = 0; mi < config.minc; mi++) {
        if (gcd(config.mstart + mi, config.d) == 1) {
            valid_m.push_back(mi);
        }
    }

    prob_gap_norm.clear();
    prob_gap_norm.resize(2*config.sieve_length+1, 0);

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
        read_unknown_line(mi, unknown_file, unknown_low, unknown_high);

        if (config.verbose >= 2 && probs_seen.empty()) {
            for (size_t i = 0; i < unknown_low[i]; i += 1) {
                if (2 * unknown_low[i] > MISSING_GAPS_LOW) {
                    printf("MISSING_GAP prob pair: ~~%.2e "
                            "if both prime: ~~%.2e (~%.1f tests)\n",
                            prob_combined[2*i],
                            prob_great_nth[2*i-2],
                            1 / prob_great_nth[2*i-2]);
                    break;
                }
            }
        }

        // Note slightly different from N_log
        float log_start_prime = K_log + log(m);

        // TODO: should this actually be (1 - prob_great_nth[size] + prob_unknown_extended)^2
        double prob_seen = (1 - prob_great_nth[unknown_high.size()]) *
            (1 - prob_great_nth[unknown_low.size()]);

        double prob_record = 0;
        double prob_is_missing_gap = 0;

        size_t min_j = unknown_high.size() - 1;
        for (size_t i = 0; i < unknown_low.size(); i++) {
            uint32_t gap_low = unknown_low[i];
            while ((min_j > 0) && (gap_low + unknown_high[min_j] > min_record_gap)) {
                min_j -= 1;
            }

            size_t max_j = std::min(unknown_high.size(), prob_combined.size() - i);

            for (size_t j = 0; j < max_j; j++) {
                uint32_t gap_high = unknown_high[j];
                uint32_t gap = gap_low + gap_high;

                // Same as prob_prime_nth[i] * prob_prime_nth[j];
                float prob_this_gap = prob_combined[i + j];

                prob_gap_norm[gap] += prob_this_gap;

                if (j >= min_j && records[gap] > log_start_prime) {
                    // assert(i + j < prob_combined.size());
                    prob_record += prob_this_gap;

                    if (MISSING_GAPS_LOW <= gap && gap <= MISSING_GAPS_HIGH &&
                            records[gap] == GAP_INF) {
                        prob_is_missing_gap += prob_this_gap;
                    }
                }
            }
        }

        double e_prev = 0, e_next = 0;
        double prob_record_outer = 0;

        const float PROB_HIGH_GREATER = prob_great_nth[unknown_high.size()];
        const float PROB_LOW_GREATER  = prob_great_nth[unknown_low.size()];

        for (size_t i = 0; i < std::max(unknown_low.size(), unknown_high.size()); i++) {
            float prob_i = prob_prime_nth[i];

            // unknown[i'th] is prime, on the otherside have prime be outside of sieve.
            if (i < unknown_low.size()) {
                float conditional_prob = prob_record_extended_gap[unknown_low[i]];
                assert(conditional_prob >= 0);

                prob_record_outer += prob_i * PROB_HIGH_GREATER * conditional_prob;
                e_prev += unknown_low[i] * prob_i;
            }
            if (i < unknown_high.size()) {
                float conditional_prob = prob_record_extended_gap[unknown_high[i]];
                assert(conditional_prob >= 0);

                prob_record_outer += prob_i * PROB_LOW_GREATER * conditional_prob;
                e_next += unknown_high[i] * prob_i;
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

        if (prob_record_combined > max_p_record) {
            max_p_record = prob_record_combined;
            if (config.verbose >= 1) {
                printf("RECORD :%-6ld line %-5ld  unknowns: %4ld, %4ld "
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
                printf("MISSING:%-6ld line %-5ld  unknowns: %4ld, %4ld "
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
        printf("prob record inside sieve: %.5f   prob outside: %.5f\n",
                sum_prob_inner, sum_prob_outer);
        printf("\tsum(prob any gap): %.5f\n", average_v(prob_gap_norm) * prob_gap_norm.size());
        printf("\tavg seen prob:   %.7f\n", average_v(probs_seen));
    }
}


void prime_gap_stats(const struct Config config) {
    const unsigned int SIEVE_LENGTH = config.sieve_length;
    const unsigned int SL = SIEVE_LENGTH;
    assert( SL > 1000 );

    auto  s_start_t = high_resolution_clock::now();

    // ----- Read from unknown file
    std::ifstream unknown_file;
    {
        std::string fn = gen_unknown_fn(config, ".txt");
        if (config.verbose >= 0) {
            printf("\nReading from %s'\n\n", fn.c_str());
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
    float N_log = K_log + log(config.mstart);

    if (config.verbose >= 2) {
        size_t min_gap = config.min_merit * N_log;
        printf("Min Gap ~= %ld (for merit > %.1f)\n\n", min_gap, config.min_merit);
    }

    // ----- Get Record Prime Gaps
    vector<float> records = get_record_gaps();

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

    // ----- Sieve stats
    const double PROB_PRIME = 1 / N_log - 1 / (N_log * N_log);
    const double UNKNOWNS_AFTER_SIEVE = 1 / (log(config.max_prime) * exp(GAMMA));
    const double PROB_PRIME_AFTER_SIEVE = PROB_PRIME / UNKNOWNS_AFTER_SIEVE;
    if (config.verbose >= 2) {
        printf("prob prime             : %.7f\n", PROB_PRIME);
        printf("prob prime after sieve : %.5f\n\n", PROB_PRIME_AFTER_SIEVE);
    }

    // Prob 0th, 1st, 2nd is prime
    vector<float> prob_prime_nth_sieve;
    // Prob < nth is not prime or conversely >= nth is prime.
    vector<float> prob_great_nth_sieve;
    prob_nth_prime(
        PROB_PRIME_AFTER_SIEVE,
        prob_prime_nth_sieve, prob_great_nth_sieve);


    /* prob_combined_sieve[i+j] = prime * (1 - prime)^i * (1 - prime)^j * prime */
    vector<float> prob_combined_sieve;
    prob_combined_gap(
        PROB_PRIME_AFTER_SIEVE,
        prob_combined_sieve);


    // Prob record with gap[i] and other gap > SL
    vector<float> prob_record_extended_gap(SL+1, 0.0);
    {
        size_t EXPECTED_UNKNOWN = round(UNKNOWNS_AFTER_SIEVE * SL);
        assert(EXPECTED_UNKNOWN < prob_prime_nth_sieve.size());

        /**
         * Doesn't need to be too large.
         * Used with one sided sum
         *   prob_record +=
         *      prob_nth[i] * prob_greater[unknown_j.size()] * prob_great_nth[unknown_i[i]]
         */
        /**
         * Aim for totaly accuracy 1e-10
         *      Sum(prob_nth[i]) = 1
         *      prob_greater ~= prob_greater[UNKNOWNS_AFTER_SIEVE * SL]
         *      Sum(prob_record_extended_gap) < 1.0 (based on what precent are records)
         * Set prob_greater_n accordingly
         */
        const float PROB_GREATER_CUTOFF = 1e-10 / prob_great_nth_sieve[EXPECTED_UNKNOWN];

        auto s_start_t = high_resolution_clock::now();

        prob_extended_gap(
            config, K,
            PROB_PRIME, PROB_GREATER_CUTOFF,
            prob_great_nth_sieve[EXPECTED_UNKNOWN],
            poss_record_gaps,
            prob_record_extended_gap
        );

        if (config.verbose >= 1) {
            auto s_stop_t = high_resolution_clock::now();
            double   secs = duration<double>(s_stop_t - s_start_t).count();
            printf("Extended prob records considered (%.2f seconds)\n\n", secs);
        }
    }
    mpz_clear(K);

    /* Over all m values */
    vector<float> prob_gap_norm;

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
        records, poss_record_gaps.front(),
        prob_prime_nth_sieve, prob_great_nth_sieve,
        prob_combined_sieve,
        prob_record_extended_gap,
        /* sieve input */
        unknown_file,
        /* output */
        prob_gap_norm,
        M_vals,
        expected_prev, expected_next,
        probs_seen,
        probs_record, probs_missing, probs_highmerit
    );

    // Compute sum_missing_prob, sum_record_prob @1,5,10,20,50,100%
    if (config.verbose >= 1) {
        prob_stats("RECORD", probs_record);
    }
    if (config.verbose >= 2) {
        if (average_v(probs_missing) > 0) {
            printf("\n");
            prob_stats("MISSING", probs_missing);
        }
    }

    if (config.save_unknowns) {
        auto s_stop_t = high_resolution_clock::now();
        double   secs = duration<double>(s_stop_t - s_start_t).count();

        store_stats(
            config, K_log,
            secs,
            prob_gap_norm,
            M_vals,
            expected_prev, expected_next,
            probs_seen, probs_record
        );
    }
}
