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
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <istream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <vector>
#include <unordered_map>

#include <gmp.h>
#include <sqlite3.h>

#include "gap_common.h"


using std::cout;
using std::endl;
using std::pair;
using std::tuple;
using std::vector;
using std::unordered_map;
using namespace std::chrono;


typedef const double cdouble;


// Limits the size of record list
const uint32_t MAX_GAP = 1'000'000;
const float    GAP_INF = std::numeric_limits<float>::max();   // log(starting_prime)

// Generated from https://primegap-list-project.github.io/lists/missing-gaps/
// Range of missing gap to search, values are loaded from gaps_db.
const uint32_t MISSING_GAPS_MIN  = 122290;


class ProbNth {
    public:
        /**
         * Precalculate and cache two calculations
         * See `prob_nth_prime`
         *
         * Change that the first prime is later than i'th unknown (1 indexed)
         * greater_nth[i] = (1 - prob_prime)^i
         *
         * Change that ith unknown number (in sieve) is the first prime (0 indexed)
         * prime_nth[0] = prob_prime
         * prime_nth[i] = (1 - prob_prime)^(i-1) * prob_prime = greater_nth[i-1] * prob_prime
         */
        vector<float> prime_nth;
        vector<float> greater_nth;

        /**
         * Probability that prev_prime & next_prime have X unknown composites in middle
         * prob_combined_sieve[i+j] = prime_nth[i] * prime_nth[j]
         *                          = prime * (1 - prime)^i * (1 - prime)^j * prime
         *                          = prime ^ 2 * (1 - prime)^(i+j)
         */
        vector<float> combined_nth;

        /**
         * wheel_d is gcd(D, 2*3*5*7)
         * Relevant statistics are lookup with (m % wheel_d)
         * ((m % wheel_d) * K + X) with factors of 2,3,5,7 can be eliminated
         * where these X aren't coprime to K
         */
        int wheel_d = -1;
        /**
         * Probability of gap[i] on prev side, next gap > SL, and is record.
         * Sum(prob_combined_sieve[i-1 + unknowns[side] + j-1,
         *      where unknowns[i] + extended[j] is record)
         *
         * Because this uses m, need to handle prev and next side differently
         */
        map<int, vector<float>> extended_record_next;
        /**
         * Similar to extended_record
         *
         * For(extended_gap_prev)
         *     For(extended_gap_next)
         *        record += prob gap prev * prob gap next
         *
         * Caller needs to adjust for prob both outside SL
         */
        map<int, double> extended_extended_record;

        /**
         * Average probability of gap > 2*SL (assuming gap > SL)
         * = pow(Prob(prime | coprime), len(extended coprimes))
         */
        map<int, double> prob_greater_extended;
};


class ProbM {
    public:
        /**
         * Probability of this m & partial sieved interval generating a record gap
         * Sum of record_inner, record_extended, record_extended2
         */
        float record = 0.0;
        float record_inner = 0.0;
        float record_extended = 0.0;
        float record_extended2 = 0.0;

        /**
         *  Probability of Record that was covered by the calculation
         *  Related to sieve-length
         *  Sum of direct, direct-extended, extended-extended
         */
        float seen = 0.0;

        /** Other probabilities that are tracked */
        float is_missing_gap = 0.0;
        float highmerit = 0.0;

        /** Expected gap size */
        float expected_gap_prev = 0.0;
        float expected_gap_next = 0.0;
};


void prob_record_vs_plimit(struct Config config);
void prime_gap_stats(struct Config config);
bool is_range_already_processed(const struct Config& config);

static double average_v(vector<float> &a) {
    return a.empty() ? 0.0 : std::accumulate(a.begin(), a.end(), 0.0) / a.size();
}

static void prob_stats(char const *name, vector<float> &probs) {
    vector<float> sorted = probs;
    std::sort(sorted.begin(), sorted.end(), std::greater<>());

    printf("\n");
    for (auto percent : {1, 5, 10, 20, 50, 100}) {
        size_t count = probs.size() * percent / 100;
        if (count == 0)
            continue;

        printf("\t%-7s: top %3d%% (%6ld)", name, percent, count);

        double sum_prob = std::accumulate(sorted.begin(), sorted.begin() + count, 0.0);
        if (strncmp(name, "EXPECTED", 8) != 0) {
            printf(" sum(prob) = %.2e", sum_prob);
        }
        printf(" (avg: %.2e)\n", sum_prob / count);

        if (sorted[count-1] == 0)
            break;
    }
}


int main(int argc, char* argv[]) {

    Config config = Args::argparse(argc, argv, Args::Pr::STATS);

    if (config.verbose >= 3) {
        printf("\tCompiled with GMP %d.%d.%d\n\n",
            __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);
    }

    if (config.valid == 0) {
        Args::show_usage(argv[0], Args::Pr::STATS);
        return 1;
    }

    if (!config.save_unknowns && !config.testing) {
        printf("Not saving unknowns (--save-unknowns=0)\n");
    } else if (is_range_already_processed(config)) {
        cout << "Range already processed!" << endl;
        return 1;
    }

    if (config.minc == 1 && config.mstart != 1) {
        prob_record_vs_plimit(config);
        return 0;
    }

    prime_gap_stats(config);
    return 0;
}


//---------------------------------------------------------------------------//


vector<float> get_record_gaps(const struct Config& config) {
    uint32_t sieve_interval = 2 * config.sieve_length + 1;
    const size_t records_size = std::max(MAX_GAP, sieve_interval);

    vector<float> records(records_size, GAP_INF);

    DB db(config.gaps_db.c_str());

    /* Create SQL statement */
    char sql[] = "SELECT gapsize, merit FROM gaps";
    char *zErrMsg = nullptr;

    /* Execute SQL statement */
    int rc = sqlite3_exec(db.get_db(), sql, [](void* recs, int argc, char **argv, char **azColName)->int {
        char *test;
        int64_t gap = strtol(argv[0], &test, 10);
        assert(test != argv[0] && *test == '\0');
        assert(gap > 0 && gap < 30'000'000);

        auto *recs_vec = static_cast<vector<float>*>(recs);
        if ((size_t)gap < recs_vec->size()) {
            // Recover log(startprime)
            (*recs_vec)[gap] = gap / strtof(argv[1], nullptr);
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
    // XXX: Records only have 5 sig figs so this method could possible counts records
    //      from same K but smaller m as possible.

    cdouble max_merit_region = N_log < 1400 ? 32 : 28;
    for (size_t g = 2; g < records.size(); g += 2) {
        // All records > max_merit_region are high likely records, See MAX_RECORD and max_e_c_i
        if (g / N_log > max_merit_region) {
            break;
        }

        if (records[g] > N_log) {
            poss_record_gaps.push_back(g);
        }
    }

    assert(is_sorted(poss_record_gaps.begin(), poss_record_gaps.end()));
    assert(poss_record_gaps.front() <= MISSING_GAPS_MIN);
}


bool is_range_already_processed(const struct Config& config) {
    DB db_helper(config.search_db.c_str());
    sqlite3 *db = db_helper.get_db();

    uint64_t hash = db_helper.config_hash(config);
    char sql[200];
    sprintf(sql, "SELECT count(*) FROM range WHERE rid = %ld and time_stats > 0", hash);
    char *zErrMsg = nullptr;

    int count = 0;
    int rc = sqlite3_exec(db, sql, [](void* data, int argc, char **argv, char **azColName)->int {
        assert( argc == 1 );
        int64_t test = strtol(argv[0], nullptr, 10);
        assert(test >= 0 && test <= 1);
        *static_cast<int*>(data) = test;
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
    char *zErrMsg = nullptr;

    double time = 0;
    int rc = sqlite3_exec(db, sql, [](void* data, int argc, char **argv, char **azColName)->int {
        assert( argc == 1 );
        *static_cast<double*>(data) = strtof(argv[0], nullptr);
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
        double time_stats,
        /* Over all M values */
        vector<float>& prob_gap_norm,
        vector<float>& prob_gap_prev,
        vector<float>& prob_gap_next,
        /* Per m value */
        vector<uint32_t>& valid_mi,
        unordered_map<uint64_t,ProbM> &M_stats) {

    assert( !is_range_already_processed(config) );

    DB db_helper(config.search_db.c_str());
    sqlite3 *db = db_helper.get_db();
    // Wait up to 60s to try and commit these records (range is most important)
    sqlite3_busy_timeout(db, 60'000);

    char *zErrMsg = nullptr;
    if (sqlite3_exec(db, "BEGIN TRANSACTION", nullptr, nullptr, &zErrMsg) != SQLITE_OK) {
        printf("BEGIN TRANSACTION failed: %s\n", zErrMsg);
        exit(1);
    }

    const uint64_t rid = db_helper.config_hash(config);
    const size_t num_rows = M_stats.size();
    char sSQL[500];
    sprintf(sSQL,
        "INSERT INTO range(rid, P,D, m_start,m_inc,"
                          "sieve_length,max_prime,"
                          "min_merit,"
                          "num_m,"
                          "time_stats)"
        "VALUES(%ld,  %d,%d, %ld,%ld,"
                "%d,%ld,  %.3f,"
                "%ld, %.2f)"
        "ON CONFLICT(rid) DO UPDATE SET time_stats=%.2f",
            rid,  config.p, config.d, config.mstart, config.minc,
            config.sieve_length, config.max_prime,
            config.min_merit,
            num_rows,
            time_stats, time_stats);

    {
        int rc = sqlite3_exec(db, sSQL, nullptr, nullptr, &zErrMsg);
        if (rc != SQLITE_OK) {
            printf("\nrange INSERT/UPDATE failed %d: %s\n",
                   rc, sqlite3_errmsg(db));
            exit(1);
        }
    }

#define BIND_OR_ERROR(func, stmt, index, value)                             \
    if (func(stmt, index, value) != SQLITE_OK) {                            \
        printf("Failed to bind param %d: %s\n", index, sqlite3_errmsg(db)); \
        break;                                                              \
    }

    /* Create SQL statement to INSERT into range_stats. */
    char insert_range_stats[] = (
        "INSERT OR IGNORE INTO range_stats(rid, gap, prob_combined, prob_low_side, prob_high_side)"
        " VALUES(?,?, ?,?,?)"
    );
    sqlite3_stmt *insert_range_stmt;
    {
        /* Prepare SQL statement */
        int rc = sqlite3_prepare_v2(db, insert_range_stats, -1, &insert_range_stmt, nullptr);
        if (rc != SQLITE_OK) {
            printf("Could not prepare statement: '%s'\n", insert_range_stats);
            exit(1);
        }
    }

    assert( prob_gap_norm.size() == prob_gap_prev.size() );
    assert( prob_gap_norm.size() == prob_gap_next.size() );
    size_t skipped_gap_stats = 0;
    for (int g = (int) prob_gap_norm.size() - 1; g >= 0; g--) {
        if (prob_gap_norm[g] < 1e-10 &&
            prob_gap_prev[g]  < 1e-10 &&
            prob_gap_next[g] < 1e-10) {

            // All skipped gaps are summed at prob[0]
            prob_gap_norm[0] += prob_gap_norm[g];
            prob_gap_prev[0]  += prob_gap_prev[g];
            prob_gap_next[0] += prob_gap_next[g];
            skipped_gap_stats += 1;
            continue;
        }

        BIND_OR_ERROR(sqlite3_bind_int64, insert_range_stmt, 1, rid);

        BIND_OR_ERROR(sqlite3_bind_int,    insert_range_stmt, 2, g);
        BIND_OR_ERROR(sqlite3_bind_double, insert_range_stmt, 3, prob_gap_norm[g]);
        BIND_OR_ERROR(sqlite3_bind_double, insert_range_stmt, 4, prob_gap_prev[g]);
        BIND_OR_ERROR(sqlite3_bind_double, insert_range_stmt, 5, prob_gap_next[g]);

        int rc = sqlite3_step(insert_range_stmt);
        if (rc != SQLITE_DONE) {
            printf("\nrange_stats insert failed (%d): %d: %s\n", g, rc, sqlite3_errmsg(db));
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
        printf("Saved %ld rows to 'range_stats' table\n", prob_gap_norm.size() - skipped_gap_stats);
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
    {
        /* Prepare SQL statement */
        int rc = sqlite3_prepare_v2(db, insert_m_stats, -1, &stmt, nullptr);
        if (rc != SQLITE_OK) {
            printf("Could not prepare statement: '%s'\n", insert_m_stats);
            exit(1);
        }
    }

    if (config.verbose >= 2) {
        printf("\n");
    }

    size_t row_i = 0;
    for (uint64_t mi : valid_mi) {
        uint64_t m = config.mstart + mi;
        ProbM &stats = M_stats[mi];

        float e_next = stats.expected_gap_next;
        float e_prev = stats.expected_gap_prev;

        size_t r = ++row_i;
        if (config.verbose >= 2 && (
                    (r <= 2) ||
                    (r <= 200 && r % 100 == 0) ||
                    (r <= 2000 && r % 1000 == 0) ||
                    (r <= 20000 && r % 10000 == 0) ||
                    (r <= 200000 && r % 100000 == 0) ||
                    (r % 1000000 == 0) ||
                    (r == num_rows))) {
            printf("Saving Row: %6ld/%ld %6ld: %.1f, %.1f | R: %.1e M: %.1e HM(%.1f): %.1e\n",
                    r, num_rows, m,
                    e_next, e_prev,
                    stats.record, stats.is_missing_gap,
                    config.min_merit, stats.highmerit);
        }

        BIND_OR_ERROR(sqlite3_bind_int64, stmt, 1, rid);

        // P, D, m
        BIND_OR_ERROR(sqlite3_bind_int, stmt, 2, config.p);
        BIND_OR_ERROR(sqlite3_bind_int, stmt, 3, config.d);
        BIND_OR_ERROR(sqlite3_bind_int64, stmt, 4, m);

        // prob_record, prob_missing, prob_merit
        BIND_OR_ERROR(sqlite3_bind_double, stmt, 5, stats.record);
        BIND_OR_ERROR(sqlite3_bind_double, stmt, 6, stats.is_missing_gap);
        BIND_OR_ERROR(sqlite3_bind_double, stmt, 7, stats.highmerit);

        // e_next, e_prev
        BIND_OR_ERROR(sqlite3_bind_double, stmt, 8, e_next);
        BIND_OR_ERROR(sqlite3_bind_double, stmt, 9, e_prev);

        int rc = sqlite3_step(stmt);
        if (rc != SQLITE_DONE) {
            printf("\nm_stats insert failed %ld: (%d, %d, %ld): %d: %s\n",
                row_i, config.p, config.d, m, rc, sqlite3_errmsg(db));
            break;
        }

        if (sqlite3_reset(stmt) != SQLITE_OK) {
            printf("Failed to reset statement\n");
        }

        if (sqlite3_clear_bindings(stmt) != SQLITE_OK) {
            printf("Failed to clear bindings\n");
        }
    }

    if (sqlite3_exec(db, "END TRANSACTION", nullptr, nullptr, &zErrMsg) != SQLITE_OK) {
        printf("END TRANSACTION failed: %s\n", zErrMsg);
        exit(1);
    }

    if (config.verbose >= 0) {
        printf("Saved %ld rows to 'm_stats' table\n",
                num_rows);
    }
}


cdouble NTH_PRIME_CUTOFF = 1e-12;
// Chosen randomly (not sure how to improve)
cdouble COMBINED_CUTOFF = 1e-14;

float nth_prob_or_zero(const vector<float>& prob_nth, size_t nth) {
    return nth < prob_nth.size() ? prob_nth[nth] : 0.0f;
}

void prob_nth_prime(
        double prob_prime,
        vector<float>& prob_prime_nth,
        vector<float>& prob_great_nth) {
    assert(prob_prime_nth.empty() && prob_great_nth.empty());

    double prob_still_prime = 1.0;
    while (prob_still_prime > NTH_PRIME_CUTOFF) {
        prob_prime_nth.push_back(prob_still_prime * prob_prime);
        prob_great_nth.push_back(prob_still_prime);
        prob_still_prime *= 1 - prob_prime;
    }

}


void prob_combined_gap(
        double prob_prime,
        vector<float>& prob_combined) {
    double prob = prob_prime * prob_prime;
    while (prob > COMBINED_CUTOFF) {
        prob_combined.push_back(prob);
        prob *= 1 - prob_prime;
    }
}


void prob_extended_gap(
        const struct Config& config,
        cdouble PROB_PRIME,
        const vector<float>& records,
        const vector<uint32_t>& poss_record_gaps,
        ProbNth &gap_probs) {

    cdouble N_log = calc_log_K(config) + log(config.mstart);

    const unsigned int SL = config.sieve_length;
    const unsigned int MIN_RECORD = poss_record_gaps.front();
    // Gaps larger than this are assumed to be record
    const unsigned int MAX_RECORD = poss_record_gaps.back();

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

    // Extended size should be up to MAX_RECORD (~30 merit)
    // but not more than 100'000 or 2*SL
    const size_t EXT_SIZE = std::min(MAX_RECORD, std::max(2*SL, 100'000u));

    if (config.verbose >= 2) {
        printf("Extended size: %ld (%.1f merit)\n", EXT_SIZE, EXT_SIZE / N_log);
    }
    if (config.sieve_length >= EXT_SIZE) {
        // Empty extended.
        gap_probs.wheel_d = 1;
        gap_probs.prob_greater_extended[0] = 1;
        gap_probs.prob_greater_extended[1] = 1;

        vector<float> prob_extended_record(SL+1, 0.0);
        gap_probs.extended_record_next[0] = prob_extended_record;
        gap_probs.extended_record_next[1] = prob_extended_record;
        gap_probs.extended_extended_record[0] = 0.0;
        gap_probs.extended_extended_record[1] = 0.0;
        return;
    }

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

    vector<uint32_t> wheel_primes;
    map<uint32_t, uint32_t> k_mod_p;

    uint32_t wheel = 1;
    for (uint32_t p : {2u, 3u, 5u, 7u}) {
        if (config.d % p == 0) {
            wheel *= p;
            prob_prime_coprime /= (1 - 1.0 / p);

            uint32_t k_mod = 1;
            for (auto k : K_primes) {
                if (config.d % k != 0)
                    k_mod = (k_mod * k) % p;
            }

            wheel_primes.push_back(p);
            k_mod_p[p] = k_mod;
        }
    }
    gap_probs.wheel_d = wheel;

    // Same as prob_prime_nth but not considering sieving (because outside of interval).
    vector<float> prob_prime_nth_out;
    vector<float> prob_greater_nth_out;
    prob_nth_prime(prob_prime_coprime, prob_prime_nth_out, prob_greater_nth_out);

    // For each wheel mark off a divisors of small primes in d.
    map<int, vector<char>> coprime_ms;

    { // Calculate <m, boolean array of coprime[0 < X < 2*SL]>
        float average_inner_coprime = 0;
        float average_extended_coprime = 0;

        for (uint64_t m = 0; m < wheel; m++) {
            // Hack to make `prob_record_vs_plimit` much faster
            if (config.minc == 1) {
                uint32_t temp = config.mstart % wheel;
                if (temp != m && temp != (wheel - m))
                    continue;
            }

            // Don't need any where m is coprime with d
            if (gcd(m, wheel) > 1) continue;

            // Copy is_coprime
            vector<char> is_coprime_m(is_coprime);

            // Mark off multiples of d primes
            for (const auto p : wheel_primes) {
                assert(config.d % p == 0);

                // (m * K) % p;
                uint64_t first = ((__int128) m * k_mod_p[p]) % p;

                // first multiple on the positive side: -(m*K) % p
                for (size_t i = p - first; i < EXT_SIZE; i += p) {
                    is_coprime_m[i] = false;
                }
            }
            coprime_ms[m] = is_coprime_m;

            {
                auto middle = is_coprime_m.begin() + SL;
                long inner_coprime    = std::count(is_coprime_m.begin(), middle, true);
                long extended_coprime = std::count(middle, is_coprime_m.end(), true);
                average_inner_coprime += inner_coprime;
                average_extended_coprime += extended_coprime;

                cdouble prob_outer = nth_prob_or_zero(prob_greater_nth_out, extended_coprime);
                gap_probs.prob_greater_extended[m] = prob_outer;

                if (config.verbose >= 3) {
                    printf("\tWheel: %-3ld %ld/%d inner, %ld/%d extended coprime, prob last: %.3g\n",
                        m, inner_coprime, SL, extended_coprime, SL, prob_outer);
                }
            }
        }

        average_inner_coprime /= coprime_ms.size();
        average_extended_coprime /= coprime_ms.size();

        cdouble prob_inner = nth_prob_or_zero(prob_greater_nth_out, average_inner_coprime);
        cdouble prob_outer = nth_prob_or_zero(prob_greater_nth_out, average_extended_coprime);

        if (config.verbose >= 2) {
            printf("Using Wheel: %d for extended probs\n", wheel);
            printf("\tAverage %5.0f inner    coprimes => %.3g%% prob_greater\n",
                average_inner_coprime,     100 * prob_inner);
            printf("\tAverage %5.0f extended coprimes => %.3g%% prob_greater\n",
                average_extended_coprime, 100 * prob_outer);
        }
    }

    for (uint32_t m = 0; m < wheel; m++) {
        // Hack to make `prob_record_vs_plimit` much faster
        if (config.minc == 1) {
            uint32_t temp = config.mstart % wheel;
            if (temp != m && temp != (wheel - m))
                continue;
        }
        if (gcd(m, wheel) > 1) continue;

        vector<char> &is_coprime_m = coprime_ms.at(m);
        // -m % wheel => wheel - m => (wheel - m) % wheel to handle wheel = 1
        vector<char> &is_coprime_m_prev = coprime_ms.at((wheel - m) % wheel);

        // List of coprime X, SL < X < EXT_SIZE
        vector<uint32_t> extended_coprime;

        // Probability of prev <= SL, next > SL (extended)
        {
            // partial_sum of is_coprime_m starting at SL+1
            vector<uint32_t> count_coprime_m(EXT_SIZE, 0);
            {
                for (size_t x = SL+1; x < EXT_SIZE; x++) {
                    if (is_coprime_m[x])
                        extended_coprime.push_back(x);
                    count_coprime_m[x] = extended_coprime.size();
                }
            }

            vector<float> prob_extended_record(SL+1, 0.0);

            // Essentially a double loop over gap_prev, gap_next but the inner
            // gap_next loop is replaced with a smaller loop over record_gap.
            for (size_t gap_prev = 1; gap_prev <= SL; gap_prev++) {
                // only needed for values that can be coprime with K
                if (!is_coprime_m_prev[gap_prev]) {
                    prob_extended_record[gap_prev] = std::nan("");
                    continue;
                }

                double prob_record = 0;
                for (uint32_t record_gap : poss_record_gaps ) {
                    uint32_t dist = record_gap - gap_prev;
                    if (dist <= SL) continue;

                    if (dist >= EXT_SIZE) break;

                    // dist can never be prime.
                    if (!is_coprime_m[dist]) continue;

                    // This is the nth possible prime after SL
                    uint32_t num_coprime = count_coprime_m[dist];
                    if (num_coprime >= prob_prime_nth_out.size()) break;

                    // chance of dist being first prime.
                    prob_record += prob_prime_nth_out[num_coprime];
                }

                // if prev + next > MAX_RECORD add the remainder of prob_next_larger
                uint32_t gap_max_record = MAX_RECORD - gap_prev;
                if (SL < gap_max_record && gap_max_record < EXT_SIZE) {
                    size_t coprimes_max_record = count_coprime_m[gap_max_record];
                    // sum(prob_prime_nth[0..i]) + prob_greater_nth_out[i+1] == 1.0
                    prob_record += nth_prob_or_zero(prob_greater_nth_out, coprimes_max_record+1);
                }

                // Prob record gap, with 1 <= gap_prev <= SL, SL <= gap_next
                assert(0 <= prob_record && prob_record <= 1);
                prob_extended_record[gap_prev] = prob_record;
            }
            gap_probs.extended_record_next[m] = prob_extended_record;
        }

        // Probability of prev, next > SL (extended^2)
        {
            // XXX: this seems to overestimate by 5-20%,
            // The best explanation I have is that this relates to how factors
            // aren't uniformly distributed. One piece of evidence is that
            // reducing wheel_primes causes worse estimation. I suspect that
            // This relates to the SL_factors_of_d and how primes are
            // over-represented then under-represented after

            // gap_prev + extended_coprime[min_i] <= MIN_RECORD
            size_t min_e_c_i = extended_coprime.size();
            // gap_prev + extended_coprime[max_i] <= MAX_RECORD
            size_t max_e_c_i = extended_coprime.size() - 1;

            double prob_e2_record = 0;
            size_t extended_coprimes_prev = 0;
            for (size_t gap_prev = SL + 1; gap_prev < EXT_SIZE; gap_prev++) {
                if (!is_coprime_m_prev[gap_prev]) {
                    continue;
                }
                // gap_prev is a coprime
                extended_coprimes_prev += 1;

                // If to many coprimes far any reasonable chance.
                if (extended_coprimes_prev >= prob_prime_nth_out.size()) {
                    break;
                }

                // NOTE: It would probably be faster to loop over coprimes (vs records)
                // This loops handles [2*SL, 4*SL+] which is generally >20-40 merit

                while (max_e_c_i && (gap_prev + extended_coprime[max_e_c_i] > MAX_RECORD)) {
                    max_e_c_i -= 1;
                }

                while (min_e_c_i && (gap_prev + extended_coprime[min_e_c_i - 1] >= MIN_RECORD)) {
                    min_e_c_i -= 1;
                }

                if (max_e_c_i == 0) {
                    assert(min_e_c_i == 0);
                    // Every gap_prev + extended_coprime[i] > MAX_RECORD
                    assert(extended_coprimes_prev >= 1);
                    assert(extended_coprimes_prev <= prob_greater_nth_out.size());
                    prob_e2_record += prob_greater_nth_out[extended_coprimes_prev - 1];
                    break;
                }

                assert( min_e_c_i == extended_coprime.size() ||
                        gap_prev + extended_coprime[min_e_c_i] >= MIN_RECORD );
                assert( gap_prev + extended_coprime[max_e_c_i] <= MAX_RECORD );
                assert( (min_e_c_i == 0) || (gap_prev + extended_coprime[min_e_c_i-1] < MIN_RECORD) );
                assert( (max_e_c_i == extended_coprime.size() - 1) ||
                        (gap_prev + extended_coprime[max_e_c_i+1] > MAX_RECORD) );

                float prob_e_e = 0;

                // When we get to an arbitrary large merit assume all things larger are record
                size_t max_i = std::min(max_e_c_i, prob_prime_nth_out.size() - extended_coprimes_prev);
                for (size_t i = min_e_c_i; i < max_i; i++) {
                    size_t gap = gap_prev + extended_coprime[i];
                    if (records[gap] > N_log) {
                        prob_e_e += prob_prime_nth_out[i];
                    }
                }

                // Everything >= max_i is assumed record (or very small prob)
                prob_e_e += nth_prob_or_zero(prob_greater_nth_out, max_i);
                prob_e2_record += prob_e_e * nth_prob_or_zero(prob_prime_nth_out, extended_coprimes_prev);
            }
            gap_probs.extended_extended_record[m] = prob_e2_record;
        }
    }
}


void setup_prob_nth(
        const struct Config &config,
        const vector<float> &records,
        const vector<uint32_t> &poss_record_gaps,
        ProbNth &gap_probs) {
    cdouble N_log = calc_log_K(config) + log(config.mstart);

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
        gap_probs.prime_nth, gap_probs.greater_nth);

    prob_combined_gap(
        PROB_PRIME_AFTER_SIEVE,
        gap_probs.combined_nth);

    // Prob record with gap[i] and other gap > SL
    {
        auto s_start_t = high_resolution_clock::now();

        prob_extended_gap(
            config,
            PROB_PRIME,
            records,
            poss_record_gaps,
            gap_probs
        );

        if (config.verbose >= 1) {
            auto s_stop_t = high_resolution_clock::now();
            double   secs = duration<double>(s_stop_t - s_start_t).count();
            printf("Extended prob records setup (%.2f seconds)\n\n", secs);
        }
    }
}

ProbM calculate_probm(
        const uint64_t &m, float log_M,
        const vector<int32_t> &unknown_prev,
        const vector<int32_t> &unknown_next,
        const Config &config,
        const vector<float> &records,
        const uint32_t &min_record_gap,
        const uint32_t &min_gap_min_merit,
        const ProbNth &gap_probs,
        vector<float> &prob_gap_norm, vector<float> &prob_gap_prev, vector<float> &prob_gap_next) {

    int m_wheel_next = m % gap_probs.wheel_d;
    // want -m % wheel_d => wheel_d - m
    int m_wheel_prev = (gap_probs.wheel_d - m_wheel_next) % gap_probs.wheel_d;

    /**
     * Directly examined (1 - PROB_PREV_GREATER) * (1 - PROB_NEXT_GREATER)
     * +
     * Extended examined (1 - PROB_PREV_GREATER) * (PROB_NEXT_GREATER * (1 - prob_extended))
     * +
     * Extended^2        (PROB_PREV_GREATER * PROB_NEXT_GREATER) * (1 - prob_extended)^2
     * =
     *
     * ^
     * | ??????????????????????
     * |                      ?
     * |-------------------   ?
     * e extended | extend |  ?
     * x record   | extend |  ?
     * t prev     | record |  ?
     * | here     |        |  ?
     * |----------.--------.  ?
     * p DIRECT   | extend |  ?
     * r COMPUT   | record |  ?
     * e ATION    | high   |  ?
     * v HERE     |        |  ?
     * *- next(SL)|EXT_SIZE| >EXT_SIZE--->
     */

    cdouble PROB_PREV_GREATER = nth_prob_or_zero(gap_probs.greater_nth, unknown_prev.size());
    cdouble PROB_NEXT_GREATER = nth_prob_or_zero(gap_probs.greater_nth, unknown_next.size());
    cdouble prob_prev_extended = gap_probs.prob_greater_extended.at(m_wheel_prev);
    cdouble prob_next_extended = gap_probs.prob_greater_extended.at(m_wheel_next);

    ProbM result;
    result.seen = (1 - PROB_PREV_GREATER * prob_prev_extended) *
                  (1 - PROB_NEXT_GREATER * prob_next_extended);

    result.record_extended2 = PROB_NEXT_GREATER * PROB_PREV_GREATER *
                              gap_probs.extended_extended_record.at(m_wheel_next);

    bool save_gap_probs = !prob_gap_norm.empty();

    { // Direct probability (both primes <= SL)
        uint32_t min_interesting_gap = std::min(min_gap_min_merit, min_record_gap);
        size_t max_i = std::min(unknown_prev.size(), gap_probs.combined_nth.size());
        size_t min_j = unknown_next.size();
        for (size_t i = 0; i < max_i; i++) {
            int32_t gap_prev = unknown_prev[i];
            while ((min_j > 0) && ((uint32_t) gap_prev + unknown_next[min_j-1] >= min_interesting_gap)) {
                min_j -= 1;
            }

            size_t max_j = std::min(unknown_next.size(), gap_probs.combined_nth.size() - i);

            // Starting at min_j causes some `prob_this_gap` to be skipped,
            // but is a sizeable speedup for large gaps.
            size_t j = (save_gap_probs || config.sieve_length < 100'000) ? 0 : min_j;
            if (j > 0) {
                // sum of gap_probs.prime_nth[0..j-1]
                float sum_prob_jth = 1.0 - gap_probs.greater_nth[j];
                prob_gap_norm[0] += gap_probs.prime_nth[i] * sum_prob_jth;
            }

            for (; j < max_j; j++) {
                int32_t gap_next = unknown_next[j];
                uint32_t gap = gap_prev + gap_next;

                // Same as prob_prime_nth[i] * prob_prime_nth[j];
                float prob_this_gap = gap_probs.combined_nth[i + j];

                if (save_gap_probs) {
                    // Used in gap_test_plotting.py, has performance impact to save.
                    prob_gap_norm[gap] += prob_this_gap;
                }

                if (gap >= min_gap_min_merit) {
                    result.highmerit += prob_this_gap;
                }

                if (gap >= min_record_gap && records[gap] > log_M) {
                    result.record_inner += prob_this_gap;

                    if (MISSING_GAPS_MIN <= gap && records[gap] == GAP_INF) {
                        result.is_missing_gap += prob_this_gap;
                    }
                }
            }
        }
    }

    { // Extended gap (one prime <= SL, one prime > SL)
        // See `prob_extended_gap`
        const vector<float> &extended_record_next = gap_probs.extended_record_next.at(m_wheel_next);
        const vector<float> &extended_record_prev  = gap_probs.extended_record_next.at(m_wheel_prev);

        // Smallest side to have reached --min-merit when other side is extended
        const int32_t min_side_for_extended_min_merit = min_gap_min_merit - config.sieve_length;

        size_t max_i = std::max(unknown_prev.size(), unknown_next.size());
        // i >= prime_nth.size() has tiny probability (see DOUBLE_NTH_PRIME_CUTOFF)
        for (size_t i = 0; i < std::min(max_i, gap_probs.prime_nth.size()); i++) {
            float prob_i = gap_probs.prime_nth[i]; // 0-indexed
            assert(0 <= prob_i && prob_i <= 1.0);

            // unknown[i'th] is prime, on the other side have prime be outside of sieve.
            if (i < unknown_prev.size()) {
                float conditional_prob = extended_record_next[unknown_prev[i]];
                assert(0 <= conditional_prob && conditional_prob <= 1);

                result.record_extended += prob_i * PROB_NEXT_GREATER * conditional_prob;

                int32_t gap_prev = unknown_prev[i];
                result.expected_gap_prev += gap_prev * prob_i;

                if (save_gap_probs)
                    prob_gap_prev[gap_prev] += prob_i;

                if (gap_prev >= min_side_for_extended_min_merit)
                    result.highmerit += prob_i * PROB_NEXT_GREATER;
            }

            if (i < unknown_next.size()) {
                float conditional_prob = extended_record_prev[unknown_next[i]];
                assert(0 <= conditional_prob && conditional_prob <= 1);

                result.record_extended += prob_i * PROB_PREV_GREATER * conditional_prob;

                int32_t gap_next = unknown_next[i];
                result.expected_gap_next += gap_next * prob_i;

                if (save_gap_probs)
                    prob_gap_next[gap_next] += prob_i;

                if (gap_next >= min_side_for_extended_min_merit)
                    result.highmerit += prob_i * PROB_PREV_GREATER;
            }
        }
    }

    result.record = result.record_inner + result.record_extended + result.record_extended2;
    return result;
}


/**
 * Calculate prob record at various max prime values
 * Takes unknown_file contains [(prime1, X1), (prime2, X2), ...]
 *
 * prime1 is a factor of X1
 * prime2 is the next smallest factor and divides X2
 * ...
 *
 * unknown_file is created with primegapverify print_factors branch
 * time primegapverify/large_sieve 73 1511 2190 -15000 30000 2000000000 > 1511_2190_73_1_s15000_l2000M.txt
 */
void prob_record_vs_plimit(struct Config config) {
    const unsigned int SL = config.sieve_length;
    const unsigned int SIEVE_INTERVAL = 2 * SL + 1;
    assert( SL > 1000 );
    assert(config.minc == 1);

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
    mpz_t N;
    mpz_t test;
    mpz_init(test);

    int64_t START_PRIME = std::max(100'000u, 4 * SL);

    double K_log;
    K_stats(config, N, nullptr, &K_log);
    cdouble N_log = K_log + log(config.mstart);

    mpz_mul_ui(N, N, config.mstart);

    // ----- Get Record Prime Gaps
    vector<float> records = get_record_gaps(config);

    // gap that would be a record with m*P#/d
    vector<uint32_t> poss_record_gaps;
    load_possible_records(N_log, records, poss_record_gaps);
    assert( poss_record_gaps.size() >= 2);
    const uint32_t min_record_gap = poss_record_gaps.front();

    // Passing an empty vector to prob_gap_norm causes those values to be skipped.
    vector<float> empty_vector;

    vector<bool> composite(SIEVE_INTERVAL, false);

    int64_t prime = 0, offset = -1;
    char delim = -1;
    while (unknown_file.good()) {
        unknown_file >> prime;
        assert(prime >= 2 && prime <= 10'000'000'000'000LL);
        config.max_prime = prime;

        unknown_file >> delim;
        assert(delim > 0);
        assert(delim == ',');

        unknown_file >> offset;

        assert(offset >= 0 && offset <= SIEVE_INTERVAL);
        {
            mpz_set(test, N);
            mpz_sub_ui(test, test, SL);
            mpz_add_ui(test, test, offset);
            assert( 0 == mpz_fdiv_ui(test, prime) );
        }

        composite[offset] = true;

        if (prime < START_PRIME) {
            continue;
        }

        vector<int32_t> unknown_prev, unknown_next;
        for (size_t x = 1; x <= SL; x++) {
            if (!composite[SL - x]) { unknown_prev.push_back(x); }
            if (!composite[SL + x]) { unknown_next.push_back(x); }
        }

        // Suppress printing in setup_prob_nth.
        config.verbose = 0;
        ProbNth gap_probs;
        setup_prob_nth(config, records, poss_record_gaps, gap_probs);

        ProbM probm = calculate_probm(
                config.mstart, N_log, unknown_prev, unknown_next,
                config, records, min_record_gap, /* min_gap_min_merit */ min_record_gap,
                gap_probs, /* prob_gap_norm */ empty_vector, empty_vector, empty_vector);

        if (config.verbose >= 2) {
            // Breakdown of prob inner, extended, extended^2
            printf("%7ld, %.7f = %0.3g + %0.3g + %0.3g (unaccounted: %.7f)\n",
                    config.max_prime, probm.record,
                    probm.record_inner, probm.record_extended, probm.record_extended2,
                    1 - probm.seen);
        }
        cout << config.max_prime << ", " << probm.record << endl;
    }
    mpz_clear(test);
    mpz_clear(N);
}


void run_gap_file(
        /* input */
        const struct Config& config,
        const mpz_t &K,
        const float K_log,
        const vector<float>& records,
        const uint32_t min_record_gap,
        const uint32_t min_gap_min_merit,
        const ProbNth &gap_probs,
        vector<uint32_t>& valid_mi,
        std::ifstream& unknown_file,
        /* output */
        vector<float>& prob_gap_norm,
        vector<float>& prob_gap_prev,
        vector<float>& prob_gap_next,
        unordered_map<uint64_t,ProbM> &M_stats) {

    auto s_start_t = high_resolution_clock::now();

    const uint32_t SIEVE_INTERVAL = 2 * config.sieve_length + 1;


    prob_gap_norm.clear();
    prob_gap_prev.clear();
    prob_gap_next.clear();

    // NOTE: prob_gap_prev uses values <= SL, all being the same size helps in store_stats
    prob_gap_norm.resize(SIEVE_INTERVAL, 0);
    prob_gap_prev .resize(SIEVE_INTERVAL, 0);
    prob_gap_next.resize(SIEVE_INTERVAL, 0);

    // Keep sum & max of several records
    ProbM sum = {};
    ProbM max_r = {};
    max_r.record = max_r.highmerit = max_r.is_missing_gap= 1e-10;

    if (config.verbose >= 1) {
        printf("\n%ld tests M_start(%ld) + mi(%d to %d)\n\n",
                valid_mi.size(), config.mstart,
                valid_mi.front(), valid_mi.back());
    }

    BitArrayHelper helper(config, K);

    #pragma omp parallel num_threads(config.threads)
    {
        // Declare thread local prob_gap norm/low/high (later they get sum'ed)
        vector<float> local_p_gap_norm(SIEVE_INTERVAL, 0);
        vector<float> local_p_gap_prev(SIEVE_INTERVAL, 0);
        vector<float> local_p_gap_next(SIEVE_INTERVAL, 0);

        #pragma omp for ordered schedule(static, 1)
        //for (uint32_t mi : valid_m) {
        for (size_t i = 0; i < valid_mi.size(); i++) {
            uint32_t mi = valid_mi[i];

            vector <int32_t> unknown_prev, unknown_next;

            /**
             * Note: it would be nice to use ordered here but omp doesn't end
             * the ordered section, omp critical could be used with
             * parse_unknown_line returning <m> but still not optimal.
             * https://stackoverflow.com/questions/43540605/must-ordered
             */

            std::string line;
            #pragma omp critical
            {
                std::getline(unknown_file, line);
            }
            std::istringstream iss_line(line);

            /**
             * XXX: For P=1511 With compression == 2, >50% of runtime is reading unknown_line
             * Already optimized:
             *      unbundled file read from parsing
             *      special cased D % 2
             */
            uint64_t mtest = config.mstart + mi;
            uint64_t m = parse_unknown_line(
                config, helper,
                mtest, iss_line, unknown_prev, unknown_next);

            ProbM probm;
            { // This section in parallel
                // Note slightly different from N_log
                float log_M = K_log + log(m);
                probm = calculate_probm(m, log_M, unknown_prev, unknown_next,
                                        config, records, min_record_gap, min_gap_min_merit,
                                        gap_probs,
                                        local_p_gap_norm, local_p_gap_prev, local_p_gap_next);

                sum.record += probm.record;
                sum.record_inner += probm.record_inner;
                sum.record_extended += probm.record_extended;
                sum.record_extended2 += probm.record_extended2;

                sum.record_extended2 += line.size();
            }

            #pragma omp critical
            { // Handle saving and stats 1 thread at a time
                // Becuase of note above we validate nothing saved for this mi before.
                assert(M_stats.count(mi) == 0);
                M_stats[mi] = probm;

                if (config.verbose >= 1) {
                    if (probm.record > max_r.record) {
                        max_r.record = probm.record;

                        printf("RECORD :%-6ld line %-6ld  unknowns: %3ld, %3ld\t| "
                               "prob record: %.2e   (%.2e + %.2e + %.2e)\n",
                               m, M_stats.size(), unknown_prev.size(), unknown_next.size(),
                               probm.record, probm.record_inner,
                               probm.record_extended, probm.record_extended2);
                    }

                    if (probm.highmerit > max_r.highmerit) {
                        max_r.highmerit = probm.highmerit;
                        printf("MERIT  :%-6ld line %-6ld  unknowns: %3ld, %3ld\t| "
                               "      merit: %.2g\n",
                               m, M_stats.size(), unknown_prev.size(), unknown_next.size(),
                               probm.highmerit);
                    }
                }

                if (config.verbose >= 2) {
                    if (probm.is_missing_gap > max_r.is_missing_gap) {
                        max_r.is_missing_gap= probm.is_missing_gap;
                        printf("MISSING:%-6ld line %-6ld  unknowns: %3ld, %3ld\t| "
                               "prob record: %.2e   missing: %.4e\n",
                               m, M_stats.size(), unknown_prev.size(), unknown_next.size(),
                               probm.record, probm.is_missing_gap);
                    }
                }
            }
        }

        #pragma omp critical
        for (size_t i = 0; i < local_p_gap_norm.size(); i++) {
            prob_gap_norm[i] += local_p_gap_norm[i];
            prob_gap_prev[i] += local_p_gap_prev[i];
            prob_gap_next[i] += local_p_gap_next[i];
        }
    }

    /**
     * Workaround the openmp ordered limitations.
     * Double check that every valid_mi was stored once.
     */
    assert(valid_mi.size() == M_stats.size());
    for (uint64_t mi : valid_mi) {
        assert(M_stats.count(mi) == 1);
    }

        // Normalize the probability of gap (across all m) to per m
    for (size_t i = 0; i < prob_gap_norm.size(); i++) {
        prob_gap_norm[i] /= valid_mi.size();
        prob_gap_prev[i]  /= valid_mi.size();
        prob_gap_next[i] /= valid_mi.size();
    }

    if (config.verbose >= 0) {
        long  s_tests = M_stats.size();
        auto s_stop_t = high_resolution_clock::now();
        double   secs = duration<double>(s_stop_t - s_start_t).count();
        printf("%ld m's processed in %.2f seconds (%.2f/sec)\n",
                s_tests, secs, s_tests / secs);

        cout << endl;
        printf("prob record total: %.5f  direct: %.5f  extended: %.5f  extended^2: %.5f\n",
                sum.record, sum.record_inner,
                sum.record_extended, sum.record_extended2);
    }
    if (config.verbose >= 1) {
        vector<float> probs_seen;
        for (auto &m_stat : M_stats) {
            probs_seen.push_back(m_stat.second.seen);
        }
        cout << endl;
        printf("\tsum(prob(gap[X])): %.5f\n", average_v(prob_gap_norm) * prob_gap_norm.size());
        printf("\tavg missing prob : %.7f\n", 1 - average_v(probs_seen));
        cout << endl;
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

    if (config.verbose >= 2) {
        printf("\n");
        printf("%sieve time: %.0f seconds (%.2f hours)\n",
            exact ? "S" : "Estimated s", combined_time, combined_time / 3600);
        printf("Estimated time/m: 2 * (%.1f PRP/m / %.1f PRP/s) = %.2f seconds\n",
            estimated_prp_per_m, 1 / prp_time_est, 2 * time_per_side);
        printf("\n");
    }

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
    if (!exact) {
        printf("\tUsed estimate for combined_sieve timing.\n");
    }
    printf("\tEstimate of optimal printed with *\n");
    printf("\n");

    // Both sides & One sided at 1% (assume 80% of prob also)
    for (int side_percent : {100, 1}) {

        double sum_prob = 0.0;
        double time = combined_time;

        if (side_percent != 100) {
            printf("\tAssuming %d%% of next_prime(...) are skipped\n", 100 - side_percent);
        }

        bool max_happened = false;
        for (size_t i = 1; i <= sorted.size(); i++) {
            // Print at print_points (1, 5, 10, ... 100% and also at "optimal" percent)
            double sides_tested = 1 + side_percent / 100.0;
            double add_t = sides_tested * time_per_side;
            double add_p = sorted[i-1] * (side_percent == 100 ? 1.0 : 0.8);

            time += add_t;
            sum_prob += add_p;

            double avg = sum_prob / time;
            bool is_below_avg = ((add_p / add_t) < avg) || i == sorted.size();
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

        config.compression = Args::guess_compression(config, unknown_file);
    }

    // ----- Merit Stuff
    mpz_t K;

    int K_digits;
    double K_log;
    K_stats(config, K, &K_digits, &K_log);
    double N_log = K_log + log(config.mstart);

    uint32_t min_gap_min_merit = std::ceil(config.min_merit * N_log);
    if (config.verbose >= 2) {
        printf("Min Gap ~= %d (for merit > %.1f)\n\n", min_gap_min_merit, config.min_merit);
    }

    // ----- Get Record Prime Gaps
    const vector<float> records = get_record_gaps(config);

    // gap that would be a record with m*P#/d
    vector<uint32_t> poss_record_gaps;
    load_possible_records(N_log, records, poss_record_gaps);
    assert( poss_record_gaps.size() >= 2);
    {
        if (config.verbose >= 2) {
            printf("Found %ld possible record gaps (%d to %d)\n",
                poss_record_gaps.size(), poss_record_gaps.front(), poss_record_gaps.back());

            for (int gi = 0; gi <= 2; gi++) {
                size_t g = poss_record_gaps[gi];
                printf("\tIf found Gap: %ld (current: %.2f) would improve to %.3f\n",
                    g, g / records[g], g / N_log);
            }
        }

        if (poss_record_gaps.front() > 3 * SIEVE_LENGTH) {
            printf("\n\n\n");
            printf("\tHard to determine record prob, 3 * sieve_length < min_record_gap\n");
            printf("\tPlease set a larger sieve_length in the future\n");
            printf("\n\n\n");
        }

        if (config.verbose >= 2) {
            cout << endl;
        }
    }

    ProbNth gap_probs;
    setup_prob_nth(config, records, poss_record_gaps, gap_probs);

    vector<uint32_t> valid_mi = is_coprime_and_valid_m(config).second;

    /* Over all m values */
    vector<float> prob_gap_norm;
    vector<float> prob_gap_prev;
    vector<float> prob_gap_next;

    /* Per m stats & probabilities */
    unordered_map<uint64_t,ProbM> M_stats;

    // ----- Main calculation
    run_gap_file(
        /* Input */
        config, K, K_log,
        records, poss_record_gaps.front(), min_gap_min_merit,
        gap_probs,
        valid_mi,
        /* sieve input */
        unknown_file,
        /* output */
        prob_gap_norm, prob_gap_prev, prob_gap_next,
        M_stats
    );

    mpz_clear(K);

    vector<float> expected_gap;
    vector<float> probs_record;
    vector<float> probs_missing;
    vector<float> probs_highmerit;

    for (auto m_stat : M_stats) {
        auto stats = m_stat.second;
        expected_gap.push_back(stats.expected_gap_prev + stats.expected_gap_next);
        probs_record.push_back(stats.record);
        probs_missing.push_back(stats.is_missing_gap);
        probs_highmerit.push_back(stats.highmerit);
    }

    // Compute sum_missing_prob, sum_record_prob @1,5,10,20,50,100%
    if (config.verbose >= 2) {
        prob_stats("EXPECTED GAP", expected_gap);
        prob_stats("RECORD", probs_record);

        double avg_missing = average_v(probs_missing);
        double avg_record  = average_v(probs_record);
        // missing mostly includes > 3 * SL, which is likely to be a record.
        double uncertainty = avg_missing / (avg_missing + avg_record);

        if (uncertainty > 1e-5) {
            printf("\tRECORD : avg: %.2e | missing: %.2e | uncertainty: %.4f%% \n",
                avg_record, avg_missing, 100 * uncertainty);
        }

        if (config.verbose >= 2) {
            double avg_merit = average_v(probs_highmerit);
            if (avg_merit > 1e-5) {
                prob_stats("MERIT", probs_highmerit);
            }

            if (avg_missing > 1e-5) {
                prob_stats("MISSING", probs_missing);
            }
        }
        printf("\n");
    }

    if (config.save_unknowns && !config.testing) {
        auto s_stop_t = high_resolution_clock::now();
        double   secs = duration<double>(s_stop_t - s_start_t).count();
        // Account for use of multiple threads
        secs *= config.threads;

        store_stats(
            config,
            secs,
            prob_gap_norm, prob_gap_prev, prob_gap_next,
            valid_mi, M_stats
        );
    }

    if (config.verbose >= 1) {
        calculate_prp_top_percent(config, valid_mi.size(), N_log, probs_record);
    }
}
