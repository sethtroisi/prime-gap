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
#include <queue>
#include <tuple>
#include <utility>
#include <vector>

#include <gmp.h>
#include <sqlite3.h>

#include "gap_common.h"

using std::cout;
using std::endl;
using std::pair;
using std::priority_queue;
using std::tuple;
using std::vector;
using namespace std::chrono;

typedef tuple<float,uint64_t,vector<pair<uint32_t,uint32_t>>> missing_gap_record;
typedef priority_queue<missing_gap_record,
                       vector<missing_gap_record>,
                       std::greater<missing_gap_record>> min_queue_gaps;

const char records_db[] = "gaps.db";
const char gaps_db[]    = "prime-gaps.db";

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


int main(int argc, char* argv[]) {
    printf("\tCompiled with GMP %d.%d.%d\n\n",
        __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);

    Config config = argparse(argc, argv);

    if (config.run_prp == 1) {
        cout << "Must set --sieve-only for gap_search" << endl;
        return 1;
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


uint64_t config_hash(const struct Config& config) {
    // Hash the config to a
    uint64_t hash =    config.mstart;
    hash = hash * 31 + config.minc;
    hash = hash * 31 + config.p;
    hash = hash * 31 + config.d;
    hash = hash * 31 + config.sieve_length;
    hash = hash * 31 + config.sieve_range;
    return hash;
}


bool is_range_already_processed(const struct Config& config) {
    uint64_t hash = config_hash(config);
    char sql[200];
    sprintf(sql, "SELECT count(*) FROM to_process_range WHERE id = %ld", hash);
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
        printf("\nto_process_range SELECT failed %s | %d: %s\n",
            zErrMsg, rc, sqlite3_errmsg(db));
        exit(1);
    }
    return count > 0;
}


void store_stats(
        const struct Config& config,
        float K_log,
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

    char *zErrMsg = 0;
    if (sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &zErrMsg) != SQLITE_OK) {
        printf("BEGIN TRANSACTION failed: %s\n", zErrMsg);
        exit(1);
    }

    /* Create SQL statement */
    char sql[] = "INSERT INTO partial_result VALUES(NULL, ?, ?, ?, 0, 0, ?, ?, ?)";

    sqlite3_stmt *stmt;
    /* Prepare SQL statement */
    int rc = sqlite3_prepare_v2(db, sql, -1, &stmt, 0);
    if (rc != SQLITE_OK) {
        printf("Could not prepare statement: %s\n", sql);
        exit(1);
    }

    const size_t num_rows = M_vals.size();
    for (size_t i = 0; i < num_rows; i++) {
        uint64_t m = M_vals[i];
        int e_next = expected_next[i];
        int e_prev = expected_prev[i];
        float e_merit = (e_next + e_prev) / (log(m) + K_log);

        size_t r = i + 1;
        if ((r < 5) || (r < 500 && r % 100 == 0) || (r % 2000 == 0) || r == num_rows) {
            printf("Row: %6ld/%ld %ld, %d, %d, %.2f\n",
                r, num_rows,
                m, e_next, e_prev, e_merit);
        }

        int binded = 0;
        if (binded == 0 && sqlite3_bind_int(stmt, 1, m) != SQLITE_OK)
                binded = 1;

        if (binded == 0 && sqlite3_bind_int(stmt, 2, config.p) != SQLITE_OK)
            binded = 2;

        if (binded == 0 && sqlite3_bind_int(stmt, 3, config.d) != SQLITE_OK)
            binded = 3;

        if (binded == 0 && sqlite3_bind_int(stmt, 4, e_next) != SQLITE_OK)
            binded = 4;
        if (binded == 0 && sqlite3_bind_int(stmt, 5, e_prev) != SQLITE_OK)
            binded = 5;

        if (binded == 0 && sqlite3_bind_double(stmt, 6, e_merit) != SQLITE_OK)
            binded = 6;

        if (binded != 0) {
            printf("Failed to bind param %d: %s\n", binded, sqlite3_errmsg(db));
            break;
        }

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

    uint64_t hash = config_hash(config);
    char sSQL[200];
    sprintf(sSQL, "INSERT INTO to_process_range VALUES("
                  "%ld, %ld, %ld, %d, %d, %d, %ld, %ld)",
            hash,
            config.mstart, config.minc,
            config.p, config.d,
            config.sieve_length, config.sieve_range,
            num_rows);

    rc = sqlite3_exec(db, sSQL, NULL, NULL, &zErrMsg);
    if (rc != SQLITE_OK) {
        printf("\nto_process_range INSERT failed %d: %s\n",
            rc, sqlite3_errmsg(db));
        exit(1);
    }

    if (sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &zErrMsg) != SQLITE_OK) {
        printf("END TRANSACTION failed: %s\n", zErrMsg);
        exit(1);
    }

    sqlite3_close(db);
}


void K_stats(
        const struct Config& config,
        const mpz_t &K, int *K_digits, float *K_log) {
    *K_digits = mpz_sizeinbase(K, 10);

    long exp;
    double mantis = mpz_get_d_2exp(&exp, K);
    *K_log = log(mantis) + log(2) * exp;

    double m_log = log(config.mstart);
    int K_bits   = mpz_sizeinbase(K, 2);

    cout << endl;
    printf("K = %d bits, %d digits, log(K) = %.2f\n",
        K_bits, *K_digits, *K_log);
    printf("Min Gap ~= %d (for merit > %.1f)\n\n",
        (int) (config.min_merit * (*K_log + m_log)), config.min_merit);
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
void prob_nth_prime(
        const double prob_prime,
        vector<float>& prob_prime_nth,
        vector<float>& prob_great_nth) {
    double prob_still_prime = 1.0;
    for (; prob_still_prime > 1e-13;) {
        prob_prime_nth.push_back(prob_still_prime * prob_prime);
        prob_great_nth.push_back(prob_still_prime);
        prob_still_prime *= 1 - prob_prime;
    }
}


void prob_combined_gap(
        const double prob_prime,
        size_t expected_unknown,
        vector<float>& prob_combined) {
    double prob = prob_prime * prob_prime;
    // Want error < 1e-9 | unknown_i * unkown_j * 1e-15 ~= 2000 * 2000 * 2.5e-6 = 1e-9
    for (; prob > 2.5e-16;) {
        prob_combined.push_back(prob);
        prob *= 1 - prob_prime;
    }
}


void read_unknown_line(
        uint64_t mi,
        std::ifstream& unknown_file,
        vector<uint32_t>& unknown_high,
        vector<uint32_t>& unknown_low) {

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
        const struct Config config,
        const vector<float>& records,
        const vector<uint32_t> poss_record_gaps,
        const uint64_t M_start,
        const uint64_t M_inc,
        const uint64_t D,
        float K_log,
        const vector<float>& prob_prime_nth,
        const vector<float>& prob_great_nth,
        const vector<float>& prob_combined,
        const vector<float>& prob_record_extended_gap,
        std::ifstream& unknown_file,
        vector<uint64_t>& M_vals,
        vector<float>& expected_prev,
        vector<float>& expected_next,
        vector<float>& probs_seen,
        vector<float>& probs_record,
        min_queue_gaps& missing_gaps_search) {

    auto  s_start_t = high_resolution_clock::now();

    const uint32_t min_record_gap = poss_record_gaps.front();

    // Store max(100, MISSING_GAP_SAVE_PERCENT);
    size_t valid_m = 0;
    for (uint64_t mi = 0; mi < M_inc; mi++) {
        if (gcd(M_start + mi, D) == 1) {
            valid_m += 1;
        }
    }
    uint32_t temp = valid_m * MISSING_GAP_SAVE_PERCENT;
    const uint32_t MAX_MISSING_GAPS = temp <= 0 ? 0 : std::max(100U, temp);
    printf("\tSaving %d best missing gap tests\n", MAX_MISSING_GAPS);

    // sum prob_record_inside sieve
    // sum prob_record_outer (extended)
    float sum_prob_inner = 0.0;
    float sum_prob_outer = 0.0;

    float max_p_record = 0;
    float max_m_record = 0;
    for (uint64_t mi = 0; mi < M_inc; mi++) {
        uint64_t m = M_start + mi;
        if (gcd(m, D) != 1) {
            continue;
        }

        vector<uint32_t> unknown_low, unknown_high;
        read_unknown_line(mi, unknown_file, unknown_low, unknown_high);

        // Note slightly different from smallest_log
        float log_start_prime = K_log + log(m);

        // TODO: should this actually be (1 - prob_great_nth[size] + prob_unknown_extended)^2
        double prob_seen = (1 - prob_great_nth[unknown_high.size()]) *
                           (1 - prob_great_nth[unknown_low.size()]);

        double prob_record = 0;
        double prob_is_missing_gap = 0.0;

        vector<pair<uint32_t, uint32_t>> missing_pairs;

        /**
         * XXX: could possible look only at record_gaps
         * replace 'records[gap] > log_start_prime' with 'gap == records_gap[gi]' or something
         */
        size_t min_j = unknown_high.size() - 1;
        for (size_t i = 0; i < unknown_low.size(); i++) {
            uint32_t gap_low = unknown_low[i];
            while ((min_j > 0) && (gap_low + unknown_high[min_j] > min_record_gap)) {
                min_j -= 1;
            }

            for (size_t j = min_j; j < unknown_high.size(); j++) {
                uint32_t gap_high = unknown_high[j];
                uint32_t gap = gap_low + gap_high;

                // TODO: Determine overhead if this is turned off (make make this a compile time DEFINE?)
                if (MAX_MISSING_GAPS > 0) {
                    if (MISSING_GAPS_LOW <= gap && gap <= MISSING_GAPS_HIGH && records[gap] == GAP_INF) {
                        // TODO log min/max i+j, print approx odds given BOTH PRIME

                        // Same as  prob_prime_nth[i] * prob_prime_nth[j];
                        prob_is_missing_gap += prob_combined[i+j];
                        missing_pairs.emplace_back(std::make_pair(gap_low, gap_high));
                    }
                }

                assert(i + j < prob_combined.size());   // TODO: can break at this point
                if (records[gap] > log_start_prime) {
                    prob_record += prob_combined[i+j];
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

        double p_record = prob_record + prob_record_outer;

        sum_prob_inner += prob_record;
        sum_prob_outer += prob_record_outer;

        M_vals.push_back(m);
        expected_prev.push_back(e_prev);
        expected_next.push_back(e_next);
        probs_seen.push_back(prob_seen);
        probs_record.push_back(p_record);

        if (p_record > max_p_record) {
            max_p_record = p_record;
            printf("RECORD :%-6ld line %-5ld  unknowns: %4ld, %4ld "
                   "| e: %.1f, %.1f\t| prob record: %.2e (%.2e + %.2e)\t| %.7f\n",
                m, M_vals.size(),
                unknown_low.size(), unknown_high.size(),
                e_prev, e_next,
                p_record, prob_record, prob_record_outer,
                prob_seen);
        }

        if (MAX_MISSING_GAPS > 0 && prob_is_missing_gap > 0.0) {
            if (prob_is_missing_gap > max_m_record) {
                max_m_record = prob_is_missing_gap;
                printf("MISSING:%-6ld line %-5ld  unknowns: %4ld, %4ld "
                       "|\t\t\t| prob record: %.2e  missing: %.2e\t| %.7f\n",
                    m, M_vals.size(),
                    unknown_low.size(), unknown_high.size(),
                    p_record, prob_is_missing_gap, prob_seen);
            }

            // priority_queue has the smallest (least likely) record on top.
            if ((missing_gaps_search.size() < MAX_MISSING_GAPS) ||
                    (std::get<0>(missing_gaps_search.top()) < prob_is_missing_gap)) {
                missing_gaps_search.push({prob_is_missing_gap, m, missing_pairs});

                if (missing_gaps_search.size() > MAX_MISSING_GAPS) {
                    missing_gaps_search.pop();
                }
            }
        }
    }

    {
        long  s_tests = probs_seen.size();
        auto s_stop_t = high_resolution_clock::now();
        double   secs = duration<double>(s_stop_t - s_start_t).count();
        // Stats!
        printf("\t%ld tests %.2f seconds (%.2f/sec)\n",
            s_tests, secs, s_tests / secs);
    }
    {
        cout << endl;
        printf("avg seen prob:   %.7f\n",
            std::accumulate(probs_seen.begin(), probs_seen.end(), 0.0) / probs_seen.size());
        printf("avg record prob: %.2e (max: %.3e)\n",
            std::accumulate(probs_record.begin(), probs_record.end(), 0.0) / probs_record.size(),
            *std::max_element(probs_record.begin(), probs_record.end()));
        printf("prob record inner: %.5f   prob record extended: %.5f\n",
            sum_prob_inner, sum_prob_outer);

    }
}


void prime_gap_stats(const struct Config config) {
    const uint64_t M_start = config.mstart;
    const uint64_t M_inc = config.minc;
    const uint64_t P = config.p;
    const uint64_t D = config.d;

    const unsigned int SIEVE_LENGTH = config.sieve_length;
    const unsigned int SL = SIEVE_LENGTH;
    assert( SL > 0 );


    // ----- Read from unknown file
    std::ifstream unknown_file;
    {
        std::string fn = gen_unknown_fn(config, ".txt");
        printf("\tReading from %s'\n", fn.c_str());
        unknown_file.open(fn, std::ios::in);
        assert( unknown_file.is_open() ); // Can't open save_unknowns file
        assert( unknown_file.good() );    // Can't open save_unknowns file
    }

    // ----- Merit Stuff
    mpz_t K;
    mpz_init(K);
    mpz_primorial_ui(K, P);
    assert( 0 == mpz_tdiv_q_ui(K, K, D) );
    assert( mpz_cmp_ui(K, 1) > 0); // K <= 1 ?!?

    int K_digits;
    float K_log;
    K_stats(config, K, &K_digits, &K_log);

    // ----- Get Record Prime Gaps
    vector<float> records = get_record_gaps();

    // gap that would be a record with m*P#/d
    vector<uint32_t> poss_record_gaps;
    float smallest_log = K_log + log(M_start);
    {
        for (size_t g = 2; g < MAX_GAP; g += 2) {
            // Ignore the infintesimal odds of finding >merit 30 gap.
            if (g / smallest_log > 30) {
                break;
            }

            if (records[g] > smallest_log) {
                poss_record_gaps.push_back(g);
                if (poss_record_gaps.size() <= 2) {
                    printf("If found Gap: %ld (current: %.2f) would improve to %.3f\n",
                        g, g / records[g], g / smallest_log);
                }
            }
        }
        assert( poss_record_gaps.size() );
        printf("\tFound %ld possible record gaps (%d to %d)",
            poss_record_gaps.size(), poss_record_gaps.front(), poss_record_gaps.back());

        if (poss_record_gaps.front() > 2 * SIEVE_LENGTH) {
            printf("\tHard to determine record prob, 2 * sieve_length < min_record_gap");
        }
        cout << endl;

        assert(is_sorted(poss_record_gaps.begin(), poss_record_gaps.end()));
        assert(poss_record_gaps.front() <= MISSING_GAPS_LOW);
    }


    // ----- Generate primes for P
    vector<uint32_t> P_primes = get_sieve_primes(P);
    assert( P_primes.back() == P);

    // ----- Sieve stats
    const double PROB_PRIME = 1 / smallest_log;
    const double UNKNOWNS_AFTER_SIEVE = 1 / (log(config.sieve_range) * exp(GAMMA));
    const double PROB_PRIME_AFTER_SIEVE = PROB_PRIME / UNKNOWNS_AFTER_SIEVE;
    {
        cout << endl;
        printf("prob prime             : %.7f\n", PROB_PRIME);
        printf("prob prime after sieve : %.5f\n\n", PROB_PRIME_AFTER_SIEVE);
    }

    // Prob 0th, 1st, 2nd is prime
    vector<float> prob_prime_nth_sieve;
    // Prob < nth is not prime or conversely >= nth is prime.
    vector<float> prob_great_nth_sieve;
    prob_nth_prime(PROB_PRIME_AFTER_SIEVE,
        prob_prime_nth_sieve, prob_great_nth_sieve);


    size_t EXPECTED_UNKNOWN = UNKNOWNS_AFTER_SIEVE * SL;
    assert(EXPECTED_UNKNOWN < prob_prime_nth_sieve.size());
    /* prob_combined_sieve[i+j] = prime * (1 - prime)^i * (1 - prime)^j * prime */
    vector<float> prob_combined_sieve;
    prob_combined_gap(
        PROB_PRIME_AFTER_SIEVE,
        EXPECTED_UNKNOWN,
        prob_combined_sieve);


    // Prob record with gap[i] and other gap > SL
    vector<float> prob_record_extended_gap(SL+1, 0.0);
    {
        // Same as above but not considering sieving done on the numbers.
        vector<float> prob_prime_nth;
        vector<float> prob_great_nth;
        prob_nth_prime(PROB_PRIME, prob_prime_nth, prob_great_nth);
        //printf("number of tests for prob < 1e-16, with sieve: %ld, without: %ld\n",
        //    prob_prime_nth_sieve.size(), prob_prime_nth.size());

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

        // Count of numbers [SL, SL+j] < i coprime to D
        // only (K, SL+i) == 1 can be prime
        vector<bool>     is_coprime;
        vector<uint32_t> count_coprime;
        {
            uint32_t count = 0;
            for (int i = SL; ; i++) {
                bool coprime = mpz_gcd_ui(NULL, K, i) == 1;
                count += coprime;
                is_coprime.push_back(coprime);
                count_coprime.push_back(count);
                if (prob_great_nth[count] < PROB_GREATER_CUTOFF) {
                    break;
                }
            }
            printf("Considering SL + %ld (%d coprime) after record extended gap\n",
                count_coprime.size(), count);
            printf("\tUsing prob_greater_cutoff %.2e * %.2e\n",
                PROB_GREATER_CUTOFF, prob_great_nth_sieve[EXPECTED_UNKNOWN]);
            assert( count <= prob_prime_nth.size() );
        }

        auto  s_start_t = high_resolution_clock::now();

        // TODO: can this be used to turn the double for loop in run_gap_search into
        //       two single for loops?

        ///*
        for (size_t gap_one = 1; gap_one <= SL; gap_one++) {
            // only needed for values that can be coprime with K
            if (mpz_gcd_ui(NULL, K, gap_one) > 1) {
                prob_record_extended_gap[gap_one] = std::nan("");
                continue;
            }

            double prob_record = 0;
            for (uint32_t record_gap : poss_record_gaps ) {
                uint32_t dist = record_gap - gap_one;
                // TODO SL should be included in sieve_length in future.
                if (dist < SL) continue;

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


        auto s_stop_t = high_resolution_clock::now();
        double   secs = duration<double>(s_stop_t - s_start_t).count();
        printf("Prob Records considered (%.2f seconds)\n", secs);
    }
    mpz_clear(K);

    vector<uint64_t> M_vals;
    vector<float> expected_prev;
    vector<float> expected_next;
    vector<float> probs_seen;
    vector<float> probs_record;

    // <prob, m, [<low, high>, <low, high>, ...]>
    min_queue_gaps missing_gaps_search;

    // ----- Main calculation
    printf("\n");
    run_gap_file(
        config,
        records,
        poss_record_gaps,
        M_start, M_inc, D,
        K_log,
        prob_prime_nth_sieve, prob_great_nth_sieve,
        prob_combined_sieve,
        prob_record_extended_gap,
        unknown_file,
        M_vals,
        expected_prev, expected_next,
        probs_seen, probs_record,
        missing_gaps_search
    );

    if (!missing_gaps_search.empty()) {
        // pop in reverse order and reverse.
        vector<missing_gap_record> missing_sorted;
        while (!missing_gaps_search.empty()) {
            missing_sorted.push_back(missing_gaps_search.top());
            missing_gaps_search.pop();
        }
        std::reverse(missing_sorted.begin(), missing_sorted.end());

        float sum_missing_prob = 0;
        for (auto missing_search : missing_sorted) {
            sum_missing_prob += std::get<0>(missing_search);
        }
        printf("sum missing prob: %.6e (max: %.3e)\n",
            sum_missing_prob, std::get<0>(missing_sorted.front()));

        std::string missing_fn = gen_unknown_fn(config, ".missing.txt");

        std::ofstream missing_gap_file(missing_fn, std::ios::out);
        assert( missing_gap_file.is_open() ); // Can't open missing_gap file
        assert( missing_gap_file.good() );

        size_t i = 0;
        for (auto missing_search : missing_sorted) {
            i++;

            float prob = std::get<0>(missing_search);
            uint64_t m = std::get<1>(missing_search);
            if (i <= 10 || (i <= 100 && i % 10 == 0) ||
                    (i % 100 == 0) || (i == missing_sorted.size())) {
                printf("MISSING TESTS %4ld:%-6ld => %.2e | pairs: %4ld\n",
                    i, m, prob, std::get<2>(missing_search).size());
            }

            missing_gap_file << prob << " : " << m << "*" << config.p << "#/" << config.d << " :";
            for (auto pair : std::get<2>(missing_search)) {
                missing_gap_file << " (" << pair.first << "," << pair.second << ")";
            }
            missing_gap_file << endl;
        }
        missing_gap_file.close();
        printf("\tSaved %ld missing-gaps tests to '%s'\n",
            missing_sorted.size(), missing_fn.c_str());
    }

    if (config.save_unknowns) {
        store_stats(
            config, K_log,
            M_vals,
            expected_prev, expected_next,
            probs_seen, probs_record
        );
    }
}
