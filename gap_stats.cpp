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
#include <fstream>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

#include <gmp.h>
#include <sqlite3.h>

#include "gap_common.h"

using std::cout;
using std::endl;
using std::pair;
using std::vector;
using namespace std::chrono;


const char records_db[] = "gaps.db";
const char gaps_db[]    = "prime-gaps.db";

// Limits the size of record list
const uint32_t MAX_GAP = 1'000'000;

// Generated from https://primegap-list-project.github.io/lists/missing-gaps/
const vector<uint32_t> MISSING_GAPS = {
    113326, 115694, 116254, 117238, 117242,
    119222, 119584, 120154, 121138, 121174,
    121366, 121832, 122290, 122666, 122686,
    123230, 123238, 123242, 123598, 123662,
    123782, 124106, 124258, 124346, 124534,
    124792, 125024, 125318, 125974, 126134,
    126206, 126236, 126298, 126376, 126394,
    126538, 126554, 126814, 127346, 127544,
    127622, 127732, 127906, 128114, 128362,
    128372, 128516, 128686, 128714, 128762,
    128872, 129298, 129406, 129538, 129698,
    129754, 129784, 130042, 130162, 130252,
    130280, 130282, 130310, 130438, 130798,
    130846, 130882, 130898, 131074, 131288,
    131378, 131402, 131446, 131530, 131536,
    131564, 131578, 131648, 131762, 131788,
    131876, 131938, 131954, 132130, 132194,
    132206, 132218, 132232, 132242, 132302,
    132314, 132446, 132506, 132548, 132598,
    132644, 132838, 132842, 132848, 132928
};

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

    // TODO move behind a flag
    assert( is_sorted(MISSING_GAPS.begin(), MISSING_GAPS.end()) );

    prime_gap_stats(config);
}


//---------------------------------------------------------------------------//


sqlite3* get_db(const char* path) {
    sqlite3 *db;
    if (sqlite3_open(path, &db) != SQLITE_OK) {
        printf("Can't open database(%s): %s\n", path, sqlite3_errmsg(db));
        exit(1);
    }
    return db;
}


vector<float> get_record_gaps() {
    // TODO technically this should be INF.
    vector<float> records(MAX_GAP, 0.0);

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
    printf("K = %d bits, %d digits, log(K) = %.2f\n",
        K_bits, *K_digits, *K_log);
    printf("Min Gap ~= %d (for merit > %.1f)\n\n",
        (int) (config.minmerit * (*K_log + m_log)), config.minmerit);
}


void prob_nth_prime(
        double prob_prime,
        vector<float>& prob_prime_nth,
        vector<float>& prob_great_nth) {
    prob_prime_nth.push_back(0.0f);
    prob_great_nth.push_back(1.0f);

    double prob_still_prime = 1.0;
    for (int i = 1; prob_still_prime > 1e-16; i++) {
        prob_prime_nth.push_back(prob_still_prime * prob_prime);
        prob_still_prime *= 1 - prob_prime;
        prob_great_nth.push_back(prob_still_prime);
    }
}

void run_gap_file(
        const struct Config config,
        const vector<float>& records,
        const vector<uint32_t> min_record_gaps,
        const uint64_t M_start,
        const uint64_t M_inc,
        const uint64_t D,
        float K_log,
        const vector<float>& prob_prime_nth,
        const vector<float>& prob_great_nth,
        const vector<float>& prob_record_extended_gap,
        std::ifstream& unknown_file,
        vector<uint64_t>& M_vals,
        vector<float>& expected_prev,
        vector<float>& expected_next,
        vector<float>& probs_seen,
        vector<float>& probs_record) {

    std::ofstream missing_gap_file;
    {
        std::string missing_fn = gen_unknown_fn(config, ".missing.txt");
        printf("\tSaving missing-gaps to '%s'\n", missing_fn.c_str());
        missing_gap_file.open(missing_fn, std::ios::out);
        assert( missing_gap_file.is_open() ); // Can't open missing_gap file
        assert( missing_gap_file.good() );
    }

    const uint32_t min_record_gap = min_record_gaps.front();
    float max_p_record = 0;
    for (uint64_t mi = 0; mi < M_inc; mi++) {
        uint64_t m = M_start + mi;
        if (gcd(m, D) != 1) {
            continue;
        }
        // Reset sieve array to unknown.
        vector<uint32_t> unknown_high, unknown_low;

        int unknown_l = 0;
        int unknown_u = 0;

        // Read a line from the file
        {
            int mtest;
            unknown_file >> mtest;
            assert( mtest >= 0 );
            if ((unsigned) mtest != mi ) {
                cout << "Mismatched mi " << mtest << " vs " << mi << endl;
            }
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

        // Note slightly different from smallest_log
        float log_merit = K_log + log(m);

        double prob_record = 0;
        double prob_seen = 0;

        // When prime is > SL (prob ~1e-4) how to figure out prob_record?
        //      for i, determine avg gap between

        prob_seen = (1 - prob_great_nth[unknown_high.size()]) *
                    (1 - prob_great_nth[unknown_low.size()]);

        float prob_is_missing_gap = 0.0;
        vector<pair<uint32_t, uint32_t>> missing_pairs;

        /**
         * TODO could possible look only at record_gaps
         * replace 'records[gap] > log_merit' with 'gap == records_gap[gi]' or something
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

                // TODO use something based on merit
                if (gap >= MISSING_GAPS.front() &&
                        gap <= MISSING_GAPS.back() &&
                        records[gap] == 0.0) {
                    float prob_joint = prob_prime_nth[i+1] * prob_prime_nth[j+1];
                    prob_is_missing_gap += prob_joint;
                    missing_pairs.emplace_back(std::make_pair(gap_low, gap_high));
                }

                if (records[gap] > log_merit) {
                    float prob_joint = prob_prime_nth[i+1] * prob_prime_nth[j+1];
                    prob_record += prob_joint;
                }
            }
        }

        double e_prev = 0, e_next = 0;
        double prob_record_estimate = 0;
        for (size_t i = 0; i < std::max(unknown_low.size(), unknown_high.size()); i++) {
            // i'th is prime, other size is >= SL
            if (i < unknown_low.size()) {
                float prob_i = prob_prime_nth[i+1];
                float prob_j_greater = prob_great_nth[unknown_high.size()];
                float conditional_prob = prob_record_extended_gap[unknown_low[i]];
                prob_record_estimate += prob_i * prob_j_greater * conditional_prob;
                e_prev += unknown_low[i] * prob_i;
            }
            if (i < unknown_high.size()) {
                float prob_i = prob_prime_nth[i+1];
                float prob_j_greater = prob_great_nth[unknown_low.size()];
                float conditional_prob = prob_record_extended_gap[unknown_high[i]];
                prob_record_estimate += prob_i * prob_j_greater * conditional_prob;
                e_next += unknown_high[i] * prob_i;
            }
        }

        M_vals.push_back(m);
        expected_prev.push_back(e_prev);
        expected_next.push_back(e_next);
        probs_seen.push_back(prob_seen);
        probs_record.push_back(prob_record);

        double p_record = prob_record + prob_record_estimate;
        if (p_record > max_p_record) {
            max_p_record = p_record;
            printf("M:%-6ld (line %ld)\tunknowns: %4ld, %4ld | e: %.1f, %.1f | prob record: %.2e (%.2e + %.2e) | %.7f\n",
                m, M_vals.size(),
                unknown_low.size(), unknown_high.size(),
                e_prev, e_next,
                p_record, prob_record, prob_record_estimate,
                prob_seen);
        }

        if (prob_is_missing_gap > 1.2e-4) {
            printf("MISSING TESTS:%-6ld => %.2e | unknowns: %4ld, %4ld | missing tests: %4ld | \n",
                m, prob_is_missing_gap,
                unknown_low.size(), unknown_high.size(),
                missing_pairs.size());

            // TODO Multiple same lower: (5, 102), (5, 110), (5, 200), (7, 104), ...
            missing_gap_file << prob_is_missing_gap << " : " << m << "*" << config.p << "#/" << config.d << " :";
            for (auto pair : missing_pairs) {
                missing_gap_file << " (" << pair.first << "," << pair.second << ")";
            }
            missing_gap_file << endl;
        }
    }

    missing_gap_file.close();
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
    // Smallest gap that would be a record with m*P#/d
    vector<uint32_t> min_record_gaps;
    float smallest_log = K_log + log(M_start);
    {
        for (size_t g = 2; g < MAX_GAP; g += 2) {
            if (records[g] > smallest_log) {
                min_record_gaps.push_back(g);
                if (min_record_gaps.size() <= 5) {
                    printf("If found Gap: %ld (current: %.2f) would improve to %.3f\n",
                        g, g / records[g], g / smallest_log);
                }
            }
        }
        printf("\tFound %ld possible record gaps", min_record_gaps.size());
        assert( min_record_gaps.size() );
        if (min_record_gaps.front() > 2 * SIEVE_LENGTH) {
            printf("\tCan't determine record prob, 2 * sieve_length < min_record_gap");
        }
        cout << endl;

        assert( is_sorted(min_record_gaps.begin(), min_record_gaps.end()) );
    }


    // ----- Generate primes for P
    vector<uint32_t> P_primes = get_sieve_primes(P);
    assert( P_primes.back() == P);

    // ----- Sieve stats
    double prob_prime = 1 / smallest_log;
    double prob_prime_after_sieve;
    {
        double unknowns_after_sieve = 1 / (log(config.sieve_range) * exp(GAMMA));
        prob_prime_after_sieve = prob_prime / unknowns_after_sieve;
        printf("prob prime            : %.7f\n", prob_prime);
        printf("prob prime after sieve: %.5f\n\n", prob_prime_after_sieve);
    }

    // Prob 0th, 1st, 2nd is prime
    vector<float> prob_prime_nth_sieve;
    // Prob < nth is not prime or conversely >= nth is prime.
    vector<float> prob_great_nth_sieve;
    prob_nth_prime(prob_prime_after_sieve,
        prob_prime_nth_sieve, prob_great_nth_sieve);

    // Prob record with gap[i] and other gap > SL
    vector<float> prob_record_extended_gap(SL+1, 0.0);
    {
        // Same as above but not considering sieving done on the numbers.
        vector<float> prob_prime_nth;
        vector<float> prob_great_nth;
        prob_nth_prime(prob_prime, prob_prime_nth, prob_great_nth);
        printf("number of tests for prob < 1e-16, with sieve: %ld, without: %ld\n",
            prob_prime_nth_sieve.size(), prob_prime_nth.size());

        // Count of numbers [SL, SL+j] < i coprime to D
        // only (D, SL+i) == 1 can be prime
        vector<bool>     is_coprime;
        vector<uint32_t> count_coprime;
        {
            uint32_t count = 0;
            for (int i = SL; ; i++) {
                bool coprime = mpz_gcd_ui(NULL, K, i) == 1;
                count += coprime;
                is_coprime.push_back(coprime);
                count_coprime.push_back(count);
                if (prob_great_nth[count] < 1e-12) {
                    break;
                }
            }
            printf("Considering %ld (%d coprime) after sieve for record_gap purpose\n",
                count_coprime.size(), count);
            assert( count <= prob_prime_nth.size() );
        }

        for (size_t gap_one = 1; gap_one <= SL; gap_one++) {
            double prob_record = 0;
            size_t count_records = 0;

            for (uint32_t record_gap : min_record_gaps ) {
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
                count_records += 1;
            }
            prob_record_extended_gap[gap_one] = prob_record;
        }
    }

    // Used for various stats
    auto  s_start_t = high_resolution_clock::now();
    //long  s_total_unknown = 0;

    vector<uint64_t> M_vals;
    vector<float> expected_prev;
    vector<float> expected_next;
    vector<float> probs_seen;
    vector<float> probs_record;

    // ----- Main calculation
    printf("\n");
    run_gap_file(
        config,
        records,
        min_record_gaps,
        M_start, M_inc, D,
        K_log,
        prob_prime_nth_sieve, prob_great_nth_sieve,
        prob_record_extended_gap,
        unknown_file,
        M_vals,
        expected_prev, expected_next,
        probs_seen, probs_record
    );

    {
        long  s_tests = probs_seen.size();
        auto s_stop_t = high_resolution_clock::now();
        double   secs = duration<double>(s_stop_t - s_start_t).count();
        // Stats!
        printf("\t%ld tests %.2f seconds (%.2f/sec)\n",
            s_tests, secs, s_tests / secs);
    }
    {
        printf("\n");
        printf("avg seen prob:   %.7f\n",
            std::accumulate(probs_seen.begin(), probs_seen.end(), 0.0) / probs_seen.size());
        printf("avg record prob: %.3e (max: %.2e)\n",
            std::accumulate(probs_record.begin(), probs_record.end(), 0.0) / probs_record.size(),
            *std::max_element(probs_record.begin(), probs_record.end()));
    }

    mpz_clear(K);

    if (config.save_unknowns) {
        store_stats(
            config, K_log,
            M_vals,
            expected_prev, expected_next,
            probs_seen, probs_record
        );
    }
}
