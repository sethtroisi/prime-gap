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
using std::vector;
using namespace std::chrono;

// Limits the size of record list
const uint32_t MAX_GAP = 1'000'000;

void prime_gap_stats(const struct Config config);

int main(int argc, char* argv[]) {
    printf("\tCompiled with GMP %d.%d.%d\n\n",
        __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);

    Config config = argparse(argc, argv);

    if (config.save_unknowns == 0) {
        cout << "Must set --save-unknowns" << endl;
        return 1;
    }
    if (config.run_prp == 1) {
        cout << "Must set --sieve-only for gap_search" << endl;
        return 1;
    }

    if (config.valid == 0) {
        show_usage(argv[0]);
        return 1;
    }

    prime_gap_stats(config);
}


vector<double> get_record_gaps() {
    vector<double> records(MAX_GAP, 0.0);

    // TODO accept db file as param.

    sqlite3 *db;
    int rc = sqlite3_open("gaps.db", &db);
    if( rc ) {
        printf("Can't open database: %s\n", sqlite3_errmsg(db));
        exit(1);
    }

    /* Create SQL statement */
    char sql[] = "SELECT gapsize, merit FROM gaps";
    char *zErrMsg = 0;

    /* Execute SQL statement */
    rc = sqlite3_exec(db, sql, [](void* recs, int argc, char **argv, char **azColName)->int {
        uint64_t gap = atol(argv[0]);
        if (gap < MAX_GAP) {
            // Recover log(startprime)
            (*static_cast<vector<double>*>(recs))[gap] = gap / atof(argv[1]);
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


std::string gen_unknown_fn(const struct Config& config, std::string suffix) {
    return std::to_string(config.mstart) + "_" +
           std::to_string(config.p) + "_" +
           std::to_string(config.d) + "_" +
           std::to_string(config.minc) + "_s" +
           std::to_string(config.sieve_length) + "_l" +
           std::to_string(config.sieve_range / 1'000'000) + "M" +
           (config.method2 ? ".m2" : "") +
           suffix;
}


void K_stats(
        const struct Config& config,
        mpz_t &K, int *K_digits, double *K_log) {
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


void prime_gap_stats(const struct Config config) {
    const uint64_t M_start = config.mstart;
    const uint64_t M_inc = config.minc;
    const uint64_t P = config.p;
    const uint64_t D = config.d;

    const unsigned int SIEVE_LENGTH = config.sieve_length;
    const unsigned int SL = SIEVE_LENGTH;
    assert( SL > 0 );


    // ----- Save Output file
    std::ifstream unknown_file;
    std::string fn = gen_unknown_fn(config, ".txt");
    printf("\tReading from '%s'\n", fn.c_str());
    unknown_file.open(fn, std::ios::in);
    assert( unknown_file.is_open() ); // Can't open save_unknowns file
    assert( unknown_file.good() );    // Can't open save_unknowns file

    // ----- Merit Stuff
    mpz_t K;
    mpz_init(K);
    mpz_primorial_ui(K, P);
    assert( 0 == mpz_tdiv_q_ui(K, K, D) );
    assert( mpz_cmp_ui(K, 1) > 0); // K <= 1 ?!?

    int K_digits;
    double K_log;
    K_stats(config, K, &K_digits, &K_log);

    // ----- Get Record Prime Gaps
    vector<double> records = get_record_gaps();
    // Smallest gap that would be a record with m*P#/d
    uint32_t min_record_gap = 0;
    double smallest_log = K_log + log(M_start);
    {
        for (size_t g = 2; min_record_gap == 0 && g < MAX_GAP; g += 2) {
            if (records[g] > smallest_log) {
                min_record_gap = g;
            }
        }
        assert( min_record_gap > 0 );
        printf("Min Record Gap: %d merit %.2f would improve to %.3f\n",
            min_record_gap,
            min_record_gap / records[min_record_gap],
            min_record_gap / smallest_log);

        if (min_record_gap > 2 * SIEVE_LENGTH) {
            printf("\tCan't determine record prob, 2 * sieve_length < min_record_gap");
        }
        cout << endl;
    }


    // ----- Generate primes for P
    vector<uint32_t> P_primes = get_sieve_primes(P);
    assert( P_primes.back() == P);

    // ----- Sieve stats
    double prob_prime_after_sieve;
    {
        double prob_prime = 1 / smallest_log;
        double unknowns_after_sieve = 1 / (log(config.sieve_range) * exp(GAMMA));
        prob_prime_after_sieve = prob_prime / unknowns_after_sieve;
        cout << "prob_prime_after_sieve: " << prob_prime_after_sieve << endl;
    }

    // Prob 0th, 1st, 2nd is prime
    vector<float> prob_prime_nth = {0};
    // Prob nth or greater.
    vector<float> prob_great_nth = {1.0};
    {
        double prob_still_prime = 1.0;
        for (int i = 1; i < 1000; i++) {
            // large enough idea.
            prob_prime_nth.push_back(prob_still_prime * prob_prime_after_sieve);
            prob_great_nth.push_back(prob_still_prime);

            prob_still_prime *= 1 - prob_prime_after_sieve;
        }
    }

    // Used for various stats
    auto  s_start_t = high_resolution_clock::now();
    //long  s_total_unknown = 0;
    long  s_tests = 0;

    vector<double> probs_seen;
    vector<double> probs_record;
//    vector<float> expected_prev;
//    vector<float> expected_next;
//    vector<float> expected_gap;

    for (uint64_t mi = 0; mi < M_inc; mi++) {
        uint64_t m = M_start + mi;
        if (gcd(m, D) != 1) {
            continue;
        }
        s_tests += 1;

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
        double log_merit = K_log + log(m);

        double prob_record = 0;
        double prob_seen = 0;
        // TODO if 2*SL < min_record_gap just skip this.
        for (size_t i = 0; i < unknown_low.size(); i++) {
            for (size_t j = 0; j < unknown_high.size(); j++) {
                double prob_joint = prob_prime_nth[i+1] * prob_prime_nth[j+1];
                prob_seen += prob_joint;

                uint32_t gap = abs(unknown_low[i]) + abs(unknown_high[j]);
                if ((records[gap] - 0.01) > log_merit) {
//                    printf("\tgap: %d curr: %5.2f, would be %5.3f (%.6f%%)\n",
//                        gap, gap / records[gap], gap / log_merit, 100.0 * prob_joint);
                    prob_record += prob_joint;
                }
            }
        }
        // Do something with i > size(), j > size();
        probs_seen.push_back(prob_seen);
        probs_record.push_back(prob_record);

        if (prob_record > 0) {
//            printf("%5ld (%ld) | %.5f %.6f\n", m, s_tests, prob_seen, prob_record);
        }
    }


    {
        auto s_stop_t = high_resolution_clock::now();
        double   secs = duration<double>(s_stop_t - s_start_t).count();
        // Stats!
        printf("\t%ld tests %.2f seconds (%.2f/sec)\n",
            s_tests, secs, s_tests / secs);
    }

    printf("avg seen prob:   %.4f\n",
        std::accumulate(probs_seen.begin(), probs_seen.end(), 0.0) / probs_seen.size());
    printf("avg record prob: %.6f (max: %.6f)\n",
        std::accumulate(probs_record.begin(), probs_record.end(), 0.0) / probs_record.size(),
        *std::max_element(probs_record.begin(), probs_record.end()));

    mpz_clear(K);
}
