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


vector<float> get_record_gaps() {
    vector<float> records(MAX_GAP, 0.0);

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
        vector<float>& probs_seen,
        vector<float>& probs_record) {

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

        //float prob_record_upper_bound = 0;

        // TODO if 2*SL < min_record_gap just skip this.
        for (size_t i = 0; i < unknown_low.size(); i++) {
            float prob_i = prob_prime_nth[i+1];
            if (unknown_low[i] + unknown_high.back() < min_record_gap) {
                float prob_any_j = 1 - prob_great_nth[unknown_high.size()];
                prob_seen += prob_i * prob_any_j;
                continue;
            }

            for (size_t j = 0; j < unknown_high.size(); j++) {
                float prob_joint = prob_i * prob_prime_nth[j+1];
                prob_seen += prob_joint;


                uint32_t gap = abs(unknown_low[i]) + abs(unknown_high[j]);
                if ((records[gap] - 0.01) > log_merit) {
//                    printf("\tgap: %d curr: %5.2f, would be %5.3f (%.6f%%)\n",
//                        gap, gap / records[gap], gap / log_merit, 100.0 * prob_joint);
                    prob_record += prob_joint;
                }
            }
        }

        double prob_record_estimate = prob_record;
        for (size_t i = 0; i < std::max(unknown_low.size(), unknown_high.size()); i++) {
            // i'th is prime, other size is >= SL
            if (i < unknown_low.size()) {
                float prob_i = prob_prime_nth[i+1];
                float prob_j_greater = prob_great_nth[unknown_high.size()];
                float conditional_prob = prob_record_extended_gap[unknown_low[i]];
                prob_record_estimate += prob_i * prob_j_greater * conditional_prob;
            }
            if (i < unknown_high.size()) {
                float prob_i = prob_prime_nth[i+1];
                float prob_j_greater = prob_great_nth[unknown_low.size()];
                float conditional_prob = prob_record_extended_gap[unknown_high[i]];
                prob_record_estimate += prob_i * prob_j_greater * conditional_prob;
            }
        }

        double p_record = prob_record + prob_record_estimate;
        if (p_record > max_p_record) {
            max_p_record = p_record;
            printf("%5ld unknowns: %3ld, %3ld | prob record: %.2e (%.2e + %.2e)\n",
                m, unknown_low.size(), unknown_high.size(),
                p_record, prob_record, prob_record_estimate);
        }

        probs_seen.push_back(prob_seen);
        probs_record.push_back(prob_record);
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
    gmp_printf("%Zd\n", K);

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
        assert( min_record_gaps.size() );
        if (min_record_gaps.front() > 2 * SIEVE_LENGTH) {
            printf("\tCan't determine record prob, 2 * sieve_length < min_record_gap");
        }
        cout << endl;
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
        printf("tests for prob < 1e-16, with sieve: %ld, without: %ld\n",
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

    vector<float> probs_seen;
    vector<float> probs_record;
//    vector<float> expected_prev;
//    vector<float> expected_next;
//    vector<float> expected_gap;

    // ----- Main calculation
    printf("\n");
    run_gap_file(
        records,
        min_record_gaps,
        M_start, M_inc, D,
        K_log,
        prob_prime_nth_sieve, prob_great_nth_sieve,
        prob_record_extended_gap,
        unknown_file,
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
        printf("avg seen prob:   %.6f\n",
            std::accumulate(probs_seen.begin(), probs_seen.end(), 0.0) / probs_seen.size());
        printf("avg record prob: %.3e (max: %.2e)\n",
            std::accumulate(probs_record.begin(), probs_record.end(), 0.0) / probs_record.size(),
            *std::max_element(probs_record.begin(), probs_record.end()));
    }

    mpz_clear(K);
}
