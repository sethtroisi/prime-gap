/* Copyright 2020 Seth Troisi

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

PRAGMA foreign_keys = ON;

/* Used in various stages of combined_sieve / gap_test */

/*
    Flow is:
        ./combined_sieve -> <params>.txt
        ./gap_stats
                INSERT row into 'range'
                INSERT range stats into 'range_stats'
                INSERT stats per m into 'm_stats'
                (if using -DSEARCH_MISSING_GAPS=1)
                       stats per m into 'm_missing_stats'

        [OPTIONAL]
                ./missing_gap_test.py
                        CHECK  'm_stats'
                        INSERT 'm_missing_stats'
                        UPDATE 'm_stats'
                            Set next_p/prev_p (with negative to indicate bounds)
                            Set counts on prp

        ./gap_test
                CHECK  'm_stats'
                UPDATE 'm_stats'
                INSERT into 'result'
*/

CREATE TABLE IF NOT EXISTS range(
        /* range id */
        rid      INTEGER PRIMARY KEY,

        m_start INTEGER,
        m_inc   INTEGER,
        P INTEGER,
        D INTEGER,

        sieve_length INTEGER,
        max_prime    INTEGER,

        /* for prob_merit in m_stats */
        min_merit DOUBLE,

        /* count of m with (m, D) = 1 */
        num_m INTEGER,

        /* number of entries in m_stats processed by gap_tests */
        num_processed INTEGER DEFAULT 0,

        /* number left to process in m_stats
         * because of --top-x-percent
         * num_m != num_processed + num_remaining
         */
        num_remaining INTEGER,

        /* Time for various tasks */
        time_sieve REAL DEFAULT 0,
        time_stats REAL DEFAULT 0,
        time_tests REAL DEFAULT 0
);

CREATE TABLE IF NOT EXISTS range_stats (
        /**
         * produced by gap_stats.py
         *
         * should be Expected values over ALL m
         */

        /* range id (foreign key) */
        rid INTEGER REFERENCES range(rid) ON UPDATE CASCADE ON DELETE CASCADE,

        /* map of <gap, prob> over all <low,high> pairs over all m's (in range) */
        gap INTEGER,
        prob_combined FLOAT,
        /* map of <gap, prob> over all <low> values over all m's (in range) */
        prob_low_side FLOAT,
        /* map of <gap, prob> over all <high> values over all m's (in range) */
        prob_high_side FLOAT,

        PRIMARY KEY (rid, gap)
);

CREATE TABLE IF NOT EXISTS m_stats (
        rid INTEGER REFERENCES range(rid) ON UPDATE CASCADE
                                          ON DELETE SET NULL,

        m INTEGER,
        P INTEGER,
        D INTEGER,

        /* next_p_interval, prev_p_interval (both positive) */
        /**
         * positive => distance to next/prev is X
         * 0        => ???
         * negative => X is prime but haven't checked all values less
         */
        next_p INTEGER,
        prev_p INTEGER,

        /* (next_p_i + prev_p_i) / log(N) */
        merit REAL,

        /* apriori probability of P(record gap), P(missing gap), P(merit > rid.min_merit) */
        prob_record REAL,
        prob_missing REAL,
        prob_merit  REAL,

        /* expected values useful for a number of printouts */
        e_gap_next REAL,
        e_gap_prev REAL,

        /* updated during gap_test / missing_gap_test */
        prp_next INTEGER,
        prp_prev INTEGER,

        test_time REAL,

        PRIMARY KEY(m, P, D)
);

CREATE TABLE IF NOT EXISTS result (
        m INTEGER,
        P INTEGER,
        D INTEGER,

        /* next_p_interval, prev_p_interval */
        next_p_i INTEGER,
        prev_p_i INTEGER,

        merit REAL,

        PRIMARY KEY(m, P, D)
);

CREATE INDEX IF NOT EXISTS    r_m_p_d ON             result(m,P,D);
CREATE INDEX IF NOT EXISTS   ms_m_p_d ON            m_stats(m,P,D);
