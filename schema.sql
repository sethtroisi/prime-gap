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

/* Gap and Credit tables to match prime-gap-list */
CREATE TABLE IF NOT EXISTS gaps(
        gapsize INTEGER,
        ismax BOOLEAN,
        primecat TEXT,
        isfirst TEXT,
        primecert TEXT,
        discoverer TEXT,
        year INTEGER,
        merit REAL,
        primedigits INTEGER,
        startprime BLOB
);

CREATE TABLE IF NOT EXISTS credits (
        abbreviation TEXT,
        name TEXT,
        ack TEXT,
        display TEXT
);

/* Used in various stages of gap_search / gap_test */

/*
    Flow is:
        Add to `to_process_range`
            find high / low
        Add to `partial_result`
            find other
        Add to `result`
*/

CREATE TABLE IF NOT EXISTS to_process_range(
        id      INTEGER PRIMARY KEY,

        m_start INTEGER,
        m_int   INTEGER,
        P INTEGER,
        D INTEGER,

        sieve_length INTEGER,
        sieve_range INTEGER,

        /* rough integer of things left to process */
        num_to_processed INTEGER
);

CREATE TABLE IF NOT EXISTS partial_result (
        prid     INTEGER PRIMARY KEY,

        m INTEGER,
        P INTEGER,
        D INTEGER,
        next_p_i INTEGER,
        prev_p_i INTEGER,

        next_expected INTEGER,
        prev_expected INTEGER,

        expected_merit REAL
);

CREATE TABLE IF NOT EXISTS missing_gap_result (
        mgrid     INTEGER PRIMARY KEY,

        m INTEGER,
        P INTEGER,
        D INTEGER,

        /*
          Some fractions of unknown numbers were tested.
          found_{next,prev}_p were the first prime found.
          0 means none of next_p_test tests were prime.
          See missing_gap_test.py
        */
        found_prev_p_i INTEGER,
        found_next_p_i INTEGER,

        prev_p_tests INTEGER,
        next_p_tests INTEGER
);

CREATE TABLE IF NOT EXISTS result (
        rid     INTEGER PRIMARY KEY,

        m INTEGER,
        P INTEGER,
        D INTEGER,
        next_p_i INTEGER,
        prev_p_i INTEGER,
        merit REAL
);

CREATE INDEX IF NOT EXISTS  pres_m_p_d ON     partial_result(m,P,D);
CREATE INDEX IF NOT EXISTS mgapr_m_p_d ON missing_gap_result(m,P,D);
CREATE INDEX IF NOT EXISTS   res_m_p_d ON             result(m,P,D);
