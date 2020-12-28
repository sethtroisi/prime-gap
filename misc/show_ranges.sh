#!/bin/bash
#
# Copyright 2020 Seth Troisi
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

set -eux

# TODO add an -u update flag
sqlite3 "prime-gap-search.db" <<EOL
.timeout 20
#UPDATE range AS ra SET
#    num_processed =
#        (SELECT count(*) FROM result r
#         WHERE r.p = ra.p AND r.d = ra.d AND
#               r.m BETWEEN ra.m_start AND ra.m_start + ra.m_inc - 1 ),
#    time_tests =
#        (SELECT COALESCE(SUM(test_time), 0.0) FROM m_stats m
#         WHERE m.p = ra.p AND m.d = ra.d AND
#               m.m BETWEEN ra.m_start AND ra.m_start + ra.m_inc - 1 );
SELECT changes() || " Range Updates";

UPDATE range SET num_remaining = num_m - num_processed;
SELECT changes() || " Count Updates";
SELECT "";

#UPDATE m_stats AS m
#SET rid = NULL
#WHERE m.rid not in (SELECT rid from range);
#SELECT changes() || " m_stats.rid Updates";
#SELECT "";

# Happens from missing_gap, not sure quite what to do.
#INSERT INTO result (P,D,m, next_p_i,prev_p_i, merit)
#  SELECT P,D,m, next_p,prev_p, merit
#  FROM m_stats WHERE
#        m NOT IN (SELECT r.m FROM result r WHERE r.p = m_stats.p AND r.d = m_stats.d)
#        AND (next_p != 0 OR prev_p != 0);
#SELECT changes() || " m_stats results added to result";
#SELECT "";

# Happens when gap_stats is rerun
#UPDATE m_stats SET (next_p, prev_p) =
#       (SELECT next_p_i, prev_p_i FROM result r WHERE m_stats.p = r.p AND m_stats.d = r.d AND m_stats.m = r.m)
#WHERE (next_p = 0 OR prev_p = 0) AND p = 3001 AND d = 2190;
EOL

sqlite3 -readonly "prime-gap-search.db" <<EOL
.timeout 20
SELECT "Table 'range'";
SELECT PRINTF("  P=%-5d D=%-9d M(%-7d) = %9d + [0, %9d] ",
              p, d, num_m, m_start, m_inc),
       PRINTF(" time(sieve(%4dB),stats,tests)= %7.1f %5.1f %9.1f (%5.1f%%) | %.4f %7.3f | %5.4f / %d",
              max_prime / 1000000000,
              time_sieve, time_stats, time_tests, 100 * time_tests / (time_sieve + time_stats + time_tests),
              time_sieve / num_m, time_tests / num_processed,
              (time_sieve + time_stats + time_tests) / num_processed,
              num_processed)
FROM range ORDER BY p, d, m_start, m_inc;
SELECT "";

SELECT "Table 'results'/'m_stats'";
SELECT PRINTF("  P=%-5d D=%-9d M(%-7d) = %-9d to %-9d | last result: %9d, primes: %7d, PRPs: %9d, prob: %6.1f, merit: %5.2f, time: %6.1f hours",
              p, d, count(*), min(m), max(m), max(m * (primes != 0)),
              sum(primes), sum(prp_tests), sum(prob_record), max(merit), sum(test_time) / 3600)
FROM (
  SELECT p,d,m,
         (next_p!=0) + (prev_p!=0) as primes,
         prp_prev + prp_next as prp_tests,
         prob_record,
         merit, test_time
  FROM m_stats
)
GROUP BY p, d ORDER BY p, d;
SELECT"";
EOL

sqlite3 -readonly "prime-gap-search.db" <<EOL
.timeout 20
SELECT "Table 'result': " || COUNT(*) FROM result;
SELECT "Table 'm_stats': " || COUNT(*) FROM m_stats;
SELECT "Table 'range_stats': " || COUNT(*) FROM range_stats;
#SELECT P,count(*) FROM m_stats group by P;
EOL
