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

# sqlite3 prime-gap-search.db <<EOL
# SELECT rid, r.p, r.d, num_m, num_processed, num_remaining, count(*)
# FROM range r
# INNER JOIN result ON
#         r.p = result.p AND r.d = result.d AND
#         result.m BETWEEN m_start AND m_start + m_inc
# GROUP BY rid;
# EOL

sqlite3 prime-gap-search.db <<EOL
UPDATE range AS ra SET
    num_processed =
        (SELECT count(*) FROM result r
         WHERE r.p = ra.p AND r.d = ra.d AND
               r.m BETWEEN ra.m_start AND ra.m_start + ra.m_inc - 1 ),
    time_tests =
        (SELECT COALESCE(SUM(test_time), 0.0) FROM m_stats m
         WHERE m.p = ra.p AND m.d = ra.d AND
               m.m BETWEEN ra.m_start AND ra.m_start + ra.m_inc - 1 );
SELECT changes() || " Range Updates";

UPDATE range SET num_remaining = num_m - num_processed;
SELECT changes() || " Count Updates";
SELECT "";

UPDATE m_stats AS m
SET rid = NULL
WHERE m.rid not in (SELECT rid from range);
SELECT changes() || " m_stats.rid Updates";
SELECT "";

SELECT "Table 'range'";
SELECT PRINTF("  P=%-5d D=%-5d M(%-5d) = %6d +[0,%6d) processed=%-5d ",
              p, d, num_m, m_start, m_inc, num_processed),
       PRINTF(" time(sieve,stats,tests,p)=%5.1f %5.1f %8.1f | per m: %5.2f   %d",
              time_sieve, time_stats, time_tests,
              (time_sieve + time_stats + time_tests) / num_processed, rid)
FROM range ORDER BY p, d, m_start, m_inc;
SELECT "";

SELECT "Table 'results'/'m_stats'";
SELECT PRINTF("  P=%-5d D=%-5d M(%-5d) = %-6d to %-6d (primes: %5d, PRPs: %7d, merit: %5.2f, time: %.0fs)",
              p, d, count(*), min(m), max(m), sum(primes), sum(prp_tests), max(merit), sum(test_time))
FROM (
  SELECT m,p,d, (next_p!=0)+(prev_p!=0) as primes, prp_prev+prp_next as prp_tests, merit, test_time FROM m_stats
    UNION
  SELECT m,p,d, 2 as primes, 2 as prp_tests, merit, 0 as test_time FROM result
)
GROUP BY p, d ORDER BY p, d;
SELECT"";

SELECT "Table 'result': " || COUNT(*) FROM result;
SELECT "Table 'm_stats': " || COUNT(*) FROM m_stats;
SELECT "Table 'range_stats': " || COUNT(*) FROM range_stats;

EOL
