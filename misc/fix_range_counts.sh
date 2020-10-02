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
        (SELECT SUM(test_time) FROM m_stats m
         WHERE m.p = ra.p AND m.d = ra.d AND
               m.m BETWEEN ra.m_start AND ra.m_start + ra.m_inc - 1 );
SELECT total_changes() || " Timing Updates";

UPDATE range SET num_remaining = num_m - num_processed;
SELECT total_changes() || " Count Updates";

SELECT PRINTF("P=%-5d D=%-5d M(%-5d) = %6d +[0,%6d) processed=%-5d ",
              P, D, num_m, m_start, m_inc, num_processed),
       PRINTF(" time(sieve,stats,tests,p)=%5.1f %5.1f %8.1f | per m: %.2f",
              time_sieve, time_stats, time_tests, (time_sieve + time_stats + time_tests) / num_processed)
FROM range;
EOL
