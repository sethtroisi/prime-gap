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

set -eu


sqlite3 prime-gap-search.db <<EOL
SELECT PRINTF("P=%-5d D=%-5d M(%-5d) = %6d +[0,%6d) time(sieve,stats)=%5.1f %5.1f",
              p, d, num_m, m_start, m_inc, time_sieve, time_stats)
FROM range
WHERE
    (num_processed = 0 AND time_tests <= 1 AND finalized = 0 AND
     0 = (SELECT COUNT(*) FROM m_stats m
          WHERE m.rid = range.rid AND (test_time > 0 OR next_p != 0 OR prev_p != 0 OR prp_next + prp_prev > 0)))
EOL

while true; do
    read -p "Drop (and cleanup) these ranges? [yN]: " yn
    case $yn in
        [Yy]* ) break;;
        * ) exit;;
    esac
done

echo "Dropping empty ranges";


sqlite3 prime-gap-search.db <<EOL
PRAGMA foreign_keys=1;
DELETE FROM range
WHERE
    (num_processed = 0 AND time_tests <= 1 AND finalized = 0 AND
     0 = (SELECT COUNT(*) FROM m_stats m
          WHERE m.rid = range.rid AND (test_time > 0 OR next_p != 0 OR prev_p != 0 OR prp_next + prp_prev > 0)));
SELECT "DELETED " || total_changes() || " range/range_stats/m_stats";
EOL
