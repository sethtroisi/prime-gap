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

make combined_sieve gap_stats RLE=1

#for m in 293609 811207 183047; do
#    for sl in 15 20 25; do
#        time ./gap_stats --unknown-filename unknowns/1511_312270_${m}_1_s${sl}000_l10000M.txt \
#                | tee data/${m}_${sl}_test2.txt
#        echo -e "\n\n"
#    done
#done


P=1511
D=312270

P=1667
D=210

for sl in {1500..45000..1500}; do
    FN="${P}_${D}_1_1000000_s${sl}_l10000M.txt"
    DB="data/test_${P}_${sl}.db"
#    rm -f "$DB";
#    sqlite3 "$DB" < schema.sql;
#    time ./combined_sieve --save-unknowns --unknown-filename $FN
#    time ./gap_stats --save-unknowns --search-db "$DB" --unknown-filename $FN
#    time ./gap_stats --unknown-filename $FN -q -q

    sqlite3 "$DB" <<EOL
        SELECT sieve_length,COUNT(*),SUM(prob_record),mean,AVG((prob_record-mean)*(prob_record-mean)) as variance FROM m_stats,
                (SELECT AVG(prob_record) as mean FROM m_stats),
                (SELECT sieve_length FROM range)
EOL

done
