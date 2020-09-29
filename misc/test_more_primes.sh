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

P=1511
D=2190
MINC=20000
SL=20000


echo "Cleaning up Sqlite"
sqlite3 prime-gap-search.db "PRAGMA foreign_keys=1; DELETE FROM range WHERE p=$P AND d=$D; SELECT total_changes()"
sqlite3 prime-gap-search.db "DELETE FROM m_stats WHERE p=$P AND d=$D AND m BETWEEN 1 AND $MINC; SELECT total_changes()"

for mp in {100,500,1000,2000,10000}; do
    time ./combined_sieve --save-unknowns -q -p $P -d $D --mstart 1 --minc $MINC --sieve-length $SL --max-prime ${mp} || true;
done

fn_base="1_${P}_${D}_${MINC}_s${SL}_l"
for fn in `ls -rt "${fn_base}"*`; do
    ./gap_stats --save-unknowns --unknown-filename $fn || true;
done

python gap_test.py --run-prp --unknown-filename ${fn_base}10000M.txt --plots
