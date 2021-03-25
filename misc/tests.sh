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

shopt -s expand_aliases
alias trace_on='set -x'
alias trace_off='{ set +x; } 2>/dev/null'

trace_on

PARAMS="-p 907 -d 2190 --mstart 1 --minc 200 --max-prime 100 --sieve-length 11000"
FN1="907_2190_1_200_s11000_l100M.txt"
FN1M1="${FN1/.txt/.m1.txt}"
FN2="953_1_1_10000_s8000_l20M.txt"
TEST_DB="local_tests.db"
COMMIT="49CC07DA"
TEST_GAPS_DB="test_gaps_$COMMIT.db"

#### MISC ####

# make gap_stats && rm -f test.db && sqlite3 test.db < schema.sql
# TMP=unknowns/1511_211890_1_1000000_s18000_l100000M.txt
# ./gap_stats --save-unknowns --search-db test.db --unknown-filename TMP
# ./gap_test.py --search-db test.db --unknown-filename TMP --no-one-side-skip --num-plots 3

#### SETUP ####

make --quiet clean
make --quiet combined_sieve gap_stats gap_test_simple -j3 VALIDATE_FACTORS=1

trace_off
mkdir -p unknowns/

rm -f "$TEST_DB" "unknowns/907_2190_1_200_s11000_l100M."{txt,m1.txt} "unknowns/$FN2"
sqlite3 $TEST_DB < schema.sql

if [ ! -f "$TEST_GAPS_DB" ]; then  # create gaps.db from last commit of 2020
    trace_on
    echo "Creating stable $TEST_GAPS_DB from last commit of 2020"
    (cd prime-gap-list; git show $COMMIT:allgaps.sql) | sqlite3 "$TEST_GAPS_DB"
    trace_off
fi
trace_on


#### COMBINED SIEVE ####

./combined_sieve --method1 -qqq --save $PARAMS --search-db $TEST_DB
./combined_sieve           -qqq --save -u $FN1 --search-db $TEST_DB

# Verify md5sum unknowns/907_2190_1_200_s11000_l100M.{txt,m1.txt}
md5sum -c <(echo "15a5cbff7301262caf047028c05f0525  unknowns/$FN1")
md5sum -c <(echo "15a5cbff7301262caf047028c05f0525  unknowns/$FN1M1")

# Tests MIDDLE_THRESHOLD
./combined_sieve           -qqq --save -u $FN2 --search-db $TEST_DB
md5sum -c <(echo "a8d85c965ddf2d9601e93c9920a9aa12  unknowns/$FN2")

# Test multithreaded
rm unknowns/907_2190_1_200_s11000_l100M.txt
./combined_sieve           -qqq --save -u $FN1 --search-db $TEST_DB -t 4
md5sum -c <(echo "15a5cbff7301262caf047028c05f0525  unknowns/907_2190_1_200_s11000_l100M.txt")

# Test --rle
./misc/convert_rle.py -u $FN1M1
md5sum -c <(echo "2ed7b0b9621dab602f204865360a3c56  unknowns/907_2190_1_200_s11000_l100M.m1.txt2")
rm unknowns/907_2190_1_200_s11000_l100M.m1.txt2

rm "unknowns/$FN1"
./combined_sieve           -qqq --save -u $FN1 --search-db $TEST_DB --rle
md5sum -c <(echo "2ed7b0b9621dab602f204865360a3c56  unknowns/907_2190_1_200_s11000_l100M.txt")

# Test --bitcompress
./misc/convert_rle.py -u $FN1M1 --bitcompress
md5sum -c <(echo "6edd5f50c4f588890a168461a7240c47  unknowns/907_2190_1_200_s11000_l100M.m1.txt2")
rm unknowns/907_2190_1_200_s11000_l100M.m1.txt2

rm "unknowns/$FN1"
./combined_sieve           -qqq --save -u $FN1 --search-db $TEST_DB --bitcompress
md5sum -c <(echo "6edd5f50c4f588890a168461a7240c47  unknowns/907_2190_1_200_s11000_l100M.txt")


#### GAP STATS ####

trace_off
DBS="--search-db $TEST_DB --prime-gaps-db $TEST_GAPS_DB"
trace_on

./gap_stats --save -u $FN1 $DBS --min-merit 8 | tee temp_tests.log

grep -q 'avg missing prob : 0.0000000' temp_tests.log
grep -q 'RECORD : top  50% (    26) sum(prob) = 1.50e-05 (avg: 5.77e-07)' temp_tests.log
grep -q 'RECORD : top 100% (    53) sum(prob) = 2.15e-05 (avg: 4.06e-07)' temp_tests.log

./gap_stats --save -u $FN2 $DBS -q -q

#### MISC ####
# python misc/double_check.py --seed=123 --unknown-filename $FN1 -c 5 --B1 1000


#### GAP_TEST ####

./gap_test_simple  -u $FN1M1 --min-merit 8 -q | tee temp_tests.log
# erase (XX.YY/sec)  Z seconds elapsed
sed -E -i -e 's#[0-9.]+( ?)(tests)?/sec\)#XX.YY\1\2/sec)#' \
          -e 's#[0-9]+ sec#Z sec#' -e 's#\.m1.txt#.txt#' temp_tests.log

md5sum -c <(echo "42b30c4e7850aa9dbdda509df20a9e92  temp_tests.log")

python gap_test.py -u $FN1 $DBS
diff <(echo "2") <(sqlite3 $TEST_DB 'SELECT COUNT(*) FROM result WHERE merit > 8')
diff <(echo "1215") <(sqlite3 $TEST_DB 'SELECT SUM(prp_next+prp_prev) FROM m_stats')

# TODO
#  - num-plots
#  - save-logs
#  - verify one sided skips


#### FINALIZE ####

trace_off
rm $TEST_DB temp_tests.log
rm -f "unknowns/$FN1" "unknowns/$FN2"

green=`tput setaf 2`
reset=`tput sgr0`

echo "${green}Success!${reset}"
