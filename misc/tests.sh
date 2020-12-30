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


PARAMS="-p 907 -d 2190 --mstart 1 --minc 200 --max-prime 100 --sieve-length 11000"
TEST_DB="local_tests.db"


#### SETUP ####

make --quiet clean
make --quiet combined_sieve gap_stats gap_test_simple -j3

rm -f $TEST_DB unknowns/907_2190_1_200_s11000_l100M.{txt,m1.txt}
sqlite3 $TEST_DB < schema.sql


#### COMBINED SIEVE ####

./combined_sieve --method1 -qqq --save-unknowns $PARAMS --search-db $TEST_DB
./combined_sieve           -qqq --save-unknowns $PARAMS --search-db $TEST_DB

# Verify md5sum unknowns/907_2190_1_200_s11000_l100M.{txt,m1.txt}
md5sum -c <(echo "080309453b4310e0310a4fb4d1779ffe  unknowns/907_2190_1_200_s11000_l100M.txt")
md5sum -c <(echo "080309453b4310e0310a4fb4d1779ffe  unknowns/907_2190_1_200_s11000_l100M.m1.txt")


#### GAP STATS ####

./gap_stats --save-unknowns --unknown-filename 907_2190_1_200_s11000_l100M.txt --search-db $TEST_DB | tee temp_tests.log

grep -q 'avg missing prob : 0.0000000' temp_tests.log
grep -q 'RECORD : top  50% (    26) sum(prob) = 1.50e-05 (avg: 5.77e-07)' temp_tests.log
grep -q 'RECORD : top 100% (    53) sum(prob) = 2.15e-05 (avg: 4.06e-07)' temp_tests.log


#### MISC ####
#python misc/double_check.py --seed=123 --unknown-filename unknowns/907_2190_1_200_s11000_l100M.txt -c 5 --B1 1000


#### GAP_TEST ####

./gap_test_simple  --unknown-filename 907_2190_1_200_s11000_l100M.txt --min-merit 8 -q | tee temp_tests.log
# erase (XX.YY/sec)  Z seconds elapsed
sed -E -i -e 's#[0-9.]+( ?)(tests)?/sec\)#XX.YY\1\2/sec)#' -e 's#[0-9]+ sec#Z sec#' temp_tests.log

md5sum -c <(echo "42b30c4e7850aa9dbdda509df20a9e92  temp_tests.log")

python gap_test.py --unknown-filename 907_2190_1_200_s11000_l100M.txt --min-merit 8 --search-db $TEST_DB
diff <(echo "2") <(sqlite3 local_tests.db 'SELECT COUNT(*) FROM result WHERE merit > 8')
diff <(echo "1215") <(sqlite3 local_tests.db 'SELECT SUM(prp_next+prp_prev) FROM m_stats')


#### FINALIZE ####

set +x;
rm $TEST_DB

green=`tput setaf 2`
reset=`tput sgr0`

echo "${green}Success!${reset}"
