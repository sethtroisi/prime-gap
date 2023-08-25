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

echo

sqlite3 gaps.db -separator $'\t' 'SELECT merit,year,discoverer,startprime FROM gaps ORDER BY merit desc LIMIT 10'

echo -e "\n"

sqlite3 gaps.db -header -separator $'\t' 'SELECT year,COUNT(*),SUM(merit > 34),ROUND(AVG(merit),1) FROM gaps WHERE year > 2010 GROUP BY 1'

echo -e "\n"

sqlite3 gaps.db 'SELECT (125000/2 - COUNT(*)) ||" Missing Below 125,000" FROM gaps WHERE gapsize <= 125000'
sqlite3 gaps.db 'SELECT (150000/2 - COUNT(*)) ||" Missing Below 150,000" FROM gaps WHERE gapsize <= 150000'

echo -e "\n\n"

sqlite3 gaps.db 'SELECT * FROM gaps WHERE gapsize in (42752, 98002) ORDER by year DESC LIMIT 5'
echo "..."
sqlite3 gaps.db 'SELECT COUNT(*)||" Remaining" FROM gaps WHERE gapsize in (42752, 98002) AND year < 2020'

echo -e "\n\n"

MANY="$(sqlite3 gaps.db 'SELECT startprime FROM gaps' | sed -n 's!^[^/]*[^0-9/]\([0-9]\+#\).*!\1!p' | sort | uniq -c | sort -nr | sed 's/#//')"

#echo "$MANY" | awk '{ ("sqlite3 gaps.db \"SELECT discoverer,count(*),round(max(merit),2),round(avg(merit),2),min(gapsize),max(gapsize) FROM gaps WHERE startprime GLOB \\\"*[ *]" $2 "*\\\" GROUP BY 1 ORDER BY 2 DESC LIMIT 1\"") | getline freq; print($0 "\t" freq); }'

#echo -e "\n\n"

echo -e "P#\tdiscoverer\tcount\tmax\tavg\tmin\tmin\tmax\t5th missing\tno lesser merit"
echo -e "\t\t\tby disc\tmerit\tmerit\tmerit\tgap\tgap\tgap"
echo -e "------------------------------------------------------------------------------------"

for P in $(echo "$MANY" | head -n 30 | awk '{print $2}'); do
        sqlite3 -separator $'\t' gaps.db "SELECT \"$P#\",*,fifth_missing-min_gap,ROUND(1.0*(fifth_missing-min_gap)/$P,2) FROM (SELECT *,(SELECT gapsize FROM gaps WHERE (gapsize/merit-15)>(min_gap/min_merit) LIMIT 1 OFFSET 4) as fifth_missing FROM (SELECT discoverer,count(*),round(max(merit),2),round(avg(merit),2),min(merit) as min_merit,min(gapsize) as min_gap,max(gapsize) FROM gaps WHERE startprime GLOB \"*[ *]$P#*\" GROUP BY 1 ORDER BY 2 DESC LIMIT 1))"
done


# Find large (33+) improved merit
# See https://www.mersenneforum.org/showpost.php?p=582856&postcount=47
