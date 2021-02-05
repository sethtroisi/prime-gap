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

usage() {
  echo "Usage: $0 [OPTION]... [UNKNOWN_FN]"
  echo ""
  echo -e "\t-u UNKNOWN_FN"
  echo -e "\t-t THREADS"
  echo -e "\t-h Print this help message"
}

UNKNOWN_FN=""
THREADS=1

while getopts u:t:h flag; do
    case "${flag}" in
        u) UNKNOWN_FN=${OPTARG};;
        t) THREADS=${OPTARG};;
        h) usage; exit 0;;
    esac
done

if [ -z "$UNKNOWN_FN" ]; then
  UNKNOWN_FN="${@:$OPTIND:1}"
fi

if [ -z "$UNKNOWN_FN" ]; then
  echo "Missing -u"
  echo
  usage
  exit 1;
fi

set -x

echo "Running --threads $THREADS --unknown-filename $UNKNOWN_FN"

make clean combined_sieve gap_stats

time ./combined_sieve -t $THREADS --save -u "$UNKNOWN_FN"
time ./gap_stats      -t $THREADS --save -u "$UNKNOWN_FN"

time ./gap_test.py -t $THREADS -u "$UNKNOWN_FN"
