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

UNKNOWN_FN="$1"

echo "Running --unknown-filename $UNKNOWN_FN"

make clean
make combined_sieve gap_stats

time ./combined_sieve --save-unknowns --unknown-filename "$UNKNOWN_FN"

time ./gap_stats --save-unknowns --unknown-filename "$UNKNOWN_FN"

echo "Run when thread are free"
echo "time python missing_gap_test.py -t <THREADS> --unknown-filename $UNKNOWN_FN"
