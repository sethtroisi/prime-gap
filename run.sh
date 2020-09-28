set -eux

UNKNOWN_FN="$1"

echo "Running --unknown-filename $UNKNOWN_FN"

make clean
make combined_sieve gap_stats

time ./combined_sieve --save-unknowns --unknown-filename "$UNKNOWN_FN"

time ./gap_stats --save-unknowns --unknown-filename "$UNKNOWN_FN"

echo "Run when thread are free"
echo "time python missing_gap_test.py -t <THREADS> --unknown-filename $UNKNOWN_FN"
