# Combined Sieve - a new program to find prime gaps.

This is the combination of a couple of ideas I had while working on gmp mpz\_prevprime.

### TODO

* [ ] Make =SL included in sieve (e.g. change < SL to <= SL)
* gap\_test.py
  * [ ] Store ALL results to sql
  * [ ] Starting at m > mstart
  * [ ] Plot average tests count
  * [ ] Sort by expected gap and PRP only top X%
  * [x] Autoscale printing to every X seconds
  * [x] Describe distribution
  * [x] Generate expected length
* gap\_stats.cpp
  * [ ] Calc record chance
  * [ ] Save to prime-gap.db
  * [ ] drop directory from `unknown-filename`
  * [ ] Investigate only testing record gap combos
  * [x] load merit from gap.db
  * [x] Load from unknown\_fn
* gap\_search.cpp
  * [x] Verify `sieve_length` math with d > 1
    * [x] Calculate `sieve_length` for all (m % d)
  * [ ] Option to output m with gcd(m, d) != 1
  * [x] (Method2) keeping all composites in memory (and immediately marking mi off)
    * No large startup cost.
    * [ ] Do all primes for small ms (to get better memory access patterns)
    * [ ] Ctrl-c then just writes out the results at that point.
    * [x] (didn't help) check GCD with P# to avoid writting to main memory.
    * [x] Printing at each 1B primes tells you PRP tests saved / time taken.
  * [ ] Allow for long mi by using bucketed `large_prime_queue`
    * [x] Store remainder and prime in same array
    * [x] don't store pi for large primes (just pass around the pair)
  * [x] Don't save to large\_prime\_queue[next\_mi] with (next\_mi, d) > 1
  * [x] Only store prime/remainder for primes that divide ANY mi.
  * [x] `sieve_range` > 4B
  * [x] Dynamic `sieve_length`
  * [x] Dynamic `sieve_range`
* record\_check.py
  * [x] Read from sql db
* double\_check.py
  * [ ] run ecm on random unknowns and verify factors found > sieve limit

## gap\_search && gap\_test
On a i7-2600k single threaded.

###

| Pn   | P#    | M/second  | PRP/second |
|------|-------|-----------|------------|
| 96   | 503   | 86        | 2081       |
| 169  | 1009  | 12        | 560        |
| 303  | 1999  | 1.03      | 92.7       |
| 670  | 5003  | 24s/test  | 9.40       |
| 1230 | 10007 | 154s/test | 2.81       |

### Memory use

| P#    | SL   | M\_inc   | Memory(MB) | Time(s) | PRP/prime |
|-------|------|----------|------------|---------|-----------|
| 10007 |  1e9 | 80,000   | 1149       | 1040    | 268       |
| 10007 |  2e9 | 80,000   | 2075       | 1176    | 260       |
| 10007 |  4e9 | 80,000   | 3350       | 1333    | 252       |
| 10007 |  8e9 | 80,000   | 4913       | 1488    | 244       |
| 10007 | 16e9 | 80,000   | 6507       | 1813    | 236       |
| 10007 | 32e9 | 80,000   | 8544       | 2327    | 230       |

## Pgsurround.pl benchmark
On a i7-2600k single threaded.

### Pgsurround.pl

| Pn   | P#    | M/second      | Estimated PRP/second (TO VERIFY) |
|------|-------|---------------|----------------------------------|
| 96   | 503   | 77.26         | 3575                             |
| 169  | 1009  | 8.68          | 690                              |
| 303  | 1999  | 0.759         | 103.5                            |
| 670  | 5003  | 0.0655        | 19.08                            |

### Just sieving

|  P#   | Width (merit 20)  | Depth     | avg\_remaining | QPS  |
|-------|-------------------|-----------|----------------|------|
| 503   | 9792              | 116490    | 703            | 1693 |
| 1009  | 19456             | 811588    | 1123           | 223  |
| 1999  | 38976             | 8660119   | 1812           | 20.9 |
| 5003  | 98624             | 173202211 | 3665           | .895 |


## Commands

```bash
$ g++ -Wall -Werror -O3 gap_search.cpp gap_common.cpp -lgmp -o gap_search
$ time ./gap_search -p 503 -d 30030 --mstart 1000000 --minc 25000 --save-unknowns --sieve-only --sieve-range 10000
$ sqlite3 prime-gaps.db < schema.sql
$ time python3 gap_test.py --save-logs --min-merit=11 --plot --unknown-filename 1000000_503_30030_25000_s3181_l10000M.txt
```

```bash
$ g++ -Wall -Werror -O3 gap_search.cpp gap_common.cpp -lgmp -o gap_search
# ./gap_search <params>
$ time ./gap_search -p 907 -d 210 --mstart 21400000 --minc 5000 --save-unknowns --sieve-only --sieve-range 1000
AUTO SET: sieve length (coprime: 1521, prob_gap longer 0.80%): 6151
...
	21400003   81 <- unknowns -> 92
	21400037  100 <- unknowns -> 95
	    tests     10         (313.17/sec)  0 seconds elapsed
	    unknowns  1794       (avg: 179.40), 98.54% composite  50.56 <- % -> 49.44%
	    large prime remaining: 3687887 (avg/test: 5820)
	21400433   94 <- unknowns -> 84
	    tests     100        (320.48/sec)  0 seconds elapsed
	    unknowns  17738      (avg: 177.38), 98.56% composite  49.84 <- % -> 50.16%
	    large prime remaining: 3458886 (avg/test: 5838)
	21402181   89 <- unknowns -> 82
	    tests     500        (331.56/sec)  2 seconds elapsed
	    unknowns  89008      (avg: 178.02), 98.55% composite  49.79 <- % -> 50.21%
	    large prime remaining: 2365877 (avg/test: 5831)
	21404371   88 <- unknowns -> 88
	    tests     1000       (355.21/sec)  3 seconds elapsed
	    unknowns  178476     (avg: 178.48), 98.55% composite  49.86 <- % -> 50.14%
	    large prime remaining: 683290 (avg/test: 5829)

real	0m21.572s
user	0m21.504s
sys	0m0.068s
```
