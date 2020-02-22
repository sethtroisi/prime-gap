# Combined Sieve - a new program to find prime gaps.

This is the combination of a couple of ideas I've had bouncing around in my
head while working on gmp mpz_prevprime.

### TODO

* gap_test.py
  * [ ] Store ALL results to sql
  * [ ] Starting at m > mstart
  * [ ] Print average to test count
  * [ ] Autoscale printing to every X seconds
  * [ ] Sort by expected gap and PRP only top X%
  * [ ] Try and describe distribution
  * [x] Generate expected length
* gap_search.cpp
  * [ ] Verify `sieve_length` math with d > 1
  * [ ] Allow for long mi by using bucketed `large_prime_queue`
    * [ ] Store remainder and prime in same array (and don't store pi for large primes)
  * [ ] Option to only output m with gcd(m, d) != 1
  * [x] Only store prime/remainder for primes that divide ANY mi.
  * [x] `sieve_range` > 4B
  * [x] Dynamic `sieve_length`
  * [x] Dynamic `sieve_range`


## Pgsurround.pl benchmark
On a i7-2600k single threaded.

### Pgsurround.pl

| Pn   | P#    | QPS (delta=0) |
|------|-------|---------------|
| 96   | 503   | 77.26         |
| 169  | 1009  | 8.68          |
| 303  | 1999  | 0.759         |
| 670  | 5003  | 0.0655        |

### Just sieving

|  P#   | Width (merit 20)  | Depth     | avg\_remaining | QPS  |
|-------|-------------------|-----------|---------------|------|
| 503   | 9792              | 116490    | 703           | 1693 |
| 1009  | 19456             | 811588    | 1123          | 223  |
| 1999  | 38976             | 8660119   | 1812          | 20.9 |
| 5003  | 98624             | 173202211 | 3665          | .895 |


## Commands

```bash
$ g++ -Wall -Werror -O3 gap_search.cpp gap_common.cpp -lgmp -o gap_search
# ./gap_search <params>
$ time ./gap_search -p 907 -d 210 --mstart 21400000 --minc 5000 --save-unknowns --sieve-only
AUTO SET: sieve length (coprime: 1521, prob_gap longer 0.80%): 6151
AUTO SET: sieve range (log(t) = ~879): 100000000
...
	21400003   99 <- unknowns -> 115
	21400037  110 <- unknowns -> 102
	    tests     10         (136.20/sec)  0 seconds elapsed
	    unknowns  2034       (avg: 203.40), 98.35% composite  49.46 <- % -> 50.54%
	    large prime remaining: 4343186 (avg/m: 13453)
	21400433  100 <- unknowns -> 98
	    tests     100        (121.27/sec)  1 seconds elapsed
	    unknowns  19945      (avg: 199.45), 98.38% composite  49.82 <- % -> 50.18%
	    large prime remaining: 4219113 (avg/m: 15368)
	21402181  101 <- unknowns -> 93
	    tests     500        (120.07/sec)  4 seconds elapsed
	    unknowns  100114     (avg: 200.23), 98.37% composite  49.77 <- % -> 50.23%
	    large prime remaining: 3447123 (avg/m: 15437)
	21404371   99 <- unknowns -> 102
	    tests     1000       (120.09/sec)  8 seconds elapsed
	    unknowns  200572     (avg: 200.57), 98.37% composite  49.84 <- % -> 50.16%
	    large prime remaining: 1404021 (avg/m: 15474)

real	0m11.539s
user	0m11.469s
sys	0m0.069s
...


```

### Average unknowns in sieve
```
# Lower side
$ g++ -Wall -Werror -O3 gap_search.cpp gap_common.cpp -lgmp -o gap_search
$ time ./gap_search -p 503 -d 1 --mstart 50000 --minc 10000 --save-unknowns --sieve-only --sieve-range 200
...
	prob_prime_coprime: 0.02275
AUTO SET: sieve length (coprime: 210, prob_gap longer 0.80%): 2011
...
	2.937% of sieve should be unknown (200M) ~= 69
	0.203% of 209 digit numbers are prime
	6.925% of tests should be prime (34.0x speedup)
	~2x14.4 = 28.9 PRP tests per m
	sieve_length=2011 is insufficient ~0.82% of time
...
	59999   74 <- unknowns -> 83
	    tests     10000      (734.73/sec)  14 seconds elapsed
	    unknowns  1373564    (avg: 137.36), 96.58% composite  50.01 <- % -> 49.99%
	    large prime remaining: 0 (avg/m: 1304)
```


## Benchmark

```bash
# Record timing and stats from sieving first 100,000 multiples after 100M
$ time ./gap_search -p 953 -d 210 --mstart 100000000 --minc 100000 --save-unknowns --sieve-only --sieve-length 8192

AUTO SET: sieve range (log(t) = ~928): 100000000

Testing m * 953#/210, m = 100000000 + [0, 100000)

sieve_length: 2x8192
sieve_range:  100000000

K = 1313 bits, 396 digits, log(K) = 909.57
Min Gap ~= 11135 (for merit > 12.0)

	PrimePi(100,000,000) = 5,761,455 (2 ... 99,999,989)
	Using 78,498 primes for SIEVE_SMALL(1,000,000)
...
	100000001  139 <- unknowns -> 130
	100000039  140 <- unknowns -> 135
	    tests     10         (109.13/sec)  0 seconds elapsed
	    unknowns  2799       (avg: 279.90), 98.29% composite  49.59 <- % -> 50.41%
	    large prime remaining: 5647957 (avg/m: 18862)
	100000429  144 <- unknowns -> 133
	    tests     100        (100.78/sec)  1 seconds elapsed
	    unknowns  27838      (avg: 278.38), 98.30% composite  49.49 <- % -> 50.51%
	    large prime remaining: 5647829 (avg/m: 20249)
	100002181  125 <- unknowns -> 129
	    tests     500        (99.27/sec)  5 seconds elapsed
	    unknowns  139125     (avg: 278.25), 98.30% composite  49.89 <- % -> 50.11%
	    large prime remaining: 5647199 (avg/m: 20551)
	100004369  149 <- unknowns -> 147
	    tests     1000       (99.11/sec)  10 seconds elapsed
	    unknowns  278194     (avg: 278.19), 98.30% composite  49.95 <- % -> 50.05%
	    large prime remaining: 5646336 (avg/m: 20590)
	100021871  139 <- unknowns -> 127
	    tests     5000       (98.79/sec)  51 seconds elapsed
	    unknowns  1391074    (avg: 278.21), 98.30% composite  49.99 <- % -> 50.01%
	    large prime remaining: 5637996 (avg/m: 20614)

real	3m51.650s
user	3m51.412s
sys	0m0.205s
```
