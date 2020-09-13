# Table of Contents

- [Benchmarks](#benchmarks)
  * [Microbenchmarks](#microbenchmarks)
  * [`gap_search --method2`](#gap_search-method2)
  * [`gap_test` (old numbers)](#gap_test-old-numbers)
  * [`Pgsurround.pl` benchmark](#pgsurroundpl-benchmark)
    + [`Pgsurround.pl`](#pgsurroundpl)
    + [Just sieving](#just-sieving)

# Benchmarks

## Microbenchmarks

Microbenchmarks of `modulo_search_*` are included in `benchmark`.
80%+ of `gap_search` is spent in `modulo_search_*` so optimizing this family of functions is vital.

**TODO**: Updated after `max_a` fix in benchmark
```bash
$ make benchmark
$ ./benchmark 100000
$ ./benchmark 100000 "# mod"
```
|  bits x count   | method\_name                    | found    | total    | time(s) | ns/iter | cycles/limb |
|-----------------|---------------------------------|----------|----------|---------|---------|-------------|
|   694 x  100000 | 503# mod <40 bit>p              | 100000   | 6088     | 0.0045  |      45 | 14.0        |
|  1390 x  100000 | 1009# mod <40 bit>p             | 100000   | 3292     | 0.0061  |      61 | 9.3         |
|  2799 x  100000 | 1999# mod <40 bit>p             | 100000   | 8627     | 0.0086  |      86 | 6.6         |
|  7099 x  100000 | 5003# mod <40 bit>p             | 100000   | 3010     | 0.0146  |     146 | 4.5         |
| 14291 x  100000 | 10007# mod <40 bit>p            | 100000   | 2119     | 0.0258  |     258 | 3.9         |
| 28588 x  100000 | 20011# mod <40 bit>p            | 100000   | 2854     | 0.0477  |     477 | 3.6         |

```bash
$ ./benchmark 100000 modulo_search
```

|  bits x count   | method\_name                     | found    | total    | time(s) | ns/iter | cycles/iter |
|-----------------|----------------------------------|----------|----------|---------|---------|-------------|
|    25 x  100000 | `modulo_search_one_op`           | 1033     | 100000   | 0.0020  |      20 | 68.8        |
|    25 x  100000 | `modulo_search_brute`            | 100000   | 100000   | 0.0316  |     316 | 1071.8      |
|    25 x  100000 | `modulo_search_euclid_small`     | 100000   | 100000   | 0.0074  |      74 | 249.9       |
|    25 x  100000 | `modulo_search_verify`           | 100000   | 100000   | 0.0253  |     253 | 857.9       |
|    25 x  100000 | `modulo_search_euclid`           | 100000   | 100000   | 0.0084  |      84 | 283.4       |
|    25 x  100000 | `modulo_search_euclid_gcd`       | 100000   | 100000   | 0.0114  |     114 | 386.1       |
|    25 x  100000 | `modulo_search_euclid_gcd2`      | 100000   | 100000   | 0.0110  |     110 | 373.2       |
|    25 x  100000 | `modulo_search_euclid_all_small` | 86763    | 292490   | 0.0277  |      95 | 321.5       |
|    25 x  100000 | `modulo_search_euclid_all_large` | 86763    | 292490   | 0.0295  |     101 | 342.5       |

|  bits x count   | method\_name                     | found    | total    | time(s) | ns/iter | cycles/iter |
|-----------------|----------------------------------|----------|----------|---------|---------|-------------|
|    40 x  100000 | `modulo_search_verify`           | 100000   | 100000   | 0.0289  |     289 | 981.2       |
|    40 x  100000 | `modulo_search_euclid`           | 100000   | 100000   | 0.0089  |      89 | 301.3       |
|    40 x  100000 | `modulo_search_euclid_gcd`       | 100000   | 100000   | 0.0108  |     108 | 365.1       |
|    40 x  100000 | `modulo_search_euclid_gcd2`      | 100000   | 100000   | 0.0093  |      93 | 315.8       |
|    40 x  100000 | `modulo_search_euclid_all_small` | 10       | 102484   | 0.0110  |     108 | 365.7       |
|    40 x  100000 | `modulo_search_euclid_all_large` | 10       | 102484   | 0.0110  |     107 | 364.6       |

## `gap_search` Method2

Testing performed on a i7-2600k compiled with `gcc-7`.

```bash
$ make gap_search
$ PARAMS="-d 210 --mstart 1 --save-unknowns --method2"
$ /usr/bin/time -v ./gap_search $PARAMS -p <P> --sieve-range <SL/1M> --minc <M_INC>
```

| P#     | `--sieve-range` | `--minc` | Memory(MB) | Time(s) | time/coprime m |
|--------|-----------------|----------|------------|---------|----------------|
|    997 |             1e9 |   10,000 | 6          | 30      | 0.013          |
|    997 |             1e9 |  100,000 | 18         | 165     | 0.0072         |
|--------|-----------------|----------|------------|---------|----------------|
|  2,003 |             1e9 |   10,000 | 8          | 44      | 0.019          |
|  2,003 |            10e9 |   10,000 | 8          | 224     | 0.098          |
|  2,003 |           100e9 |   10,000 | 8          | 2022    | 0.890          |
|--------|-----------------|----------|------------|---------|----------------|
|  4,007 |             1e9 |   10,000 | 10         | 62      | 0.027          |
|  4,007 |             1e9 |  100,000 | 56         | 490     | 0.021          |
|  4,007 |             1e9 |  500,000 | 260        | 3005    | 0.026          |
|--------|-----------------|----------|------------|---------|----------------|
| 10,007 |             1e9 |   10,000 | 18         | 133     | 0.057          |
| 10,007 |             1e9 |  100,000 | 127        | 1139    | 0.049          |

This hides a number of complexity.

* `--minc` scaling
  * linearly for small primes and sublinearly for large primes (>50x sieve-length).
  * larger values suffer from worse cache performance (or fail to allocate)
  * optimal value seem to be 50,000 - 400,000 based on P#.
* `--sieve-range`
  * Personal preference but 1-50B is reasonable.
  * Look at output to determine optimal (based on time per PRP test)
    * Stop when skipped PRP ~ 0.5-2x the rate you can test primes.
  * `gap_search` estimates skipped PRP/seconds AKA tests avoided from larger `--sieve-range`.
    * `gap_stats` allows testing of only the best `m`'s from the search (I test the top 5-20%)
    * `gap_stats` more accurately calculates probability of record with higher `--sieve-range`

Method2 output

```
2          (primes 1/1)	(seconds: 10.71/10.7 | per m: 9.4e-05)
100,003    (primes 9,592/9,593)	(seconds: 112.17/122.9 | per m: 0.0011)
500,009    (primes 7,678/41,539)	(seconds: 47.32/357.2 | per m: 0.0031)
1,000,003  (primes 36,960/78,499)	(seconds: 183.94/541.1 | per m: 0.0047)
	~ 2x 159.74 PRP/m		(9341133 new composites ~= 964316.7 skipped PRP => 5242.4 PRP/seconds)
...
10,000,019 (primes 316,066/664,580)	(seconds: 186.75/1228.8 | per m: 0.011)
50,000,017 (primes 567,480/3,001,135)	(seconds: 79.89/1757.0 | per m: 0.015)
100,000,007 (primes 2,760,321/5,761,456)	(seconds: 256.33/2013.4 | per m: 0.018)
	factors  10,805,266,428 		(interval: 77,586,357, avg m/large_prime interval: 28.1)
	unknowns 132,846,039/114287	(avg/m: 1162.39) (composite: 98.56% +0.056%)
	~ 2x 119.81 PRP/m		(5191959 new composites ~= 535370.7 skipped PRP => 2088.6 PRP/seconds)
...
```

TODO new `gap_tast` numbers.
TODO `gap_stat` numbers.
TODO `missing_gap_verify.py`

## `gap_test` (old numbers)


| Pn   | P#    | M/second  | PRP/second |
|------|-------|-----------|------------|
| 96   | 503   | 86        | 2081       |
| 169  | 1009  | 12        | 560        |
| 303  | 1999  | 1.03      | 92.7       |
| 670  | 5003  | 24s/test  | 9.40       |
| 1230 | 10007 | 154s/test | 2.81       |


## `Pgsurround.pl` benchmark
On the same i7-2600k system, running only a single thread

### `Pgsurround.pl`

| Pn   | P#    | M/second      | Estimated PRP/second (TO VERIFY) |
|------|-------|---------------|----------------------------------|
| 96   | 503   | 77.26         | 3575                             |
| 169  | 1009  | 8.68          | 690                              |
| 303  | 1999  | 0.759         | 103.5                            |
| 670  | 5003  | 0.0655        | 19.08                            |

### Just sieving

|  P#   | Width (merit 20)  | Depth       | avg\_remaining | QPS  |
|-------|-------------------|-------------|----------------|------|
| 503   | 9792              | 116,490     | 703            | 1693 |
| 1009  | 19456             | 811,588     | 1123           | 223  |
| 1999  | 38976             | 8,660,119   | 1812           | 20.9 |
| 5003  | 98624             | 173,202,211 | 3665           | .895 |

