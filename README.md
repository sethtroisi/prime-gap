# Combined Sieve - a new program to find prime gaps.

A fast prime gap searching suite.

# Table of Contents

- [Combined Sieve - a new program to find prime gaps.](#combined-sieve---a-new-program-to-find-prime-gaps)
  * [Tools](#tools)
    + [Combined Sieve (Sieve many `m * p#/d`)](#combined-sieve-sieve-many-m--pd)
      - [Method1 vs Method2 (Default)](#method1-vs-method2-default)
    + [Gap Stats](#gap-stats)
    + [Gap Test](#gap-test)
    + [Benchmark](#benchmark)
  * [Flow](#flow)
    + [Database](#database)
  * [Benchmarks](#benchmarks) | [BENCHMARKS.md](BENCHMARKS.md)
  * [Theory](#theory) | [THEORY.md](THEORY.md)
  * [Setup](#setup)
  * [Notes](#notes)
    + [Quick test of all functions](#quick-test-of-all-functions)
    + [TODO](#todo)
    + [TODONE](#todone)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>


## Tools

### Combined Sieve (Sieve many `m * p#/d`)

#### Method1 vs Method2 (Default)

`combined_sieve` parameter `--method1` changes how the sieve is performed.

```bash
$ make combined_sieve
$ time ./combined_sieve [--method1] -p <p> -d <d> --mstart <m_start> --minc <m_inc> --max-prime <LIMIT>M [--sieve-length <SR>] --save-unknowns
```

**Method1** performs a large initial setup phase to find the first `mi` that each `prime` divides before starting any sieving.
Each `mi` is then processed in order, looking at small primes and only those large primes that divide a number in the sieve interval.

A large percent of time will be spent watching dots cross this screen.
```
        Calculating first m each prime divides
        .
(later)
        .........
(eventually)
	......................................
```

**Method2** inverts this process and keeps all sieve intervals in memory at the same time.
Primes are iterated in order and all sieve intervals are handled simultaneously.

In practice Method2 is better. Method1 is only used to validate results.
* Status output is significantly nicer in Method2
* Stats are computed at regular intervals giving better insight into `--max-prime`. A
* Method2 allows for early stopping (set a very large `--max-prime` and press Ctrl+C once when you want to stop)
* Method2 is slightly faster because it use `modulo_search_euclid_all(...)` then a lookup for precomputed gcd.
  *  Method1 uses `modulo_search_euclid_gcd(...)` which calls `gcd(m, D)` which is significantly slower than lookup.

```bash
$ make combined_sieve
$ time ./combined_sieve -p 907 -d 210 --mstart 1 --minc 1000 --max-prime 1000 --save-unknowns --method1 -q
```

<details>
<summary>combined_sieve --method1 output</summary>
<p>

```bash
AUTO SET: sieve length: 7554 (coprime: 340, prob_gap longer 0.79%)

Testing m * 907#/210, m = 1 + [0, 1,000)
	Using 33860 primes for SIEVE_SMALL(400000)

	Calculating first m each prime divides
	......................................
	Sum of m1: 563865384

	Setup took 8.8 seconds

Saving unknowns to '907_210_1_1000_s7554_l1000M.m1.txt'
	1  104 <- unknowns -> 123
	    intervals 1          (1623.41/sec, with setup per m: 8.8)  0 seconds elapsed
	    unknowns  227        (avg: 227.00), 98.50% composite  45.81 <- % -> 54.19%
	    large prime remaining: 1220238 (avg/test: 7115)
	41  118 <- unknowns -> 113
	    intervals 10         (677.19/sec, with setup per m: 0.88)  0 seconds elapsed
	    unknowns  2264       (avg: 226.40), 98.50% composite  50.27 <- % -> 49.73%
	    large prime remaining: 1180663 (avg/test: 7127)
	437  118 <- unknowns -> 112
	    intervals 100        (621.68/sec, with setup per m: 0.089)  0 seconds elapsed
	    unknowns  22819      (avg: 228.19), 98.49% composite  49.61 <- % -> 50.39%
	    large prime remaining: 740948 (avg/test: 7146)
	997  120 <- unknowns -> 114
	    intervals 228        (747.02/sec, with setup per m: 0.04)  0 seconds elapsed
	    unknowns  52104      (avg: 228.53), 98.49% composite  49.78 <- % -> 50.22%
	    large prime remaining: 0 (avg/test: 7155)

real	0m9.107s
```

</p></details>


```bash
$ make combined_sieve
$ time ./combined_sieve -p 907 -d 210 --mstart 1 --minc 1000 --max-prime 1000 --save-unknowns
```

<details>
<summary>combined_sieve --method2 output</summary>
<p>

```
Testing m * 907#/210, m = 1 + [0, 1,000)
sieve_length: 2x 7,554
max-prime:       1,000,000,000   small_threshold:  137,790 (18.2 x SL)
coprime m    228/1000,  coprime i 1991/7554,  ~0MB


2          (primes 1/1)	(seconds: 0.00/0.0 | per m: 9e-06)
	factors       1,722,312 	(interval: 1,722,312 avg m/large_prime interval: 0.0)
	unknowns   363,888/228  	(avg/m: 1596.00) (composite: 89.44% +89.437% +3,080,964)
	~ 2x 698.79 PRP/m		(~ 74738.0 skipped PRP => 36330491.0 PRP/seconds)

100,003    (primes 9,592/9,593)	(seconds: 0.03/0.0 | per m: 0.00015)
	factors       3,784,530 	(interval: 2,062,218 avg m/large_prime interval: 0.0)
	unknowns    94,001/228  	(avg/m: 412.29) (composite: 97.27% +7.835% +269,887)
	~ 2x 42.07 PRP/m		(~ 299465.8 skipped PRP => 9528126.8 PRP/seconds)
...
1,000,003  (primes 36,960/78,499)	(seconds: 0.03/0.1 | per m: 0.00055)
	factors       3,971,493 	(interval: 46,806 avg m/large_prime interval: 1.3)
	unknowns    78,446/228  	(avg/m: 344.06) (composite: 97.72% +0.121% +4,167)
	~ 2x 35.06 PRP/m		(~ 844.5 skipped PRP => 24296.4 PRP/seconds)
...
100,000,007 (primes 2,760,321/5,761,456)	(seconds: 0.48/1.3 | per m: 0.0055)
	factors       4,232,762 	(interval: 35,036 avg m/large_prime interval: 0.0)
	unknowns    58,659/228  	(avg/m: 257.28) (composite: 98.30% +0.068% +2,340)
	~ 2x 26.29 PRP/m		(~ 468.8 skipped PRP => 977.0 PRP/seconds)
...
999,999,937 (primes 24,491,666/50,847,534)	(seconds: 4.60/9.6 | per m: 0.042)
	factors       4,339,924 	(interval: 30,726 avg m/large_prime interval: 0.0)
	unknowns    52,104/228  	(avg/m: 228.53) (composite: 98.49% +0.052% +1,779)
	~ 2x 23.37 PRP/m		(~ 368.8 skipped PRP => 80.2 PRP/seconds)
		Estimated ~31.0x faster to just run PRP now (CTRL+C to stop sieving)

Saving unknowns to '907_210_1_1000_s7554_l1000M.txt'

real	0m9.258s
```

</p></details>

Verify outputs are the same with

```bash
$ md5sum unknowns/907_210_1_1000_s7554_l1000M{.m1,}.txt
bfdb820237b40de0de2af71e8160ee66  907_210_1_1000_s7554_l1000M.m1.txt
bfdb820237b40de0de2af71e8160ee66  907_210_1_1000_s7554_l1000M.txt

$ diff 907_210_1_1000_s7554_l1000M.m1.txt 907_210_1_1000_s7554_l1000M.txt
<empty>

```

### Gap Stats

```bash
$ make combined_sieve gap_stats
# Use a larger `--sieve-length` for better record probability
$ time ./combined_sieve -p 907 -d 210 --mstart 1 --minc 1000 --max-prime 1000 --save-unknowns --sieve-length 15000 -qq
```

<details>
<summary>gap_stats output</summary>
<p>

```bashA
Testing m * 907#/210, m = 1 + [0, 1,000)
999,999,937 (primes 50,847,534/142,065,175)	(seconds: 9.56/9.6 | per m: 0.042)
	factors     177,931,017 	(interval: 8,887,445 avg m/large_prime interval: 0.0)
	unknowns   123,678/228  	(avg/m: 542.45) (composite: 98.19% +98.192% +6,716,550)
	~ 2x 23.37 PRP/m		(~ 382730.2 skipped PRP => 40026.5 PRP/seconds)


Saving unknowns to '907_210_1_1000_s15000_l1000M.txt'

$ time ./gap_stats --save-unknowns --unknown-filename 907_210_1_1000_s15000_l1000M.txt

Reading from 1_907_210_1000_s15000_l1000M.txt'

K = 1244 bits, 375 digits, log(K) = 861.69
Min Gap ~= 10340 (for merit > 12.0)

Found 3484 possible record gaps (20818 to 30158)
	If found Gap: 20818 (current: 23.66) would improve to 24.159
	If found Gap: 20822 (current: 24.07) would improve to 24.164
	If found Gap: 21164 (current: 24.47) would improve to 24.561

prob prime             : 0.0011592
prob prime after sieve : 0.04278

...

228 m's processed in 0.03 seconds (7475.76/sec)

	RECORD : top   1% (    2) => sum(prob) = 3.48e-06 (avg: 1.74e-06)
	RECORD : top   5% (   11) => sum(prob) = 1.35e-05 (avg: 1.23e-06)
	RECORD : top  10% (   22) => sum(prob) = 2.29e-05 (avg: 1.04e-06)
	RECORD : top  20% (   45) => sum(prob) = 3.85e-05 (avg: 8.57e-07)
	RECORD : top  50% (  114) => sum(prob) = 6.96e-05 (avg: 6.10e-07)
	RECORD : top 100% (  228) => sum(prob) = 9.47e-05 (avg: 4.15e-07)
Saved 12339 rows to 'gap_stats' table
Saved 228 rows to 'm_stats' table

real	0m0.122s
```

</p></details>

### Gap Test

Quick note, `python gap_test.py ...` and `./gap_test_simple` should be roughly equivilant.
`python gap_test.py` supports `--stats` to SLOWLY double check `./gap_stats` and `--plots` to print fun plots.


```bash
$ make gap_test_simple
$ time ./gap_test_simple --unknown-filename 907_210_1_1000_s15000_l1000M.txt -q
```

<details>
<summary>gap_test_simple output</summary>
<p>

```bash
Testing m * 907#/210, m = 1 + [0, 1,000)
Min Gap ~= 10340 (for merit > 12.0)

Reading unknowns from '907_210_1_1000_s15000_l1000M.txt'

228 tests M_start(1) + mi(0 to 996)

	m=1  244 <- unknowns -> 282 	1050 <- gap -> 600
	    tests     1          (49.43/sec)  0 seconds elapsed
	    unknowns  526        (avg: 526.00), 98.25% composite  46.39% <- % -> 53.61%
	    prp tests 39         (avg: 39.00) (1927.7 tests/sec)
	    best merit this interval: 1.915 (at m=1)
	m=41  282 <- unknowns -> 278 	 336 <- gap -> 1176
	    tests     10         (40.57/sec)  0 seconds elapsed
	    unknowns  5364       (avg: 536.40), 98.21% composite  49.53% <- % -> 50.47%
	    prp tests 537        (avg: 53.70) (2178.5 tests/sec)
	    best merit this interval: 9.038 (at m=29)
10920 12.6002  143 * 907#/210 -6796 to +4124
12164 14.0190  397 * 907#/210 -9644 to +2520
	m=437  272 <- unknowns -> 279 	1954 <- gap -> 6
	    tests     100        (42.54/sec)  2 seconds elapsed
	    unknowns  54318      (avg: 543.18), 98.19% composite  49.75% <- % -> 50.25%
	    prp tests 5089       (avg: 50.89) (2165.0 tests/sec)
	    best merit this interval: 14.019 (at m=397)
10644 12.2646  479 * 907#/210 -8038 to +2606
	m=997  257 <- unknowns -> 278 	2582 <- gap -> 250
	    tests     228        (45.33/sec)  5 seconds elapsed
	    unknowns  123678     (avg: 542.45), 98.19% composite  49.95% <- % -> 50.05%
	    prp tests 10740      (avg: 47.11) (2135.5 tests/sec)
	    best merit this interval: 12.265 (at m=479)

real	0m5.275s
```

</p></details>

or slightly more quiet

```bash
$ time ./gap_test_simple --unknown-filename 907_210_1_1000_s15000_l1000M.txt -qq
Testing m * 907#/210, m = 1 + [0, 1,000)
10920 12.6002  143 * 907#/210 -6796 to +4124
12164 14.0190  397 * 907#/210 -9644 to +2520
10644 12.2646  479 * 907#/210 -8038 to +2606
	m=997  257 <- unknowns -> 278 	2582 <- gap -> 250
	    tests     228        (44.98/sec)  5 seconds elapsed
	    unknowns  123678     (avg: 542.45), 98.19% composite  49.95% <- % -> 50.05%
	    prp tests 10740      (avg: 47.11) (2118.6 tests/sec)
```

### Benchmark

80%+ of the time is spent in `modulo_search_euclid` so it's important to optimize this function.

See [BENCHMARKS.md](BENCHMARKS.md) for more details

<details>
<summary>benchmark output</summary>
<p>

```bash
$ make benchmark
$ ./benchmark 100000
...
	|  bits x count   | method_name                    | found    | total    | time(s) | ns/iter | cycles/limb |
	|   694 x  100000 | 503# mod <40 bit>p             | 100000   | 6088     | 0.0045  |      45 | 14.0        |
	|  1390 x  100000 | 1009# mod <40 bit>p            | 100000   | 3292     | 0.0061  |      61 | 9.3         |
	|  2799 x  100000 | 1999# mod <40 bit>p            | 100000   | 8627     | 0.0086  |      86 | 6.6         |
	|  7099 x  100000 | 5003# mod <40 bit>p            | 100000   | 3010     | 0.0146  |     146 | 4.5         |
	| 14291 x  100000 | 10007# mod <40 bit>p           | 100000   | 2119     | 0.0258  |     258 | 3.9         |
	| 28588 x  100000 | 20011# mod <40 bit>p           | 100000   | 2854     | 0.0477  |     477 | 3.6         |

$ ./benchmark 100000 modulo_search
...
	|  bits x count   | method_name                    | found    | total    | time(s) | ns/iter | cycles/iter |
	|    25 x  100000 | single_mod_op                  | 1033     | 100000   | 0.0020  |      20 | 68.8        |
	|    25 x  100000 | modulo_search_brute            | 100000   | 100000   | 0.0316  |     316 | 1071.8      |
	|    25 x  100000 | modulo_search_euclid_small     | 100000   | 100000   | 0.0074  |      74 | 249.9       |
	|    25 x  100000 | modulo_search_verify           | 100000   | 100000   | 0.0253  |     253 | 857.9       |
	|    25 x  100000 | modulo_search_euclid           | 100000   | 100000   | 0.0084  |      84 | 283.4       |
	|    25 x  100000 | modulo_search_euclid_gcd       | 100000   | 100000   | 0.0114  |     114 | 386.1       |
	|    25 x  100000 | modulo_search_euclid_gcd2      | 100000   | 100000   | 0.0110  |     110 | 373.2       |
	|    25 x  100000 | modulo_search_euclid_all_small | 86763    | 292490   | 0.0277  |      95 | 321.5       |
	|    25 x  100000 | modulo_search_euclid_all_large | 86763    | 292490   | 0.0295  |     101 | 342.5       |
...

$ ./benchmark 100000 _all
	|  bits x count   | method_name                    | found    | total    | time(s) | ns/iter | cycles/iter |
	|-----------------|--------------------------------|----------|----------|---------|---------|-------------|
	|    25 x  100000 | modulo_search_euclid_all_small | 86755    | 292525   | 0.0282  |      96 | 327.1       |
	|    25 x  100000 | modulo_search_euclid_all_large | 86755    | 292525   | 0.0308  |     105 | 357.5       |

	|  bits x count   | method_name                    | found    | total    | time(s) | ns/iter | cycles/iter |
	|-----------------|--------------------------------|----------|----------|---------|---------|-------------|
	|    30 x  100000 | modulo_search_euclid_all_small | 18777    | 129541   | 0.0211  |     163 | 552.7       |
	|    30 x  100000 | modulo_search_euclid_all_large | 18777    | 129541   | 0.0217  |     168 | 569.0       |

	|  bits x count   | method_name                    | found    | total    | time(s) | ns/iter | cycles/iter |
	|-----------------|--------------------------------|----------|----------|---------|---------|-------------|
	|    31 x  100000 | modulo_search_euclid_all_small | 9755     | 115239   | 0.0206  |     179 | 606.3       |
	|    31 x  100000 | modulo_search_euclid_all_large | 9755     | 115239   | 0.0204  |     177 | 601.4       |

	|  bits x count   | method_name                    | found    | total    | time(s) | ns/iter | cycles/iter |
	|-----------------|--------------------------------|----------|----------|---------|---------|-------------|
	|    32 x  100000 | modulo_search_euclid_all_small | 4831     | 107701   | 0.0217  |     202 | 684.2       |
	|    32 x  100000 | modulo_search_euclid_all_large | 4831     | 107701   | 0.0222  |     206 | 698.6       |

	|  bits x count   | method_name                    | found    | total    | time(s) | ns/iter | cycles/iter |
	|-----------------|--------------------------------|----------|----------|---------|---------|-------------|
	|    35 x  100000 | modulo_search_euclid_all_small | 571      | 100811   | 0.0233  |     231 | 784.5       |
	|    35 x  100000 | modulo_search_euclid_all_large | 571      | 100811   | 0.0238  |     236 | 799.3       |

	|  bits x count   | method_name                    | found    | total    | time(s) | ns/iter | cycles/iter |
	|-----------------|--------------------------------|----------|----------|---------|---------|-------------|
	|    40 x  100000 | modulo_search_euclid_all_small | 18       | 100021   | 0.0240  |     240 | 813.5       |
	|    40 x  100000 | modulo_search_euclid_all_large | 18       | 100021   | 0.0242  |     242 | 821.5       |
```

</p></details>

## Flow

1. `combined_sieve`
   * Takes a set of parameters (`[m\_start, m\_start + m\_inc), P#, d, sieve-length, max-prime`)
   * Performs combined sieve algorithm
   * Produces `<PARAMS>.txt`, list of all numbers near `m * P#/d` with no small factors.
     * each row is `m: -count_low, +count_high | -18 -50 -80 ... | +6 +10 +12 ...`
   * Writes initial `range` entry (with `time_sieve`) into `range` table.
1. `gap_stats`
   * Calculates some statistics over each range in the interval
     * Expected gap in both directions
     * Probability of being a record gap, of being a missing gap, of being gap > `--min-merit`
   * Store statistics in `m_stats` table
   * Compute prob(gap) over all m and save in `range_stats` table
1. `gap_test_simple` / `gap_test.py` / `missing_gap_test` `--unknown-filename <PARAM>.txt`
   * Runs PRP tests on most likely m's from `<PARAM>.txt`
   * Stores results (as they are calculated in `result` & `m_stats`
     * Allows for saving of progress and graceful restart
   * `gap_test.py` can produce graphs of distribution (via `--plots`)


### Database

Often I run `combined_sieve` and `gap_stats` for a BUNCH of ranges not I never end up testing.

`misc/show_ranges.sh` updates and shows ranges from `prime-gap-search.db`.

`misc/drop_range.sh` deletes ranges when none of them have ever been tested

`misc/record_check.py` searches through `prime-gap-search.db` for records (over `gaps.db`)

`misc/finaliz.py -p <P> -d <D>` dumps `m_stats` to `<UNKNOWN>.csv` and `<UNKNOWN>.stats.txt` files.
It also DELETES DATA from `prime-gap-search.db` use it when you are done searching that `P#/d`

---

`misc/run.sh <UNKNOWN_FILENAME>` runs `combined_sieve`, `gap_stats`, and prints the `gap_test.py` command

`misc/status.py` to check if a filename has been processed 100% (and can be deleted)


## Benchmarks

See [BENCHMARKS.md](BENCHMARKS.md)

## Theory

Theory and justifaction for some calculations in present in [THEORY.md](THEORY.md)

## Setup

```bash
# Should probably build and install gmp head for faster next_prime
$ sudo apt install libgmp10

# slower but acceptable code is present up to commit 093a3b2b
$ sudo apt install libprimesieve-dev

$ sudo apt install libgmp-dev libmpfr-dev libmpc-dev
$ sudo apt install sqlite3 libsqlite3-dev
$ pip install --user gmpy2=2.1.0b5 tqdm primegapverify

$ sqlite3 prime-gap-search.db < schema.sql

# For misc/double_check.py
$ sudo apt install gmp-ecm
```

* [prime-gap-list](https://github.com/primegap-list-project/prime-gap-list)
    ```bash
    $ git clone https://github.com/primegap-list-project/prime-gap-list.git
    $ cd prime-gap-list
    $ mkdir -p logs/ unknowns/
    $ sqlite3 gaps.db < allgaps.sql
    $ cd <prime-gaps>
    $ ln -s ../prime-gap-list/gaps.db .
    ```

## Notes

* Anecdotally `combined_sieve` is ~8% faster when using `CC=clang++-11`
* There is some support for using OpenPFGW [PrimeWiki](https://www.rieselprime.de/ziki/PFGW) [SourceForge](https://sourceforge.net/projects/openpfgw/)
  * Unpack somewhere and link binary as `pfgw64`
  * may also need to `chmod a+x pfgw64`
  * modify `is_prime` in gap\_utils.py
* Multiple layers of verification of `combined_sieve`
  * Can compare `--method1` output with `--method2`
  * Can add use `make DEFINES="-DGMP_VALIDATE_FACTORS" combined_sieve`
  * `misc/double_check.py` double checks using `ecm`, to check that small factor aren't found for numbers in unknown-file.txt
    * `python misc/double_check.py --unknown-filename <unknown_file.txt> -c 10`
  * `skipped PRP/s` is checked in [THEORY.md](THEORY.md#skipped-prp-tests)
* Dev GMPlib | GMP 6.2.99
  * GMP 6.2.0 hasn't yet accepted my `mpz_prevprime` patch
    * `hg apply <patch>` from https://gmplib.org/list-archives/gmp-devel/2020-August/005851.html
    * If you are a developer consider asking telling them that `mpz_prevprime` would be useful

### Quick test of all functions

<details>
<summary>quick test commands and output</summary>
<p>

```bash
$ PARAMS="-p 907 -d 2190 --mstart 1 --minc 200 --max-prime 100 --sieve-length 11000"
$ make combined_sieve gap_stats gap_test_simple
$ time ./combined_sieve --method1 -qqq --save-unknowns $PARAMS
$ time ./combined_sieve           -qqq --save-unknowns $PARAMS
$ md5sum 907_2190_1_200_s11000_l100M.{txt,m1.txt}
080309453b4310e0310a4fb4d1779ffe  907_2190_1_200_s11000_l100M.txt
080309453b4310e0310a4fb4d1779ffe  907_2190_1_200_s11000_l100M.m1.txt

$ ./gap_stats --unknown-filename 907_2190_1_200_s11000_l100M.txt
# Verify RECORD @ 1,7,13,101,137
# Verify "avg seen prob:   0.9993981
# Verify "RECORD : top  50% (   26) => sum(prob) = 7.61e-07 (avg: 2.93e-08)"
# Verify "RECORD : top 100% (   53) => sum(prob) = 1.05e-06 (avg: 1.98e-08)"

# Optionaly validate combined_sieve results with ecm
$ python misc/double_check.py --seed=123 --unknown-filename 907_2190_1_200_s11000_l100M.txt -c 5 --B1 1000
1 * 907# / 2190 	unknowns 218 + 212 = 430
			+1423  had trivial factor of 2 skipping
	+6570 	should have small factor
		ecm: 2099 (1*907#/2190+6570)/2099
			+1373  had trivial factor of 2 skipping
			-1439  had trivial factor of 2 skipping
			+7789  had trivial factor of 2 skipping
			-2306  had trivial factor of 3 skipping
			+4377  had trivial factor of 2 skipping
	-4222 	should have small factor
		ecm: 3325199 (1*907#/2190+-4222)/3325199
		factor 3325199: 1319 2521
...
17 * 907# / 2190 	unknowns 216 + 208 = 424
	+2722 	shouldn't have small factor
		ecm: 17*907#/2190+2722
...
4150 trivially composite, 86 unknowns, 179 known composites
ecm found 172 composites: known 168 composite, 4 were unknown


$ ./gap_test_simple --unknown-filename 907_2190_1_200_s11000_l100M.txt -q --min-merit 8
Testing m * 907#/2190, m = 1 + [0, 200)
Min Gap ~= 7734 (for merit > 9.0)

Reading unknowns from '907_2190_1_200_s11000_l100M.txt'

53 tests M_start(1) + mi(0 to 198)

7750  9.0185  1 * 907#/2190 -450 to +7300
	m=1  218 <- unknowns -> 212 	 450 <- gap -> 7300
	    tests     1          (15.29/sec)  0 seconds elapsed
	    unknowns  430        (avg: 430.00), 98.05% composite  50.70% <- % -> 49.30%
	    prp tests 144        (avg: 144.00) (2201.5 tests/sec)
7462  8.6548  17 * 907#/2190 -6982 to +480
6938  8.0460  19 * 907#/2190 -6722 to +216
	m=37  214 <- unknowns -> 234 	1296 <- gap -> 250
	    tests     10         (34.35/sec)  0 seconds elapsed
	    unknowns  4262       (avg: 426.20), 98.06% composite  48.97% <- % -> 51.03%
	    prp tests 643        (avg: 64.30) (2208.7 tests/sec)
8520  9.8689  53 * 907#/2190 -4084 to +4436
	m=113  207 <- unknowns -> 215 	  40 <- gap -> 972
	    tests     30         (41.87/sec)  1 seconds elapsed
	    unknowns  12676      (avg: 422.53), 98.08% composite  49.31% <- % -> 50.69%
	    prp tests 1550       (avg: 51.67) (2163.0 tests/sec)
8382  9.6979  143 * 907#/2190 -1800 to +6582
7678  8.8820  163 * 907#/2190 -432 to +7246
	m=199  216 <- unknowns -> 205 	 192 <- gap -> 30
	    tests     53         (43.88/sec)  1 seconds elapsed
	    unknowns  22377      (avg: 422.21), 98.08% composite  49.97% <- % -> 50.03%
	    prp tests 2590       (avg: 48.87) (2144.5 tests/sec)


$ python gap_test.py --unknown-filename 907_2190_1_200_s11000_l100M.txt --min-merit 8
...
7750  9.0185  1 * 907#/2190 -7300 to +450
	  1  218 <- unknowns ->  212	7300 <- gap ->  450
7462  8.6548  17 * 907#/2190 -480 to +6982
6938  8.0460  19 * 907#/2190 -216 to +6722
	 37  214 <- unknowns ->  234	 250 <- gap -> 1296
	    tests     10         (35.7/sec)  0 seconds elapsed
	    unknowns  4262       (avg: 426.20), 98.06% composite  48.97% <- % -> 51.03%
	    prp tests 643        (avg: 64.30) (2294.854 tests/sec)
	    best merit this interval: 9.018 (at m=1)
8520  9.8689  53 * 907#/2190 -4436 to +4084
	113  207 <- unknowns ->  215	 972 <- gap ->   40
	    tests     30         (43.1/sec)  1 seconds elapsed
	    unknowns  12676      (avg: 422.53), 98.08% composite  49.31% <- % -> 50.69%
	    prp tests 1550       (avg: 51.67) (2226.529 tests/sec)
	    best merit this interval: 9.869 (at m=53)
8382  9.6979  143 * 907#/2190 -6582 to +1800
7678  8.8820  163 * 907#/2190 -7246 to +432
	199  216 <- unknowns ->  205	  30 <- gap ->  192
	    tests     53         (45/sec)  1 seconds elapsed
	    unknowns  22377      (avg: 422.21), 98.08% composite  49.97% <- % -> 50.03%
	    prp tests 2590       (avg: 48.87) (2199.432 tests/sec)
	    best merit this interval: 9.698 (at m=143)
```

</p></details>

### TODO

* Flow
* README.md
* THEORY.md
  * [ ] Add some theory for only doing one side test.
  * [ ] Add prob_record via greater on both sides
* Project
  * [ ] Records / day in status.py or record_check.py
  * [ ] All 63.8 bit factors?
    * [ ] Benchmark GMP_VALIDATE_LARGE_FACTORS

* combined\_sieve.cpp
* gap\_stats.cpp
* gap\_test.py
  * [ ] Sort by expected gap and PRP only top X%
  * [ ] `--top-x-percent` (see THEORY.md)
  * [ ] Read `time_sieve` and `time_stats` print optimal to restart search point
    * [ ] Leave XXX note for restart
* gap\_test.cpp
* schema.sql
  * [ ] next_p_i => next_p (fix in finalize and other places)
* benchmarking
  * [ ] Add instructions to verify `modulo\_search` is >80% of the time.

### Low Priority TODOs

* New script to update DB (and maybe delete uninteresting results
  * [ ] Finalize script
  * [ ] Set rid if null (but reasonable range exists)
* combined\_sieve.cpp
  * [ ] Option to output m with gcd(m, d) != 1
  * [ ] optimize D helper
  * [ ] (M_inc * max_prime < 63)
    * The trick is to compute (M_start + M_inc) * base_r as the combination of two additions.
    * While M_start * base_r can overflow, break it into
      * (max_m * (M_start / max_m) + (M_start % max_m)) * base_r
      * (d * ((max_m * base_r) % p) + (r * base_r)) % p
* gap\_stats.cpp
  * [ ] Produce P(record) / day estimate
  * [ ] Check if higher prob is related to unique (mi % d)
  * [ ] Option to starting at m > mstart
* gap\_test.py
  * [ ] Plan for minimize memory when using large `--top-x-percent`
  * [ ] Plot average tests count
* missing\_gap\_test.py && missing\_gap\_verify.py
  * [ ] Save to result | figure out plan for gap_test.py to reuse
* gap\_common
  * [ ] Sieve Interval up to 100M instead of using next_prime
  * [ ] Compute smaller PRP and use that to computer larger (slower) PRP estimate
* benchmarking
  * [ ] Validate `modulo\_search` is >80% execution time for `combined_sieve`.

### TODONE

<details>
<summary>Tracker for finished TODO items</summary>
<p>
* README.md
  * [x] record\_check.py guide to use
  * [x] Clarify gap\_test.py vs gap\_test.cpp
  * [x] Add Flow section based on comments in schema.sql
  * [x] Add commands for benchmarking
  * [x] Fill out gap test section
  * [x] Split out some benchmarking
* Project level
  * [x] Move UNKNOWN.txt to unknowns/
  * [x] Support --search-db everywhere
  * [x] change `ms_P_D_minc` to `P_D_ms_minc_`
  * [x] Make =SL included in sieve (e.g. change \< SL to \<= SL)
  * [x] Rename `gap_search` to `combined_sieve`
  * [x] Rename prime-gap.db to prime-gap-search.db
  * [x] Make method2 the default
  * [x] config.verbose in gap\_search, gap\_stats, gap\_test
* combined\_sieve.cpp
  * [x] Integrate faster prime iterator (libprimesieve)
  * [x] ETA for `combined_sieve` timing.
  * [x] Test number to be marked combosite see if already marked by 2 / 3 and skip (e.g. try to avoid memory lookup)
  * [x] "expect %ld left" doesn't match
  * [x] Write time_sieve into DB
  * Method2 (all composites in memory)
    * [x] `vector<char>` doesn't show improvement over `vector<bool>`
    * [x] Verify skipped PRP by testing at X and 2\*X
    * [x] Estimate PRP/s and include in status.
    * [X] Ctrl-C to quit early (but writes out the results at that point).
    * [x] Look at Method1 vs Method2 writeup and understand why outputs seem different
    * [x] Use reindex\_m trick to reduce number size of composites
    * [x] Do all primes for small ms (to get better memory access patterns)
    * [x] check GCD with P# to avoid writting to main memory. (didn't work)
    * [x] Printing at each 1B primes tells you PRP tests saved / time taken.
    * [x] Write up of Method1 vs Method2 (from memory M2 is faster but uses more memory)
  * Method1
    * [x] Store remainder and prime in same array
    * [x] don't store pi for large primes (just pass around the pair)
  * [x] Verify `sieve_length` math with d > 1
    * [x] Calculate `sieve_length` for all (m % d)
  * [x] Don't save to large\_prime\_queue[next\_mi] with (next\_mi, d) > 1
  * [x] Only store prime/remainder for primes that divide ANY mi.
  * [x] `max_prime` > 4B
  * [x] Dynamic `sieve_length`
  * [x] Dynamic `max_prime`
* gap\_stats.cpp
  * [x] Write all data that `gap_test.py` consumes to DB
  * [x] Tweak logging at different verbose levels
  * [x] Move missing gaps behind a compile flag
  * [x] drop directory from `unknown-filename`
  * [x] Missing gaps prob % filter
  * [x] Don't rewrite `partial\_results`
  * [x] Calc record chance
  * [x] Save expected/prob record to sql
  * [x] load merit from gap.db
  * [x] Load from unknown\_fn
* gap\_test.py
  * [x] Correlation between Expected Gap & Gap
  * [x] Show P(record) / day
  * [x] Multithreaded
  * [x] Update `m_stats` (with `test_time`),
  * [x] Plot P(gap > min\_merit) and P(record) sorted and unsorted.
  * [x] Option to toggle between OpenPFGW and gmp (see note in benchmarking below)
  * [x] Store all results to sql
  * [x] Autoscale printing to every X seconds
  * [x] Describe distribution
  * [x] Generate expected length
* missing\_gap\_verify.py & missing\_gap\_test.py
  * [x] read and update DB directly
  * [x] save to `m_stats`
  * [x] load/save from DB file.
  * [x] grouped output of BOTH PRIME (every X entries)
  * [x] record stats about prime found (potential candidate for high merit?)
  * [x] tee logging (for preemptible machines)
  * [x] gap\_utils.py (config parsing, isPrime, teelogger)
  * [x] --ignore-gaps 130898
  * [x] store records in some .gitignore'd file
  * [x] Multiprocessing
* benchmarking
  * [x] Experimentally compute `prime_time_estimate`
  * [x] Redo prime-time with user time (gmp better)
    * There's serious overhead using subproccess, and I see 150% CPU usage.
    * I Tried passing via stdin (pfgw supports -- but quits after one) but couldn't make it work.
    * MF suggest you can use gwnum library and hack some api but seems like a like of work for ~20-200% gain on reasonable primes.
  * [x] Add int64 `modulo_search_euclid_all`
  * [x] Add benchmarks for int32/int64 `modulo_search`
  * [x] Add benchmarks for `K % p`
* prime-gap-search.db
  * [x] A plan for how to clean up [partially] finished ranges
* misc/
  * Finalize script
    * [x] dump to csv, dump stats
  * show\_ranges.sh
    * [x] Update range.`num_processed / num_to_process`
    * [x] Verify every table.`m_stat` result in table.`result`
  * record\_check.py
    * [x] Read from sql db
  * double\_check.py
    * [x] run ecm on random unknowns and verify factors found > sieve limit
</p></details>
