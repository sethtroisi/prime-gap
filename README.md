# Combined Sieve - a new program to find prime gaps.

A fast prime gap searching suite.

# Table of Contents

- [Overview](#overview)
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

## Overview

This is a suite of tools (`combined_sieve`, `gap_stats`, and `gap_test.py`) which are useful for searching for prime gaps
centered around `m * P#/d`.

The flow is to sieve many intervals with `combined_sieve --save -u <P>_<D>_<MSTART>_<MINC>_s<SIEVE_LENGTH>_l<MAX_PRIME>.txt`
Then calculate statistics about these intervals with `gap_stats --save -u <SAME_FILENAME>.txt`
Then test some portion of these intervals with `gap_test.py -u <SAME_FILENAME>.txt --prp-top-percent 25`

A known 30 merit prime gap can be easily found with (After [Setup](#setup) has been complete)

```bash
$ make
$ sqlite3 test.db < schema.sql
$ FN=809_1841070_11244000_1000_s27000_l200M.txt
$ ./combined_sieve --save --search-db test.db -u $FN
$ ./gap_stats --save --search-db test.db -u $FN --min-merit 20
$ ./gap_test.py --search-db test.db -u $FN
...
23142  30.1485  11244911 * 809#/1841070 -5788 to +17354	RECORD!
...
```

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
	Using 33860 primes for SMALL_PRIME_LIMIT(400,000)

	Calculating first m each prime divides
	......................................
	Sum of m1: 563865384
	PrimePi(1000000000) = 50847534

	Setup took 6.4 seconds

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

real	0m6.671s
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
max_prime:        1,000,000,000

coprime m    228/1000,  coprime i     1991/7554, ~0MB
                        coprime wheel 346/7554, ~0MB

907        (primes 155/155)	(seconds: 0.00/0.0 | per m: 5.39e-07)
	factors               0 	(interval: 0 avg m/large_prime interval: 0.0)
	unknowns   157,918/228  	(avg/m: 692.62) (composite: 95.42% +95.416% +3,286,934)
	~ 2x 71.12 PRP/m		(~ 360955.6 skipped PRP => 2939330029.3 PRP/seconds)

100,000    (primes 8,363/9593)	(seconds: 0.04/0.0 | per m: 0.000194)
	factors          82,189 	(interval: 35,101 avg m/large_prime interval: 296.4)
	unknowns    94,001/228  	(avg/m: 412.29) (composite: 97.27% +97.271% +3,350,851)
	~ 2x 42.07 PRP/m		(~ 374203.7 skipped PRP => 8452672.8 PRP/seconds)
...
1,000,000  (primes 36,960/78499)	(seconds: 0.03/0.1 | per m: 0.000559)
	factors         110,815 	(interval: 8,186 avg m/large_prime interval: 21.1)
	unknowns    78,446/228  	(avg/m: 344.06) (composite: 97.72% +0.121% +4,167)
	~ 2x 35.06 PRP/m		(~ 844.5 skipped PRP => 27948.8 PRP/seconds)
...
100,000,000 (primes 544,501/5761456)	(seconds: 0.09/1.2 | per m: 0.00511)
	factors         156,334 	(interval: 907 avg m/large_prime interval: 0.2)
	unknowns    58,659/228  	(avg/m: 257.28) (composite: 98.30% +0.009% +325)
	~ 2x 26.29 PRP/m		(~ 69.0 skipped PRP => 785.3 PRP/seconds)
...
999,999,937 (primes 4,838,319/50847535)	(seconds: 0.84/8.4 | per m: 0.0368)
	factors         174,997 	(interval: 812 avg m/large_prime interval: 0.0)
	unknowns    52,104/228  	(avg/m: 228.53) (composite: 98.49% +0.008% +290)
	~ 2x 23.37 PRP/m		(~ 54.5 skipped PRP => 64.6 PRP/seconds)
		Estimated ~36.4x faster to just run PRP now (CTRL+C to stop sieving)

Saving unknowns to '907_210_1_1000_s7554_l1000M.txt'

real	0m13.558s
```

</p></details>

Verify outputs are the same with

```bash
$ md5sum unknowns/907_210_1_1000_s7554_l1000M{.m1,}.txt
d74bd47c8d43eac9a2309fc7722839ac  unknowns/907_210_1_1000_s7554_l1000M.m1.txt
d74bd47c8d43eac9a2309fc7722839ac  unknowns/907_210_1_1000_s7554_l1000M.txt

$ diff unknowns/907_210_1_1000_s7554_l1000M.m1.txt unknowns/907_210_1_1000_s7554_l1000M.txt
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

```bash
$ time ./combined_sieve -qq --save -u 907_210_1_1000_s15000_l1000M.txt
Testing m * 907#/210, m = 1 + [0, 1,000)
999,999,937 (primes 50,847,535/50847535)	(seconds: 8.36/8.4 | per m: 0.0367)
	factors         412,684 	(interval: 412,684 avg m/large_prime interval: 0.4)
	unknowns   123,678/228  	(avg/m: 542.45) (composite: 98.19% +98.192% +6,716,550)
	~ 2x 23.37 PRP/m		(~ 382730.2 skipped PRP => 45793.8 PRP/seconds)


Saving unknowns to 'unknowns/907_210_1_1000_s15000_l1000M.txt'

$ time ./gap_stats --save -u 907_210_1_1000_s15000_l1000M.txt

Reading from 1_907_210_1000_s15000_l1000M.txt'

K = 1244 bits, 375 digits, log(K) = 861.69
Min Gap ~= 15511 (for merit > 18.0)

Found 3484 possible record gaps (20818 to 30158)
	If found Gap: 20818 (current: 23.66) would improve to 24.159
	If found Gap: 20822 (current: 24.07) would improve to 24.164
	If found Gap: 21164 (current: 24.47) would improve to 24.561

prob prime             : 0.0011592
prob prime coprime     : 0.0140599
prob prime after sieve : 0.04278

Extended size: 27574 (32.0 merit)
Using Wheel: 210 for extended probs
	Average   817 inner    coprimes => 0.000898% prob_greater
	Average   939 extended coprimes => 0.000161% prob_greater
Extended prob records setup (0.11 seconds)

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

Quick note, `python gap_test.py ...` and `./gap_test_simple` should be roughly equivalent.
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
     * Use `combined_sieve --mstart 1 --minc 100 --max-prime 1 -d 4 -p <P>` for some insight on choose of `d`
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
1. `gap_test_simple` / `gap_test.py` `--unknown-filename <PARAM>.txt`
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

Theory and justification for some calculations in present in [THEORY.md](THEORY.md)

## Setup

In general this is going to be easy under Ubuntu 20.04 as you need **Python 3.7**
and **sqlite3 >= 3.24** both of which are harder to install in previous versions
of Ubuntu.

```bash
#$ sudo apt install libgmp10 libgmp-dev
#$ sudo apt install mercurial build-essential automake autoconf bison make libtool texinfo m4
# Building gmp from source is required until 6.3
$ hg clone https://gmplib.org/repo/gmp/
$ cd gmp
$ ./.bootstrap
$ mkdir build
$ cd build
$ ../configure
$ make
$ make check
$ make install
```

```
# Need sqlite >= 3.24 for upsert "on conflict" statement.
$ sudo apt install sqlite3

$ sudo apt install libmpfr-dev libmpc-dev libomp-dev libsqlite3-dev libbenchmark-dev

# Make sure this is >= python3.7 and not python2
$ python -m pip install --user gmpy2==2.1.0b5 primegapverify

# For misc/double_check.py
$ sudo apt install gmp-ecm
$ sudo apt install sympy
```

```
$ git clone https://github.com/sethtroisi/prime-gap.git
$ cd prime-gap

$ mkdir -p logs/ unknowns/
$ sqlite3 prime-gap-search.db < schema.sql

# There are a handful of constants (MODULE_SEARCH_SECS, PRIME_RANGE_SEC, ...)
# in gap_common.cpp `combined_sieve_method2_time_estimate` that can be set for
# more accurate `combined_sieve` timing.
```

* [prime-gap-list](https://github.com/primegap-list-project/prime-gap-list)
    ```bash
    $ git clone --depth 10 https://github.com/primegap-list-project/prime-gap-list.git
    $ sqlite3 gaps.db < prime-gap-list/allgaps.sql
    ```

## Notes

* Clang (`CC=clang++-12`) might speed up `combined_sieve` slightly.
* There is some support for using OpenPFGW [PrimeWiki](https://www.rieselprime.de/ziki/PFGW) [SourceForge](https://sourceforge.net/projects/openpfgw/)
  * Unpack somewhere and link binary as `pfgw64`
  * may also need to `chmod a+x pfgw64`
  * modify `is_prime` in gap\_utils.py
* Multiple layers of verification of `combined_sieve`
  * Can compare `--method1` output with `--method2`
  * Can add use `make clean combined_sieve VALIDATE_FACTORS=1` or `make clean combined_sieve VALIDATE_LARGE=1`
  * `misc/double_check.py` double checks using `ecm`, to check that small factor aren't found for numbers in unknown-file.txt
    * `python misc/double_check.py --unknown-filename <unknown_file.txt> -c 10`
  * `skipped PRP/s` is checked in [THEORY.md](THEORY.md#skipped-prp-tests)
* Compression
  * method1 doesn't support compression at this time so it also makes verifying `combined_sieve` slightly harder.
    * `misc/convert_rle.py` can convert existing files to RLE if you want to double check / shrink old files.
  * Run Length Encoding
    * `can be forced with `--rle` saves ~60% of space
  * `--bitcompress`
    * Saves ~95% of space, but makes code very hard
    * Stores 7 bits per byte (to avoid ascii control characters) of is_composite array
    * Reduces size of is_composite by removing any number coprime to K# or d


### Quick test of all functions

This is somewhat automated with `./misc/tests.sh` which should print a green `Success!`

<details>
<summary>quick test commands and output</summary>
<p>

```bash
$ PARAMS="-p 907 -d 2190 --mstart 1 --minc 200 --max-prime 100 --sieve-length 11000"
$ make combined_sieve gap_stats gap_test_simple
$ time ./combined_sieve --method1 -qqq --save-unknowns $PARAMS
$ time ./combined_sieve           -qqq --save-unknowns $PARAMS
$ md5sum unknowns/907_2190_1_200_s11000_l100M.{txt,m1.txt}
15a5cbff7301262caf047028c05f0525  unknowns/907_2190_1_200_s11000_l100M.txt
15a5cbff7301262caf047028c05f0525  unknowns/907_2190_1_200_s11000_l100M.m1.txt

$ ./gap_stats --unknown-filename 907_2190_1_200_s11000_l100M.txt
# Verify RECORD @ 1,7,13,101,137
# Verify "avg missing prob : 0.0000000"
# Verify "RECORD : top  50% (    26) sum(prob) = 1.50e-05 (avg: 5.77e-07)"
# Verify "RECORD : top 100% (    53) sum(prob) = 2.15e-05 (avg: 4.06e-07)"

# Optionally validate combined_sieve results with ecm
$ python misc/double_check.py --seed=123 --unknown-filename unknowns/907_2190_1_200_s11000_l100M.txt -c 5 --B1 1000
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
ecm found 178 composites: known 171 composite, 7 were unknown


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
* Project
  * [x] --update flag for `misc/show_ranges.sh`
    * [x] Finalize field
      * [ ] Sum and set test-time in finalize? (IGNORE very large values which are likely pause)
      * [ ] Finalize script should update ranges being finalized
    * [ ] Faster finalize.py
  * [ ] can modulo_search (equation 2/3) be modified so that it only returns odd multiples?
  * [ ] Check that correct indexes exist in sql db
  * [ ] Consider new names for prp-top-percent, no-one-side-skip, sieve-length
* combined\_sieve.cpp
  * [ ] Better sharing work for medium_primes (or skip unknown counting)
    * [ ] Unknown counting could be faster with boost or `bit_vector`
  * [ ] ncurse or similiar status line with all large_prime intervals & eta
    * Idea that made sense was to printf then reset position to start of status (so that future prints override them)
* gap\_stats.cpp
* gap\_test\_simple
  * [ ] Do in process for very small P
  * [ ] records to DB what range was tested
* gap\_test.py
  * [ ] Save prp-percentage finalized for faster skipping of complete files
  * [ ] Validate the unknown-file matches expectations (unknowns / line, ...)
  * [ ] Ctrl+C sometimes hangs waiting for one more result
* gap\_test\_simple.cpp
* schema.sql
* benchmarking
  * [ ] Add instructions to verify `modulo\_search` is >80% of the time.
  * [ ] Try out multi-threaded combined_sieve

### Low Priority TODOs

* Project
  * [ ] --rle could be 1/2 the size if I used coprime_x index instead of x
  * [ ] Figure out how to load (in c & python) and set config occasionally
  * [ ] Records / day in status.py or record_check.py
    * [ ] Test if sum(prob_record) matches with --no-one-side-skip
* THEORY.md
* New script to update DB (and maybe delete uninteresting results
  * [ ] Set rid if null (but reasonable range exists)
  * [ ] Check if range_stats needs cleanup after deleting range?
* combined\_sieve.cpp
  * [ ] Record final PRP/thread-second rate
  * [ ] Improve optimize D helper
    * [ ] Document D helper
  * [ ] Benchmark GMP_VALIDATE_LARGE_FACTORS
* gap\_stats.cpp
  * [ ] Produce P(record) / day / core estimate
  * [ ] Check if higher prob is related to unique (mi % d)
* gap\_test.py
  * [ ] Test least dense side (or portion) first (See MF Prime Gap Theory #125)
  * [ ] MEGAGAPS (what's remaining here?)
    * [ ] Always test next_p = -1 results (regardless of prp-top-percent)
  * [ ] Clean up fallback prev_prime
  * [ ] Do PRP top X% in order (for more dynamic cutoff)
  * [ ] Plot average tests count
* gap\_common
  * [ ] Load time estimate values from config file
  * [ ] Compute smaller PRP and use that to computer larger (slower) PRP estimate
* benchmarking
  * [ ] re-benchmark `combined_sieve`.

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
  * [x] convert unknown to start with m, not mi
  * [x] Change min-merit to 15, 18, or 20
  * [x] Run length encoding to reduce filesize
  * [x] Move UNKNOWN.txt to unknowns/
  * [x] Support --search-db everywhere
  * [x] change `ms_P_D_minc` to `P_D_ms_minc_`
  * [x] Make =SL included in sieve (e.g. change \< SL to \<= SL)
  * [x] Rename `gap_search` to `combined_sieve`
  * [x] Rename prime-gap.db to prime-gap-search.db
  * [x] Make method2 the default
  * [x] config.verbose in gap\_search, gap\_stats, gap\_test
  * [x] Make sure that next_p = 0, = -1, is handled correctly in places.
* combined\_sieve.cpp
  * [x] `M_inc * max_prime > 2^64` use `modulo_search_euclid_all_large`
  * [x] Improve error messages (gap_stats missing, SL size)
  * [X] Calculating coprime [L, R] * K^-1 mod p for medium p
  * [x] Wheel for coprime i (saves 50-80% of space)
  * [x] Add some theory for only doing one side test.
  * [x] Optimize D helper
  * [x] Integrate faster prime iterator (libprimesieve)
  * [x] ETA for `combined_sieve` timing.
  * [x] Test number to be marked composite see if already marked by 2 / 3 and skip (e.g. try to avoid memory lookup)
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
    * [x] check GCD with P# to avoid writing to main memory. (didn't work)
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
  * [x] Accumulate missing / skipped gap_prob at 0
  * [x] P(record) is missing double extended
  * [x] Benchmark reindex_m_wheel @ 6, 30, 210
  * [x] Print (1 - seen) / prob as "unknown"
  * [x] Wheel for extended range
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
  * [x] benchmark not calling commit() on every save result
  * [x] Load min-merit from db
  * [x] Read min-merit from DB in gap_test.py
  * [x] Reduce memory if not --stats / --num-plots
  * [x] First Ctrl+C stops new results, 2nd breaks.
  * [x] Use primegapverify.sieve for better prev_prime
  * [x] Dynamic one sided test threshold and logging.
  * [x] Extended gap for one sided tests?
  * [x] `--top-x-percent` (see THEORY.md)
  * [x] Correlation between Expected Gap & Gap
  * [x] Show P(record) / day
  * [x] Multi-threaded
  * [x] Update `m_stats` (with `test_time`),
  * [x] Plot P(gap > min\_merit) and P(record) sorted and unsorted.
  * [x] Option to toggle between OpenPFGW and gmp (see note in benchmarking below)
  * [x] Store all results to sql
  * [x] Autoscale printing to every X seconds
  * [x] Describe distribution
  * [x] Generate expected length
* benchmarking
  * [x] Experimentally compute `prime_time_estimate`
  * [x] Redo prime-time with user time (gmp better)
    * There's serious overhead using subprocess, and I see 150% CPU usage.
    * I Tried passing via stdin (pfgw supports -- but quits after one) but couldn't make it work.
    * MF suggest you can use gwnum library and hack some api but seems like a like of work for ~20-200% gain on reasonable primes.
  * [x] Add int64 `modulo_search_euclid_all`
  * [x] Add benchmarks for int32/int64 `modulo_search`
  * [x] Add benchmarks for `K % p`
* prime-gap-search.db
  * [x] A plan for how to clean up [partially] finished ranges
* misc/
  * tests.sh
    * [x] misc/tests.sh breaks because records improve (so chance goes down) used fixed gaps.db
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
