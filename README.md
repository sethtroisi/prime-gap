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
  * [Theory](#theory) | [THEORY.md](THEORY.md)
  * [Setup](#setup)
  * [Notes](#notes)
    + [Quick test of all functions](#quick-test-of-all-functions)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>

## Overview

This is a suite of tools (`combined_sieve`, `gap_stats`, `gap_test.py` and `gap_test_gpu`) which are useful for searching for prime gaps
centered around `m * P#/d`.

The flow is to sieve many intervals with `combined_sieve --save -u <filename>`

Then calculate statistics about these intervals with `gap_stats --save -u <filename>`

Then test some portion of these intervals with `gap_test.py -u <filename>.txt --prp-top-percent 25`

To quickly test a known 30 merit prime gap ([Setup](#setup) must already be complete)

```bash
$ make
$ sqlite3 test.db < schema.sql
$ FN=809_1841070_11244000_1000_s15000_l200M.txt
$ ./combined_sieve --save --search-db test.db -u $FN
$ ./gap_stats --save --search-db test.db -u $FN --min-merit 20
$ ./gap_test.py --search-db test.db -u $FN
...
23142  30.1485  11244911 * 809#/1841070 -5788 to +17354	RECORD!
...
```

see [EXTENDED.md](for labourious notes and details)

## Tools

See detailed instructions on each `combined_sieve`, `gap_stats`, `gap_test`, and more in
[EXTENDED.md tools](EXTENDED.md#tools).

## Flow

1. `combined_sieve`
   * Takes a set of parameters (`[m_start, m_start + m_inc), P#, d, sieve-length, max-prime`)
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

See [EXTENDED.md benchmark](EXTENDED.md#benchmark) and slightly out of date [BENCHMARKS.md](BENCHMARKS.md)

## Theory

Theory and justification for some calculations in present in [THEORY.md](THEORY.md)

## Setup

```bash
$ sudo apt install libgmp10 libgmp-dev
$ sudo apt install build-essential automake autoconf bison make libtool texinfo m4

# Building gmp from source is required until 6.3 (not needed in Ubuntu 24.04)
# Mercurial is not working in Ubuntu 24.04 so from source
# $ hg clone https://gmplib.org/repo/gmp/
# $ cd gmp

# $ wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
# $ tar -xJf gmp-6.3.0.tar.xz
# $ cd gmp-6.3.0

# $ ./.bootstrap
# $ mkdir build
# $ cd build
# $ ../configure
# $ make
# $ make check
# $ make install
```

```
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

* `gap_test_gpu` is experimental using https://github.com/NVlabs/CGBN
  * It's **much** faster for <2,048 bits.
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
