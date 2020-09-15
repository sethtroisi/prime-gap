# Combined Sieve - a new program to find prime gaps.

A fast prime gap searching suite

# Table of Contents

- [Combined Sieve - a new program to find prime gaps.](#combined-sieve---a-new-program-to-find-prime-gaps)
  * [Tools](#tools)
    + [Gap Search (Sieve many `m * p#/d`)](#gap-search-sieve-many-m--pd)
      - [Method1 vs Method2 (Default)](#method1-vs-method2-default)
    + [Gap Stats](#gap-stats)
    + [Gap Test](#gap-test)
    + [Benchmark](#benchmark)
  * [Benchmarks](#benchmarks) | [BENCHMARKS.md](BENCHMARKS.md)
  * [Theory](#theory) | [THEORY.md](THEORY.md)
  * [Setup](#setup)
  * [Notes](#notes)
    + [Quick test of all functions](#quick-test-of-all-functions)
    + [TODO](#todo)
    + [TODONE](#todone)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>


## Tools

### Gap Search (Sieve many `m * p#/d`)

#### Method1 vs Method2 (Default)

`gap_search` parameter `--method1` changes how the sieve is performed.

```bash
$ make gap_search
$ time ./gap_search [--method1] -p <p> -d <d> --mstart <m_start> --minc <m_inc> --sieve-range <SR>M [--sieve-length <SR>] --save-unknowns
```

**Method1** performs a large initial setup phase to find the first `mi` that each `prime` divides before starting any sieving.
Each `mi` is then processed in order, looking at small primes and only those large primes that divide a number in the sieve interval.

Half the time will be spent watching dots cross this screen.
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

Status output is significantly nicer in Method2. Stats are computed at regular intervals giving better insight into `--sieve-range`.

In practice Method2 is better. Method1 is only used to validate results.

```bash
$ make gap_search
$ time ./gap_search -p 907 -d 210 --mstart 21400000 --minc 1000 --sieve-range 1000 --save-unknowns --method1
AUTO SET: sieve length (coprime: 353, prob_gap longer 0.79%): 7746

	Calculating first m each prime divides
	......................................
	PrimePi(1000000000) = 50847534 guessed 50185140
	Setup took 16.8 seconds
...
	21400003  102 <- unknowns -> 122
	    intervals 1          (326.87/sec, with setup per m: 17)  0 seconds elapsed
	    unknowns  224        (avg: 224.00), 98.54% composite  45.54 <- % -> 54.46%
	    large prime remaining: 1226493 (avg/test: 7470)
	21400037  127 <- unknowns -> 119
	    intervals 10         (332.96/sec, with setup per m: 1.7)  0 seconds elapsed
	    unknowns  2328       (avg: 232.80), 98.48% composite  50.47 <- % -> 49.53%
	    large prime remaining: 1186700 (avg/test: 7272)
	21400433  119 <- unknowns -> 108
	    intervals 100        (349.03/sec, with setup per m: 0.17)  0 seconds elapsed
	    unknowns  23193      (avg: 231.93), 98.49% composite  49.90 <- % -> 50.10%
	    large prime remaining: 762206 (avg/test: 7279)
	21400999  140 <- unknowns -> 107
	    intervals 230        (398.46/sec, with setup per m: 0.074)  1 seconds elapsed
	    unknowns  53383      (avg: 232.10), 98.49% composite  49.92 <- % -> 50.08%
	    large prime remaining: 0 (avg/test: 7280)
```

```bash
$ make gap_search
$ time ./gap_search -p 907 -d 210 --mstart 21400000 --minc 1000 --sieve-range 1000 --save-unknowns
AUTO SET: sieve length (coprime: 353, prob_gap longer 0.79%): 7746

Testing m * 907#/210, m = 21400000 + [0, 1,000)

coprime m:   230/1000, coprime i 2058/7746
sieve_range:  1,000,000,000   small_threshold:  387,300

2          (primes 1/1)	(seconds: 0.00/0.0 | per m: 1.6e-05)
	factors  1,781,580 		(interval: 1,781,580, avg m/large_prime interval: 0.0)
	unknowns   379,040/230  	(avg/m: 1648.00) (composite: 89.36% +89.362%)
	~ 2x 711.66 PRP/m		(3183890 new composites ~= 38390.6 skipped PRP => 10536271.7 PRP/seconds)

100,003    (primes 9,592/9,593)	(seconds: 0.06/0.1 | per m: 0.00028)
	factors  3,913,349 		(interval: 2,131,769, avg m/large_prime interval: 0.0)
	unknowns    97,222/230  	(avg/m: 422.70) (composite: 97.27% +7.910%)
	~ 2x 42.85 PRP/m		(281818 new composites ~= 153826.3 skipped PRP => 2497147.8 PRP/seconds)

1,000,003  (primes 36,960/78,499)	(seconds: 0.07/0.3 | per m: 0.0013)
	factors  4,178,373 		(interval: 48,639, avg m/large_prime interval: 1.3)
	unknowns    81,125/230  	(avg/m: 352.72) (composite: 97.72% +0.119%)
	~ 2x 35.70 PRP/m		(4243 new composites ~= 433.8 skipped PRP => 6632.5 PRP/seconds)

...

999,999,937 (primes 24,491,666/50,847,534)	(seconds: 8.40/17.6 | per m: 0.076)
	factors  4,562,430 		(interval: 32,262, avg m/large_prime interval: 0.0)
	unknowns    54,176/230  	(avg/m: 235.55) (composite: 98.48% +0.052%)
	~ 2x 23.80 PRP/m		(1856 new composites ~= 189.5 skipped PRP => 22.6 PRP/seconds)

```

```bash
$ diff 21400000_907_210_1000_s7674_l1000M.txt 21400000_907_210_1000_s7674_l1000M.m1.txt
<empty>
```

### Gap Stats

```bash
$ make gap_stats
$ time ./gap_stats --save-unknowns --unknown-filename "<M>_<P>_<D>_<MINC>_s<SL>_l<SR>M.txt"
```

### Gap Test

TODO

### Benchmark

80%+ of the time is spent in `modulo_search_euclid` so it's important to optimize this function.

More benchmarks are present in [BENCHMARKS.md](BENCHMARKS.md)

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
	|    25 x  100000 | modulo_search_one_op           | 1033     | 100000   | 0.0020  |      20 | 68.8        |
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

## Benchmarks

See [BENCHMARKS.md](BENCHMARKS.md)

## Theory

Theory and justifaction for some calculations in present in [THEORY.md](THEORY.md)

## Setup

```bash
# Should probably build and install gmp head for faster next_prime
$ sudo apt install libgmp10

$ sudo apt install libgmp-dev libmpfr-dev libmpc-dev
$ sudo apt install sqlite3 libsqlite3-dev
$ pip install --user gmpy2=2.1.0b5 tqdm

$ sqlite3 prime-gap-search.db < schema.sql

```

* [prime-gap-list](https://github.com/primegap-list-project/prime-gap-list)
        ```bash
        $ git clone https://github.com/primegap-list-project/prime-gap-list.git
        $ cd prime-gap-list
        $ sqlite3 gaps.db < allgaps.sql
        $ cd <prime-gaps>
        $ ln -s ../prime-gap-list/gaps.db .
        ```

## Notes

* Anecdotally `gap_search` is ~8% faster when using `CC=clang++-11`
* There is some support for using OpenPFGW [PrimeWiki](https://www.rieselprime.de/ziki/PFGW) [SourceForge](https://sourceforge.net/projects/openpfgw/)
  * Unpack somewhere and link binary as `pfgw64`
  * may also need to `chmod a+x pfgw64`
  * modify `is_prime` in gap\_utils.py
* Multiple layers of verification of `gap_search`
  * Can compare `--method1` result with `--method2`
  * Can add `-DGMP_VALIDATE_FACTORS=1` to `CFLAGS` in Makefile
  * `misc/double_check.py` double checks using `ecm`, if small factor found number shouldn't appear in unknown-file.txt
    * `python misc/double_check.py --unknown-filename <unknown_filen> -c 10`
  * skipped PRP is checked with
```bash
$ make gap_search DEFINES=-DSAVE_INCREMENTS
# Saves a 100m, 200m, 300m, 400m, 500m, 1B, 2B, 3B, 4B, 5B
$ time ./gap_search --unknown-filename 1_1999_210_2000_s4000_l5000M.txt --save-unknowns
100,000,007 (primes 2,760,321/5,761,456)	(seconds: 0.90/2.4 | per m: 0.0053)
	~ 2x 58.96 PRP/m		(~ 1055.9 skipped PRP => 2357.2 PRP/seconds)

200,000,033 (primes 5,317,482/11,078,938)	(seconds: 1.73/4.2 | per m: 0.0091)
	~ 2x 56.82 PRP/m		(~ 979.3 skipped PRP => 1130.0 PRP/seconds)

300,000,007 (primes 5,173,388/16,252,326)	(seconds: 1.70/5.9 | per m: 0.013)
	~ 2x 55.64 PRP/m		(~ 540.6 skipped PRP => 637.2 PRP/seconds)

400,000,009 (primes 5,084,001/21,336,327)	(seconds: 1.69/7.6 | per m: 0.017)
	~ 2x 54.84 PRP/m		(~ 370.2 skipped PRP => 437.2 PRP/seconds)

500,000,003 (primes 5,019,541/26,355,868)	(seconds: 1.67/9.2 | per m: 0.02)
	~ 2x 54.23 PRP/m		(~ 279.8 skipped PRP => 334.3 PRP/seconds)

1,000,000,007 (primes 24,491,667/50,847,535)	(seconds: 8.40/17.6 | per m: 0.039)
	~ 2x 52.41 PRP/m		(~ 830.7 skipped PRP => 197.7 PRP/seconds)

2,000,000,011 (primes 47,374,753/98,222,288)	(seconds: 16.85/34.5 | per m: 0.075)
	~ 2x 50.72 PRP/m		(~ 776.9 skipped PRP => 92.2 PRP/seconds)

3,000,000,019 (primes 46,227,250/144,449,538)	(seconds: 17.60/52.1 | per m: 0.11)
	~ 2x 49.77 PRP/m		(~ 431.6 skipped PRP => 49.1 PRP/seconds)

4,000,000,007 (primes 45,512,275/189,961,813)	(seconds: 17.93/70.0 | per m: 0.15)
	~ 2x 49.13 PRP/m		(~ 296.6 skipped PRP => 33.1 PRP/seconds)

4,999,999,937 (primes 44,992,410/234,954,223)	(seconds: 18.00/88.0 | per m: 0.19)
	~ 2x 48.63 PRP/m		(~ 224.8 skipped PRP => 25.0 PRP/seconds)

real	1m28.037s
user	1m28.030s
sys	0m0.008s

# Then testing each interval seperatly
$ for fn in `ls -tr 1_1999*`; do echo -e "\n\nProcessing $fn"; /usr/bin/time -f "\nReal\t%E" ./gap_test --run-prp --unknown-filename "$fn" -qq; done

Processing 1_1999_210_2000_s4000_l30M.txt
	m=1999   51 <- unknowns -> 57  	 240 <- gap -> 9072
	    tests     458        (1.03 secs/test)  472 seconds elapsed
	    unknowns  48104      (avg: 105.03), 98.69% composite  50.16% <- % -> 49.84%
	    prp tests 32674      (avg: 71.34) (69.2 tests/sec)

Real	7:51.93

Processing 1_1999_210_2000_s4000_l100M.txt
	m=1999   51 <- unknowns -> 51  	 240 <- gap -> 9072
	    tests     458        (1.01/sec)  454 seconds elapsed
	    unknowns  44958      (avg: 98.16), 98.77% composite  50.18% <- % -> 49.82%
	    prp tests 30585      (avg: 66.78) (67.4 tests/sec)

Real	7:33.93

Processing 1_1999_210_2000_s4000_l200M.txt
	m=1999   49 <- unknowns -> 51  	 240 <- gap -> 9072
	    tests     458        (1.01/sec)  455 seconds elapsed
	    unknowns  43290      (avg: 94.52), 98.82% composite  50.22% <- % -> 49.78%
	    prp tests 29477      (avg: 64.36) (64.7 tests/sec)

Real	7:35.42

Processing 1_1999_210_2000_s4000_l300M.txt
	m=1999   45 <- unknowns -> 48  	 240 <- gap -> 9072
	    tests     458        (1.02/sec)  450 seconds elapsed
	    unknowns  42362      (avg: 92.49), 98.84% composite  50.18% <- % -> 49.82%
	    prp tests 28858      (avg: 63.01) (64.1 tests/sec)

Real	7:30.48

Processing 1_1999_210_2000_s4000_l400M.txt
	m=1999   45 <- unknowns -> 47  	 240 <- gap -> 9072
	    tests     458        (1.03/sec)  445 seconds elapsed
	    unknowns  41709      (avg: 91.07), 98.86% composite  50.20% <- % -> 49.80%
	    prp tests 28417      (avg: 62.05) (63.8 tests/sec)

Real	7:25.21

Processing 1_1999_210_2000_s4000_l1000M.txt
	m=1999   42 <- unknowns -> 44  	 240 <- gap -> 9072
	    tests     458        (1.05/sec)  434 seconds elapsed
	    unknowns  39892      (avg: 87.10), 98.91% composite  50.07% <- % -> 49.93%
	    prp tests 27199      (avg: 59.39) (62.6 tests/sec)

Real	7:14.16

Processing 1_1999_210_2000_s4000_l2000M.txt
	m=1999   41 <- unknowns -> 42  	 240 <- gap -> 9072
	    tests     458        (1.07/sec)  426 seconds elapsed
	    unknowns  38597      (avg: 84.27), 98.95% composite  50.07% <- % -> 49.93%
	    prp tests 26310      (avg: 57.45) (61.7 tests/sec)

Real	7:06.31

Processing 1_1999_210_2000_s4000_l3000M.txt
	m=1999   40 <- unknowns -> 42  	 240 <- gap -> 9072
	    tests     458        (1.08/sec)  423 seconds elapsed
	    unknowns  37857      (avg: 82.66), 98.97% composite  50.12% <- % -> 49.88%
	    prp tests 25794      (avg: 56.32) (61.0 tests/sec)

Real	7:03.12

Processing 1_1999_210_2000_s4000_l4000M.txt
	m=1999   40 <- unknowns -> 41  	 240 <- gap -> 9072
	    tests     458        (1.08/sec)  423 seconds elapsed
	    unknowns  37359      (avg: 81.57), 98.98% composite  50.11% <- % -> 49.89%
	    prp tests 25459      (avg: 55.59) (60.2 tests/sec)

Real	7:02.94

Processing 1_1999_210_2000_s4000_l5000M.txt
	m=1999   40 <- unknowns -> 40  	 240 <- gap -> 9072
	    tests     458        (1.08/sec)  424 seconds elapsed
	    unknowns  36992      (avg: 80.77), 98.99% composite  50.12% <- % -> 49.88%
	    prp tests 25202      (avg: 55.03) (59.5 tests/sec)

Real	7:03.66

```


* Dev GMPlib | GMP 6.2.99
  * GMP 6.2.0 hasn't yet accepted my `mpz_prevprime` patch
    * `hg apply <patch>` from https://gmplib.org/list-archives/gmp-devel/2020-August/005851.html
    * If you are a developer consider asking telling them that `mpz_prevprime` would be useful

### Quick test of all functions
```bash
$ PARAMS="-p 907 -d 2190 --mstart 1 --minc 200 --sieve-range 100 --sieve-length 11000"
$ make gap_search gap_stats gap_test
$ time ./gap_search --method1 --save-unknowns $PARAMS
$ time ./gap_search           --save-unknowns $PARAMS
$ md5sum 1_907_2190_200_s11000_l100M.{txt,m1.txt}
080309453b4310e0310a4fb4d1779ffe  1_907_2190_200_s11000_l100M.txt
080309453b4310e0310a4fb4d1779ffe  1_907_2190_200_s11000_l100M.m1.txt

$ ./gap_stats --unknown-filename 1_907_2190_200_s11000_l100M.txt
# Verify RECORD @ 1,7,13,101,137
# Verify "avg seen prob:   0.9994038
# Verify "avg record prob: 4.54e-06 (max: 8.123e-06)"

# TODO gap_test.cpp, gap_test.py
```

### TODO

* [ ] Make =SL included in sieve (e.g. change < SL to <= SL)
* [ ] rename `next_p_i` to `next_p` in schema.sql
* [ ] Rename `gap_search` to `combined_sieve`
* gap\_utils.py
* README.md
  * [ ] gap\_test.py vs gap\_test.cpp
* gap\_search.cpp
  * [ ] Option to output m with gcd(m, d) != 1
  * (Method2)
    * [ ] test dynamically choosing `vector<bool>` vs `vector<char>
* gap\_test.py
  * [ ] Rename composite = [[], []] to unknowns
  * [ ] Option to starting at m > mstart
  * [ ] Plot Prob(record)
  * [ ] Plot P(gap > min\_merit) and P(record) sorted and unsorted.
  * [ ] Plot average tests count
  * [ ] Sort by expected gap and PRP only top X%
  * [ ] Option to toggle between OpenPFGW and gmp
* benchmarking
  * [ ] Add instructions to verify `modulo\_search` is >80% of the time.
  * [ ] Redo prime-time with user time (because pfgw uses more than 1 thread)
* gap\_stats.cpp
  * [ ] Tweak logging at different verbose levels
* prime-gap-search.db
  * A plan for how to clean up [partially] finished ranges

### TODONE

* README.md
  * [x] Add commands for benchmarking
  * [x] Fill out gap test section
  * [x] Split out some benchmarking
* Project level
  * [x] Rename prime-gap.db to prime-gap-search.db
  * [x] Make method2 the default
  * [x] config.verbose in gap\_search, gap\_stats, gap\_test
* gap\_search.cpp
  * Method2 (all composites in memory)
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
  * [x] `sieve_range` > 4B
  * [x] Dynamic `sieve_length`
  * [x] Dynamic `sieve_range`
* gap\_stats.cpp
  * [x] Move missing gaps behind a compile flag
  * [x] drop directory from `unknown-filename`
  * [x] Missing gaps prob % filter
  * [x] Don't rewrite `partial\_results`
  * [x] Calc record chance
  * [x] Save expected/prob record to sql
  * [x] load merit from gap.db
  * [x] Load from unknown\_fn
* gap\_test.py
  * [x] Store all results to sql
  * [x] Autoscale printing to every X seconds
  * [x] Describe distribution
  * [x] Generate expected length
* missing\_gap\_verify.py & missing\_gap\_test.py
  * [x] load/save from DB file.
  * [x] grouped output of BOTH PRIME (every X entries)
  * [x] record stats about prime found (potential candidate for high merit?)
  * [x] tee logging (for preemptible machines)
  * [x] gap\_utils.py (config parsing, isPrime, teelogger)
  * [x] --ignore-gaps 130898
  * [x] store records in some .gitignore'd file
  * [x] Multiprocessing
* benchmarking
  * [x] Add int64 `modulo_search_euclid_all`
  * [x] Add benchmarks for int32/int64 `modulo_search`
  * [x] Add benchmarks for `K % p`
* record\_check.py
  * [x] Read from sql db
* double\_check.py
  * [x] run ecm on random unknowns and verify factors found > sieve limit
