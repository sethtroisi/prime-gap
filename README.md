# Combined Sieve - a new program to find prime gaps.

A fast prime gap searching suite

# Table of Contents

- [Combined Sieve - a new program to find prime gaps.](#combined-sieve---a-new-program-to-find-prime-gaps)
  * [Tools](#tools)
    + [Gap Search (Sieve many `m * p#/d`)](#gap-search-sieve-many-m--pd)
      - [Method1 vs Method2](#method1-vs-method2)
    + [Gap Stats](#gap-stats)
    + [Gap Test](#gap-test)
    + [Benchmark](#benchmark)
  * [Benchmarks](#benchmarks)
  * [Setup](#setup)
  * [Notes](#notes)
    + [Quick test of all functions](#quick-test-of-all-functions)
    + [TODO](#todo)
    + [TODONE](#todone)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>


## Tools

### Gap Search (Sieve many `m * p#/d`)

#### Method1 vs Method2

`gap_search` parameter `--method2` changes how the sieve is performed.

```bash
$ make gap_search
$ time ./gap_search [--method2] -p <p> -d <d> --mstart <m_start> --minc <m_inc> --sieve-range <SR>M [--sieve-length <SR>] --save-unknowns
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
$ time ./gap_search -p 907 -d 210 --mstart 21400000 --minc 1000 --sieve-range 1000 --save-unknowns
AUTO SET: sieve length (coprime: 353, prob_gap longer 0.79%): 7746

	Calculating first m each prime divides
	......................................
	PrimePi(1000000000) = 50847534 guessed 50185140
	Setup took 16.8 seconds
...
	21400037  127 <- unknowns -> 120
	    valid m   10         (316.95/sec, with setup per m: 1.7)  0 seconds elapsed
	    unknowns  2353       (avg: 235.30), 98.48% composite  50.36 <- % -> 49.64%
	    large prime remaining: 1196058 (avg/test: 7338)
	21400433  120 <- unknowns -> 110
	    valid m   100        (333.38/sec, with setup per m: 0.17)  0 seconds elapsed
	    unknowns  23530      (avg: 235.30), 98.48% composite  49.88 <- % -> 50.12%
	    large prime remaining: 768403 (avg/test: 7347)
	21400999  142 <- unknowns -> 110
	    valid m   230        (369.18/sec, with setup per m: 0.076)  1 seconds elapsed
	    unknowns  54176      (avg: 235.55), 98.48% composite  49.91 <- % -> 50.09%
	    large prime remaining: 0 (avg/test: 7349)
```

```bash
$ make gap_search
$ time ./gap_search -p 907 -d 210 --mstart 21400000 --minc 1000 --sieve-range 1000 --save-unknowns --method2
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
$ diff 21400000_907_210_1000_s7746_l1000M.txt 21400000_907_210_1000_s7746_l1000M.m2.txt
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

**TODO**: Updated after `max_a` fix in benchmark
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
```

## Benchmarks

See [BENCHMARKS.md](BENCHMARKS.md)

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
  * `misc/double_check.py` double checks using `ecm`, if factor found number shouldn't appear in unknown-file.txt
    * `python misc/double_check.py --unknown-filename <unknown_filen> -c 10`

### Quick test of all functions
```bash
$ PARAMS="-p 907 -d 2190 --mstart 1 --minc 200 --sieve-range 100 --sieve-length 11000"
$ make gap_search gap_stats gap_test
$ time ./gap_search --method2 --save-unknowns $PARAMS
$ time ./gap_search           --save-unknowns $PARAMS
$ md5sum 1_907_2190_200_s11000_l100M.{txt,m2.txt}
080309453b4310e0310a4fb4d1779ffe  1_907_2190_200_s11000_l100M.txt
080309453b4310e0310a4fb4d1779ffe  1_907_2190_200_s11000_l100M.m2.txt

$ ./gap_stats --unknown-filename 1_907_2190_200_s11000_l100M.m2.txt
# Verify RECORD @ 1,7,13,101,137
# Verify "avg seen prob:   0.9994038
# Verify "avg record prob: 4.54e-06 (max: 8.123e-06)"

# TODO gap_test.cpp, gap_test.py
```

### TODO

* [ ] Make =SL included in sieve (e.g. change < SL to <= SL)
* [ ] rename `next_p_i` to `next_p` in schema.sql
* [ ] Rename `gap_search` to `combined_sieve`
* [ ] config.verbose in gap\_search, gap\_stats, gap\_test
* gap\_utils.py
  * [ ] Option for gmp vs OpenPFGW
* README.md
  * [ ] gap\_test.py vs gap\_test.cpp
* benchmarking
  * [ ] Try to avoid int128 in `modulo_search_euclid` (used by `_all`)
  * [ ] Make sure `max_a` is working as intended.
  * [ ] Add instructions to verify `modulo_search` is >80% of the time.
  * [ ] Redo prime-time with user time (because pfgw uses more than 1 thread)
* gap\_test.py
  * [ ] Rename composite = [[], []] to unknowns
  * [ ] Option to starting at m > mstart
  * [ ] Plot Prob(record)
  * [ ] Plot P(gap > min\_merit) and P(record) sorted and unsorted.
  * [ ] Plot average tests count
  * [ ] Sort by expected gap and PRP only top X%
  * [ ] Option to toggle between OpenPFGW and gmp
* gap\_stats.cpp
  * [ ] Move missing gaps behind a compile flag
    * [ ] Create two seperate Makefile targets
* gap\_search.cpp
  * [ ] Option to output m with gcd(m, d) != 1
  * (Method2)
    * [ ] Make method2 the default
    * [ ] Estimate PRP/s and include in status.
    * [ ] Ctrl-c then just writes out the results at that point.
* prime-gap-search.db
  * A plan for how to clean up [partially] finished ranges

### TODONE

* README.md
  * [x] Add commands for benchmarking
  * [x] Fill out gap test section
  * [x] Split out some benchmarking
* Project level
  * [x] Rename prime-gap.db to prime-gap-search.db
* gap\_search.cpp
  * Method2 (all composites in memory)
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
