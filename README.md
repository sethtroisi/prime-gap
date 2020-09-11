# Combined Sieve - a new program to find prime gaps.

This is the combination of a couple of ideas I had while working on gmp mpz\_prevprime.

### TODO

* [ ] python util.py (config parsing, isPrime, teelogger)
* [ ] Make =SL included in sieve (e.g. change < SL to <= SL)
* [ ] Rename prime-gap.db
* [ ] missing\_gap\_verify.py & missing\_gap\_test.py
  * [ ] --ignore-gaps 130898
  * [ ] --skip-till-m mi (for resuming)
  * [ ] grouped output (every X entries)
  * [ ] tee logging by default (for premptible machines)
  * [ ] record when no primes found (candidate for high merit?)
  * [ ] store records in some .gitignore'd file
  * [x] Multiprocessing
* benchmarking
  * [ ] Try to avoid int128 in `modulo_search_euclid` (used by `_all`)
  * [ ] Make sure `max_a` is working as intended.
  * [ ] Add instructions to verify `modulo_search` is >80% of the time.
  * [x] Add int64 `modulo_search_euclid_all`
  * [x] Add benchmarks for int32/int64 `modulo_search`
  * [x] Add benchmarks for `K % p`
* gap\_test.py
  * [ ] Starting at m > mstart
  * [ ] Plot average tests count
  * [ ] Sort by expected gap and PRP only top X%
  * [x] Store all results to sql
  * [x] Autoscale printing to every X seconds
  * [x] Describe distribution
  * [x] Generate expected length
* gap\_stats.cpp
  * [ ] Move missing gaps behind a compile flag
    * [ ] Create two seperate Makefile targets
  * [x] drop directory from `unknown-filename`
  * [x] Missing gaps prob % filter
  * [x] Don't rewrite `partial\_results`
  * [x] Calc record chance
  * [x] Save expected/prob record to sql
  * [x] load merit from gap.db
  * [x] Load from unknown\_fn
* gap\_search.cpp
  * [ ] Option to output m with gcd(m, d) != 1
  * [x] (Method2) keeping all composites in memory (and immediately marking mi off)
    * [ ] Look at Method1 vs Method2 writeup and understand why outputs seem different
    * [ ] Make method2 the default
    * [ ] Consider calculating skipped PRP based on index (earlier is 1.0, end of sieve is 0.005)
    * [ ] Ctrl-c then just writes out the results at that point.
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


## Prerequisites

```bash
$ sudo apt install libgmp-dev libmpfr-dev libmpc-dev
$ pip install --user gmpy2=2.1.0b5
$ sudo apt install sqlite3 libsqlite3-dev

# Should probably build and install from head for faster next_prime
```

* OpenPFGW [PrimeWiki](https://www.rieselprime.de/ziki/PFGW) [SourceForge](https://sourceforge.net/projects/openpfgw/)
  * Unpack somewhere and link binary as `pfgw64`
  * may also need to `chmod a+x pfgw64`

* [prime-gap-list](https://github.com/primegap-list-project/prime-gap-list)
        ```bash
        $ git clone https://github.com/primegap-list-project/prime-gap-list.git
        $ cd prime-gap-list
        $ sqlite3 gaps.db < allgaps.sql
        $ cd <prime-gaps>
        $ ln -s ../prime-gap-list/gaps.db .
        ```

## Commands


### Gap Search (Sieve many `m * p#/d`)

#### Method1 vs Method2

`gap_search` parameter `--method2` which changes how the sieve is performed.

```bash
$ make gap_search
$ time ./gap_search [--method2] -p 503 -d 30030 --mstart 1000000 --minc 1000 --save-unknowns --sieve-only --sieve-range 1000
```

**Method1** performs a large initial setup phase which finds the first `mi` that each prime divides before starting any sieving.
Each `mi` is then processed in order, looking at small primes and only the large primes that divide a number in the sieve interval.

**Method2** inverts this process and keeps all sieve intervals in memory at the same time.
Primes are iterated in order and all sieve intervals are handled simultaneously.
Status output is significantly nicer in Method2. Stats are computed at regular intervals giving better insight into `--sieve-range`.

In practice I always suggest Method2, the only downside seems to be slightly more memory. Method1 is mostly used to validate results.

```bash
$ make gap_search
$ time ./gap_search -p 907 -d 210 --mstart 21400000 --minc 5000 --save-unknowns --sieve-only --sieve-range 1000
AUTO SET: sieve length (coprime: 353, prob_gap longer 0.79%): 7746
...
        21400037  127 <- unknowns -> 120
            tests     10         (257.62/sec)  0 seconds elapsed
            unknowns  2353       (avg: 235.30), 98.48% composite  50.36 <- % -> 49.64%
            large prime remaining: 4400229 (avg/test: 7338)
        21400433  120 <- unknowns -> 110
            tests     100        (261.09/sec)  0 seconds elapsed
            unknowns  23530      (avg: 235.30), 98.48% composite  49.88 <- % -> 50.12%
            large prime remaining: 4130914 (avg/test: 7347)
        21402181  121 <- unknowns -> 107
            tests     500        (270.10/sec)  2 seconds elapsed
            unknowns  117692     (avg: 235.38), 98.48% composite  49.81 <- % -> 50.19%
            large prime remaining: 2839408 (avg/test: 7342)
        21404371  119 <- unknowns -> 117
            tests     1000       (288.24/sec)  3 seconds elapsed
            unknowns  235709     (avg: 235.71), 98.48% composite  49.84 <- % -> 50.16%
            large prime remaining: 831741 (avg/test: 7340)

real    0m22.424s
user    0m22.356s
sys     0m0.068s
$
```

```bash
$ make gap_search
$ time ./gap_search -p 907 -d 210 --mstart 21400000 --minc 5000 --save-unknowns --sieve-only --sieve-range 1000 --method2

100,003    (prime_i: 9,592/9,593) (seconds: 0.27/0.3 | per m: 0.00026)
	factors found interval: 19,450,638, total: 28,304,316, avg m/large_prime interval: 0.0
	unknowns 484,174  	(avg/m: 423.60) (delta: 2.734%)
	~ 2x 42.85 PRP/m	(1399490 new composites ~= 764449.9 skipped PRP => 2784403.7 PRP/seconds)

1,000,003  (prime_i: 36,960/78,499) (seconds: 0.29/1.3 | per m: 0.0011)
	factors found interval: 20,674,342, total: 130,158,153, avg m/large_prime interval: 6.5
	unknowns 404,012  	(avg/m: 353.47) (delta: 2.282%)
	~ 2x 35.70 PRP/m	(21207 new composites ~= 2155.6 skipped PRP => 7315.1 PRP/seconds)

999,999,937 (prime_i: 24,491,666/50,847,534) (seconds: 8.26/20.9 | per m: 0.018)
	factors found interval: 20,592,889, total: 438,561,307, avg m/large_prime interval: 0.0
	unknowns 269,510  	(avg/m: 235.79) (delta: 1.522%)
	~ 2x 23.80 PRP/m	(9367 new composites ~= 941.5 skipped PRP => 114.0 PRP/seconds)
```

```
$ diff 21400000_907_210_5000_s7746_l1000M.m2.txt 21400000_907_210_5000_s7746_l1000M.txt
<empty>
```

### Gap Stats

```bash
$ make gap_stats
$ time ./gap_stats --save-unknowns --sieve-only --unknown-filename <M_P_D_MINC_sX_lSL.txt>
```


### Benchmark

80%+ of the time is spent in `modulo_search_euclid` so it's important to optimize this function.

**TODO**: Updated after `max_a` fix in benchmark
```bash
$ make benchmark
$ ./benchmark 100000
$ ./benchmark 100000 "# mod"
...
	|  bits x count   | method_name                    | found    | total    | time(s) | ns/iter | cycles/limb |
	|   694 x  100000 | 503# mod <40 bit>p             | 100000   | 6088     | 0.0045  |      45 | 14.0        |
	|  1390 x  100000 | 1009# mod <40 bit>p            | 100000   | 3292     | 0.0061  |      61 | 9.3         |
	|  2799 x  100000 | 1999# mod <40 bit>p            | 100000   | 8627     | 0.0086  |      86 | 6.6         |
	|  7099 x  100000 | 5003# mod <40 bit>p            | 100000   | 3010     | 0.0146  |     146 | 4.5         |
	| 14291 x  100000 | 10007# mod <40 bit>p           | 100000   | 2119     | 0.0258  |     258 | 3.9         |
	| 28588 x  100000 | 20011# mod <40 bit>p           | 100000   | 2854     | 0.0477  |     477 | 3.6         |

$ ./benchmark 100000 modulo_search
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
	|  bits x count   | method_name                    | found    | total    | time(s) | ns/iter | cycles/iter |
	|    40 x  100000 | modulo_search_verify           | 100000   | 100000   | 0.0289  |     289 | 981.2       |
	|    40 x  100000 | modulo_search_euclid           | 100000   | 100000   | 0.0089  |      89 | 301.3       |
	|    40 x  100000 | modulo_search_euclid_gcd       | 100000   | 100000   | 0.0108  |     108 | 365.1       |
	|    40 x  100000 | modulo_search_euclid_gcd2      | 100000   | 100000   | 0.0093  |      93 | 315.8       |
	|    40 x  100000 | modulo_search_euclid_all_small | 10       | 102484   | 0.0110  |     108 | 365.7       |
	|    40 x  100000 | modulo_search_euclid_all_large | 10       | 102484   | 0.0110  |     107 | 364.6       |
```

### Quick test of all functions
```bash
$ PARAMS="-p 907 -d 2190 --mstart 1 --minc 100 --sieve-range 100 --sieve-length 11000"
$ make gap_search gap_stats gap_test
$ time ./gap_search --method2 --save-unknowns --sieve-only $PARAMS
$ time ./gap_search           --save-unknowns --sieve-only $PARAMS
$ diff 1_907_2190_100_s11000_l100M.m2.txt 1_907_2190_100_s11000_l100M.txt > /dev/null && echo "Good" || echo "Bad!";
$ ./gap_stats --sieve-only --unknown-filename 1_907_2190_100_s11000_l100M.m2.txt
# Verify RECORD @ 1,7,13
# Verify "avg seen prob:   0.9994067
# Verify "avg record prob: 4.44e-06 (max: 7.599e-06)"

# TODO gap_test.cpp, gap_test.py
```
