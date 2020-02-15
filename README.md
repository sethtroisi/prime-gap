### A new program to find prime gaps.
---

This is the combination of a couple of ideas I've had bouncing around in my
head while working on gmp mpz_prevprime.

### TODO

* gap_test.py
  * [ ] starting at m > mstart.
  * [ ] sort and PRP only top expected_length
  * [ ] try and describe distribution
  * [x] generate expected length
* gap_search.cpp
  * [ ] `sieve_range` > 4B
  * [ ] Verify `sieve_length` math with d > 1
  * [x] Dynamic `sieve_length`
  * [x] Dynamic `sieve_range`


## Commands

```bash
$ g++ -Wall -Werror -O3 gap_search.cpp gap_common.cpp -lgmp -o gap_search
# ./gap_search <params>
$ time ./gap_search -p 907 -d 210 --mstart 21400000 --minc 5000 --min-merit 16
...
	21400500  185 <- unknowns -> 212 	2341 <- gap -> 1087
	    tests     501        (15.15/sec)  33 seconds elapsed
	    unknowns  123083     (avg: 245.67), 97.93% composite  49.87 <- % -> 50.13%
	    prp tests 26675      (avg: 53.24)
	    fallback prev_gap 7 (1.4%), next_gap 9 (1.8%)
	    best merit this interval: 14.25 (at m=21400357)
	    large prime remaining: 4212388 (avg/test: 5812)
15618  17.7766  21400751 * 907#/210 -12018 to +3600
	21401000  160 <- unknowns -> 158 	  81 <- gap -> 1093
	    tests     1001       (14.97/sec)  67 seconds elapsed
	    unknowns  246000     (avg: 245.75), 97.93% composite  49.95 <- % -> 50.05%
	    prp tests 54232      (avg: 54.18)
	    fallback prev_gap 16 (1.6%), next_gap 18 (1.8%)
	    best merit this interval: 17.78 (at m=21400751)
	    large prime remaining: 4027538 (avg/test: 5809)
22182  25.2478  21402323 * 907#/210 -4262 to +17920
	21404999   94 <- unknowns -> 69  	3134 <- gap -> 864
	    tests     5000       (15.05/sec)  332 seconds elapsed
	    unknowns  1232988    (avg: 246.60), 97.92% composite  49.97 <- % -> 50.03%
	    prp tests 266060     (avg: 53.21)
	    fallback prev_gap 69 (1.4%), next_gap 85 (1.7%)
	    best merit this interval: 25.25 (at m=21402323)
	    large prime remaining: 0 (avg/test: 5809)

real	5m40.449s
user	5m40.369s
sys	0m0.057s
...


```

### Average unknowns in sieve
```
# Lower side
$ g++ -Wall -Werror -O3 gap_search.cpp gap_common.cpp -lgmp -o gap_search
$ time ./gap_search -p 701 -d 1 --mstart 50000 --minc 10000 --save-unknowns --sieve-range 200
...
	prob_prime_coprime: 0.01718
AUTO SET: sieve length (coprime: 266, prob_gap longer 0.99%): 2690
...
	2.937% of sieve should be unknown (200M) ~= 92
	0.146% of 293 digit numbers are prime
	4.978% of tests should be prime (34.0x speedup)
	~2x20.1=40.2 PRP tests per m
	sieve_length=2690 is insufficient ~0.99% of time
...
	59999  110 <- unknowns -> 103 	 919 <- gap -> 797
	    tests     10000      (38.26/sec)  261 seconds elapsed
	    unknowns  1836020    (avg: 183.60), 96.59% composite  50.02 <- % -> 49.98%
	    prp tests 397833     (avg: 39.78)
	    fallback prev_gap 101 (1.0%), next_gap 99 (1.0%)
	    best merit this interval: 9.94 (at m=57983)
	    large prime remaining: 0 (avg/test: 2829)

$ cat 50000_701_1_10000_s2690_l200M.txt | awk ' { print $7 }' | ministat
    N           Min           Max        Median           Avg        Stddev
x 10000          -120           -61           -92      -91.8368     7.7448411

$ cat 50000_701_1_10000_s2690_l200M.txt | awk ' { print $8 }' | ministat
    N           Min           Max        Median           Avg        Stddev
x 10000            66           122            92       91.7652      7.728379

$ cat 50000_701_1_10000_s2690_l200M.txt | cut -d' ' -f 10- | sed 's/ /\n/g' | sort | uniq -c
   3316 +1619
   3324 -1811
   3339 -2141
   3343 -937
   3346 +881
   3351 -1777
   3554 +2659
   3557 -1367
   3558 -1621
   3561 +1597
   3577 -2531
   3581 +2153
```


## Benchmark

```bash
# Record timing and stats from sieving first 100,000 multiples after 100M
$ time ./gap_search -p 953 -d 210 --mstart 100000000 --minc 100000 --sieve-only --sieve-length 8192

AUTO SET: sieve range (log(t) = ~928): 100000000

Testing m * 953#/210, m = 100000000 + [0, 100000)

sieve_length: 2x8192
sieve_range:  100000000

K = 1313 bits, 396 digits, log(K) = 909.57
Min Gap ~= 9279 (for merit > 10.0)

	PrimePi(100000000) = 5761455 (2 ... 99999989)
	Using 7837 primes for SIEVE_SMALL(80000)
...
	Calculating remainder and inverse for each prime
	Calculating first m each large prime divides
	............................................
	Sum of m1: 19437733456
	Setup took 3.8 seconds

	Starting m=100000000

	100000000  217 <- unknowns -> 231 	   0 <- gap -> 0
	100000001  139 <- unknowns -> 130 	   0 <- gap -> 0
	100000010  216 <- unknowns -> 239 	   0 <- gap -> 0
	100000100  205 <- unknowns -> 221 	   0 <- gap -> 0
	    tests     101        (623.55/sec)  0 seconds elapsed
	    unknowns  37572      (avg: 372.00), 97.73% composite  49.77 <- % -> 50.23%
	    large prime remaining: 5718594 (avg/test: 8009)
	100000500  265 <- unknowns -> 286 	   0 <- gap -> 0
	    tests     501        (612.98/sec)  1 seconds elapsed
	    unknowns  185859     (avg: 370.98), 97.74% composite  49.84 <- % -> 50.16%
	    large prime remaining: 5718465 (avg/test: 8011)
...
	100005000  286 <- unknowns -> 273 	   0 <- gap -> 0
	    tests     5001       (603.42/sec)  8 seconds elapsed
	    unknowns  1848860    (avg: 369.70), 97.74% composite  49.99 <- % -> 50.01%
	    large prime remaining: 5716742 (avg/test: 8014)
...
	100095000  288 <- unknowns -> 284 	   0 <- gap -> 0
	    tests     95001      (607.34/sec)  156 seconds elapsed
	    unknowns  35104716   (avg: 369.52), 97.74% composite  50.00 <- % -> 50.00%
	    large prime remaining: 4806679 (avg/test: 8016)
	100099999  128 <- unknowns -> 148 	   0 <- gap -> 0
	    tests     100000     (608.59/sec)  164 seconds elapsed
	    unknowns  36950332   (avg: 369.50), 97.74% composite  50.00 <- % -> 50.00%
	    large prime remaining: 0 (avg/test: 8016)

real	2m48.625s
user	2m48.561s
sys	0m0.045s
```
