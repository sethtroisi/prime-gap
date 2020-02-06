### A new program to find prime gaps.
---

This is the combination of a couple of ideas I've had bouncing around in my
head while working on gmp mpz_prevprime.

## Commands

```bash
$ g++ -Wall -Werror -O3 gap_search.cpp -lgmp -o gap_search
# ./gap_search <params>
$ time ./gap_search 200000000 1000000 887 210 20

```

## Benchmark

```bash
# Record timing and stats from first 100,000 multiples after 100M
$ time ./gap_search 100000000 100000 953 210 20
Testing m * 953#/210, m = 100000000 + [0, 100000)
K = 1313 bits, 396 digits, log(K) = 909.57
Min Gap ~= 18559 (for merit > 20.0)

	PrimePi(30000000) = 1857859 (2 ... 29999999)
	Using 5133 primes for SIEVE_SMALL(50000)

	Calculating modulos and inverses
	Calculating prime steps
	.......................
	Sum of m1:   2283262310

	100000000  225 <- unknowns ->  248	   0 <- gap ->    0
	    tests     1          (358.02/sec)  0 seconds elapsed
	    unknowns  473        (avg: 473.00), 97.11% composite  47.57 <- % -> 52.43%
	    prp tests 0          (avg: 0.00)
	    fallback prev_gap 0, next_gap 0
	    best merit this interval: 0.00 (at m=0)
	    large prime queue size: 1849351 (avg/test: 7567)
...
	100001000  231 <- unknowns ->  253	   0 <- gap ->    0
	    tests     1001       (373.86/sec)  3 seconds elapsed
	    unknowns  396912     (avg: 396.52), 97.58% composite  49.89 <- % -> 50.11%
	    prp tests 0          (avg: 0.00)
	    fallback prev_gap 0, next_gap 0
	    best merit this interval: 0.00 (at m=-1)
	    large prime queue size: 1849321 (avg/test: 7600)
...
	100099999  133 <- unknowns ->  156	   0 <- gap ->    0
	    tests     100000     (368.73/sec)  271 seconds elapsed
	    unknowns  39534663   (avg: 395.35), 97.59% composite  50.00 <- % -> 50.00%
	    prp tests 0          (avg: 0.00)
	    fallback prev_gap 0, next_gap 0
	    best merit this interval: 0.00 (at m=-1)
	    large prime queue size: 0 (avg/test: 7604)
```
