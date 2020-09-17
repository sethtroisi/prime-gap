# Table of Contents

- [Theory](#Theory)
  * [`gap_search`](#gap_search)
    * [`--sieve_length`](#sieve_length)
    * [Skipped PRP Tests](#skipped-prp-tests)
  * [`gap_stats`](#gap_stats)
  * [Out of order testing](#out-of-order-testing)


# Theory

Information and experimental validation that doesn't fit in [README.md](README.md)

## `gap_search`

### `--sieve_length`

`gap_search` autosets `--sieve_length` (SL) if not passed.

This section abbreviates gcd(a, b) = c as (a, b) = c

1. Finds SL such that `prob_prime ** |{i | (K, i) = 1, -SL < i < SL}|` is < 0.008
  * `prob_prime ** goal_coprime < 0.008` ⇔ `goal_coprime > log(0.008, prob_prime)`
1. `prob_prime` is improved be adjusted for sharing no factor with K
  1. `prob_prime ≈ (1 / log(m * K)) / Prod(1 - 1/p, p <= P)`
1. Over all `0 <= mi < M_inc` find
  ```python
      min(len(i for i in range(-SL, SL) if math.gcd(m * K + i, K) == 1)
          for m in range(M_start, M_start+M_inc) if math.gcd(m, D) == 1)
  ```
  1. Because D often includes many small factors can improve accuracy by also handling factor of D
    1. **Approxiamated** by replacing `m` with `m % D`
    1. Check `math.gcd((m % D) * K + i, K*D) == 1`
1. Check if `min(...) >= goal_coprime`
  1. Expensive to calculate so use these optimizations
    * Cache `len(i for i in range(-SL+1, SL) if (m * K + i, K*D) == 1)`
    * `((m % D) * K + i, K*D) > 1` => `(i, K) > 1 or ((m % D) * (K % D) + i, D) > 1`
    * if `(i, K) > 1` advance to next (all `m * K + i` will be coprime)
    * All (m % D) have same count, so only consider m with (m, d) == 1


### Skipped PRP Tests

Given `m * K` and a prime limit how many PRP tests to find the next/previous prime?

* `prob_prime(X) ≈ 1 / log(X)`
  * This result is due to the [Prime Number Theorum](https://en.wikipedia.org/wiki/Prime_number_theorem)
  * This approximates `prob_prime(random.randint(2, X))`, not `prob_prime(X)`
* A Better estimate is given by `prob_prime(X) ≈ (log(X) - 1) / log(X)^2`
  * `prob_prime(X) ≈ 1/log(X) - 1/log(X)^2`

Accounting for no small prime factors

* Let `prob_p = 1/log(m * K) - 1/log(m * K)^2` be the probability with no information
* Adjust for all the primes we know don't divide it `prob_p / Prod(1 - 1/p, p <= P_limit)`
  * Mertens' third theorem gives very tight approx of `Prod(1 - 1/p, p <= P_limit)`
* `prob_prime ≈ prob_p / (1  / (log(P_limit) * exp(GAMMA)))`
* `prob_prime ≈ prob_p * log(P_limit * exp(GAMMA))`

How many expected tests to find the next prime? one over the probability of prime.

* `E(test) = 1 / prob_prime`

As `P_limit` increases calculate how many fewer expected tests will be performed.

```python
def ExpectedTests(test, m * K, old_limit, new_limit):
        GAMMA = 0.57721566
        prob_p = 1 / log(m * K) - 1 / log(m * K) ** 2
        # 1 / (prob_p * log(old_limit) * exp(GAMMA)) - 1 / (prob_p * log(new_limit) * exp(GAMMA)
        return 1 / (prob_p * math.exp(GAMMA)) * (1/math.log(old_limit) - 1/math.log(new_limit))
```


#### Experimental evidence

Sieve a range to multiple depths (100M, 200M, ... 3B, 4B)

```bash
$ make gap_search DEFINES=-DSAVE_INCREMENTS
$ time ./gap_search --unknown-filename 1_1009_210_2000_s7000_l4000M.txt --save-unknowns

100,000,007 (primes 2,760,321/5,761,456)	(seconds: 0.99/3.1 | per m: 0.0067)
	~ 2x 29.22 PRP/m		(~ 1046.7 skipped PRP => 1013.7 PRP/seconds)

200,000,033 (primes 5,317,482/11,078,938)	(seconds: 1.81/4.9 | per m: 0.011)
	~ 2x 28.16 PRP/m		(~ 970.8 skipped PRP => 512.7 PRP/seconds)

300,000,007 (primes 5,173,388/16,252,326)	(seconds: 1.86/7.0 | per m: 0.015)
	~ 2x 27.58 PRP/m		(~ 535.9 skipped PRP => 287.9 PRP/seconds)

400,000,009 (primes 5,084,001/21,336,327)	(seconds: 1.82/8.8 | per m: 0.019)
	~ 2x 27.18 PRP/m		(~ 366.9 skipped PRP => 201.3 PRP/seconds)

500,000,003 (primes 5,019,541/26,355,868)	(seconds: 1.83/10.6 | per m: 0.023)
	~ 2x 26.88 PRP/m		(~ 277.4 skipped PRP => 151.9 PRP/seconds)

1,000,000,007 (primes 24,491,667/50,847,535)	(seconds: 9.14/19.8 | per m: 0.043)
	~ 2x 25.98 PRP/m		(~ 823.4 skipped PRP => 90.1 PRP/seconds)

2,000,000,011 (primes 47,374,753/98,222,288)	(seconds: 18.06/37.8 | per m: 0.083)
	~ 2x 25.14 PRP/m		(~ 770.1 skipped PRP => 42.7 PRP/seconds)

3,000,000,019 (primes 46,227,250/144,449,538)	(seconds: 18.53/56.4 | per m: 0.12)
	~ 2x 24.67 PRP/m		(~ 427.8 skipped PRP => 23.1 PRP/seconds)

3,999,999,979 (primes 45,512,274/189,961,812)	(seconds: 19.18/75.5 | per m: 0.16)
	~ 2x 24.35 PRP/m		(~ 294.0 skipped PRP => 15.3 PRP/seconds)
```

Then testing each interval seperatly

```bash
for fn in `ls -tr 1_1009*`; do
  echo -e "\n\nProcessing $fn";
  /usr/bin/time -f "\nReal\t%E" ./gap_test --run-prp --unknown-filename "$fn" -qq;
done
```

Can compute how many real PRP tests were skipped by each increase in `--sieve-range`

| `--sieve-range` | PRP needed | delta from last | Skipped PRP estimate | error |
|-----------------|------------|-----------------|----------------------|-------|
| 100M            | 25699      | N/A             | N/A                  | N/A   |
| 200M            | 24775      | 924             | 971                  |  5.1% |
| 300M            | 24262      | 513             | 536                  |  4.5% |
| 400M            | 23904      | 358             | 367                  |  2.5% |
| 500M            | 23642      | 262             | 277                  |  5.7% |
| 1000M           | 22856      | 786             | 823                  |  4.7% |
| 2000M           | 22050      | 806             | 770                  | -4.5% |
| 3000M           | 21645      | 405             | 427                  |  5.4% |
| 4000M           | 21388      | 257             | 294                  | 14.4% |


## `gap_stats`

TODO


## Out of order testing

Some thoughts:
        Should store all primes into db? file?L
        For slow tests I should store composites too

        For in progress maybe store in db?
                when would I come back and delete?

        Maybe mark row with tested-all-prob-record-gaps
                tested X pairs that could have been a record gap non were
                but later which are record gaps could have changed
