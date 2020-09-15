# Table of Contents

- [Theory](#Theory)
  * [`gap_search`](#gap_search)
    * [`--sieve_length`](#sieve_length)
    * [Skipped PRP Tests](#skipped-prp-tests)
  * [`gap_stats`](#gap_stats)

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

