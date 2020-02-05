### A new program to find prime gaps.
---

This is the combination of a couple of ideas I've had bouncing around in my
head while working on gmp mpz_prevprime.

## Commands

```bash
$ g++ -Wall -O3 gap_search.cpp -lgmp -o gap_search
# ./gap_search <params>
$ time ./gap_search 200000000 1000000 887 210 20

```

