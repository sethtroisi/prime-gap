// Copyright 2020 Seth Troisi
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <queue>
#include <utility>
#include <vector>

#include <gmp.h>

using std::cout;
using std::endl;
using std::vector;


// 13 = 8192, 14 = 16384
#define SIEVE_BITS      14
#define SIEVE_LENGTH    (1 << SIEVE_BITS)

// TODO determine which is fastest
//#define SIEVE_RANGE   300'000'000
#define SIEVE_RANGE   100'000'000
//#define SIEVE_RANGE    30'000'000
//#define SIEVE_RANGE     1'000'000
#define SIEVE_SMALL        40'000

#define MAX_INT     ((1L << 32) - 1)


void prime_gap_search(long M, long M_inc, int P, int D, float min_merit);
void show_usage(char* name);


int main(int argc, char* argv[]) {
    if (argc != 6) {
        show_usage(argv[0]);
        return 1;
    }

    long m     = atol(argv[1]);
    long m_inc = atol(argv[2]);
    long p     = atol(argv[3]);
    long d     = atol(argv[4]);
    float min_merit = atof(argv[5]);

    printf("GMP %d.%d.%d\n",
        __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);
    printf("SIEVE_LENGTH: 2x%d, SIEVE_RANGE: %d\n\n", SIEVE_LENGTH, SIEVE_RANGE);
    printf("showing results with merit > %.1f\n\n", min_merit);
    printf("Testing m * %ld#/%ld, m = %ld + [0, %ld)\n\n", p, d, m, m_inc);

    // Some mod logic below demands this for now.
    assert( (m + m_inc) < MAX_INT );

    prime_gap_search(m, m_inc, p, d,   min_merit);
}

void show_usage(char* name) {
    cout << "Usage: " << name << " m  m_inc  p  d  min_merit" << endl
         << "\tcalculates prime_gaps for m * p#/d, with merit > min_merit" << endl;
}


vector<int> get_sieve_primes() {
    vector<int> primes = {2};
    vector<bool> is_prime(SIEVE_RANGE+1, true);
    for (int p = 3; p <= SIEVE_RANGE; p += 2) {
        if (is_prime[p]) {
            primes.push_back(p);
            int p2 = p * p;
            if (p2 > SIEVE_RANGE) break;

            for (int m = p2; m <= SIEVE_RANGE; m += p)
                is_prime[m] = false;
        }
    }
    for (int p = primes.back() + 2; p <= SIEVE_RANGE; p += 2) {
        if (is_prime[p])
            primes.push_back(p);
    }
    return primes;
}


inline void sieve_small_primes(
        const long m,  const int SIEVE_SMALL_PRIME_PI,
        const vector<int>& primes, int *remainder, bool is_composite[2][SIEVE_LENGTH]) {

    // For small primes that we don't do trick things with.
    for (int pi = 0; pi < SIEVE_SMALL_PRIME_PI; pi++) {
        int prime = primes[pi];
        long modulo = (remainder[pi] * m) % prime;

//        if (0) {
//            mpz_t test; mpz_init(test); mpz_mul_ui(test, K, m);
//            int mod = mpz_fdiv_ui(test, prime);
//            if (mod != modulo) {
//                cout << m << " " << prime << "\t" << mod << " vs " << modulo << endl;
//                assert( false );
//            }
//            mpz_clear(test);
//        }

        for (int d = modulo; d < SIEVE_LENGTH; d += prime) {
            is_composite[0][d] = false;
        }
        int first_negative = modulo == 0 ? 0 : -(modulo - prime);
        assert( 0 <= first_negative && first_negative < prime );
        for (int d = first_negative; d < SIEVE_LENGTH; d += prime) {
            is_composite[1][d] = false;
        }
    }
}


void prime_gap_search(long M, long M_inc, int P, int D, float min_merit) {
    // ----- Generate primes under SIEVE_RANGE.
    const vector<int> primes = get_sieve_primes();
    printf("PrimePi(%d) = %ld (2 ... %d)\n",
        SIEVE_RANGE, primes.size(), primes.back());

    const int P_pi = std::distance(primes.begin(),
        std::find(primes.begin(), primes.end(), P));
    if (primes[P_pi] != P) {
        cout << "p# not prime (p=" << P << ")" << endl;
        exit(1);
    }

    // SIEVE_SMALL deals with all primes can mark of two items in SIEVE_LENGTH.
    assert( SIEVE_SMALL > 2 * SIEVE_LENGTH );
    const long SIEVE_SMALL_PRIME_PI = std::distance(primes.begin(),
        std::lower_bound(primes.begin(), primes.end(), SIEVE_SMALL));
    printf("\tUsing %ld primes for SIEVE_SMALL(%d)\n\n",
        SIEVE_SMALL_PRIME_PI, SIEVE_SMALL);
    assert( primes[SIEVE_SMALL_PRIME_PI] > SIEVE_SMALL );

    // ----- Setup helper
    mpz_t K;
    mpz_init(K);
    mpz_primorial_ui(K, P);
    assert( 0 == mpz_tdiv_q_ui(K, K, D) );
    int K_bits   = mpz_sizeinbase(K, 2);
    int K_digits = mpz_sizeinbase(K, 10);
    float K_log;
    float m_log = log(M);
    {
        long exp;
        double mantis = mpz_get_d_2exp(&exp, K);
        K_log = log(mantis) + log(2) * exp;
    }
    printf("K = %d bits, %d digits, log(K) = %.2f\n",
        K_bits, K_digits, K_log);
    printf("min gap (for merit > %.1f) ~= %d\n\n",
        min_merit, (int) (min_merit * (K_log + m_log)));

    // ----- Allocate memory for a handful of utility functions.

    // Remainders of (p#/d) % prime
    int* remainder = (int*) malloc(sizeof(int) * primes.size());

    for (size_t pi = 0; pi < primes.size(); pi++) {
        int prime = primes[pi];

        // Genius idea is to only do this once.
        long mod = mpz_fdiv_ui(K, prime);
        assert( 0 <= mod && mod < prime );
        remainder[pi] = mod;
    }
    cout << "\tCalculated all modulos\n\n" << endl;

    bool is_composite[2][SIEVE_LENGTH];

    // pair<next_m, prime> for large primes that only rarely divide a sieve
    typedef std::pair<int, int> mpair;
    std::priority_queue<mpair, vector<mpair>, std::greater<mpair>> next_m;
    {
        // Add everything to the queue for eval on mi=0 to avoid duplicate logic
        for (size_t pi = SIEVE_SMALL_PRIME_PI; pi < primes.size(); pi++) {
            next_m.push(std::make_pair(0, pi));
        }
    }

    // ----- Main sieve loop.

    long total_unknown = 0;
    for (int mi = 0; mi < M_inc; mi++) {
        long m = M + mi;

        // Reset last sieve to True;
        std::fill_n(*is_composite, 2*SIEVE_LENGTH, 1);

        sieve_small_primes(
            m, SIEVE_SMALL_PRIME_PI,
            primes, remainder, is_composite);

        int unknown_small_u = std::count(is_composite[1], is_composite[1]+SIEVE_LENGTH, true);
        int unknown_small_l = std::count(is_composite[0], is_composite[0]+SIEVE_LENGTH, true);

        int tested = 0, valid = 0;
        // /*
        while (next_m.size() && next_m.top().first == mi) {
            // Check if prime divides this, otherwise push to somewhere
            int pi = next_m.top().second;
            next_m.pop();
            int prime   = primes[pi];
            int base_r  = remainder[pi];
            long modulo = (base_r * m) % prime;
            tested += 1;

            if (0) {
                mpz_t test; mpz_init(test); mpz_mul_ui(test, K, m);
                int mod = mpz_fdiv_ui(test, prime);
                if (mod != modulo) {
                    cout << m << " " << prime << "\t" << mod << " vs " << modulo << endl;
                    assert( false );
                }
                mpz_clear(test);
            }

            if (modulo < SIEVE_LENGTH) {
                valid += 1;
                // Just past a multiple
                is_composite[0][modulo] = false;
            } else {
                // Don't have to deal with 0 case anymore.
                int first_negative = -(modulo - prime);
                if (first_negative < SIEVE_LENGTH) {
                    valid += 1;
                    // Just before a multiple
                    is_composite[1][first_negative] = false;
                }
            }

            // Find next m this will divide
            // solve (modulo + base_r * i) % prime < SIEVE_LENGTH
            //  or   (modulo + base_r * i) + SIEVE_LENGTH % prime < SIEVE_LENGTH
            // =>
            // solve (modulo + SIEVE_LENGTH + base_r * i) % prime < 2 * SIEVE_LENGTH
            //
            // Do the easy thing for now FIXME, TODO
            // Some sort of big/small step?
            {
                int temp = (modulo + SIEVE_LENGTH + base_r);
                if (temp >= prime) temp -= prime;

                for (int m2 = mi + 1; m2 < M_inc; m2++) {
                    if (temp >= prime) temp -= prime;
                    if (temp < (2*SIEVE_LENGTH)) {
                        next_m.push(std::make_pair(m2, pi));
                        break;
                    }
                    temp += base_r;
                }
            }
        }
        // */

        int unknown_u = std::count(is_composite[1], is_composite[1]+SIEVE_LENGTH, true);
        int unknown_l = std::count(is_composite[0], is_composite[0]+SIEVE_LENGTH, true);
        total_unknown += unknown_u + unknown_l;

        if ( (m % 1000 == 0) || ((mi+1) == M_inc) ) {

            printf("%ld\t unknown: %4d, %4d  | (total: %ld, avg: %.1f) (pqueue: %ld)\n",
                m,
                unknown_l, unknown_u,
                total_unknown, total_unknown / (float) (mi+1),
                next_m.size()
            );
            // TODO verbose flag
            if (0) {
                printf("\tfast_sieve %4d, %4d | large_sieve tested: %d, %d\n",
                    unknown_small_l, unknown_small_u,
                    tested, valid
                );
            }
        }

        if (1) {
            mpz_t center, ptest;
            mpz_init(center); mpz_init(ptest);
            mpz_mul_ui(center, K, m);

            int prev_p_i = 0;
            int next_p_i = 0;
            for (int i = 1; (next_p_i == 0 || prev_p_i == 0) && i < SIEVE_LENGTH; i++) {
                if (prev_p_i == 0 && is_composite[0][i]) {
                    mpz_sub_ui(ptest, center, i);
                    if (mpz_probab_prime_p(ptest, 25)) {
                        prev_p_i = i;
                    }
                }
                if (next_p_i == 0 && is_composite[1][i]) {
                    mpz_add_ui(ptest, center, i);
                    if (mpz_probab_prime_p(ptest, 25)) {
                        next_p_i = i;
                    }
                }
            }

            if (next_p_i == 0) {
                // Using fallback to slower gmp routine
                mpz_add_ui(ptest, center, SIEVE_LENGTH-1);
                mpz_nextprime(ptest, ptest);
                mpz_sub(ptest, ptest, center);
                next_p_i = mpz_get_ui(ptest);
                cout << "\tfalling back to mpz_nextprime" << endl;
            }

            if (prev_p_i == 0) {
                // REALLY UGLY FALLBACK
                cout << "\tUGLY prevprime hack" << endl;
                mpz_sub_ui(ptest, center, 2*SIEVE_LENGTH-1);
                mpz_nextprime(ptest, ptest);
                if (mpz_cmp(ptest, center) > 1) {
                    cout << m << "What!" << endl;
                    exit(1);
                }

                while (mpz_cmp(ptest, center) < 0) {
                    // save distance
                    mpz_sub(center, center, ptest);
                    prev_p_i = mpz_get_ui(center);
                    mpz_nextprime(ptest, ptest);
                }
            }

            int gap = next_p_i + prev_p_i;
            float merit = gap / (K_log + m_log);
            // TODO parameter or merit.
            if (merit > min_merit)  {
                printf("m: %ld, gap: %5d = %5d + %5d, %.3f\n",
                    m, gap, prev_p_i, next_p_i, merit);
            }

            mpz_clear(center); mpz_clear(ptest);
        }
    }


    // ----- cleanup

    free(remainder);
    mpz_clear(K);
}


/*

1M
31009999	 unknown:  168,  174  | (total: 4661232, avg: 466.1) (pqueue: 0)
real	0m11.390s

30M
31009999	 unknown:  143,  142  | (total: 3740846, avg: 374.1) (pqueue: 0)
real	1m3.345s

100M
31009999	 unknown:  132,  135  | (total: 3496592, avg: 349.7) (pqueue: 0)
real	2m23.030s

300M
31009999	 unknown:  127,  125  | (total: 3299835, avg: 330.0) (pqueue: 0)
real	5m36.869s



*/

// Old code for large prime search ~20x slower.
        // SIEVE_SMALL has dealt with small primes.
        /*
        for (int pi = SIEVE_SMALL_PRIME_PI; pi < primes.size(); pi++) {
            int prime = primes[pi];
            long modulo = (remainder[pi] * m) % prime;

            if (modulo < SIEVE_LENGTH) {
                is_composite[1][modulo] = false;
            } else {
                // Don't have to deal with 0 case anymore.
                int first_negative = -(modulo - prime);
                if (first_negative < SIEVE_LENGTH) {
                    is_composite[0][first_negative] = false;
                }
            }

            if (0) {
                mpz_t test; mpz_init(test); mpz_mul_ui(test, K, m);
                int mod = mpz_fdiv_ui(test, prime);
                if (mod != modulo) {
                    cout << m << " " << prime << "\t" << mod << " vs " << modulo << endl;
                    assert( false );
                }
                mpz_clear(test);
            }
        }
        */

