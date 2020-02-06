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
#include <chrono>
#include <cstdio>
#include <iostream>
#include <queue>
#include <utility>
#include <vector>

#include <gmp.h>

using std::cout;
using std::endl;
using std::vector;
using namespace std::chrono;


#define SKIP_PRP        0

// Aim for 99% of gaps smaller?
#define SIEVE_LENGTH    8'192

// TODO determine which is fastest
// Dynamically set smaller if M_inc is tiny
//#define SIEVE_RANGE   300'000'000
//#define SIEVE_RANGE   100'000'000
#define SIEVE_RANGE    30'000'000
//#define SIEVE_RANGE    10'000'000
//#define SIEVE_RANGE     3'000'000
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

    printf("Compiled with GMP %d.%d.%d\n\n",
        __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);
    printf("SIEVE_LENGTH: 2x%d, SIEVE_RANGE: %d\n\n", SIEVE_LENGTH, SIEVE_RANGE);
    printf("Testing m * %ld#/%ld, m = %ld + [0, %ld)\n", p, d, m, m_inc);

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
        const vector<int>& primes, int *remainder, bool composite[2][SIEVE_LENGTH]) {

    // For small primes that we don't do trick things with.
    for (int pi = 0; pi < SIEVE_SMALL_PRIME_PI; pi++) {
        const int prime = primes[pi];
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
            composite[0][d] = false;
        }
        int first_negative = modulo == 0 ? 0 : -(modulo - prime);
        assert( 0 <= first_negative && first_negative < prime );
        for (int d = first_negative; d < SIEVE_LENGTH; d += prime) {
            composite[1][d] = false;
        }
    }
}


void prime_gap_search(long M, long M_inc, int P, int D, float min_merit) {
    // ----- Merit STuff
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
    printf("Min Gap ~= %d (for merit > %.1f)\n\n",
        (int) (min_merit * (K_log + m_log)), min_merit);


    // ----- Generate primes under SIEVE_RANGE.
    const vector<int> primes = get_sieve_primes();
    printf("\tPrimePi(%d) = %ld (2 ... %d)\n",
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


    // ----- Allocate memory for a handful of utility functions.

    // Remainders of (p#/d) % prime
    int* remainder   = (int*) malloc(sizeof(int) * primes.size());
    int* rem_inverse = (int*) malloc(sizeof(int) * primes.size());

    {
        cout << "\tCalculating modulos and inverses" << endl;

        mpz_t temp, m_prime;
        mpz_init(temp);
        mpz_init(m_prime);

        for (size_t pi = 0; pi < primes.size(); pi++) {
            const int prime = primes[pi];

            // Big improvement over surround_prime is reusing this for each m.
            long mod = mpz_fdiv_ui(K, prime);
            assert( 0 <= mod && mod < prime );
            remainder[pi] = mod;

            // Measure vs powermod(mod, prime-2, prime);
            if (mod == 0) {
                rem_inverse[pi] = 0;
                assert( prime <= P );
            } else {
                mpz_set_ui(temp, mod);
                mpz_set_ui(m_prime, prime);
                mpz_invert(temp, temp, m_prime);
                long inverse = mpz_get_ui(temp);
                assert( inverse * mod % prime == 1 );
                rem_inverse[pi] = inverse;
            }
        }

        mpz_clear(temp);
        mpz_clear(m_prime);
    }

    // Big improvement over surround_prime is avoiding checking each large prime.
    // pair<next_m, prime> for large primes that only rarely divide a sieve
    typedef std::pair<int, int> mpair;
    std::priority_queue<mpair, vector<mpair>, std::greater<mpair>> next_m;
    {
        // Find next m this will divide
        // solve (base_r * i) % prime < SIEVE_LENGTH
        //  or   (base_r * i) + SIEVE_LENGTH % prime < SIEVE_LENGTH
        // =>
        // solve (SIEVE_LENGTH + base_r * i) % prime < 2 * SIEVE_LENGTH
        //
        // find inverse_r such that inverse_r * base_r = 1 mod prime
        //   distance * inverse * base_r = distance mod prime
        // find distance in [-SIEVE_LENGTH, +SIEVE_LENGTH]
        //   that minimizes (distance * inverse) % prime

        // is it faster to just look search one-by-one for next range?
        //   expected_searches = prime/(2*SL)
        //   (searches are addition and 2 * conditional)
        // otherwise
        //   do exactly 2*SL checks
        //   (check involves modulo
        const long SL = SIEVE_LENGTH;
        const int brute_speedup = 3;
        long search_threshold = brute_speedup * 4L * SL * SL;
        // Potentially can be optimized if M_inc is small
        if (M_inc <= brute_speedup * 2 * SL) {
            search_threshold = SIEVE_RANGE;
        }
        printf("\tCalculating prime steps\n");
        if (search_threshold < SIEVE_RANGE) {
            printf("\tThreshold: %12ld\n", search_threshold);
        }
        // Print "."s during, equal in length to 'Calculat...'
        unsigned int print_dots = 24;

        long first_m_sum = 0;
        cout << "\t";
        for (size_t pi = SIEVE_SMALL_PRIME_PI; pi < primes.size(); pi++) {
            if ((pi * print_dots) % primes.size() < print_dots) {
                cout << "." << std::flush;
            }

            const int prime = primes[pi];
            const int base_r = remainder[pi];
            const int modulo = (base_r * M) % prime;
            if ( (modulo < SL) || (modulo + SL) > prime) {
                next_m.push(std::make_pair(0, pi));
                assert( (modulo + SL) % prime < 2*SL );
                continue;
            }

            if (prime < search_threshold) {
                // just look for next M if prime is small.
                int temp = (modulo + SL - 1) - base_r;
                for (int mi = 0; mi < M_inc; mi++) {
                    temp += base_r;
                    if (temp >= prime) temp -= prime;
                    assert( temp < prime );
                    if (temp < (2*SL-1)) {
                        first_m_sum += mi;

                        assert( (base_r * (M + mi) + (SL - 1)) % prime < 2*SL );
                        next_m.push(std::make_pair(mi, pi));
                        break;
                    }
                }
                continue;
            }

            // This method is slow but takes constant time (w.r.t. prime size)
            const long inverse = rem_inverse[pi];
            long first_m = M_inc+1;

            long distance = prime - (modulo - SL + 1);
            assert (0 < distance && distance < prime);
            long test = ((distance * inverse) % prime) + inverse;
            for (int k = 0; k < 2 * SL - 1; k++) {
                // Need negative distance but C mod negative is weird.
                test -= inverse;
                if (test < 0) test += prime;

                if (test < first_m) {
                    first_m = test;
                    assert( ( base_r * (M + test) + SL) % prime <= 2*SL );
                }
            }

            if (first_m < M_inc) {
                first_m_sum += first_m;
                next_m.push(std::make_pair(first_m, pi));
            }
        }
        cout << endl;
        printf("\tSum of m1: %12ld\n", first_m_sum);
    }


    // ----- Main sieve loop.
    cout << "\n\n\tStarting m=" << M << "\n\n" << endl;
    bool composite[2][SIEVE_LENGTH];

    // Used for various stats
    auto  s_start_t = high_resolution_clock::now();
    long  s_total_unknown = 0;
    long  s_total_prp_tests = 0;
    long  s_gap_out_of_sieve_prev = 0;
    long  s_gap_out_of_sieve_next = 0;
    float s_best_merit_interval = 0;
    long  s_best_merit_interval_m = 0;
    long  s_bad_next_m = 0;

    for (int mi = 0; mi < M_inc; mi++) {
        long m = M + mi;
        // TODO if gcd(m, d) != 1 continue?

        // Reset last sieve to True;
        std::fill_n(*composite, 2*SIEVE_LENGTH, 1);

        sieve_small_primes(
            m, SIEVE_SMALL_PRIME_PI,
            primes, remainder, composite);

        // Maybe useful for some stats later.
        // int unknown_small_u = std::count(composite[1], composite[1]+SIEVE_LENGTH, true);
        // int unknown_small_l = std::count(composite[0], composite[0]+SIEVE_LENGTH, true);

        int s_large_primes_tested = 0;
        // /*
        while (next_m.size() && next_m.top().first == mi) {
            s_large_primes_tested += 1;

            // Large prime should divide some number in SIEVE for this m
            // When done push to next mi.
            int pi = next_m.top().second;
            next_m.pop();
            const int prime   = primes[pi];
            int base_r  = remainder[pi];
            long modulo = (base_r * m) % prime;

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
                // Just past a multiple
                composite[0][modulo] = false;
            } else {
                // Don't have to deal with 0 case anymore.
                int first_negative = -(modulo - prime);
                if (first_negative < SIEVE_LENGTH) {
                    // Just before a multiple
                    composite[1][first_negative] = false;
                } else {
                    s_bad_next_m += 1;
                }
            }

            // TODO this COULD be improved by caching next_m,
            // See 'This method is slow' section.
            // The trick is to record a series of record near m2
            // that chance modulo by <= 2*SIEVE_LENGTH.
            {
                int temp = (modulo + SIEVE_LENGTH - 1 + base_r);
                if (temp >= prime) temp -= prime;
                for (int m2 = mi + 1; m2 < M_inc; m2++) {
                    if (temp >= prime) temp -= prime;
                    if (temp < (2*SIEVE_LENGTH-1)) {
                        next_m.push(std::make_pair(m2, pi));
                        break;
                    }
                    temp += base_r;
                }
            }
        }
        // */

        int unknown_u = std::count(composite[1], composite[1]+SIEVE_LENGTH, true);
        int unknown_l = std::count(composite[0], composite[0]+SIEVE_LENGTH, true);
        s_total_unknown += unknown_u + unknown_l;

        // TODO break out to function, also count tests.
        int prev_p_i = 0;
        int next_p_i = 0;
        if (SKIP_PRP) {
            mpz_t center, ptest;
            mpz_init(center); mpz_init(ptest);
            mpz_mul_ui(center, K, m);

            for (int i = 1; (next_p_i == 0 || prev_p_i == 0) && i < SIEVE_LENGTH; i++) {
                if (prev_p_i == 0 && composite[0][i]) {
                    s_total_prp_tests += 1;

                    mpz_sub_ui(ptest, center, i);
                    if (mpz_probab_prime_p(ptest, 25)) {
                        prev_p_i = i;
                    }
                }
                if (next_p_i == 0 && composite[1][i]) {
                    s_total_prp_tests += 1;

                    mpz_add_ui(ptest, center, i);
                    if (mpz_probab_prime_p(ptest, 25)) {
                        next_p_i = i;
                    }
                }
            }

            if (next_p_i == 0) {
                s_gap_out_of_sieve_next += 1;
                // Using fallback to slower gmp routine
                cout << "\tfalling back to mpz_nextprime" << endl;
                mpz_add_ui(ptest, center, SIEVE_LENGTH-1);
                mpz_nextprime(ptest, ptest);
                mpz_sub(ptest, ptest, center);
                next_p_i = mpz_get_ui(ptest);
            }

            if (prev_p_i == 0) {
                s_gap_out_of_sieve_prev += 1;
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
                    mpz_add(center, center, ptest);
                    mpz_nextprime(ptest, ptest);
                }
            }

            int gap = next_p_i + prev_p_i;
            float merit = gap / (K_log + m_log);
            // TODO parameter or merit.
            if (merit > min_merit)  {
                // TODO write to file.
                printf("%d  %.4f  %ld * %d#/%d - %d to +%d\n",
                    gap, merit, m, P, D, prev_p_i, next_p_i);
            }
            if (merit > s_best_merit_interval) {
                s_best_merit_interval = merit;
                s_best_merit_interval_m = m;
            }

            mpz_clear(center); mpz_clear(ptest);
        }


        if ( (mi == 1 || mi == 100 || mi == 1000) || (m % 5000 == 0) || ((mi+1) == M_inc) ) {
            auto s_stop_t = high_resolution_clock::now();
            double   secs = duration<double>(s_stop_t - s_start_t).count();

            printf("\t%ld %4d, %4d unknown\t%4d, %4d gap\n",
                m,
                unknown_l, unknown_u,
                prev_p_i, next_p_i);
            printf("\t    tests     %-10d (%.1f/sec)  %.0f elapsed\n",
                mi+1, (mi+1) / secs, secs);
            printf("\t    unknowns  %-10ld (avg: %.2f)\n"
                   "\t    prp tests %-10ld (avg: %.2f)\n",
                s_total_unknown, s_total_unknown / (float) (mi+1),
                s_total_prp_tests, s_total_prp_tests / (float) (mi+1));
            printf("\t    fallback prev_gap %ld, next_gap %ld\n",
                s_gap_out_of_sieve_prev, s_gap_out_of_sieve_next);
            printf("\t    best merit (interval): %.2f (at m=%ld)\n",
                s_best_merit_interval, s_best_merit_interval_m);
            printf("\t    large prime queue size: %ld (%d used here) (bad next_m %ld)\n",
                next_m.size(), s_large_primes_tested, s_bad_next_m);

            s_best_merit_interval = 0;
            s_best_merit_interval_m = -1;
        }
    }


    // ----- cleanup

    free(remainder);
    free(rem_inverse);
    mpz_clear(K);
}


/*

// Old code for large prime search ( ~20x slower ).
        // SIEVE_SMALL has dealt with small primes.
        for (int pi = SIEVE_SMALL_PRIME_PI; pi < primes.size(); pi++) {
            int prime = primes[pi];
            long modulo = (remainder[pi] * m) % prime;

            if (modulo < SIEVE_LENGTH) {
                composite[1][modulo] = false;
            } else {
                // Don't have to deal with 0 case anymore.
                int first_negative = -(modulo - prime);
                if (first_negative < SIEVE_LENGTH) {
                    composite[0][first_negative] = false;
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

// Old code for finding first range with mult of large p
// This is 15-20x slower because of many uses of modulo.
            const long inverse = rem_inverse[pi];
            long first_m = M_inc+1;
            // Have to do a little extra work because M doesn't start at zero.
            for (int distance = modulo - SL + 1; distance < modulo + SL; distance++) {
                // Need negative distance but C mod negative is weird.
                long test = ((-distance % prime) + prime) % prime;
                test = (test * inverse) % prime;
                if (test < first_m) {
                    first_m = test;
                    assert( ( base_r * (M + test) + SL) % prime <= 2*SL );
                }
            }

*/

