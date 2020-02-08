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
#include <getopt.h>
#include <iostream>
#include <vector>

#include <gmp.h>

using std::cout;
using std::endl;
using std::vector;
using namespace std::chrono;


#define SKIP_PRP        1

// Aim for ~98% of gaps short
//#define SIEVE_LENGTH    24'676
//#define SIEVE_LENGTH    16'384
#define SIEVE_LENGTH    8'192
//#define SIEVE_LENGTH    7'040
//#define SIEVE_LENGTH    1'040

// For SIEVE_LENGTH=8192
//  SIEVE_RANGE = 1000M, SIEVE_SMALL=60K seems to work best

// TODO determine which is fastest
// TODO dynamically set (or accept from CLI)
//#define SIEVE_RANGE   2'000'000'000
#define SIEVE_RANGE   1'000'000'000
//#define SIEVE_RANGE     300'000'000
//#define SIEVE_RANGE     100'000'000
//#define SIEVE_RANGE      30'000'000
//#define SIEVE_RANGE      20'000'000
//#define SIEVE_RANGE      10'000'000
//#define SIEVE_RANGE       3'000'000
//#define SIEVE_RANGE       1'000'000

#define SIEVE_SMALL       80'000

#define MAX_INT     ((1L << 32) - 1)


void show_usage(char* name);
struct Config {
    int valid;
    int mstart;
    int minc;
    int p;
    int d;
    float minmerit;
};
Config argparse(int argc, char* argv[]);
void prime_gap_search(const struct Config config);

Config argparse(int argc, char* argv[]) {
    // TODO add NOPRP option.
    // TODO except SIEVE_LENGTH, SIEVE_RANGE

    static struct option long_options[] = {
        {"help",     no_argument,       0,  'h' },
        {"mstart",   required_argument, 0,   1  },
        {"minc",     required_argument, 0,   2  },
        {"p",        required_argument, 0,  'p' },
        {"d",        required_argument, 0,  'd' },
        {"minmerit", required_argument, 0,   3  },
        {0,         0,                 0,  0 }
    };

    Config config;
    config.valid = 1;

    int option_index = 0;
    char c;
    while ((c = getopt_long(argc, argv, "hp:d:", long_options, &option_index)) >= 0) {
        switch (c) {
            case 'h':
                show_usage(argv[0]);
                exit(0);
            case 'p':
                config.p = atoi(optarg);
                break;
            case 'd':
                config.d = atoi(optarg);
                break;
            case 1:
                config.mstart = atoi(optarg);
                break;
            case 2:
                config.minc = atoi(optarg);
                break;
            case 3:
                config.minmerit = atof(optarg);
                break;

            case 0:
                printf("option %s arg %s\n", long_options[option_index].name, optarg);
                config.valid = 0;
                break;
            case '?':
                break;
            default:
                printf("getopt returned \"%d\"\n", c);
        }
    }

    if (optind < argc) {
        config.valid = 0;
        printf("unknown positional arguements: ");
        while (optind < argc) {
            printf("%s ", argv[optind++]);
        }
        printf("\n");
    }

    // ----- Validation

    if (((long) config.mstart + config.minc) >= MAX_INT) {
        config.valid = 0;
        cout << "mstart + minc must be < 1e9" << endl;
    }

    if (config.minc >= 50'000'000) {
        config.valid = 0;
        cout << "minc > 50M will use to much memory" << endl;
    }

    {
        mpz_t ptest;
        mpz_init_set_ui(ptest, config.p);
        bool valid = 0 != mpz_probab_prime_p(ptest, 25);
        mpz_clear(ptest);
        if (!valid) {
            config.valid = 0;
            cout << "p# not prime (p=" << config.p << ")" << endl;
        }
    }

    if (config.d <= 0) {
        config.valid = 0;
        cout << "d must be greater than 0: " << config.d << endl;
    }

    if (config.minc <= 0) {
        config.valid = 0;
        cout << "minc must be greater than 0: " << config.minc << endl;
    }

    return config;
}


int main(int argc, char* argv[]) {
    Config config = argparse(argc, argv);
    if (config.valid == 0) {
        show_usage(argv[0]);
        return 1;
    }

    printf("Compiled with GMP %d.%d.%d\n\n",
        __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);
    printf("SIEVE_LENGTH: 2x%d, SIEVE_RANGE: %d\n\n", SIEVE_LENGTH, SIEVE_RANGE);

    printf("Testing m * %d#/%d, m = %d + [0, %d)\n",
        config.p, config.d, config.mstart, config.minc);

    // Some mod logic below demands this for now.
    assert( SIEVE_RANGE < MAX_INT );

    prime_gap_search(config);
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

            for (int m = p2; m <= SIEVE_RANGE; m += 2*p)
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
        int first_negative = modulo == 0 ? 0 : prime - modulo;
//        assert( 0 <= first_negative && first_negative < prime );
        for (int d = first_negative; d < SIEVE_LENGTH; d += prime) {
            composite[1][d] = false;
        }
    }
}

int modulo_search_euclid(int p, int A, int L, int R);
int modulo_search(int p, int A, int L, int R) {
    assert( R - L == 2 * SIEVE_LENGTH - 2 );

    // if expect a hit within 16 just brute force.
    if (16 * SIEVE_LENGTH > p) {
        int temp = 0;
        for (int i = 1; i <= 20; i++) {
            temp += A;
            if (temp >= p) temp -= p;
            if (L <= A && A <= R) {
                return i;
            }
        }
    }

    return modulo_search_euclid(p, A, L, R);
}


int modulo_search_euclid(int p, int A, int L, int R) {
    // min i : L <= (A * i) % p <= L
/*
    assert( 0 <= A && A < p );
    assert( 0 <= L && L < p );
    assert( 0 <= R && R < p );
    assert(      L <= R     );
*/

    if (L == 0) return 0;

    if (2L * A > p) {
        std::swap(L, R);
        A = p - A;
        L = p - L;
        R = p - R;
    }

    // check if small i works
    if (A <= R) {
        long mult = ((long) L + A-1L) / A;
        long test = mult * A;
        if (test <= R) {
            // assert( mult >= 0 );
            return mult;
        }
    }

    // Reduce to simplier problem
    long new_a = ((long) A) + ((-p) % A);
    assert( 0 <= new_a && new_a < A );
    long k = modulo_search_euclid(A, new_a, L % A, R % A);

    long tL = L + p * k;
    long mult = (tL + A-1L) / A;
/*
    long tR = R + p * k;
    long test = mult * A;
    assert(       test <= tR );
    assert( tL <= test );
*/
    return mult;

    if (0) {
        long temp = 0;
        for (int i = 1; i < p; i++) {
            temp += A;
            if (temp > p) temp -= p;
            if (L <= temp && temp <= R) {
                cout << p << " " << A << " " << L << " " << R << " | " << mult << " " << i << endl;
                assert( mult == i );
                return i;
            }
        }
        assert( false );
    }
}


void prime_gap_search(const struct Config config) {
    const long M = config.mstart;
    const long M_inc = config.minc;
    const long P = config.p;
    const long D = config.d;
    const float min_merit = config.minmerit;

    // ----- Merit STuff
    mpz_t K;
    mpz_init(K);
    mpz_primorial_ui(K, P);
    assert( 0 == mpz_tdiv_q_ui(K, K, D) );

    int K_digits = mpz_sizeinbase(K, 10);
    float K_log;
    {
        long exp;
        double mantis = mpz_get_d_2exp(&exp, K);
        K_log = log(mantis) + log(2) * exp;
        float m_log = log(M);
        int K_bits   = mpz_sizeinbase(K, 2);

        printf("K = %d bits, %d digits, log(K) = %.2f\n",
            K_bits, K_digits, K_log);
        printf("Min Gap ~= %d (for merit > %.1f)\n\n",
            (int) (min_merit * (K_log + m_log)), min_merit);
    }


    // ----- Generate primes under SIEVE_RANGE.
    vector<int> const primes = get_sieve_primes();
    printf("\tPrimePi(%d) = %ld (2 ... %d)\n",
        SIEVE_RANGE, primes.size(), primes.back());

    // SIEVE_SMALL deals with all primes can mark of two items in SIEVE_LENGTH.
    assert( SIEVE_SMALL > 2 * SIEVE_LENGTH );
    const long SIEVE_SMALL_PRIME_PI = std::distance(primes.begin(),
        std::lower_bound(primes.begin(), primes.end(), SIEVE_SMALL));
    printf("\tUsing %ld primes for SIEVE_SMALL(%d)\n\n",
        SIEVE_SMALL_PRIME_PI, SIEVE_SMALL);
    assert( primes[SIEVE_SMALL_PRIME_PI] > SIEVE_SMALL );


    // ----- Sieve stats
    {
        double unknowns_after_sieve = 1;
        for (long prime : primes) unknowns_after_sieve *= (prime - 1.0) / prime;

        double prob_prime = 1 / K_log;
        double prob_prime_after_sieve = prob_prime / unknowns_after_sieve;

        double prob_gap_shorter_hypothetical =
            1 - pow(1 - prob_prime, SIEVE_LENGTH);

        // In pratice we get slightly better sieving than expected.
        // TODO refine this constant
        double prob_gap_shorter_experimental =
            1 - pow(1 - prob_prime, 0.37 * SIEVE_LENGTH);

        printf("\t%.3f%% of sieve should be unknown (%dM)\n",
            100 * unknowns_after_sieve, SIEVE_RANGE/1'000'000);
        printf("\t%.3f%% of %d digit numbers are prime\n",
            100 * prob_prime, K_digits);
        printf("\t%.3f%% of tests should be prime (%.1fx speedup)\n",
            100 * prob_prime_after_sieve, 1 / unknowns_after_sieve);
        printf("\t~%.1f PRP tests per m\n",
            2 / prob_prime_after_sieve);
        printf("\tSIEVE_LENGTH=%d is sufficient %.2f%%, experimentally: ~%.2f%%\n",
            SIEVE_LENGTH,
            100 * prob_gap_shorter_hypothetical, 100 * prob_gap_shorter_experimental);
        cout << endl;
    }


    // ----- Allocate memory for a handful of utility functions.

    // Remainders of (p#/d) mod prime
    int *remainder   = (int*) malloc(sizeof(int) * primes.size());
    // Inverse of Remainder mod prime
    int *rem_inverse = (int*) malloc(sizeof(int) * primes.size());
    {
        cout << "\tCalculating remainder and inverse for each prime" << endl;

        mpz_t temp, m_prime;
        mpz_init(temp);
        mpz_init(m_prime);

        for (size_t pi = 0; pi < primes.size(); pi++) {
            const long prime = primes[pi];

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
    // vector<m, vector<pi>> for large primes that only rarely divide a sieve
    int s_large_primes_rem = 0;
    std::vector<int> *large_prime_queue = new vector<int>[M_inc];
    {
        // Find next m this will divide
        // solve (base_r * i) % prime < SIEVE_LENGTH
        //  or   (base_r * i) + SIEVE_LENGTH % prime < SIEVE_LENGTH
        // =>
        // solve (SIEVE_LENGTH + base_r * i) % prime < 2 * SIEVE_LENGTH
        //
        const long SL = SIEVE_LENGTH;
        printf("\tCalculating first m each large prime divides\n");

        // Print "."s during, equal in length to 'Calculat...'
        unsigned int print_dots = 44 + 1;

        long first_m_sum = 0;
        cout << "\t";
        for (size_t pi = SIEVE_SMALL_PRIME_PI; pi < primes.size(); pi++) {
            if ((pi * print_dots) % primes.size() < print_dots) {
                cout << "." << std::flush;
            }

            const long prime = primes[pi];
            const long base_r = remainder[pi];
            const long modulo = (base_r * M) % prime;
            if ( (modulo < SL) || (modulo + SL) > prime) {
                large_prime_queue[0].push_back(pi);
                s_large_primes_rem += 1;
                assert( (modulo + SL) % prime < 2*SL );
                continue;
            }

            // solve base_r * (M + mi) + (SL - 1)) % prime < 2 * SL
            //   0 <= (base_r * M + SL - 1) + base_r * mi < 2 * SL mod prime
            //
            // let shift = (base_r * M + SL - 1) % prime
            //   0 <= shift + base_r * mi < 2 * SL mod prime
            // add (prime - shift) to all three
            //
            //  (prime - shift) <= base_r * mi < (prime - shift) + 2 * SL mod prime

            long shift = modulo + (SL - 1);
            assert( 0 <= shift && shift < prime );
            long low  = (prime - shift);
            long high = low + (2*SL-2);
            assert( 0 <= low && high < prime );

            long mi = modulo_search(prime, base_r, low, high);
            assert( low <= (mi * base_r) % prime );
            assert(        (mi * base_r) % prime <= high );

            assert( (base_r * (M + mi) + (SL - 1)) % prime < (2*SL-1) );
            if (mi < M_inc) {
                large_prime_queue[mi].push_back(pi);
                s_large_primes_rem += 1;
                first_m_sum += mi;
            }

            if (0) {
                // Brute force doublecheck
                long temp = (modulo + SL - 1) - base_r;
                int mi2 = - 1;
                for (mi2 = 0; mi2 < M_inc; mi2++) {
                    temp += base_r;
                    if (temp >= prime) temp -= prime;
                    assert( temp < prime );
                    if (temp < (2*SL-1)) {
                        first_m_sum += mi2;

                        assert( (base_r * (M + mi2) + (SL - 1)) % prime < 2*SL );
                        cout << prime << " " << base_r << " " << modulo << " | " << mi << " " << mi2 << endl;
                        assert( mi2 == mi );
                        break;
                    }
                }
            }
        }
        cout << endl;
        printf("\tSum of m1: %ld\n", first_m_sum);
    }


    // ----- Main sieve loop.
    cout << "\n\tStarting m=" << M << "\n" << endl;
    bool composite[2][SIEVE_LENGTH];

    // Used for various stats
    auto  s_start_t = high_resolution_clock::now();
    long  s_total_unknown = 0;
    long  s_t_unk_low = 0;
    long  s_t_unk_hgh = 0;
    long  s_total_prp_tests = 0;
    long  s_gap_out_of_sieve_prev = 0;
    long  s_gap_out_of_sieve_next = 0;
    float s_best_merit_interval = 0;
    long  s_best_merit_interval_m = 0;
    long  s_large_primes_tested = 0;

    for (int mi = 0; mi < M_inc; mi++) {
        long m = M + mi;
        // TODO if gcd(m, d) != 1 continue?

        // Reset last sieve to True;
        std::fill_n(*composite, 2*SIEVE_LENGTH, 1);

        sieve_small_primes(
            m, SIEVE_SMALL_PRIME_PI,
            primes, remainder, composite);

        // Maybe useful for some stats later.
        // int unknown_small_l = std::count(composite[0], composite[0]+SIEVE_LENGTH, true);
        // int unknown_small_u = std::count(composite[1], composite[1]+SIEVE_LENGTH, true);

        for (int pi : large_prime_queue[mi]) {
            s_large_primes_tested += 1;
            s_large_primes_rem -= 1;

            // Large prime should divide some number in SIEVE for this m
            // When done find next mi prime divides.
            const long prime   = primes[pi];
            long base_r  = remainder[pi];
            long modulo = (base_r * m) % prime;

            if (0) {
                mpz_t test; mpz_init(test); mpz_mul_ui(test, K, m);
                long mod = mpz_fdiv_ui(test, prime);
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
                long first_negative = prime - modulo;

                assert( first_negative < SIEVE_LENGTH); // Bad next m!
                //cout << "Bad next m: " << m << " " << prime << " mod: " << modulo << endl;

                // Just before a multiple
                composite[1][first_negative] = false;
            }

            // Find next mi that primes divides part of SIEVE
            {
                const long SL = SIEVE_LENGTH;
                // next modulo otherwise m = 0;
                long shift = (modulo + base_r) + (SL - 1);
                if (shift >= prime) shift -= prime;
                if (shift >= prime) shift -= prime;

                assert( 0 <= shift && shift < prime );

                long low  = (prime - shift);
                long high = low + (2*SL-2);

                int next_mi = M_inc;
                if (high >= prime) {
                    next_mi = mi + 1;
                } else {
                    assert( 0 <= low && high < prime );
                    long m2 = modulo_search(prime, base_r, low, high);
                    assert( low <= (m2 * base_r) % prime );
                    assert(        (m2 * base_r) % prime <= high );
                    next_mi = mi + 1 + m2;
                }

                assert( (base_r * (M + next_mi) + (SL - 1)) % prime < (2*SL-1) );
                if (next_mi < M_inc) {
                    large_prime_queue[next_mi].push_back(pi);
                    s_large_primes_rem += 1;
                }

                if (0) {
                    long temp = (modulo + base_r) + SL - 1;
                    if (temp >= prime) temp -= prime;
                    for (int m3 = m + 1; m3 < M_inc; m3++) {
                        if (temp >= prime) temp -= prime;
                        if (temp < (2*SL-1)) {
                            assert( next_mi == m3 );
                            break;
                        }
                        temp += base_r;
                    }
                }
            }
        }
        large_prime_queue[mi].clear();
        large_prime_queue[mi].shrink_to_fit();

        int unknown_l = std::count(composite[0], composite[0]+SIEVE_LENGTH, true);
        int unknown_u = std::count(composite[1], composite[1]+SIEVE_LENGTH, true);
        s_total_unknown += unknown_l + unknown_u;
        s_t_unk_low += unknown_l;
        s_t_unk_hgh += unknown_u;

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
                //cout << "\tfalling back to mpz_nextprime" << endl;
                mpz_add_ui(ptest, center, SIEVE_LENGTH-1);
                mpz_nextprime(ptest, ptest);
                mpz_sub(ptest, ptest, center);
                next_p_i = mpz_get_ui(ptest);
            }

            if (prev_p_i == 0) {
                s_gap_out_of_sieve_prev += 1;
                /*
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
                // */

                // /*
                // Medium ugly fallback.
                for (int i = SIEVE_LENGTH+1; ; i++) {
                    bool composite = false;
                    for (int pi = 0; pi < 2000; pi++) {
                        const long prime = primes[pi];
                        long modulo = (remainder[pi] * m) % prime;
                        if (i % prime == modulo) {
                            composite = true;
                            break;
                        }
                    }
                    if (!composite) {
                        mpz_sub_ui(ptest, center, i);
                        if (mpz_probab_prime_p(ptest, 25)) {
                            prev_p_i = i;
                            break;
                        }
                    }
                }
                // */
            }

            int gap = next_p_i + prev_p_i;
            float merit = gap / (K_log + log(m));
            // TODO parameter or merit.
            if (merit > min_merit)  {
                // TODO write to file.
                printf("%d  %.4f  %ld * %ld#/%ld -%d to +%d\n",
                    gap, merit, m, P, D, prev_p_i, next_p_i);
            }
            if (merit > s_best_merit_interval) {
                s_best_merit_interval = merit;
                s_best_merit_interval_m = m;
            }

            mpz_clear(center); mpz_clear(ptest);
        }


        if ( (mi == 1 || mi == 10 || mi == 100 || mi == 1000) ||
             (m % 5000 == 0) || ((mi+1) == M_inc) ) {
            auto s_stop_t = high_resolution_clock::now();
            double   secs = duration<double>(s_stop_t - s_start_t).count();

            printf("\t%ld %4d <- unknowns -> %-4d\t%4d <- gap -> %-4d\n",
                m,
                unknown_l, unknown_u,
                prev_p_i, next_p_i);
            if (mi <= 10) continue;
            // Stats!
            printf("\t    tests     %-10d (%.2f/sec)  %.0f seconds elapsed\n",
                mi+1, (mi+1) / secs, secs);
            printf("\t    unknowns  %-10ld (avg: %.2f), %.2f%% composite  %.2f <- %% -> %.2f%%\n",
                s_total_unknown, s_total_unknown / (float) (mi+1),
                100.0 * (1 - s_total_unknown / (2.0 * (SIEVE_LENGTH - 1) * (mi+1))),
                100.0 * s_t_unk_low / s_total_unknown,
                100.0 * s_t_unk_hgh / s_total_unknown);
            printf("\t    prp tests %-10ld (avg: %.2f)\n",
                s_total_prp_tests, s_total_prp_tests / (float) (mi+1));
            printf("\t    fallback prev_gap %ld, next_gap %ld\n",
                s_gap_out_of_sieve_prev, s_gap_out_of_sieve_next);
            printf("\t    best merit this interval: %.2f (at m=%ld)\n",
                s_best_merit_interval, s_best_merit_interval_m);
            printf("\t    large prime remaining: %d (avg/test: %ld)\n",
                s_large_primes_rem, s_large_primes_tested / (mi+1));

            s_best_merit_interval = 0;
            s_best_merit_interval_m = -1;
        }
    }

    // Should be cleaning up after self.
    for(int mi = 0; mi < M_inc; mi++)  {
        assert( large_prime_queue[mi].size() == 0 );
    }

    // ----- cleanup

    delete[] large_prime_queue;
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

