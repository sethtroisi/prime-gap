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
#include <chrono>
#include <clocale>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include <gmp.h>

#include "gap_common.h"

using std::cout;
using std::endl;
using std::pair;
using std::map;
using std::vector;
using namespace std::chrono;

// Tweaking this doesn't seem to method1 much.
// method2 seems very sensative and is controlled elsewhere.
#define SIEVE_SMALL       200'000

void set_defaults(struct Config& config);
void prime_gap_search(const struct Config config);
void prime_gap_parallel(const struct Config config);


int main(int argc, char* argv[]) {
    printf("\tCompiled with GMP %d.%d.%d\n\n",
        __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);

    Config config = argparse(argc, argv);
    set_defaults(config);

    if (config.save_unknowns == 0) {
        cout << "Must set --save-unknowns" << endl;
        return 1;
    }
    if (config.run_prp == 1) {
        cout << "Must set --sieve-only for gap_search" << endl;
        return 1;
    }

    if (config.valid == 0) {
        show_usage(argv[0]);
        return 1;
    }

    printf("\n");
    printf("Testing m * %u#/%u, m = %ld + [0, %ld)\n",
        config.p, config.d, config.mstart, config.minc);

    printf("\n");
    printf("sieve_length: 2x%d\n", config.sieve_length);
    printf("sieve_range:  %ld\n", config.sieve_range);
    printf("\n");

    if (config.method2) {
        prime_gap_parallel(config);
    } else {
        prime_gap_search(config);
    }
}


void set_defaults(struct Config& config) {
    if (config.valid == 0) {
        // Don't do anything if argparse didn't work.
        return;
    }

    float logK;
    mpz_t K;

    {
        mpz_init(K);
        mpz_primorial_ui(K, config.p);
        assert( 0 == mpz_tdiv_q_ui(K, K, config.d) );
        long exp;
        double mantis = mpz_get_d_2exp(&exp, K);
        logK = log(config.mstart) + log(mantis) + log(2) * exp;
    }

    if (config.sieve_length == 0) {
        double prob_prime_coprime = 1 / (logK + log(config.mstart));
        // factors of K = p#/d
        vector<uint32_t> K_primes = get_sieve_primes(config.p);
        {
            for (size_t pi = K_primes.size()-1; ; pi--) {
                prob_prime_coprime /= 1 - 1.0 / K_primes[pi];
                if (config.d % K_primes[pi] == 0) {
                    K_primes.erase(K_primes.begin() + pi);
                }
                if (pi == 0) break;
            }
        }

        // K = #p/d
        // only numbers K+i has no factor <= p
        //      => (K+i, i) == 1
        //      => only relatively prime i
        //
        // factors of d are hard because they depend on m*K
        //  some of these m are worse than others so use worst m


        assert( config.p >= 503 );

        // Search till chance of shorter gap is small.
        {
            assert( config.d <= 30030 );

            uint32_t base = mpz_fdiv_ui(K, config.d);

            // count of (m*K) % d over all m
            vector<uint32_t> count_by_mod_d(config.d, 0);
            {
                for (uint64_t mi = 0; mi < config.minc; mi++) {
                    uint64_t m = config.mstart + mi;
                    if (gcd(m, config.d) == 1) {
                        uint32_t center = (m * base) % config.d;
                        uint32_t center_down = (config.d - center) % config.d;

                        // distance heading up
                        count_by_mod_d[ center ] += 1;
                        // distance heading up
                        count_by_mod_d[ center_down ] += 1;
                    }
                }
            }

            // Note: By averaging over counts prob_larger could be improve here.
            map<uint32_t, uint32_t> coprime_by_mod_d;
            for (size_t i = 0; i < config.d; i++) {
                if (count_by_mod_d[i] > 0) {
                    coprime_by_mod_d[i] = 0;
                }
            }

            for (size_t tSL = 1; ; tSL += 1) {
                bool all_divisible = false;
                for (int prime : K_primes) {
                    if ((tSL % prime) == 0) {
                        all_divisible = true;
                        break;
                    }
                }
                if (all_divisible) continue;

                // check if tSL is divisible for all center mods
                for (auto& coprime_counts : coprime_by_mod_d) {
                    const auto center = coprime_counts.first;
                    // if some factor in d will mark this off don't count it
                    if (gcd(center + tSL, config.d) == 1) {
                        coprime_counts.second += 1;
                    }
                }

                // Find the smallest number of coprimes
                uint32_t min_coprime = tSL;
                for (auto& coprime_counts : coprime_by_mod_d) {
                    min_coprime = std::min(min_coprime, coprime_counts.second);
                }

                // Assume each coprime is independent
                double prob_gap_shorter = pow(1 - prob_prime_coprime, min_coprime);

                // This seems to balance PRP fallback and sieve_size
                if (prob_gap_shorter <= 0.008) {
                    config.sieve_length = tSL;
                    printf("AUTO SET: sieve length (coprime: %d, prob_gap longer %.2f%%): %ld\n",
                        min_coprime, 100 * prob_gap_shorter, tSL);
                    break;
                }
            }
        }
        assert( config.sieve_length > 100 ); // Something went wrong above.
    }

    if (config.sieve_range == 0) {
        // each additional numbers removes unknowns / prime
        // and takes log2(prime / sieve_length) time

        // TODO improve this.
        // Potentially based on:
        //  sieve_length
        //  min_merit
        if (logK >= 1500) {
            // Largest supported number right now
            config.sieve_range = 4'000'000'000;
            // 2020-02-09 tuning notes
            //  P 1627 SL=8192
            //      log(t) ~= 1600, ~= 4 primes in SL
            //       100M => 620/s | 293 unknowns
            //       200M => 549/s | 282 unknowns
            //       400M => 480/s | 273 unknowns
            //          81.6 PRP / test => 1.70s/test
            //       800M => 440/s | 264 unknowns
            //          78.6 PRP / test => 1.75s/test
            //      1600M => 440/s | 254 unknowns
            //          76.2 PRP / test => 1.78s/test
        } else {
            config.sieve_range =   1'000'000'000;
        }

        printf("AUTO SET: sieve range (log(t) = ~%.0f): %ld\n",
            logK, config.sieve_range);
    }

    mpz_clear(K);
}


uint32_t modulo_search_brute(uint32_t p, uint32_t A, uint32_t L, uint32_t R) {
    // A + p must not overflow.
    assert( A + p > p );

    uint32_t temp = 0;
    for (int i = 1; ; i++) {
        temp += A;
        if (temp >= p) temp -= p;
        if (L <= temp && temp <= R) {
            return i;
        }
    }
}


uint32_t modulo_search_euclid_small(uint32_t p, uint32_t a, uint32_t l, uint32_t r) {
    // min i : l <= (a * i) % p <= l
    if (l == 0) return 0;

    if (a > (p >> 1)) {
        std::swap(l, r);
        a = p - a; l = p - l; r = p - r;
    }

    if (a <= r) {
        uint64_t mult = (l-1) / a + 1;
        uint64_t test = mult * a;
        if (test <= r) {
            return mult;
        }
    }

    // reduce to simplier problem
    uint32_t new_a = a - (p % a);
    assert( 0 <= new_a && new_a < a );
    uint64_t k = modulo_search_euclid_small(a, new_a, l % a, r % a);

    uint64_t tl = k * p + l;
    uint64_t mult = (tl - 1) / a + 1;

    assert( mult < p );
    return mult;
}


uint64_t modulo_search_euclid(uint64_t p, uint64_t a, uint64_t l, uint64_t r) {
    // min i : l <= (a * i) % p <= l
/*
    assert( 0 <= a && a < p );
    assert( 0 <= l && l < p );
    assert( 0 <= r && r < p );
    assert(      l <= r     );
// */
    if (l == 0) return 0;

    if (p < 0xFFFFFFFF) {
        return modulo_search_euclid_small(p, a, l, r);
    }

    // 2 *a > p, but avoids int max
    if (a > (p >> 1)) {
        std::swap(l, r);
        a = p - a; l = p - l; r = p - r;
    }

    // check if small i works
    if (a <= r) {
        // find next multiple of a >= l
        // check if <= r
        uint64_t mult = (l-1) / a + 1;
        uint64_t test = mult * a;
        if (test <= r) {
            return mult;
        }
    }

    // reduce to simplier problem
    uint64_t new_a = a - (p % a);
    assert( 0 <= new_a && new_a < a );
    uint64_t k = modulo_search_euclid(a, new_a, l % a, r % a);

    __int128 tl = (__int128) p * k + l;
    uint64_t mult = (tl-1) / a + 1;

    assert( mult < p );
/*
    __int128 tr = r + p * k;
    uint64_t test = mult * a;
    assert(       test <= tr );
    assert( tl <= test );
// */
    return mult;
}


uint64_t modulo_search_euclid_gcd(
        uint64_t M, uint64_t D, uint64_t max_m, uint64_t SL,
        uint64_t prime, uint64_t base_r) {
    uint64_t mi = 0;

    // TODO validate no overflow
    uint64_t modulo = (base_r * (M + mi)) % prime;
    while (mi < max_m) {
        if ( (modulo < SL) || (modulo + SL) > prime) {
            if (gcd(M + mi, D) > 1) {
                mi += 1;
                modulo += base_r;
                if (modulo >= prime) modulo -= prime;
                continue;
            }
            //assert( mi == z );
            return mi;
        }

        uint64_t shift = modulo + (SL - 1);
        assert( 0 <= shift && shift < prime );
        uint64_t low  = (prime - shift);
        uint64_t high = low + (2*SL-2);
        assert( 0 <= low && high < prime );

        uint64_t mi2 = modulo_search_euclid(prime, base_r, low, high);
        mi += mi2;

        __int128 mult = (__int128) base_r * (M + mi);
        modulo = mult % prime;

        assert( (modulo < SL) || (modulo + SL) > prime );
    }
    return max_m;
}


void modulo_search_euclid_all(
        uint64_t M, uint64_t max_m, uint64_t SL,
        uint64_t prime, uint64_t base_r,
        std::function<void (uint64_t)> lambda) {
    uint64_t mi = 0;

    // TODO validate no overflow
    uint64_t modulo = (base_r * (M + mi)) % prime;
    while (mi < max_m) {
        if ( (modulo < SL) || (modulo + SL) > prime) {
            lambda(mi);
            mi += 1;
            modulo += base_r;
            if (modulo >= prime) modulo -= prime;
            continue;
        }

        uint64_t shift = modulo + (SL - 1);
        assert( 0 <= shift && shift < prime );
        uint64_t low  = (prime - shift);
        uint64_t high = low + (2*SL-2);
        assert( 0 <= low && high < prime );

        uint64_t mi2 = modulo_search_euclid(prime, base_r, low, high);
        mi += mi2;

        // TODO can just add mi * modulo?
        __int128 mult = (__int128) base_r * (M + mi);
        modulo = mult % prime;

        assert( (modulo < SL) || (modulo + SL) > prime );
    }
}


std::string gen_unknown_fn(const struct Config& config, std::string suffix) {
    return std::to_string(config.mstart) + "_" +
           std::to_string(config.p) + "_" +
           std::to_string(config.d) + "_" +
           std::to_string(config.minc) + "_s" +
           std::to_string(config.sieve_length) + "_l" +
           std::to_string(config.sieve_range / 1'000'000) + "M" +
           suffix;
}


void K_stats(
        const struct Config& config,
        mpz_t &K, int *K_digits, double *K_log) {
    *K_digits = mpz_sizeinbase(K, 10);

    long exp;
    double mantis = mpz_get_d_2exp(&exp, K);
    *K_log = log(mantis) + log(2) * exp;

    double m_log = log(config.mstart);
    int K_bits   = mpz_sizeinbase(K, 2);
    printf("K = %d bits, %d digits, log(K) = %.2f\n",
        K_bits, *K_digits, *K_log);
    printf("Min Gap ~= %d (for merit > %.1f)\n\n",
        (int) (config.minmerit * (*K_log + m_log)), config.minmerit);
}


// Todo float => double for K_log
double prob_prime_and_stats(
        const struct Config& config, mpz_t &K, vector<uint32_t> &primes) {

    int K_digits;
    double K_log;
    K_stats(config, K, &K_digits, &K_log);

    double prob_prime = 1 / (K_log + log(config.mstart));

    // From Mertens' 3rd theorem
    double unknowns_after_sieve = 1 / (log(config.sieve_range) * exp(GAMMA));
    double prob_prime_after_sieve = prob_prime / unknowns_after_sieve;

    size_t count_coprime = 0;
    double prob_prime_coprime = 0;
    double prob_gap_hypothetical = prop_gap_larger(
        config, prob_prime, &prob_prime_coprime, &count_coprime);

    printf("\tavg %.0f left from %.3f%% of %u interval with %ldM sieve\n",
            count_coprime * (unknowns_after_sieve / prob_prime_coprime),
            100 * unknowns_after_sieve,
            config.sieve_length, config.sieve_range/1'000'000);
    printf("\t%.3f%% of %d digit numbers are prime\n",
            100 * prob_prime, K_digits);
    printf("\t%.3f%% of tests should be prime (%.1fx speedup)\n",
            100 * prob_prime_after_sieve, 1 / unknowns_after_sieve);
    printf("\t~2x%.1f = %.1f PRP tests per m\n",
            1 / prob_prime_after_sieve, 2 / prob_prime_after_sieve);
    printf("\tsieve_length=%d is insufficient ~~%.3f%% of time\n",
            config.sieve_length, 100 * prob_gap_hypothetical);
    cout << endl;

    return prob_prime;
}


void prime_gap_search(const struct Config config) {
    const uint64_t M_start = config.mstart;
    const uint64_t M_inc = config.minc;
    const uint64_t P = config.p;
    const uint64_t D = config.d;

    const unsigned int SIEVE_LENGTH = config.sieve_length;
    const unsigned int SL = SIEVE_LENGTH;

    const uint64_t SIEVE_RANGE = config.sieve_range;

    mpz_t test;
    mpz_init(test);

    // ----- Merit Stuff
    mpz_t K;
    mpz_init(K);
    mpz_primorial_ui(K, P);
    assert( 0 == mpz_tdiv_q_ui(K, K, D) );
    assert( mpz_cmp_ui(K, 1) > 0); // K <= 1 ?!?

    // ----- Save Output file
    std::ofstream unknown_file;
    if (config.save_unknowns) {
        std::string fn = gen_unknown_fn(config, ".txt");
        printf("\tSaving unknowns to '%s'\n", fn.c_str());
        unknown_file.open(fn, std::ios::out);
        assert( unknown_file.is_open() ); // Can't open save_unknowns file
    }

    // ----- Generate primes under SIEVE_RANGE.
    vector<uint32_t> small_primes = get_sieve_primes(SIEVE_SMALL);
    const size_t SIEVE_SMALL_PRIME_PI = small_primes.size();
    {
        // SIEVE_SMALL deals with all primes that can mark off two items in SIEVE_LENGTH.
        assert( SIEVE_SMALL > 2 * SIEVE_LENGTH );
        printf("\tUsing %'ld primes for SIEVE_SMALL(%'d)\n\n",
            SIEVE_SMALL_PRIME_PI, SIEVE_SMALL);
        assert( small_primes[SIEVE_SMALL_PRIME_PI-1] < SIEVE_SMALL);
        assert( small_primes[SIEVE_SMALL_PRIME_PI-1] + 200 > SIEVE_SMALL);
    }

    // ----- Sieve stats
    prob_prime_and_stats(config, K, small_primes);

    // ----- Allocate memory for a handful of utility functions.
    auto  s_setup_t = high_resolution_clock::now();

    // Remainders of (p#/d) mod prime
    typedef pair<uint64_t,uint64_t> p_and_r;
    vector<p_and_r> prime_and_remainder;

    // Big improvement over surround_prime is avoiding checking each large prime.
    // vector<m, vector<pi>> for large primes that only rarely divide a sieve
    int s_large_primes_rem = 0;

    // To save space, only save remainder for primes that divide ANY m in range.
    // This helps with memory usage when SIEVE_RANGE >> X * MINC;

    // TODO only allocate 10'000 at a time, further out mi go to a waiting vector.
    // Downside is this takes more memory (have to store pi, mi).
    std::vector<p_and_r> *large_prime_queue = new vector<p_and_r>[M_inc];
    {
        size_t pr_pi = 0;
        printf("\tCalculating first m each prime divides\n");
        // large_prime_queue size can be approximated by
        // https://en.wikipedia.org/wiki/Meissel–Mertens_constant

        // Print "."s during, equal in length to 'Calculat...'
        size_t print_dots = 38;

        size_t expected_primes = common_primepi.count(SIEVE_RANGE) ?
            common_primepi[SIEVE_RANGE] : 1.04 * SIEVE_RANGE / log(SIEVE_RANGE);

        long first_m_sum = 0;
        double expected_large_primes = 0;
        cout << "\t";
        size_t pi = 0;
        get_sieve_primes_segmented_lambda(SIEVE_RANGE, [&](const uint64_t prime) {
            pi += 1;
            if ((pi * print_dots) % expected_primes < print_dots) {
                cout << "." << std::flush;
            }

            // Big improvement over surround_prime is reusing this for each m.
            const uint64_t base_r = mpz_fdiv_ui(K, prime);

            if (prime <= SIEVE_SMALL) {
                prime_and_remainder.emplace_back(prime, base_r);
                pr_pi += 1;
                return;
            }

            expected_large_primes += (2.0 * SL - 1) / prime;

            // solve base_r * (M + mi) + (SL - 1)) % prime < 2 * SL
            //   0 <= (base_r * M + SL - 1) + base_r * mi < 2 * SL mod prime
            //
            // let shift = (base_r * M + SL - 1) % prime
            //   0 <= shift + base_r * mi < 2 * SL mod prime
            // add (prime - shift) to all three
            //
            //  (prime - shift) <= base_r * mi < (prime - shift) + 2 * SL mod prime
            uint64_t mi = modulo_search_euclid_gcd(
                    M_start, D, M_inc, SL, prime, base_r);

            // signals mi > M_inc
            if (mi == M_inc) return;

            assert (mi < M_inc);

            __int128 mult = (__int128) base_r * (M_start + mi) + (SL - 1);
            assert( mult % prime < (2*SL-1) );

            //assert ( gcd(M + mi, D) == 1 );

            large_prime_queue[mi].emplace_back(prime, base_r);
            pr_pi += 1;

            s_large_primes_rem += 1;
            first_m_sum += mi;
        });
        cout << endl;

        printf("\tSum of m1: %ld\n", first_m_sum);
        setlocale(LC_NUMERIC, "");
        if (expected_primes == pi) {
            printf("\tPrimePi(%ld) = %ld\n", SIEVE_RANGE, pi);
        } else {
            printf("\tPrimePi(%ld) = %ld guessed %ld\n", SIEVE_RANGE, pi, expected_primes);
            assert(common_primepi.count(SIEVE_RANGE) == 0);
        }

        printf("\t%ld primes not needed (%.1f%%)\n",
            (pi - SIEVE_SMALL_PRIME_PI) - pr_pi,
            100 - (100.0 * pr_pi / (pi - SIEVE_SMALL_PRIME_PI)));
        printf("\texpected large primes/m: %.1f (theoretical: %.1f)\n",
            expected_large_primes,
            (2 * SL - 1) * (log(log(SIEVE_RANGE)) - log(log(SIEVE_SMALL))));
        setlocale(LC_NUMERIC, "C");

        assert(prime_and_remainder.size() == small_primes.size());
    }
    {
        auto  s_stop_t = high_resolution_clock::now();
        double   secs = duration<double>(s_stop_t - s_setup_t).count();
        printf("\n\tSetup took %.1f seconds\n", secs);
    }

    // ----- Main sieve loop.
    cout << "\nStarting m=" << M_start << "\n" << endl;

    // vector<bool> uses bit indexing which is ~5% slower.
    vector<char> composite[2] = {
        vector<char>(SIEVE_LENGTH, 0),
        vector<char>(SIEVE_LENGTH, 0)
    };
    assert( composite[0].size() == SIEVE_LENGTH );
    assert( composite[1].size() == SIEVE_LENGTH );

    // Used for various stats
    long  s_tests = 0;
    auto  s_start_t = high_resolution_clock::now();
    long  s_total_unknown = 0;
    long  s_t_unk_low = 0;
    long  s_t_unk_hgh = 0;
    long  s_large_primes_tested = 0;

    uint64_t mi_last = M_inc - 1;
    for (; mi_last > 0 && gcd(mi_last, D) > 1; mi_last -= 1);
    assert( mi_last > 0 );

    for (uint64_t mi = 0; mi < M_inc; mi++) {
        uint64_t m = M_start + mi;
        bool good_m = gcd(m, D) == 1;

        if (!good_m) {
            assert( large_prime_queue[mi].empty() );
            continue;
        }

        // Reset sieve array to unknown.
        std::fill_n(composite[0].begin(), SIEVE_LENGTH, 0);
        std::fill_n(composite[1].begin(), SIEVE_LENGTH, 0);
        // center is always composite.
        composite[0][0] = composite[1][0] = 1;

        // For small primes that we don't do trick things with.
        for (const auto& pr : prime_and_remainder) {
            const uint64_t modulo = (pr.second * m) % pr.first;
            //            const auto& [prime, remainder] = prime_and_remainder[pi];
            //            const uint64_t modulo = (remainder * m) % prime;

            for (size_t d = modulo; d < SIEVE_LENGTH; d += pr.first) {
                composite[0][d] = true;
            }

            // Not technically correct but fine to skip modulo == 0
            int first_negative = pr.first - modulo;
            for (size_t d = first_negative; d < SIEVE_LENGTH; d += pr.first) {
                composite[1][d] = true;
            }
        }

        // Maybe useful for some stats later.
        // int unknown_small_l = std::count(composite[0].begin(), composite[0].end(), false);
        // int unknown_small_u = std::count(composite[1].begin(), composite[1].end(), false);

        for (const auto& pr : large_prime_queue[mi]) {
            s_large_primes_tested += 1;
            s_large_primes_rem -= 1;

            const auto& prime = pr.first;
            const auto& remainder = pr.second;

            // Large prime should divide some number in SIEVE for this m
            // When done find next mi where prime divides a number in SIEVE.
            const uint64_t modulo = (remainder * m) % prime;

            if (0) {
                mpz_mul_ui(test, K, m);
                uint64_t mod = mpz_fdiv_ui(test, prime);
                assert( mod == modulo );
            }

            if (good_m) {
                if (modulo < SIEVE_LENGTH) {
                    // Just past a multiple
                    composite[0][modulo] = true;
                } else {
                    // Don't have to deal with 0 case anymore.
                    int64_t first_positive = prime - modulo;
                    assert( first_positive < SIEVE_LENGTH); // Bad next m!
                    // Just before a multiple
                    composite[1][first_positive] = true;
                }
            }

            // Find next mi where primes divides part of SIEVE
            {
                uint64_t start = mi + 1;
                uint64_t next_mi = start + modulo_search_euclid_gcd(
                        M_start + start, D, M_inc - start, SL, prime, remainder);
                if (next_mi == M_inc) continue;

                //assert( (remainder * (M + next_mi) + (SL - 1)) % prime < (2*SL-1) );
                __int128 mult = (__int128) remainder * (M_start + next_mi) + (SL - 1);
                assert ( mult % prime <= (2 * SL - 1) );

                //assert ( gcd(M_start + next_mi, D) == 1 );

                large_prime_queue[next_mi].push_back(pr);
                s_large_primes_rem += 1;
            }
        }
        large_prime_queue[mi].clear();
        large_prime_queue[mi].shrink_to_fit();

        s_tests += 1;

        // 2-3% of runtime, could be optimized into save_unknowns loop..
        int unknown_l = std::count(composite[0].begin(), composite[0].end(), false);
        int unknown_u = std::count(composite[1].begin(), composite[1].end(), false);
        s_total_unknown += unknown_l + unknown_u;
        s_t_unk_low += unknown_l;
        s_t_unk_hgh += unknown_u;

        // Save unknowns
        if (config.save_unknowns) {

            unknown_file << mi;
            unknown_file << " : -" << unknown_l << " +" << unknown_u << " |";

            for (int d = 0; d <= 1; d++) {
                char prefix = "-+"[d];

                for (size_t i = 1; i < SL; i++) {
                    if (!composite[d][i]) {
                        unknown_file << " " << prefix << i;
                    }
                }
                if (d == 0) {
                    unknown_file << " |";
                }
            }
            unknown_file << "\n";
        }

        if ( (s_tests == 1 || s_tests == 10 || s_tests == 100 || s_tests == 500 || s_tests == 1000) ||
                (s_tests % 5000 == 0) || (mi == mi_last) ) {
            auto s_stop_t = high_resolution_clock::now();
            double   secs = duration<double>(s_stop_t - s_start_t).count();

            printf("\t%ld %4d <- unknowns -> %-4d\n",
                    m, unknown_l, unknown_u);

            if (mi <= 10) continue;

            // Stats!
            printf("\t    tests     %-10ld (%.2f/sec)  %.0f seconds elapsed\n",
                    s_tests, s_tests / secs, secs);
            printf("\t    unknowns  %-10ld (avg: %.2f), %.2f%% composite  %.2f <- %% -> %.2f%%\n",
                    s_total_unknown, s_total_unknown / ((double) s_tests),
                    100.0 * (1 - s_total_unknown / (2.0 * (SIEVE_LENGTH - 1) * s_tests)),
                    100.0 * s_t_unk_low / s_total_unknown,
                    100.0 * s_t_unk_hgh / s_total_unknown);
            printf("\t    large prime remaining: %d (avg/test: %ld)\n",
                    s_large_primes_rem, s_large_primes_tested / s_tests);
        }
    }

    // Should be cleaning up after self.
    for(uint32_t mi = 0; mi < M_inc; mi++)  {
        assert( large_prime_queue[mi].size() == 0 );
    }

    // ----- cleanup

    delete[] large_prime_queue;
    mpz_clear(K);
    mpz_clear(test);
}

void prime_gap_parallel(const struct Config config) {
    // Method2
    const uint64_t M_start = config.mstart;
    const uint64_t M_inc = config.minc;
    const uint64_t P = config.p;
    const uint64_t D = config.d;

    const uint32_t SIEVE_LENGTH = config.sieve_length;
    const uint32_t SL = SIEVE_LENGTH;

    const uint64_t SIEVE_RANGE = config.sieve_range;

    const uint64_t SMALL_THRESHOLD = 50 * SIEVE_LENGTH;

    // SMALL_THRESHOLD deals with all primes that can mark off two items in SIEVE_LENGTH.
    assert( SMALL_THRESHOLD > 2 * SIEVE_LENGTH );

    mpz_t test;
    mpz_init(test);

    // ----- Merit Stuff
    mpz_t K;
    mpz_init(K);
    mpz_primorial_ui(K, P);
    assert( 0 == mpz_tdiv_q_ui(K, K, D) );
    assert( mpz_cmp_ui(K, 1) > 0); // K <= 1 ?!?

    // ----- Save Output file
    std::ofstream unknown_file;
    if (config.save_unknowns) {
        std::string fn = gen_unknown_fn(config, ".m2.txt");
        printf("\tSaving unknowns to '%s'\n", fn.c_str());
        unknown_file.open(fn, std::ios::out);
        assert( unknown_file.is_open() ); // Can't open save_unknowns file
    }

    // ----- Generate primes for P
    vector<uint32_t> P_primes = get_sieve_primes(P);
    assert( P_primes.back() == P);

    // ----- Sieve stats
    double prob_prime = prob_prime_and_stats(config, K, P_primes);

    // ----- Allocate memory
    uint32_t SIEVE_INTERVAL = 2 * SIEVE_LENGTH - 1;

    vector<char> coprime_composite(SIEVE_INTERVAL, false);
    for (uint32_t prime : P_primes) {
        if (D % prime != 0) {
            int first = (SIEVE_LENGTH-1) % prime;
            assert ( ((SIEVE_LENGTH-1) - first) % prime == 0 );
            for (size_t d = first; d < SIEVE_INTERVAL; d += prime) {
                coprime_composite[d] = true;
            }
        }
    }

    vector<int32_t> m_reindex(M_inc, -1);
    size_t valid_ms = 0;
    {
        for (uint64_t mi = 0; mi < M_inc; mi++) {
            if (gcd(M_start + mi, D) == 1) {
                m_reindex[mi] = valid_ms;
                valid_ms++;
            }
        }
    }

    // <bool> is slower than <char>, but uses 1/8th the memory.
    vector<bool> composite[valid_ms];
    {
        auto  s_interval_t = high_resolution_clock::now();
        for (size_t i = 0; i < valid_ms; i++) {
            // Improve this setup.
            composite[i].resize(SIEVE_INTERVAL, true);
            for (size_t d = 0; d < SIEVE_INTERVAL; d++) {
                if (coprime_composite[d] == false) {
                    composite[i][d] = false;
                }
            }
        };
        auto   s_stop_t = high_resolution_clock::now();
        double int_secs = duration<double>(s_stop_t - s_interval_t).count();
        printf("coprime setup took %.1f (valid coprime m %ld of %ld)\n",
            int_secs, valid_ms, M_inc);
    }

    // Used for various stats
    auto  s_start_t = high_resolution_clock::now();
    auto  s_interval_t = high_resolution_clock::now();
    long  s_total_unknowns = SIEVE_INTERVAL * M_inc;
    long  s_small_primes_tested = 0;
    long  s_large_primes_tested = 0;
    uint64_t  s_next_print = 0;
    uint64_t  next_mult = SMALL_THRESHOLD <= 100000 ? 100000 : 1000000;
    double s_prp_needed = 1 / prob_prime;

    // Note: Handling small primes here had better localized memory access
    // But wasn't worth the extra code IMHO.

    size_t pi = 0;
    get_sieve_primes_segmented_lambda(SIEVE_RANGE, [&](const uint64_t prime) {
        pi += 1;
        // Big improvement over surround_prime is reusing this for each m.
        const uint64_t base_r = mpz_fdiv_ui(K, prime);

        if (prime <= SMALL_THRESHOLD) {
            // Handled by coprime_composite above
            if (D % prime != 0 && prime <= P) {
                return;
            }

            for (uint64_t mi = 0; mi < M_inc; mi++) {
                int32_t mii = m_reindex[mi];
                if (mii < 0) continue;

                uint64_t m = M_start + mi;

                uint64_t modulo = (base_r * m) % prime;

                uint32_t flip = modulo + prime - (SIEVE_LENGTH % prime);
                if (flip >= prime) flip -= prime;
                uint32_t first = prime - flip - 1;

                uint32_t shift = prime;
                if (prime > 2) {
                    //bool lowIsEven = ((SIEVE_LENGTH % 2) == 0) ^ ((D % 2) == 1) || ((m % 2) == 0);
                    //bool evenFrom = (first % 2) == 0;

                    bool lowIsEven = (D & 1) || ((m & 1) == 0);
                    bool evenFromCenter = (first & 1) == (SIEVE_LENGTH & 1);

                    if (lowIsEven ^ evenFromCenter) {
                        // Same parity (divisible by 2) move to next prime
                        assert( (first >= SIEVE_INTERVAL) || (composite[mii][first] > 0) );
                        first += prime;
                    }

                    // Don't need to count cross off even multiples.
                    shift = 2*prime;
                }

                for (size_t d = first; d < SIEVE_INTERVAL; d += shift) {
                    composite[mii][d] = true;
                    s_small_primes_tested += 1;
                }
            }
        } else {
            modulo_search_euclid_all(M_start, M_inc, SL, prime, base_r, [&](const uint64_t mi) {
                // TODO validate no overflow
                // TODO return first from lambda?
                uint64_t first = (base_r * (M_start + mi) + (SL-1)) % prime;
                assert( first < SIEVE_INTERVAL );
                first = SIEVE_INTERVAL - first - 1;
                assert( 0 <= first && first < SIEVE_INTERVAL );

                if (coprime_composite[first] == true) {
                    return;
                }
                int32_t mii = m_reindex[mi];
                if (mii < 0) {
                    return;
                }

                if (0) {
                    mpz_mul_ui(test, K, M_start + mi);
                    mpz_sub_ui(test, test, SIEVE_LENGTH-1);
                    mpz_add_ui(test, test, first);
                    uint64_t mod = mpz_fdiv_ui(test, prime);
                    assert( mod == 0 );
                }

                composite[mii][first] = true;
                s_large_primes_tested += 1;
            });
        }

        if (prime >= s_next_print) {
            if (prime >= s_next_print) {
                if (s_next_print == 5 * next_mult) {
                    next_mult = 2 * s_next_print;
                    s_next_print = 0;
                }
                s_next_print += next_mult;
                if (s_next_print >= SIEVE_RANGE) {
                    // Would be nice to have a prev prime :)
                    s_next_print = SIEVE_RANGE - 50;
                }
            }

            auto   s_stop_t = high_resolution_clock::now();
            double     secs = duration<double>(s_stop_t - s_start_t).count();
            double int_secs = duration<double>(s_stop_t - s_interval_t).count();

            double unknowns_after_sieve = 1 / (log(prime) * exp(GAMMA));
            double prob_prime_after_sieve = prob_prime / unknowns_after_sieve;

            uint64_t t_total_unknowns = 0;
            for (size_t i = 0; i < valid_ms; i++) {
                t_total_unknowns += std::count(composite[i].begin(), composite[i].end(), false);
            }
            uint64_t saved_prp = s_total_unknowns - t_total_unknowns;
            double skipped_prp = M_inc * (s_prp_needed - 1/prob_prime_after_sieve);

            printf("%10ld %5.2f/%-6.1f seconds |", prime, int_secs, secs);
            printf(" %ld primes used (avg/m: %.1f) \n",
                s_large_primes_tested,
                1.0 * s_large_primes_tested / M_inc);
            printf("\tunknowns %ld (avg: %.2f) (%.3f%%)\n",
                t_total_unknowns,
                1.0 * t_total_unknowns / valid_ms,
                100.0 * t_total_unknowns / (SIEVE_INTERVAL * M_inc));
            printf("\t~ 2x %.2f PRP/m (%ld new composites ~= %4.1f skipped PRP => %.1f PRP/seconds)\n",
                1 / prob_prime_after_sieve,
                saved_prp,
                skipped_prp,
                1.0 * skipped_prp / int_secs);

            s_total_unknowns = t_total_unknowns;
            s_interval_t = s_stop_t;
            s_prp_needed = 1 / prob_prime_after_sieve;
        }
    });

    if (config.save_unknowns) {
        for (uint64_t mi = 0; mi < M_inc; mi++) {
            if (gcd(M_start + mi, D) > 1) {
                continue;
            }
            int32_t mii = m_reindex[mi];
            assert( mii >= 0 );

            const auto& comp = composite[mii];

            size_t unknown_l = std::count(
                comp.begin(),       comp.begin() + SL, false);
            size_t unknown_u = std::count(
                comp.begin() + SL,  comp.begin() + 2*SL-1, false);

            unknown_file << mi << " : -" << unknown_l << " +" << unknown_u << " |";
            for (int d = 0; d <= 1; d++) {
                char prefix = "-+"[d];

                for (size_t i = 1; i < SL; i++) {
                    if (!comp[(SL-1) + (d == 0 ? -1 : 1) * i]) {
                        unknown_file << " " << prefix << i;
                    }
                }
                if (d == 0) {
                    unknown_file << " |";
                }
            }
            unknown_file << "\n";
        }
    }

    mpz_clear(K);
    mpz_clear(test);
}

