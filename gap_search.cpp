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
#include <fstream>
#include <iostream>
#include <vector>
#include <locale.h>

#include <gmp.h>

#include "gap_common.h"

using std::cout;
using std::endl;
using std::vector;
using namespace std::chrono;

// Tweaking this doesn't seem to change speed much.
#define SIEVE_SMALL       400'000

void set_defaults(struct Config& config);
void prime_gap_search(const struct Config config);


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
    printf("Testing m * %d#/%d, m = %d + [0, %d)\n",
        config.p, config.d, config.mstart, config.minc);

    printf("\n");
    printf("sieve_length: 2x%d\n", config.sieve_length);
    printf("sieve_range:  %ld\n", config.sieve_range);
    printf("\n");


    prime_gap_search(config);
}


void set_defaults(struct Config& config) {
    if (config.valid == 0) {
        // Don't do anything if argparse didn't work.
        return;
    }

    float logK;
    {
        mpz_t K;
        mpz_init(K);

        mpz_primorial_ui(K, config.p);
        assert( 0 == mpz_tdiv_q_ui(K, K, config.d) );
        long exp;
        double mantis = mpz_get_d_2exp(&exp, K);
        logK = log(config.mstart) + log(mantis) + log(2) * exp;
        mpz_clear(K);
    }

    if (config.sieve_length == 0) {
        // factors of K = p#/d
        vector<uint32_t> K_primes = get_sieve_primes(config.p);
        {
            int td = config.d;
            while (td > 1) {
                bool change = false;
                for (size_t pi = 0; pi <= K_primes.size(); pi++) {
                    if (td % K_primes[pi] == 0) {
                        td /= K_primes[pi];
                        K_primes.erase(K_primes.begin() + pi);
                        change = true;
                        // Indexes are messed changes so do easy thing.
                        break;
                    }
                }
                assert( change ); // d is not made up of primes <= p.
            }
        }

        // Prob prime is slightly higher in practice
        double prob_prime_coprime = 1 / logK;
        for (int prime : K_primes) {
            prob_prime_coprime /= (1 - 1.0/prime);
        }
        printf("\tprob_prime_coprime: %.5f\n", prob_prime_coprime);

        // K = #p/d
        // only numbers K+i has no factor <= p
        //      => (K+i, i) == 1
        //      => only relatively prime i
        // p is generally large enough that SL <= p*p
        //      => i is prime (not two product of two primes > p)
        //          or
        //         i is composite of factors in d
        assert( config.p > 51 );

        // Search till something counting chance of shorter gap.
        {
            size_t P = config.p;
            int count_coprime = 0;
            for (size_t tSL = 1; tSL <= P*P; tSL += 1) {
                count_coprime += 1;
                for (int prime : K_primes) {
                    if ((tSL % prime) == 0) {
                        count_coprime -= 1;
                        break;
                    }
                }

                // Assume each coprime is independent (not quite true)
                double prob_gap_shorter = pow(1 - prob_prime_coprime, count_coprime);

                // This seems to balance PRP fallback and sieve_size
                if (prob_gap_shorter <= 0.008) {
                    config.sieve_length = tSL;
                    printf("AUTO SET: sieve length (coprime: %d, prob_gap longer %.2f%%): %ld\n",
                        count_coprime, 100 * prob_gap_shorter, tSL);
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
            config.sieve_range =   100'000'000;
        }

        printf("AUTO SET: sieve range (log(t) = ~%.0f): %ld\n",
            logK, config.sieve_range);
    }

}


uint32_t gcd(uint32_t a, uint32_t b) {
    if (b == 0) return a;
    return gcd(b, a % b);
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

uint64_t modulo_search_euclid(uint64_t p, uint64_t a, uint64_t l, uint64_t r) {
    // min i : l <= (a * i) % p <= l
/*
    assert( 0 <= a && a < p );
    assert( 0 <= l && l < p );
    assert( 0 <= r && r < p );
    assert(      l <= r     );
// */

    if (l == 0) return 0;

    // 2 *a > p, but avoids int max
    if (a > (p >> 1)) {
        std::swap(l, r);
        a = p - a;
        l = p - l;
        r = p - r;
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
    // k < a

    // TODO make sure this doesn't overflow.
    //uint64_t tl = l + p * k;
    //uint64_t mult = (tl-1) / a + 1;
    __int128 tl = l + p * k;
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


void prime_gap_search(const struct Config config) {
    const uint64_t M = config.mstart;
    const uint64_t M_inc = config.minc;
    const uint64_t P = config.p;
    const uint64_t D = config.d;
    const float min_merit = config.minmerit;

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

    // ----- Save Output file
    std::ofstream unknown_file;
    if (config.save_unknowns) {
        std::string fn =
            std::to_string(M) + "_" +
            std::to_string(P) + "_" +
            std::to_string(D) + "_" +
            std::to_string(M_inc) + "_s" +
            std::to_string(config.sieve_length) + "_l" +
            std::to_string(config.sieve_range / 1'000'000) + "M" +
            ".txt";
        printf("\tSaving unknowns to '%s'\n", fn.c_str());
        unknown_file.open(fn, std::ios::out);
        assert( unknown_file.is_open() ); // Can't open save_unknowns file
    }

    // ----- Generate primes under SIEVE_RANGE.
    auto  s_primes_start_t = high_resolution_clock::now();
    vector<uint64_t> primes = get_sieve_primes_segmented(SIEVE_RANGE);
    auto  s_primes_stop_t = high_resolution_clock::now();
    const size_t SIEVE_SMALL_PRIME_PI = std::distance(primes.begin(),
        std::lower_bound(primes.begin(), primes.end(), SIEVE_SMALL));
    {
        double   secs = duration<double>(s_primes_stop_t - s_primes_start_t).count();

        setlocale(LC_NUMERIC, "");
        printf("\tPrimePi(%'ld) = %'ld (2 ... %'lu)\n",
            SIEVE_RANGE, primes.size(), primes.back());
        printf("\tSegmented prime sieve took %.1f seconds\n", secs);

        // SIEVE_SMALL deals with all primes can mark off two items in SIEVE_LENGTH.
        assert( SIEVE_SMALL > 2 * SIEVE_LENGTH );
        printf("\tUsing %'ld primes for SIEVE_SMALL(%'d)\n\n",
            SIEVE_SMALL_PRIME_PI, SIEVE_SMALL);
        assert( primes[SIEVE_SMALL_PRIME_PI] > SIEVE_SMALL );
        setlocale(LC_NUMERIC, "C");
    }

    // ----- Sieve stats
    {
        // Can also use Mertens' 3rd theorem
        double unknowns_after_sieve = 1;
        for (int64_t prime : primes) unknowns_after_sieve *= (prime - 1.0) / prime;

        double prob_prime = 1 / (K_log + log(M));
        double prob_prime_coprime = 1;
        double prob_prime_after_sieve = prob_prime / unknowns_after_sieve;

        for (size_t pi = 0; primes[pi] <= P; pi++) {
            if (D % primes[pi] != 0) {
                prob_prime_coprime *= (1 - 1.0/primes[pi]);
            }
        }

        int count_coprime = SL-1;
        for (size_t i = 1; i < SL; i++) {
            for (uint64_t prime : primes) {
                if (prime > P) break;
                if ((i % prime) == 0 && (D % prime) != 0) {
                    count_coprime -= 1;
                    break;
                }
            }
        }
        double chance_coprime_composite = 1 - prob_prime / prob_prime_coprime;
        double prob_gap_shorter_hypothetical = pow(chance_coprime_composite, count_coprime);

        // count_coprime already includes some parts of unknown_after_sieve
        printf("\t%.3f%% of sieve should be unknown (%ldM) ~= %.0f\n",
            100 * unknowns_after_sieve,
            SIEVE_RANGE/1'000'000,
            count_coprime * (unknowns_after_sieve / prob_prime_coprime));
        printf("\t%.3f%% of %d digit numbers are prime\n",
            100 * prob_prime, K_digits);
        printf("\t%.3f%% of tests should be prime (%.1fx speedup)\n",
            100 * prob_prime_after_sieve, 1 / unknowns_after_sieve);
        printf("\t~2x%.1f = %.1f PRP tests per m\n",
            1 / prob_prime_after_sieve, 2 / prob_prime_after_sieve);
        printf("\tsieve_length=%d is insufficient ~%.2f%% of time\n",
            SIEVE_LENGTH, 100 * prob_gap_shorter_hypothetical);
        cout << endl;
    }


    // ----- Allocate memory for a handful of utility functions.
    auto  s_setup_t = high_resolution_clock::now();

    // Remainders of (p#/d) mod prime
    vector<uint64_t> remainder;
    {
        cout << "\t";
        for (size_t pi = 0; pi < SIEVE_SMALL_PRIME_PI; pi++) {
            const uint64_t prime = primes[pi];

            // Big improvement over surround_prime is reusing this for each m.
            uint64_t mod = mpz_fdiv_ui(K, prime);
            //assert( 0 <= mod && mod < prime );
            remainder.push_back(mod);
        }
    }
    cout << endl;

    // Big improvement over surround_prime is avoiding checking each large prime.
    // vector<m, vector<pi>> for large primes that only rarely divide a sieve
    int s_large_primes_rem = 0;

    // To save space, only save remainder for primes that divide ANY m in range.
    // This helps with memory usage when SIEVE_RANGE > X * MINC;

    // TODO only allocate 10'000 at a time, further out mi go to a waiting vector.
    // Downside is this takes more memory (have to store pi, mi).
    std::vector<uint32_t> *large_prime_queue = new vector<uint32_t>[M_inc];
    {
        size_t new_pi = SIEVE_SMALL_PRIME_PI;
        printf("\tCalculating first m each prime divides\n");
        // large_prime_queue size can be approximated by
        // https://en.wikipedia.org/wiki/Meissel–Mertens_constant

        // Print "."s during, equal in length to 'Calculat...'
        size_t print_dots = 38;

        long first_m_sum = 0;
        double expected_large_primes = 0;
        cout << "\t."; // 0 * 38 % primes.size() < 38
        for (size_t pi = SIEVE_SMALL_PRIME_PI; pi < primes.size(); pi++) {
            if ((pi * print_dots) % primes.size() < print_dots) {
                cout << "." << std::flush;
            }

            const uint64_t prime = primes[pi];
            expected_large_primes += (2.0 * SL - 1) / prime;

            // Big improvement over surround_prime is reusing this for each m.
            uint64_t base_r = mpz_fdiv_ui(K, prime);
            const uint64_t modulo = (base_r * M) % prime;
            if ( (modulo < SL) || (modulo + SL) > prime) {
                remainder.push_back(base_r);
                large_prime_queue[0].push_back(new_pi);
                new_pi += 1;

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

            uint64_t shift = modulo + (SL - 1);
            assert( 0 <= shift && shift < prime );
            uint64_t low  = (prime - shift);
            uint64_t high = low + (2*SL-2);
            assert( 0 <= low && high < prime );

            uint64_t mi = modulo_search_euclid(prime, base_r, low, high);
            assert( low <= (mi * base_r) % prime );
            assert(        (mi * base_r) % prime <= high );

            assert( (base_r * (M + mi) + (SL - 1)) % prime < (2*SL-1) );
            if (mi < M_inc) {
                remainder.push_back(base_r);
                large_prime_queue[mi].push_back(new_pi);
                new_pi += 1;

                s_large_primes_rem += 1;
                first_m_sum += mi;
            } else {
                // Delete this prime
                primes[pi] = 0;
            }

            if (0) {
                // Brute force doublecheck
                assert( mi == modulo_search_brute(prime, base_r, low, high) );
            }
        }

        size_t pre_num_primes = primes.size();
        primes.erase(std::remove(primes.begin(), primes.end(), 0), primes.end());
        size_t post_num_primes = primes.size();
        assert( new_pi == post_num_primes );

        cout << endl;
        printf("\tSum of m1: %ld\n", first_m_sum);
        printf("\tDeleted %ld primes (%.1f%%)\n",
            pre_num_primes - post_num_primes,
            100 - (100.0 * new_pi / pre_num_primes));
        printf("\texpected large primes/m: %.1f (theoretical: %.1f)\n",
            expected_large_primes,
            (2 * SL - 1) * (log(log(SIEVE_RANGE)) - log(log(SIEVE_SMALL))));
    }
    {
        auto  s_stop_t = high_resolution_clock::now();
        double   secs = duration<double>(s_stop_t - s_setup_t).count();
        printf("\n\tSetup took %.1f seconds\n", secs);
    }


    // ----- Main sieve loop.
    cout << "\nStarting m=" << M << "\n" << endl;

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

    uint64_t mi_last = M_inc;
    for (; mi_last > 0 && gcd(mi_last, D) > 1; mi_last -= 1);
    assert( mi_last > 0 );

    for (uint64_t mi = 0; mi < M_inc; mi++) {
        uint64_t m = M + mi;

        // Reset sieve array to unknown.
        std::fill_n(composite[0].begin(), SIEVE_LENGTH, 0);
        std::fill_n(composite[1].begin(), SIEVE_LENGTH, 0);
        // center is always composite.
        composite[0][0] = composite[1][0] = 1;

        // For small primes that we don't do trick things with.
        for (size_t pi = 0; pi < SIEVE_SMALL_PRIME_PI; pi++) {
            const uint32_t prime = primes[pi];
            long modulo = (remainder[pi] * m) % prime;

            for (size_t d = modulo; d < SIEVE_LENGTH; d += prime) {
                composite[0][d] = true;
            }

            // Not technically correct but fine to skip modulo == 0
            int first_negative = prime - modulo;
            for (size_t d = first_negative; d < SIEVE_LENGTH; d += prime) {
                composite[1][d] = true;
            }
        }

        // Maybe useful for some stats later.
        // int unknown_small_l = std::count(composite[0].begin(), composite[0].end(), false);
        // int unknown_small_u = std::count(composite[1].begin(), composite[1].end(), false);

        for (uint32_t pi : large_prime_queue[mi]) {
            s_large_primes_tested += 1;
            s_large_primes_rem -= 1;

            // Large prime should divide some number in SIEVE for this m
            // When done find next mi where prime divides a number in SIEVE.
            const uint64_t prime   = primes[pi];
            const uint64_t base_r  = remainder[pi];
            const uint64_t modulo = (base_r * m) % prime;

            if (0) {
                mpz_mul_ui(test, K, m);
                uint64_t mod = mpz_fdiv_ui(test, prime);
                assert( mod == modulo );
            }

            if (modulo < SIEVE_LENGTH) {
                // Just past a multiple
                composite[0][modulo] = true;
            } else {
                // Don't have to deal with 0 case anymore.
                int64_t first_negative = prime - modulo;
                assert( first_negative < SIEVE_LENGTH); // Bad next m!
                // Just before a multiple
                composite[1][first_negative] = true;
            }

            // Find next mi that primes divides part of SIEVE
            {
                // next modulo otherwise modulo_search returns 0;
                uint64_t shift = (modulo + base_r) + (SL - 1);
                if (shift >= prime) shift -= prime;
                if (shift >= prime) shift -= prime;

                assert( 0 <= shift && shift < prime );

                uint64_t low  = (prime - shift);
                uint64_t high = low + (2*SL-2);

                uint64_t next_mi;
                if (high >= prime) {
                    next_mi = mi + 1;
                } else {
                    assert( 0 <= low && high < prime );
                    int64_t m2 = modulo_search_euclid(prime, base_r, low, high);
                    assert( low <= (m2 * base_r) % prime );
                    assert(        (m2 * base_r) % prime <= high );
                    next_mi = mi + 1 + m2;
                }

                assert( (base_r * (M + next_mi) + (SL - 1)) % prime < (2*SL-1) );
                if (next_mi < M_inc) {
                    large_prime_queue[next_mi].push_back(pi);
                    s_large_primes_rem += 1;
                }
            }
        }
        large_prime_queue[mi].clear();
        large_prime_queue[mi].shrink_to_fit();

        if (gcd(m, D) > 1) continue;

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
            printf("\t    large prime remaining: %d (avg/m: %ld)\n",
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

