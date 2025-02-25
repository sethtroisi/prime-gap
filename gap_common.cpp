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

#include "gap_common.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <limits>

/* for dirname(3) */
#include <libgen.h>

/* for primesieve::iterator */
#include <primesieve.hpp>


using std::cout;
using std::endl;
using std::pair;
using std::string;
using std::vector;
using namespace std::chrono;



string UNKNOWNS_DIR = "unknowns";

static const std::map<uint64_t,uint64_t> common_primepi = {
    {       10'000'000,        664'579},
    {      100'000'000,      5'761'455},
    {      200'000'000,     11'078'937},
    {      400'000'000,     21'336'326},
    {      800'000'000,     41'146'179},
    {    1'000'000'000,     50'847'534},
    {    2'000'000'000,     98'222'287},
    {    3'000'000'000,    144'449'537},
    {    4'000'000'000,    189'961'812},
    {    5'000'000'000,    234'954'223},
    {    6'000'000'000,    279'545'368},
    {   10'000'000'000,    455'052'511},
    {   15'000'000'000,    670'180'516},
    {   20'000'000'000,    882'206'716},
    {   25'000'000'000,  1'091'987'405},
    {   30'000'000'000,  1'300'005'926},
    {   40'000'000'000,  1'711'955'433},
    {   50'000'000'000,  2'119'654'578},
    {   60'000'000'000,  2'524'038'155},
    {  100'000'000'000,  4'118'054'813},
    {  200'000'000'000,  8'007'105'059},
    {  300'000'000'000,  11'818'439'135},
    {  400'000'000'000,  15'581'005'657},
    {  500'000'000'000,  19'308'136'142},
    {1'000'000'000'000,  37'607'912'018},
    {2'000'000'000'000,  73'301'896'139},
    {3'000'000'000'000, 108'340'298'703},
    {4'000'000'000'000, 142'966'208'126},
    {5'000'000'000'000, 177'291'661'649}
};


static void assert_file_exists(string path) {
    std::ifstream f(path);
    if (!f.good()) {
        printf("'%s' doesn't exist\n", path.c_str());
        exit(1);
    }
}

bool has_prev_prime_gmp() {
    return (
        (__GNU_MP_VERSION > 6) ||
        (__GNU_MP_VERSION == 6 && __GNU_MP_VERSION_MINOR >= 3) ||
        (__GNU_MP_VERSION == 6 && __GNU_MP_VERSION_MINOR == 2 && __GNU_MP_VERSION_PATCHLEVEL == 99)
    );
}


uint64_t gcd(uint64_t a, uint64_t b) {
    if (b == 0) return a;
    return gcd(b, a % b);
}


double _log(const mpz_t &K) {
    long exp;
    double mantis = mpz_get_d_2exp(&exp, K);
    return log(mantis) + log(2) * exp;
}


double calc_log_K(const struct Config& config) {
    mpz_t K;
    init_K(config, K);
    double log = _log(K);
    mpz_clear(K);
    return log;
}


void init_K(const struct Config& config, mpz_t &K) {
    mpz_init(K);
    mpz_primorial_ui(K, config.p);
    assert(0 == mpz_tdiv_q_ui(K, K, config.d));
    assert(mpz_cmp_ui(K, 1) > 0);  // K <= 1 ?!?
}


void K_stats(
        const struct Config& config,
        mpz_t &K, int *K_digits, double *K_log) {
    init_K(config, K);
    *K_log = _log(K);

    if (K_digits != nullptr) {
        int base10 = mpz_sizeinbase(K, 10);
        *K_digits = base10;

        if (config.verbose >= 2) {
        int K_bits   = mpz_sizeinbase(K, 2);
        printf("K = %d bits, %d digits, log(K) = %.2f\n",
                K_bits, base10, *K_log);
        }
    }
}


/**
 * Return estimated time (in seconds) to PRP test a composite with no small factor
 */
double prp_time_estimate_composite(double K_log, int verbose) {
    // TODO: For large K_log, time smaller PRP then upscale with polynomial

    // Some rough estimates at
    // https://github.com/sethtroisi/misc-scripts/tree/master/prime-time

    float K_log_2 = K_log * K_log;
    float t_estimate_poly = -1.1971e-03
        +  5.1072e-07 * K_log
        +  9.4362e-10 * K_log_2
        +  1.8757e-13 * K_log_2 * K_log
        + -1.9582e-18 * K_log_2 * K_log_2;
    float t_estimate = std::max(1e-3f, t_estimate_poly);

    if (verbose >= 2) {
        if (t_estimate > 0.3) {
            printf("Estimated secs/PRP: %.1f\n", t_estimate);
        } else {
            // Benchmark in thread

            // Create some non-trivial semi-primes.
            mpz_t n, p, q;
            mpz_inits(n, p, q, nullptr);

            size_t bits = K_log * 1.442;
            assert( bits > 50 );

            mpz_set_ui(n, 1);

            // Multiply "large" static primes (25 bits+) to get number of size N
            size_t bit_goal = bits - 24;
            while (bit_goal > 0) {
                // Important to not ever choose small p
                size_t p_size = bit_goal < 50 ? bit_goal : 25;
                assert(p_size >= 25);
                mpz_ui_pow_ui(p, 2, p_size);
                mpz_nextprime(p, p);
                mpz_mul(n, n, p);
                bit_goal -= p_size;
            }
            mpz_set(p, n);

            // Smaller prime for fast nextprime.
            // Large enough to avoid being found with trial division.
            mpz_ui_pow_ui(q, 2, 25);
            mpz_nextprime(q, q);

            double t = 0;
            size_t count = 0;
            // time a reasonable number (or for 5 seconds)
            for (; count < 15 || t < 5; count++) {
                mpz_mul(n, p, q);
                assert( mpz_sizeinbase(n, 2) >= bits );

                auto  s_start_t = high_resolution_clock::now();

                assert( mpz_probab_prime_p(n, 25) == 0 );

                t += duration<double>(high_resolution_clock::now() - s_start_t).count();
                mpz_nextprime(q, q);
            }

            printf("Estimating PRP/s: %ld / %.2f = %.1f/s vs polyfit estimate of %.1f/s\n",
                    count, t, count / t, 1 / t_estimate);
            t_estimate = t / count;

            mpz_clears(n, p, q, nullptr);
        }
    }

    return t_estimate;
}


// See misc/benchmark.cpp
static double benchmark_primorial_modulo(const mpz_t& K, size_t count) {
    auto t_start = high_resolution_clock::now();

    uint64_t z = 0;

    // Benchmark suggest this doesn't really depend on size but use 34 bits
    // As this is size of "most" of primes (and > 32)
    uint64_t p = 1LL << 34;
    for (size_t i = 0; i < count; i++) {
        z += mpz_fdiv_ui(K, p + i);
    }

    double time = duration<double>(high_resolution_clock::now() - t_start).count();
    // Keep compiler from optimizing out loop.
    double eps = 1e-100 * z;

//    printf("\tK mod / s Estimated %ld/%.2g = %.2g\n", count, time, count / time);
    return time / count + eps;
}


/**
 * Count of numbers [0, SL] coprime to K
 */
static
size_t count_coprime_sieve(const struct Config& config) {
    vector<char> is_coprime(config.sieve_length + 1, true);
    for (auto prime : get_sieve_primes(config.p)) {
        if (config.d % prime == 0)
            continue;

        for (size_t i = 0; i < is_coprime.size(); i += prime)
            is_coprime[i] = false;
    }
    return std::count(is_coprime.begin(), is_coprime.end(), true);
}


/**
 * Count of numbers coprime to d less than end; sum( gcd(m, d) == 1 for m in range(n, n+i) )
 * Uses inclusion exclusion on prime factorization of d
 */
static
uint64_t _r_count_num_m(uint64_t n, const vector<int> &factors_d, int i) {
    if (n == 0) return 0;
    if (i < 0) return n;

    return _r_count_num_m(n, factors_d, i-1) - _r_count_num_m(n / factors_d[i], factors_d, i-1);
}

/**
 * Count number of m [ms, ms + mi) coprime to d
 */
size_t count_num_m(long ms, long mi, uint64_t d) {
    if (d == 1)
        return mi;

    if (ms + mi < 10000) {
        size_t count = 0;
        for (long m = ms; m < ms + mi; m++)
            count += (gcd(m, d) == 1);
        return count;
    }

    vector<int> D_factors;
    {
        uint64_t temp = d;
        for (long p = 2; p*p <= temp; p++) {
            while (temp % p == 0) {
                D_factors.push_back(p);
                temp /= p;
            }
        }
        if (temp > 1)
            D_factors.push_back(temp);
    }

    return _r_count_num_m(ms + mi - 1, D_factors, D_factors.size()-1) -
           _r_count_num_m(ms - 1,      D_factors, D_factors.size()-1);
}


/**
 * Vector of mi, such that gcd(config.mstart + mi, config.d)
 * Returns a copy, but copy is "fast" compared to cost of computing vector
 */
pair<vector<bool>, vector<uint32_t>> is_coprime_and_valid_m(const struct Config& config) {
    const uint64_t M_start = config.mstart;
    const uint64_t M_inc = config.minc;
    assert(M_inc < std::numeric_limits<uint32_t>::max());

    const uint32_t D = config.d;
    const vector<uint32_t> P_primes = get_sieve_primes(config.p);
    assert( P_primes.back() == config.p);

    vector<uint32_t> valid_mi;
    vector<bool> is_m_coprime(M_inc, 1);

    for (uint32_t p : P_primes) {
        if (D % p == 0) {
            // mark off any m = m_start + mi that shares factor with d
            uint64_t first = (p - (M_start % p)) % p;
            assert((M_start + first) % p == 0);
            for (uint64_t mi = first; mi < M_inc; mi += p) {
                is_m_coprime[mi] = 0;
            }
        }
    }

    // Slower than dynamic bitset, but fast enough
    size_t count = std::count(is_m_coprime.begin(), is_m_coprime.end(), 1);
    valid_mi.reserve(count);

    for (uint32_t mi = 0; mi < M_inc; mi++) {
        if (is_m_coprime[mi]) {
            assert(gcd(M_start + mi, D) == 1);
            valid_mi.push_back(mi);
        }
    }

    return {is_m_coprime, valid_mi};
}


pair<uint64_t, uint64_t> calculate_thresholds_method2(
        const struct Config config,
        size_t count_coprime_sieve,
        size_t valid_ms) {
    uint32_t sieve_interval = 2 * config.sieve_length + 1;

    // (small vs modulo_search)  MULT  vs  log2(MULT) * (M_inc/valid_ms)
    float SMALL_MULT = std::max(8.0, log(8) * config.minc / valid_ms);

    // (small vs medium)         valid_m  vs  count_coprime_sieve * (M_inc / prime)
    uint64_t MEDIUM_CROSSOVER_SMALL = 1.0 * count_coprime_sieve * config.minc / valid_ms;

    // (medium vs modulo_search)  count_coprime_sieve vs M*S/P * (log2(P) - log2(SL))
    float M_PER_P_CROSSOVER = 1.0 * config.minc * sieve_interval / count_coprime_sieve;
    // correct for how much work it takes to skip to next m
    float MEDIUM_MULT = std::max(1.9, 0.65 * log2(M_PER_P_CROSSOVER / count_coprime_sieve));
    uint64_t MEDIUM_CROSSOVER_SEARCH = MEDIUM_MULT * M_PER_P_CROSSOVER;

    // XXX: What would it look like to do this more dynamically?
    // Everytime prime >= next_mult run a couple through both MEDIUM & LARGE prime and choose faster.

    uint64_t SMALL_THRESHOLD = std::min((uint64_t) SMALL_MULT * sieve_interval, MEDIUM_CROSSOVER_SMALL);
    if (SMALL_THRESHOLD < sieve_interval) {
        SMALL_THRESHOLD = sieve_interval + 1;
    }

    uint64_t MEDIUM_THRESHOLD = std::max(SMALL_THRESHOLD, MEDIUM_CROSSOVER_SEARCH);
    MEDIUM_THRESHOLD = std::min(MEDIUM_THRESHOLD, config.max_prime);

    return {SMALL_THRESHOLD, MEDIUM_THRESHOLD};
}


double combined_sieve_method2_time_estimate(
        const struct Config& config,
        const mpz_t &K,
        uint64_t valid_ms,
        double prp_time_est) {
    // XXX: pull these from config file or somewhere
    const double INVERSES_SECS = 18e-9;
    const double MODULE_SEARCH_SECS = 125e-9;
    // much less important to correctly set.
    const double COUNT_VECTOR_BOOL_PER_SEC = 6871000500;
    // ~ `primesieve -t1 500e9 --dist 1e9'
    const double PRIME_RANGE_SEC = 0.26 / 1e9;

	const size_t coprimes = 2 * count_coprime_sieve(config);
	const auto THRESHOLDS = calculate_thresholds_method2(config, coprimes, valid_ms);
	const size_t s_threshold_primes = primepi_estimate(THRESHOLDS.first);
    const size_t m_threshold_primes = primepi_estimate(THRESHOLDS.second);
    const size_t expected_primes = primepi_estimate(config.max_prime);

	// Time to compute all (primes % K)
    const double K_log = _log(K);
    const double mod_time_est = benchmark_primorial_modulo(K, 100'000 * (K_log < 2000 ? 20 : 1));
    const double k_mod_time = expected_primes * mod_time_est;

	// Time for SMALL_THRESHOLD to MEDIUM_THRESHOLD
	const size_t inverses = (m_threshold_primes - s_threshold_primes) * coprimes;
	const double inverse_time = inverses * INVERSES_SECS;

	// Time for solving module_search
    const size_t interval = 2 * config.sieve_length + 1;
    const size_t expected_m_stops =
        (log(log(config.max_prime)) - log(log(THRESHOLDS.second))) * interval * config.minc;
	const size_t solves = (expected_m_stops + (expected_primes - m_threshold_primes));
    const double m_search_time = solves * MODULE_SEARCH_SECS;

    const size_t count_prints = 5 * (log10(config.max_prime) - 4);
    const double extra_time =
        // PrimePi takes ~0.3s / billion
        config.max_prime * PRIME_RANGE_SEC +
        // 5 prints per log10 * std::count(all_unknowns)
        count_prints * 1.0 * valid_ms * coprimes / COUNT_VECTOR_BOOL_PER_SEC;
    const double total_estimate = k_mod_time + m_search_time + inverse_time + extra_time;

    // Estimate still needs to account for:
    //      small primes
    //      marking off factors (small and large)

    if (config.verbose >= 2 && config.show_timing) {
        const double N_log = K_log + log(config.mstart);
        const double prob_prime = 1 / N_log - 1 / (N_log * N_log);
        const double estimated_prp_per_m = 1 / (prob_prime * log(config.max_prime) * exp(GAMMA));
        const double test_estimate = 2 * valid_ms * estimated_prp_per_m * prp_time_est;

        printf("Estimated misc (PrimePi, count unknown, ...) time: %.0f (%.1f%% total)\n",
            extra_time, 100.0 * extra_time / total_estimate);

        printf("Estimated K mod/s: %'.0f, estimated time for all mods: %.0f (%.1f%% total)\n",
            1 / mod_time_est, k_mod_time, 100.0 * k_mod_time / total_estimate);

        printf("Estimated modulo_searches(million): %ld, time: %.0f (%.1f%% total)\n",
                (expected_m_stops + expected_primes) / 1'000'000,
                m_search_time, 100.0 * m_search_time / total_estimate);

        printf("Estimated sieve time: %.0f seconds (%.2f hours) (%.3f%%)\n",
                total_estimate, total_estimate / 3600,
                100 * total_estimate / (test_estimate + total_estimate));
        printf("Estimated test  time: %.0f hours (%.1f%%)\n",
                test_estimate / 3600,
                100 * test_estimate / (test_estimate + total_estimate));

        printf("\n");
    }

    return total_estimate;
}


/**
 * Handles approx count of divisors by d
 * See "Optimizing Choice Of D" in THEORY.md
 * Return: Counts the number of coprimes of (N, i), -sl <= i <= sl
 */
std::tuple<double, uint32_t, double, double> count_K_d(const struct Config& config) {
    uint64_t K_mod_d;
    double N_log;
    const uint64_t d = config.d;
    {
        double K_log;
        mpz_t K;
        K_stats(config, K, nullptr, &K_log);

        // Looks a little silly (P# / d) % d
        K_mod_d = mpz_fdiv_ui(K, d);
        mpz_clear(K);

        N_log = K_log + log(config.mstart);
    }

    vector<uint32_t> P_primes = get_sieve_primes(config.p);
    assert( P_primes.back() == config.p );

    // Prob prime if no factor of number less than P
    double prob_prime_adj = prob_prime_coprime(config);

    if (config.verbose >= 3) {
        printf("prob_prime: %.6f => %.6f\n",
            1 / N_log - 1 / (N_log * N_log),
            prob_prime_adj);
    }

    // Find factors of D
    vector<uint32_t> D_primes;
    for (uint32_t prime : P_primes) {
        if (d % prime == 0)
            D_primes.push_back(prime);
    }

    const size_t sl = config.sieve_length;
    const size_t length = 2 * sl + 1;

    // composite from coprime K
    char compositeK[length];
    std::fill(compositeK, compositeK + length, false);

    for (uint32_t prime : P_primes) {
        if (d % prime != 0) {
            // mark off all multiples of prime
            uint32_t first = (sl % prime);
            for (size_t m = first; m < length; m += prime) {
                compositeK[m] = true;
            }
        }
    }

    double expected_length = 0;
    size_t expected_count = 0;
    double remaining_prob = 0;

    char composite[length];

    uint64_t m = config.mstart;

    // Periodic in d, but d might be X00'000'000 so limit to 5'000
    const uint64_t intervals = std::min(d, 5'000UL);
    size_t m_count = 0;
    for (; m_count < intervals; m++) {
        if (m >= config.mstart + config.minc) break; // Tested all values.
        if (d > 1 && gcd(m, d) > 1) continue;
        m_count++;

        // Reset to composites from coprime K
        std::copy(compositeK, compositeK + length, composite);

        // Handle d primes for this m
        for (uint32_t p : D_primes) {
            // -((m * K) - SL) % p => (m * K_mod_d + p - (sl % p)) % p
            assert(K_mod_d % p != 0);
            uint64_t first = (p - (((m % p) * (K_mod_d % p) + p - (sl % p)) % p)) % p;
            for (size_t mi = first; mi < length; mi += p) {
                composite[mi] = true;
            }
        }

        if (config.verbose >= 3 && m <= 6) {
            size_t count_unknown = std::count(composite + sl, composite + 2*sl, false);
            printf("%ld * %d#/%ld | %ld | ", m, config.p, d, count_unknown);

            for (int x = 0; (size_t) x <= sl; x++)
                if (!composite[sl + x])
                    printf("%d ", x);
            printf("\n");
        }

        for (int dir = -1; dir <= 1; dir += 2) {
            double expected = 0;
            double prob = 1.0;
            for (int x = 0; (size_t) x <= sl; x++) {
                if (!composite[sl + dir * x]) {
                    expected += x * prob * prob_prime_adj;
                    prob *= 1 - prob_prime_adj;
                    expected_count += 1;
                }
            }
            expected += sl * prob;
            expected_length += expected;
            remaining_prob += prob;
        }
    }
    return {expected_length / m_count,
            expected_count / (m_count * 2),
            remaining_prob / (m_count * 2),
            prob_prime_adj
    };
}


double prob_prime_and_stats(const struct Config& config, mpz_t &K) {
    int K_digits;
    double K_log;
    K_stats(config, K, &K_digits, &K_log);

    if (config.verbose >= 2) {
        // From Mertens' 3rd theorem
        double unknowns_after_sieve = 1 / (log(config.max_prime) * exp(GAMMA));
        const double N_log = K_log + log(config.mstart);
        const double prob_prime = 1 / N_log - 1 / (N_log * N_log);
        double prob_prime_after_sieve = prob_prime / unknowns_after_sieve;

        auto stats = count_K_d(config);
        size_t count_coprime_p = std::get<1>(stats);
        double prob_prime_coprime_p = std::get<3>(stats);
        double prob_gap_hypothetical = std::get<2>(stats);

        float expected = count_coprime_p * (prob_prime_coprime_p / prob_prime_after_sieve);
        printf("\n");
        printf("\texpect %.0f left = 2 * %.0f (%.3f%%) of %u after %ldM\n",
                2 * expected, expected,  100.0 * expected / (config.sieve_length + 1),
                config.sieve_length, config.max_prime/1'000'000);
        printf("\t%.3f%% of %d digit numbers are prime\n",
                100 * prob_prime, K_digits);
        printf("\t%.3f%% of tests should be prime (%.1fx speedup)\n",
                100 * prob_prime_after_sieve, 1 / unknowns_after_sieve);
        printf("\t~2x%.1f = %.1f PRP tests per m\n",
                1 / prob_prime_after_sieve, 2 / prob_prime_after_sieve);
        printf("\tsieve_length=%d is insufficient ~~%.3f%% of time\n",
                config.sieve_length, 100 * prob_gap_hypothetical);
        printf("\n");
    }

    return K_log;
}


/**
 * Change that a number near K is prime
 * GIVEN no factor of K or D => no factor of P#
 */
double prob_prime_coprime(const struct Config& config) {
    double N_log = calc_log_K(config) + log(config.mstart);
    double prob_prime_coprime_P = 1 / N_log - 1 / (N_log * N_log);

    // Adjust for prob_prime for no primes <= P
    for (auto prime : get_sieve_primes(config.p)) {
        prob_prime_coprime_P /= (1 - 1.0 / prime);
    }

    return prob_prime_coprime_P;
}


// Small sieve of Eratosthenes.
vector<uint32_t> get_sieve_primes(uint32_t n) {
    assert(n < 1'001'000); // Use libprimesieve for larger intervals

    vector<uint32_t> primes = {2};
    uint32_t half_n = n >> 1;
    vector<bool> is_prime(half_n + 1, true);

    for (uint32_t p = 3; p <= n; p += 2) {
        if (is_prime[p >> 1]) {
            primes.push_back(p);
            uint64_t p2 = p * p;
            if (p2 > n) break;

            for (uint32_t m = p2 >> 1; m <= half_n; m += p)
                is_prime[m] = false;
        }
    }
    for (uint32_t p = primes.back() + 2; p <= n; p += 2) {
        if (is_prime[p >> 1])
            primes.push_back(p);
    }
    return primes;
}


bool is_prime_brute(uint32_t n) {
    if ((n & 1) == 0)
        return false;
    for (uint32_t p = 3; p * p <= n; p += 2)
        if (n % p == 0)
            return false;
    return true;
}


size_t primepi_estimate(uint64_t max_prime) {
    // Lookup primepi for common max_prime values.
    if (common_primepi.count(max_prime)) {
        return common_primepi.at(max_prime);
    }
    return 1.04 * max_prime / log(max_prime);

}



void Args::show_usage(char* name, Pr program) {
    Config defaults;

    cout << "Usage: " << name << endl;
    cout << "[REQUIRED]" << endl;
    cout << "  -p <p>" << endl;
    cout << "  -d <p>" << endl;
    cout << "  --mstart <start>" << endl;
    cout << "  --minc   <int>" << endl;
    cout << "OR" << endl;
    cout << "  -u, --unknown-filename <filename>" << endl;
    cout << "    parse p, d, mstart, minc, sieve-length, max-prime from filename" << endl;
    cout << endl;
    cout << "[OPTIONALLY]" << endl;
if (program == Pr::SIEVE || program == Pr::STATS) {
    cout << "  -t, --threads N" << endl;
    cout << "    Use N threads (OpenMP)" << endl;
}
    cout << "  --min-merit <min_merit>" << endl;
    cout << "    only display prime gaps with merit >= min_merit" << endl;
if (program == Pr::TEST_GPU) {
    cout << "  --mskip <start at this m>" << endl;
    cout << "    allows for partial resume of a unknown-file" << endl;
}
if (program == Pr::SIEVE) {
    cout << "  --sieve-length" << endl;
    cout << "    how large the positive/negative sieve arrays should be" << endl;
    cout << "  --max-prime" << endl;
    cout << "    use primes <= max-prime (in millions) for checking composite" << endl;
    cout << endl;
    cout << "  --save-unknowns" << endl;
    cout << "    save unknowns to a temp file where they are processed in a 2nd pass." << endl;
    cout << "  --rle" << endl;
    cout << "    save in run-length encoded format" << endl;
    cout << "  --bitcompress" << endl;
    cout << "    save in new bitcompressed format" << endl;
    cout << "  --maxmem <max memory in GB>" << endl;
    cout << "    Combined sieve will print a warning if it's likely to use more memory." << endl;
}
    cout << endl;
    cout << "[OPTIONAL]" << endl;
if (program == Pr::SIEVE || program == Pr::STATS) {
    cout << "  --search-db" << endl;
    cout << "    Database for this project (Default: '" << defaults.search_db << "')" << endl;
    cout << "  --prime-gaps-db" << endl;
    cout << "    Prime gap prime gap search db (Default: '" << defaults.gaps_db << "')" << endl;
}
    cout << endl;
    cout << "  -q, --quiet" << endl;
    cout << "    suppress some status output (twice for more suppression)" << endl;
    cout << "  -h, --help" << endl;
    cout << "    print this help message" << endl;
    cout << endl;
    cout << "calculates prime_gaps for (mstart + mi) * p#/d, mi <= minc " << endl;
}


std::string Args::gen_unknown_fn(const struct Config& config, std::string suffix) {
    if (!config.unknown_filename.empty()) {
        // dirname (unknown/ or input) handled by parse.

        // re-generating unknown_fn can cause issue (with losing dirname)
        return config.unknown_filename;
    }

    return "unknowns/" +
           std::to_string(config.p) + "_" +
           std::to_string(config.d) + "_" +
           std::to_string(config.mstart) + "_" +
           std::to_string(config.minc) + "_s" +
           std::to_string(config.sieve_length) + "_l" +
           std::to_string(config.max_prime / 1'000'000) + "M" +
           (config.method1 ? ".m1" : "") +
           suffix;
}


int Args::guess_compression(const struct Config& config, std::ifstream& unknown_file) {
    // Get current position
    int pos = unknown_file.tellg();
    assert(pos == 0);

    // 100 characters gets past <m>: -count count | <OFFSETS>

    // Check that <m> is <m> not <mi>
    {
        int64_t mtest = -1;
        unknown_file >> mtest;
        assert(mtest >= 0);

        int64_t m = config.mstart;
        for (; gcd(m, config.d) > 1; m++);

        if (m != mtest) {
            cout << endl;
            cout << "file format has changed," << endl;
            cout << "lines should start with <mstart + mi> not <mi>" << endl;
            cout << "\texpected: " << m << " found: " << mtest << endl;
            cout << "you can add <mstart> to each line, recreate, or git checkout 74241f7c" << endl;
            cout << "Sorry" << endl;
            exit(1);
        }
    }

    char t[100] = {0};
    unknown_file.read(t, sizeof(t) - 1);
    unknown_file.seekg(pos, std::ios_base::beg);

    // Compression 2 uses || seperator
    for (size_t i = 0; i < strlen(t) - 1; i++) {
        if (t[i] == '|') {
            assert(i + 1 < strlen(t));
            if (t[i + 1] == '|')
                return 2;
            break;
        }
    }

    bool has_space = false;
    bool has_high_range = false;
    for (size_t i = 50; i < strlen(t) - 1; i++) {
        has_space      |= t[i] == ' ' && t[i+1] != '|' && t[i-1] != '|';
        has_high_range |= t[i] > '9';
    }
    assert(has_space ^ has_high_range);

    return has_high_range ? 1 : 0;
}


Config Args::argparse(int argc, char* argv[], Pr program) {
    // NOTE: Remember to add to getopt_long(argc, argv, OPTIONS_STRING, ...) below
    static struct option long_options[] = {
        {"p",                required_argument, 0,  'p' },
        {"d",                required_argument, 0,  'd' },

        {"mstart",           required_argument, 0,   1  },
        {"minc",             required_argument, 0,   2  },
        {"mskip",            required_argument, 0,  16  },

        {"unknown-filename", required_argument, 0,  'u' },

        {"sieve-length",     required_argument, 0,   4  },
        {"max-prime",        required_argument, 0,   5  },

        {"threads",          required_argument, 0,  't' },

        {"min-merit",        required_argument, 0,   3  },

        {"save",             no_argument,       0,   7  },
        {"save-unknowns",    no_argument,       0,   7  },
        {"rle",              no_argument,       0,  13  },
        {"bitcompressed",    no_argument,       0,  15  },
        {"uncompressed",     no_argument,       0,  17  },
        {"max-mem",          required_argument, 0,  14  },

        {"search-db",        required_argument, 0,   9  },
        {"prime-gaps-db",    required_argument, 0,  10  },

        {"method1",          no_argument,       0,   8  },

        // Secret option
        {"hide-timing",      no_argument,       0,  11  },
        {"testing",          no_argument,       0,  12  },

        {"quiet",            no_argument,       0,  'q' },
        {"help",             no_argument,       0,  'h' },
        {0,                  0,                 0,   0  }
    };

    Config config;
    config.valid = 1;

    int option_index = 0;
    char c;
    while ((c = getopt_long(argc, argv, "qhp:d:u:t:", long_options, &option_index)) >= 0) {
        switch (c) {
            case 'h':
                show_usage(argv[0], program);
                exit(0);

            case 'q':
                config.verbose--;
                break;

            case 'p':
                config.p = atoi(optarg);
                break;
            case 'd':
                config.d = atoi(optarg);
                break;
            case 1:
                config.mstart = atoll(optarg);
                break;
            case 2:
                config.minc = atoll(optarg);
                break;

            case 16:
                config.mskip = atoll(optarg);
                break;

            case 'u':
                {
                    // Ugh, change to c++17 filesystem::path at some later point
                    char* t = strdup(optarg);
                    string dir = dirname(t);
                    free(t);

                    char* copy = strdup(optarg);
                    t = basename(optarg);
                    assert(*t != 0);
                    assert(strcmp(t, ".") != 0);

                    // Add "unknowns/" if no directory present
                    dir = (dir == ".") ? UNKNOWNS_DIR : dir;
                    config.unknown_filename = dir + "/" + t;

                    assert( std::count(t, t + strlen(t), '_')  == 5);

                    config.p = atoi(t);
                    t = std::strchr(t, '_');
                    t++;

                    config.d = atoi(t);
                    t = std::strchr(t, '_');
                    t++;

                    config.mstart = atoll(t);
                    t = std::strchr(t, '_');
                    t++;

                    config.minc = atoll(t);
                    t = std::strchr(t, '_');
                    assert( t[0] == '_' && t[1] == 's' );
                    t += 2;

                    config.sieve_length = atoi(t);
                    t = std::strchr(t, '_');
                    assert( t[0] == '_' && t[1] == 'l' );
                    t += 2;

                    config.max_prime = atol(t) * 1'000'000;
                    t = std::strchr(t, 'M');

                    config.method1 = (t[3] == '1');

                    assert( std::strcmp(t, "M.txt") == 0 || std::strcmp(t, "M.m1.txt") == 0 );
                    free(copy);
                }
                break;

            case 't':
                config.threads = atoi(optarg);
                break;

            case 3:
                config.min_merit = atof(optarg);
                break;
            case 4:
                config.sieve_length = atoi(optarg);
                break;
            case 5:
                config.max_prime = atol(optarg) * 1'000'000;
                break;

            case 7:
                config.save_unknowns = true;
                break;
            case 8:
                config.method1 = true;
                break;

            case 9:
                config.search_db = optarg;
                assert_file_exists(optarg);
                break;
            case 10:
                config.gaps_db = optarg;
                assert_file_exists(optarg);
                break;

            case 11:
                config.show_timing = false;
                break;

            case 12:
                config.testing = true;
                break;

            case 13:
                config.compression = 1;
                break;
            case 15:
                config.compression = 2;
                break;
            case 17:
                config.compression = 3;
                break;

            case 14:
                config.max_mem = atol(optarg);
                break;

            case 0:
                printf("option %s arg %s\n", long_options[option_index].name, optarg);
                config.valid = 0;
                break;
            case '?':
                config.valid = 0;
                break;
            default:
                config.valid = 0;
                printf("getopt returned \"%d\"\n", c);
        }
    }

    if (optind < argc) {
        config.valid = 0;
        printf("unknown positional argument: ");
        while (optind < argc) {
            printf("%s ", argv[optind++]);
        }
        printf("\n");
    }

    if (config.testing) {
        config.save_unknowns = false;
    }

    // ----- Validation
#ifdef RLE
#error "Don't build with RLE=1 anymore instead pass --rle or let combined_sieve auto select"
#endif

    if (config.mstart <= 0) {
        config.valid = 0;
        cout << "mstart must be greater than 0: " << config.mstart << endl;
    }

    int64_t last_m = config.mstart + config.minc;
    if (last_m <= 0 || last_m > 10'000'000'001 ) {
        config.valid = 0;
        cout << "mstart + minc must be <= 10e9" << endl;
    }

    if (config.minc <= 0) {
        config.valid = 0;
        cout << "minc must be greater than 0: " << config.minc << endl;
    }
    if (config.minc >= std::numeric_limits<int32_t>::max()) {
        config.valid = 0;
        cout << "minc must be less than 2B " << config.minc << endl;
    }

    if (config.max_prime > 40'000'000'000'000) {
        /**
         * improved primeiterator can find all primes < 1T in ~400s
         * mpz_mod is slow part, modulo_search is always fast.
         * module_search_..._large avoids overflow.
         */
        config.valid = 0;
        cout << "max_prime > 40 Trillion not supported" << endl;
    }

    /**
     * Overflow happens when base_r * (M + mi) > int64.
     * Given base_r < p, mi < max_m this happens rarely when log2(...) = 65
     * But more and more frequently after.
     * For 1-2% performance modulo_search_euclid_all_large handles these safely.
     */
    if (config.method1) {
        uint64_t max_m = std::numeric_limits<uint64_t>::max() / config.max_prime;
        if (max_m < 1000 || (max_m + 1000) <= (size_t) last_m) {
            config.valid = 0;
            printf("max_prime * last_m(%ld) would overflow int64, log2(...) = %.3f\n",
                last_m, log2(1.0 * last_m * config.max_prime));
        }
    }

    {
        // check if p is valid
        bool valid = config.p < 1'000'0000;
        for (size_t t = 2; valid && t*t <= config.p; t++) {
            valid = (config.p % t) > 0;
        }

        if (!valid) {
            config.valid = 0;
            cout << "p# not prime (p=" << config.p << ")" << endl;
            exit(1);
        }
    }

    {
        // Check that SL % 500 = 0 or (SL, K) = 1
        size_t SL = config.sieve_length;
        if (SL % 500 != 0) {
            for (size_t p = 2; p < config.p; p += 1 + (p > 2)) {
                if (is_prime_brute(p)) {
                    if (SL % p == 0 && config.d % p != 0) {
                        config.valid = 0;
                        cout << "SL=" << SL << " not coprime (p=" << p << ")" << endl;
                    }
                }
            }
        }
    }

    if (config.d <= 0) {
        config.valid = 0;
        cout << "d must be greater than 0: " << config.d << endl;
    }

    if (config.threads <= 0 || config.threads > 16) {
        config.valid = 0;
        cout << "invalid number of threads(" << config.threads << ") only 1 - 16 supported" << endl;
    }

    if (config.valid == 0) {
        cout << endl;
    }

    return config;
}


DB::DB(const char* path) {
    assert_file_exists(path);

    if (sqlite3_open(path, &db) != SQLITE_OK) {
        printf("Can't open database(%s): %s\n", path, sqlite3_errmsg(db));
        exit(1);
    }
};


uint64_t DB::config_hash(const struct Config& config) {
    // Hash the config to a uint64
    uint64_t hash =    config.mstart;
    hash = hash * 31 + config.minc;
    hash = hash * 31 + config.p;
    hash = hash * 31 + config.d;
    hash = hash * 31 + config.sieve_length;
    hash = hash * 31 + config.max_prime;
    return hash;
}


BitArrayHelper::BitArrayHelper(const struct Config& config, const mpz_t &K) {
    const unsigned int SL = config.sieve_length;
    const unsigned int SIEVE_INTERVAL = 2 * SL + 1;
    const unsigned int D = config.d;

    SL_mod_d = SL % D;
    neg_K_mod_d = mpz_cdiv_ui(K, D);
    if (D > 1) {
        assert(neg_K_mod_d != 0);
    }

    is_offset_coprime.resize(SIEVE_INTERVAL);
    std::fill(is_offset_coprime.begin(), is_offset_coprime.end(), 1);

    for (auto prime : get_sieve_primes(config.p)) {
        if (config.d % prime == 0) {
            D_primes.push_back(prime);
        } else {
            P_primes.push_back(prime);
        }
    }
    assert(D_primes.size() <= 9);  // 23# > 2^32
    assert( (config.d == 1) || (!D_primes.empty() && D_primes.front() >= 2) );

    for (uint32_t prime : P_primes) {
        uint32_t first = SL % prime;
        assert( 0 <= first && first < prime );
        assert( (SL - first) % prime == 0 );
        for (size_t x = first; x < SIEVE_INTERVAL; x += prime) {
            is_offset_coprime[x] = 0;
        }
    }

    // assume m % 2 == 1 => X % 2 == 0
    is_offset_coprime_even = is_offset_coprime;
    {
        for (size_t x = (SL + 1) % 2; x < SIEVE_INTERVAL; x += 2) {
            is_offset_coprime_even[x] = 0;
        }
    }

    for (int x = -SL; x <= (signed) SL; x++) {
        bool test = (mpz_gcd_ui(nullptr, K, abs(x)) == 1);
        bool test_array = is_offset_coprime[x + SL];
        assert(test == test_array);
        if (test) {
            coprime_X.push_back(SL + x);
            if (D % 2 == 0 && x % 2 == 0) {
                coprime_X_even.push_back(SL + x);
            }
        }
    }

    assert(coprime_X.size() == (unsigned) std::count(
        is_offset_coprime.begin(), is_offset_coprime.end(), 1));
    assert(coprime_X_even.size() == (unsigned) std::count(
        is_offset_coprime_even.begin(), is_offset_coprime_even.end(), 1));
};


// TODO: debup with load_and_verify_unknowns in gap_test_simple
/** Parse line (potentially with rle) to two positive lists */
int64_t parse_unknown_line(
        const struct Config& config,
        const BitArrayHelper& helper,
        uint64_t m_expected,
        std::istream& input_line,
        vector<int32_t>& unknown_prev,
        vector<int32_t>& unknown_next) {

    int unknown_l = 0;
    int unknown_u = 0;

    // Read a line from the file
    {
        int64_t m_test = -1;
        input_line >> m_test;
        assert( m_test >= 0 );

        if (config.threads == 1) {
            assert( (size_t) m_test == m_expected );
        }

        std::string delim = "ERROR";
        char delim_char;
        input_line >> delim;
        assert( delim == ":" );

        if (config.compression == 2) {
            const unsigned int SL = config.sieve_length;
            const unsigned int SIEVE_INTERVAL = 2 * SL + 1;

            const bool d_even = config.d % 2 == 0;
            vector<char> is_offset_fully_coprime(
                d_even ? helper.is_offset_coprime_even : helper.is_offset_coprime);

            // Logic borrowed from combined_sieve
            bool centerOdd = d_even && (m_test & 1);
            bool lowIsEven = centerOdd == (SL & 1);

            int num_coprimes = (d_even ? helper.coprime_X_even : helper.coprime_X).size();
            for (uint32_t d : helper.D_primes) {
                if (d == 2) continue;  // Handled by is_offset_coprime_even / coprime_X_even

                // First multiple = -(m * K - SL) % d = (m * -K + SL) % d
                uint64_t first = (m_test * helper.neg_K_mod_d + SL) % d;

                bool evenFromLow = (first & 1) == 0;
                bool firstIsEven = lowIsEven == evenFromLow;
                if (firstIsEven) {
                    assert( (first >= SIEVE_INTERVAL) || is_offset_fully_coprime[first] == 0 );
                    first += d;
                }
                uint32_t shift = 2 * d;

                for (uint64_t mult = first; mult < SIEVE_INTERVAL; mult += shift) {
                    if (is_offset_fully_coprime[mult]) {
                        is_offset_fully_coprime[mult] = 0;
                        num_coprimes -= 1;
                    }
                }
            }

            int unknown_total = -1;
            input_line >> unknown_total;
            assert(unknown_total > 0);

            int bytes_check = -1;
            input_line >> bytes_check;
            assert(bytes_check > 0);

            // This helps verify that num_coprimes calculation was correct
            int bytes_needed = (num_coprimes + 6) / 7;
            assert(bytes_check == bytes_needed);

            input_line >> delim;
            assert( delim == "||" );
            delim_char = input_line.get(); // get space character
            assert( delim_char == ' ');

            char buffer[bytes_needed];
            input_line.read(buffer, bytes_needed);

            // NOTE: on average these will be 50% full, 60% tries to avoid a final resize.
            unknown_prev.reserve(unknown_total * 6 / 10);
            unknown_next.reserve(unknown_total * 6 / 10);

            // Note: I tried a clever double for loop to avoid the if (bit == 7)
            // and it hurt performance as is_offset_fully_coprime[x] is called
            // the same number of times

            unsigned char b = buffer[0];
            uint8_t bits = 0;
            uint16_t index = 0;
            for (int32_t x : d_even ? helper.coprime_X_even : helper.coprime_X) {
                if (is_offset_fully_coprime[x]) {
                    if (bits == 7) {
                        b = buffer[++index];
                        bits = 0;
                    }
                    if ((b & 1) == 0) {
                        if (x <= (signed) SL) {
                            unknown_prev.push_back(-x + SL);
                        } else {
                            unknown_next.push_back(x - SL);
                        }
                    }
                    b >>= 1;
                    bits += 1;
                }
            }

            size_t unknowns_found = unknown_prev.size() + unknown_next.size();

            assert((unsigned) unknown_total == unknowns_found);
            assert(input_line.peek() == '\n' || input_line.peek() == EOF);

            // Reverse unknown_prev
            std::reverse(unknown_prev.begin(), unknown_prev.end());
            return m_test;

        } else if (config.compression == 0 || config.compression == 1) {
            input_line >> unknown_l;
            unknown_l *= -1;
            input_line >> unknown_u;

            input_line >> delim;
            assert( delim == "|" );
            delim_char = input_line.get(); // get space character
            assert( delim_char == ' ');

            unsigned char a, b;
            int c = 0;
            for (int k = 0; k < unknown_l; k++) {
                if (config.compression == 1) {
                    // Read bits in pairs (see save_unknowns_method2)
                    a = input_line.get();
                    b = input_line.get();
                    c += (a - 48) * 128 + (b - 48);
                } else {
                    input_line >> c;
                    c *= -1;
                }
                unknown_prev.push_back((unsigned) c);
            }

            input_line >> delim;
            assert( delim == "|" );
            delim_char = input_line.get(); // get space character
            assert( delim_char == ' ');

            c = 0;
            for (int k = 0; k < unknown_u; k++) {
                if (config.compression) {
                    a = input_line.get();
                    b = input_line.get();
                    c += (a - 48) * 128 + (b - 48);
                } else {
                    input_line >> c;
                }
                unknown_next.push_back((unsigned) c);
            }

            assert( unknown_l >= 0 && unknown_u >= 0 );
            assert( (size_t) unknown_l == unknown_prev.size() );
            assert( (size_t) unknown_u == unknown_next.size() );

            //assert( is_sorted(unknowns_prev.begin(), unknowns_prev.end()) );
            //assert( is_sorted(unknowns_next.begin(), unknowns_next.end()) );

            return m_test;
        }

        cout << "Unsupported config.compression(" << config.compression << ")" << endl;
        assert(false); // bad config.compression
    }
}
