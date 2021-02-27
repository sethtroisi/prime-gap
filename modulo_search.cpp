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

/**
 * Performs modulo_search (Equation 3 from paper)
 *
 * Find m such that
 *      L <= m * base_r <= R  mod P
 *
 * Note that m is always smaller than P
 *
 * base_r < P, m_end * P < 64 bits
 *      modulo_search_euclid
 *          works for P < 64 bits, because it doesn't use M_start
 */


#include "modulo_search.h"

#include <cassert>
#include <iostream>
#include <vector>

using std::vector;


static
uint32_t gcd(uint32_t a, uint32_t b) {
    if (b == 0) return a;
    return gcd(b, a % b);
}


uint32_t modulo_search_brute(uint32_t p, uint32_t A, uint32_t L, uint32_t R) {
    // A + p must not overflow.
//    uint32_t check = p + A;
//    assert( check > p );
//    assert( check > A );

/*
    uint32_t temp = 0;
    for (uint32_t i = 1; ; i++) {
        temp += A;
        if (temp >= p) temp -= p;
        if (L <= temp && temp <= R) {
            return i;
        }
    }
*/

    uint32_t goal = R - L;
    uint32_t temp = p - L;
    for (uint32_t i = 1; ; i++) {
        temp += A;
        if (temp >= p) temp -= p;
        if (temp <= goal) {
            return i;
        }
    }
}


uint32_t modulo_search_euclid_small(uint32_t p, uint32_t a, uint32_t l, uint32_t r) {
    if (a > (p >> 1)) {
        std::swap(l, r);
        a = p - a; l = p - l; r = p - r;

        // Seems like this should be 1 or 2 ops faster
        // uint32_t t = p - l
        // a = p - a; l = p - r; r = t;
    }

    // Check if l <= ceil(l/a)*a <= r
    uint32_t l_div = l / a;
    uint32_t l_mod = l - l_div * a;
    uint32_t r_div = r / a + (l_mod == 0); // Wrong but faster than (l_div < r_div) || (l_div == 0)
    uint32_t r_mod = r - r_div * a;
    if (l_div < r_div) {
        return l_div + (l_mod > 0);
    }

    uint32_t new_a = a - (p % a);
    uint64_t k = modulo_search_euclid_small(a, new_a, l_mod, r_mod);

    uint64_t tl = k * p + l;
    uint64_t mult = (tl - 1) / a + 1;

    return mult;
}


// k1 = modulo(a1, new_a, l % a, r % a)
//      k2 = modulo
//      tl_inner = k * a1 + l
//      mult = (tl_inner - 1) / new_a + 1
//      return mult
// tl = k1 * p + l
//    = ((tl_inner - 1) / new_a + 1) * p + l
// mult = (tl - 1) / new_a + 1
//      = (((tl_inner - 1) / new_a + 1) * p + l - 1) / new_a + 1


uint64_t modulo_search_euclid(uint64_t p, uint64_t a, uint64_t l, uint64_t r) {
/*
    assert( 0 <= a && a < p );
    assert( 0 <= l && l < p );
    assert( 0 <= r && r < p );
    assert(      l <= r     );
*/
    if (l == 0) return 0;

    if (p < 0xFFFFFFFF) {
        return modulo_search_euclid_small(p, a, l, r);
    }

    // 2 * a > p, but avoids int max
    if (a > (p >> 1)) {
        std::swap(l, r);
        a = p - a; l = p - l; r = p - r;
    }

    uint64_t l_div = l / a;
    uint64_t l_mod = l - l_div * a;
    uint64_t r_div = r / a + (l_mod == 0); // Wrong but faster than (l_div < r_div) || (l_div == 0)
    uint64_t r_mod = r - r_div * a;
    if (l_div < r_div) {
        return l_div + (l_mod > 0);
    }

    // reduce to simpler problem
    uint64_t new_a = a - (p % a);
//    assert( 0 <= new_a && new_a < a );
    uint64_t k = modulo_search_euclid(a, new_a, l_mod, r_mod);

    __int128 tl = (__int128) p * k + l;
    uint64_t mult = (tl-1) / a + 1;

//    assert( mult < p );
/*
    __int128 tr = r + p * k;
    uint64_t test = mult * a;
    assert(       test <= tr );
    assert( tl <= test );
// */
    return mult;
}

uint64_t modulo_search_euclid_stack(uint64_t p, uint64_t a, uint64_t l, uint64_t r) {
    if (p < 0xFFFFFFFF) {
        return modulo_search_euclid_small(p, a, l, r);
    }

    static thread_local uint64_t stack[4 * 64];

    uint64_t stack_i = 0;
    uint64_t k = 0;

    assert(l < r);
    uint64_t delta = r - l;

    while (true) {
        // Check for modulo_search_euclid_small inside loop doesn't seem to help.

        // 2 * a > p, but avoids int max
        if (a > (p >> 1)) {
            l = l + delta;
            assert(l <= p);  // technically l < p, but l = p also magically works
            a = p - a;
            l = p - l;
        }

        uint64_t l_div = (l - 1) / a;
        uint64_t l_mod = l - l_div * a;
        uint64_t r_mod = l_mod + delta;
        if (r_mod >= a) {
            k = l_div + 1;
            break;
        }

        // reduce to simpler problem
        uint64_t div = p / a;
        uint64_t rem = p - div * a; // p % a
        uint64_t new_a = a - rem; // a - (p % a);

        stack[stack_i++] = rem;
        stack[stack_i++] = div;
        stack[stack_i++] = a;
        stack[stack_i++] = l;
        p = a;
        a = new_a;
        l = l_mod;
    }

    for (; stack_i > 0; ) {
        auto l = stack[--stack_i];
        auto a = stack[--stack_i];
        auto div = stack[--stack_i]; // p / a
        auto rem = stack[--stack_i]; // p % a

        auto k_1 = k * div + 1;
        uint64_t k_2;
        //if (rem < 0x7FFFFFFF && k < 0x7FFFFFFF) {
        //    k_2 = (k * rem + l - 1) / a;
        //} else {
            k_2 = ((__int128) k * rem + l - 1) / a;
        //}
        k = k_1 + k_2;

        //k = ((__int128)k * p + l - 1) / a + 1;
    }
    return k;
}

/**
 * find (base_r * (M + mi) + SL) <= 2 * SL
 * Skips solutions where gcd((M + mi), D) > 1
 * returns:
 *      0 <= mi < max_m (valid solution)
 *      max_m (no more solutions)
 */
/* Only used by Benchmark */
uint64_t modulo_search_euclid_gcd2(
        uint64_t M, uint64_t D, uint64_t max_m, uint64_t SL,
        uint64_t prime, uint64_t base_r) {

    uint64_t modulo = ((__int128) base_r * M + SL) % prime;
    uint64_t init = modulo;
    uint32_t two_SL = SL << 1;

    uint64_t mi = 0;
    while (true) {
        if ( modulo <= two_SL ) {
            if (gcd(D, (M + mi) % D) > 1) {
                mi += 1;
                if (mi >= max_m)
                    return max_m;

                modulo += base_r;
                if (modulo >= prime) modulo -= prime;
                continue;
            }
            return mi;
        }

        //assert( 0 <= modulo && modulo < prime );
        uint64_t low  = prime - modulo;
        uint64_t high = low + two_SL;
        //assert( 0 <= low && high < prime );

        mi += modulo_search_euclid(prime, base_r, low, high);
        if (mi >= max_m)
            return max_m;

        // 0 < base_r, mi < prime.

        /*
        __int128 mult = (__int128) base_r * mi;
        modulo = mult % prime;
        // Do with uint64_t instead of with __int128
        //modulo += init;
        if (modulo >= prime) modulo -= prime;
        */

        uint64_t mult = base_r * mi + init;
        modulo = mult % prime;

//        assert( (modulo <= SL) || (modulo + SL) >= prime );
    }
}

uint64_t modulo_search_euclid_gcd(
        uint64_t M, uint64_t D, uint64_t max_m, uint64_t SL,
        uint64_t prime, uint64_t base_r) {
    uint64_t mi = 0;

    uint64_t modulo = ((__int128) base_r * M) % prime;
    while (mi < max_m) {
        if ( (modulo <= SL) || (modulo + SL) >= prime) {
            if (gcd(M + mi, D) > 1) {
                mi += 1;
                modulo += base_r;
                if (modulo >= prime) modulo -= prime;
                continue;
            }
            return mi;
        }

        uint64_t shift = modulo + SL;
        assert( 0 <= shift && shift < prime );
        uint64_t low  = (prime - shift);
        uint64_t high = low + 2*SL;
        assert( 0 <= low && high < prime );

        uint64_t mi2 = modulo_search_euclid(prime, base_r, low, high);
        mi += mi2;
        if (mi >= max_m)
            return max_m;

        __int128 mult = (__int128) base_r * (M + mi);
        modulo = mult % prime;

        assert( (modulo <= SL) || (modulo + SL) >= prime );
    }
    return max_m;
}


/* Used when (M + max_m) * p fits in uint64 */
void modulo_search_euclid_all_small(
        uint32_t M, uint32_t max_m, uint32_t SL,
        uint64_t prime, uint64_t base_r,
        std::function<void (uint32_t, uint64_t)> lambda) {

    // (M + max_m) * p fits in uint64 (see gap_common.cpp)

    uint64_t mi = 0; // mi can be incremented by value up to prime.
    uint64_t modulo = (base_r * M + SL) % prime;
    uint64_t initial_modulo = modulo;
    uint32_t two_SL = SL << 1;

    while (true) {
        if ( modulo <= two_SL ) {
            if (mi >= max_m) return;

            // Note: combined sieve computes first = (base_r * (M_start + mi) + SL) % prime;
            // this is modulo.
            lambda(mi, modulo);
            mi += 1;

            if (mi >= max_m) return;

            modulo += base_r;
            if (modulo >= prime) modulo -= prime;
            continue;
        }

        /* using gcd2 optimizations */
        uint64_t low  = prime - modulo;
        uint64_t high = low + two_SL;

        mi += modulo_search_euclid(prime, base_r, low, high);
        if (mi >= max_m) return;

        /**
         * Guaranteed not to overflow.
         * base_r < P, mi < max_m
         * P * max_m < uint64_t
         */
        modulo = base_r * mi + initial_modulo;
        modulo %= prime;

        assert( modulo <= two_SL );
    }
}


/* Used when (M + max_m) * p greater than uint64 */
void modulo_search_euclid_all_large(
        uint32_t M, uint32_t max_m, uint64_t SL,
        uint64_t prime, uint64_t base_r,
        std::function<void (uint32_t, uint64_t)> lambda) {

    uint64_t mi = 0;
    uint64_t modulo = ((__int128) base_r * M + SL) % prime;
    uint64_t initial_modulo = modulo;
    uint32_t two_SL = SL << 1;
    while (true) {
        if ( modulo <= two_SL ) {
            if (mi >= max_m)
                return;

            lambda(mi, modulo);
            mi += 1;

            if (mi >= max_m) return;

            modulo += base_r;
            if (modulo >= prime) modulo -= prime;
            continue;
        }

        uint64_t low  = prime - modulo;
        uint64_t high = low + two_SL;

        mi += modulo_search_euclid_stack(prime, base_r, low, high);
        if (mi >= max_m) return;

        __int128 mult = (__int128) base_r * mi + initial_modulo;
        modulo = mult % prime;

        //assert( (modulo <= SL) || (modulo + SL) >= prime );
        assert( modulo <= two_SL );
    }
}

