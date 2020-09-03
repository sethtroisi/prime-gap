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

#pragma once

#include <cassert>
#include <functional>


uint32_t modulo_search_brute(uint32_t p, uint32_t A, uint32_t L, uint32_t R);
uint32_t modulo_search_euclid_small(uint32_t p, uint32_t a, uint32_t l, uint32_t r);

uint64_t modulo_search_euclid(uint64_t p, uint64_t a, uint64_t l, uint64_t r);
uint64_t modulo_search_euclid_gcd(
        uint64_t M, uint64_t D, uint64_t max_m, uint64_t SL,
        uint64_t prime, uint64_t base_r);

void modulo_search_euclid_all(
        uint64_t M, uint64_t max_m, uint64_t SL,
        uint64_t prime, uint64_t base_r,
        std::function<void (uint64_t)> lambda);
