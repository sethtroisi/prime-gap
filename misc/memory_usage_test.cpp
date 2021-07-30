// Copyright 2021 Seth Troisi
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

#include <cassert>
#include <chrono>
#include <cstdio>
#include <mutex>
#include <thread>
#include <vector>

using std::vector;
using std::mutex;

int main(int argc, char **argv) {
    // Input validation
    assert(1 <= argc && argc <= 3);

    size_t y = 300;
    size_t max_x = 1'000'000;

    if (argc >= 2) {
        max_x = atol(argv[1]);
    }
    if (argc >= 3) {
        y = atol(argv[2]);
    }

    size_t MB = 1024 * 1024;
    size_t guess = max_x * y / (8 * MB);

    size_t vb_size = sizeof(vector<bool>);
    size_t mt_size = sizeof(mutex);

    printf("Testing memory usage %ldx%ld ~= %ld MB\n", max_x, y, guess);
    printf("\tsizeof(vector<bool>) = %ld\n", vb_size);
    printf("\tsizeof(mutex)        = %ld\n", mt_size);
    printf("\tvector<bool>[%ld]    = %ld MB\n", max_x, vb_size * max_x / MB);

    assert( guess <= 10 * 1024 );

    vector<bool> *composite = new vector<bool>[max_x];

    for (size_t x = 10; x <= max_x; x = 23 * x / 10) {
        guess = x * y / (8 * MB);
        printf("\t%ldx%ld => %ld MB\n", x, y, guess);

        for (int i = 0; i < x; i++) {
            // This resizes old entries multiple times, should be fine.
            composite[i].resize(y, false);
        }

        int seconds = guess < 200 ? 1 : 5;
        std::this_thread::sleep_for(std::chrono::seconds(seconds));
    }

    return 0;
}
