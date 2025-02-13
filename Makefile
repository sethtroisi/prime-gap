# Copyright 2020 Seth Troisi
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


OBJS	= gap_common.o modulo_search.o gap_test_common.o
OUT	= combined_sieve gap_stats gap_test_simple benchmark benchmark_google
CC	= g++
CFLAGS	= -Wall -Werror -g -fopenmp
NVCC	= nvcc
CUDA_FLAGS	= -Xcompiler -Wall -Xcompiler -Werror -g -Xcompiler -fopenmp
BITS    = 1024

LDFLAGS	= -lgmp -lsqlite3
# Need for local gmp / primesieve
LDFLAGS+= -L /usr/local/lib

DEFINES =
ifdef VALIDATE_FACTORS
DEFINES += -DGMP_VALIDATE_FACTORS
endif
ifdef VALIDATE_LARGE
DEFINES += -DGMP_VALIDATE_LARGE_FACTORS
endif

all: $(OUT)

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS) $(DEFINES)

%_cuda.o: %.cpp
	$(NVCC) -c -o $@ $< $(CUDA_FLAGS) $(DEFINES)

combined_sieve: combined_sieve.cpp $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) -lprimesieve $(LDFLAGS) $(DEFINES)

combined_sieve_small: combined_sieve_small.cpp gap_common.o gap_test_common.o
	$(CC) -o $@ $^ $(CFLAGS) -lprimesieve $(LDFLAGS) $(DEFINES)

gap_stats: gap_stats.cpp gap_common.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

gap_test_simple: gap_test_simple.cpp gap_common.o gap_test_common.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

gap_test_gpu: gap_test_gpu.cu miller_rabin.h gap_common_cuda.o gap_test_common_cuda.o combined_sieve_small_cuda.o
	nvcc -o $@ -I../CGBN/include \
		-DGPU_BITS=$(BITS) \
		$(filter-out miller_rabin.h, $^) \
		-arch=sm_61 $(CUDA_FLAGS) \
		$(filter-out -fopenmp, $(LDFLAGS)) -lprimesieve

benchmark: misc/benchmark.cpp modulo_search.o
	$(CC) -o $@ $^ $(CFLAGS) -lprimesieve $(LDFLAGS) -I.

benchmark_google: misc/benchmark_google.cpp modulo_search.o
	$(CC) -o $@ $^ $(CFLAGS) -std=c++14 -lprimesieve -lbenchmark -lpthread -I.


.PHONY: clean


clean:
	rm -f $(OBJS) $(OUT) gap_test_gpu gap_common_cuda.o gap_test_common_cuda.o
