/***

Copyright (c) 2018-2019, NVIDIA CORPORATION.  All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.

***/

#include <cassert>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <vector>

#include <cuda.h>
#include <gmp.h>
#include "cgbn/cgbn.h"

//#include "../utility/support.h"
void to_mpz(mpz_t r, uint32_t *x, uint32_t count) {
  mpz_import(r, count, -1, sizeof(uint32_t), 0, 0, x);
}

void from_mpz(mpz_t s, uint32_t *x, uint32_t count) {
  size_t words;

  if(mpz_sizeinbase(s, 2)>count*32) {
    fprintf(stderr, "from_mpz failed -- result does not fit\n");
    exit(1);
  }

  mpz_export(x, &words, -1, sizeof(uint32_t), 0, 0, s);
  while(words<count)
    x[words++]=0;
}

// support routines
void cuda_check(cudaError_t status, const char *action=NULL, const char *file=NULL, int32_t line=0) {
  // check for cuda errors

  if(status!=cudaSuccess) {
    printf("CUDA error occurred: %s\n", cudaGetErrorString(status));
    if(action!=NULL)
      printf("While running %s   (file %s, line %d)\n", action, file, line);
    exit(1);
  }
}
#define CUDA_CHECK(action) cuda_check(action, #action, __FILE__, __LINE__)


void cgbn_check(cgbn_error_report_t *report, const char *file=NULL, int32_t line=0) {
  // check for cgbn errors

  if(cgbn_error_report_check(report)) {
    printf("\n");
    printf("CGBN error occurred: %s\n", cgbn_error_string(report));

    if(report->_instance!=0xFFFFFFFF) {
      printf("Error reported by instance %d", report->_instance);
      if(report->_blockIdx.x!=0xFFFFFFFF || report->_threadIdx.x!=0xFFFFFFFF)
        printf(", ");
      if(report->_blockIdx.x!=0xFFFFFFFF)
      printf("blockIdx=(%d, %d, %d) ", report->_blockIdx.x, report->_blockIdx.y, report->_blockIdx.z);
      if(report->_threadIdx.x!=0xFFFFFFFF)
        printf("threadIdx=(%d, %d, %d)", report->_threadIdx.x, report->_threadIdx.y, report->_threadIdx.z);
      printf("\n");
    }
    else {
      printf("Error reported by blockIdx=(%d %d %d)", report->_blockIdx.x, report->_blockIdx.y, report->_blockIdx.z);
      printf("threadIdx=(%d %d %d)\n", report->_threadIdx.x, report->_threadIdx.y, report->_threadIdx.z);
    }
    if(file!=NULL)
      printf("file %s, line %d\n", file, line);
    exit(1);
  }
}
#define CGBN_CHECK(report) cgbn_check(report, __FILE__, __LINE__)

/////////////////////


// For this example, there are quite a few template parameters that are used to generate the actual code.
// In order to simplify passing many parameters, we use the same approach as the CGBN library, which is to
// create a container class with static constants and then pass the class.

// The CGBN context uses the following three parameters:
//   TBP             - threads per block (zero means to use the blockDim.x)
//   MAX_ROTATION    - must be small power of 2, imperically, 4 works well
//   SHM_LIMIT       - number of bytes of dynamic shared memory available to the kernel
//   CONSTANT_TIME   - require constant time algorithms (currently, constant time algorithms are not available)

// Locally it will also be helpful to have several parameters:
//   TPI             - threads per instance
//   BITS            - number of bits per instance
//   WINDOW_BITS     - number of bits to use for the windowed exponentiation

template<uint32_t tpi, uint32_t bits, uint32_t window_bits>
class mr_params_t {
  public:
  // parameters used by the CGBN context
  static const uint32_t TPB=0;                     // get TPB from blockDim.x
  static const uint32_t MAX_ROTATION=4;            // good default value
  static const uint32_t SHM_LIMIT=0;               // no shared mem available
  static const bool     CONSTANT_TIME=false;       // constant time implementations aren't available yet

  // parameters used locally in the application
  static const uint32_t TPI=tpi;                   // threads per instance
  static const uint32_t BITS=bits;                 // instance size
  static const uint32_t WINDOW_BITS=window_bits;   // window size
};


template<class params>
class miller_rabin_t {
  public:
  static const uint32_t window_bits=params::WINDOW_BITS;  // used a lot, give it an instance variable

  // Definition of instance_t type.  Note, the size of instance_t is not a multiple of 128-bytes, so array loads and stores
  // will not be 128-byte aligned.  It's ok in this example, because there are so few loads and stores compared to compute,
  // but using non-128-byte aligned types could be a performance limiter for different load/store/compute balances.

  typedef struct {
    cgbn_mem_t<params::BITS> candidate;
    uint32_t                 passed;
  } instance_t;

  typedef cgbn_context_t<params::TPI, params>    context_t;
  typedef cgbn_env_t<context_t, params::BITS>    env_t;
  typedef typename env_t::cgbn_t                 bn_t;
  typedef typename env_t::cgbn_local_t           bn_local_t;
  typedef typename env_t::cgbn_wide_t            bn_wide_t;

  context_t _context;
  env_t     _env;
  int32_t   _instance;

  __device__ __forceinline__ miller_rabin_t(cgbn_monitor_t monitor, cgbn_error_report_t *report, int32_t instance) : _context(monitor, report, (uint32_t)instance), _env(_context), _instance(instance) {
  }

  __device__ __forceinline__ void powm(bn_t &x, const bn_t &power, const bn_t &modulus) {
    bn_t       t;
    bn_local_t window[1<<window_bits];
    int32_t    index, position, offset;
    uint32_t   np0;

    // conmpute x^power mod modulus, using the fixed window algorithm
    // requires:  x<modulus,  modulus is odd

    // compute x^0 (in Montgomery space, this is just 2^BITS - modulus)
    cgbn_negate(_env, t, modulus);
    cgbn_store(_env, window+0, t);

    // convert x into Montgomery space, store into window table
    np0=cgbn_bn2mont(_env, x, x, modulus);
    cgbn_store(_env, window+1, x);
    cgbn_set(_env, t, x);

    // compute x^2, x^3, ... x^(2^window_bits-1), store into window table
    #pragma nounroll
    for(index=2;index<(1<<window_bits);index++) {
      cgbn_mont_mul(_env, x, x, t, modulus, np0);
      cgbn_store(_env, window+index, x);
    }

    // find leading high bit
    position=params::BITS - cgbn_clz(_env, power);

    // break the exponent into chunks, each window_bits in length
    // load the most significant non-zero exponent chunk
    offset=position % window_bits;
    if(offset==0)
      position=position-window_bits;
    else
      position=position-offset;
    index=cgbn_extract_bits_ui32(_env, power, position, window_bits);
    cgbn_load(_env, x, window+index);

    // process the remaining exponent chunks
    while(position>0) {
      // square the result window_bits times
      #pragma nounroll
      for(int sqr_count=0;sqr_count<window_bits;sqr_count++)
        cgbn_mont_sqr(_env, x, x, modulus, np0);

      // multiply by next exponent chunk
      position=position-window_bits;
      index=cgbn_extract_bits_ui32(_env, power, position, window_bits);
      cgbn_load(_env, t, window+index);
      cgbn_mont_mul(_env, x, x, t, modulus, np0);
    }

    // we've processed the exponent now, convert back to normal space
    cgbn_mont2bn(_env, x, x, modulus, np0);
  }

  __device__ __forceinline__ uint32_t miller_rabin(const bn_t &candidate, uint32_t *primes, uint32_t prime_count) {
    int       k, trailing, count;
    bn_t      x, power, minus_one;
    bn_wide_t w;

    cgbn_sub_ui32(_env, power, candidate, 1);
    trailing=cgbn_ctz(_env, power);
    cgbn_rotate_right(_env, power, power, trailing);

    for(k=0;k<prime_count;k++) {
      cgbn_set_ui32(_env, x, primes[k]);
      powm(x, power, candidate);
      cgbn_sub_ui32(_env, minus_one, candidate, 1);
      if(!cgbn_equals_ui32(_env, x, 1) && !cgbn_equals(_env, x, minus_one)) {
        // x is neither 1, nor candidate-1
        if(trailing==1)
          return k;

        // in the case of random data, trailing=ctz(candidate-1) is on average quite small.  If trailing is large,
        // then you might want to do the reduction with Barrett remainders or in Montgomery space.
        count=trailing;
        while(true) {
          cgbn_sqr_wide(_env, w, x);
          cgbn_rem_wide(_env, x, w, candidate);
          if(cgbn_equals(_env, x, minus_one))
            break;
          if(--count==0 || cgbn_equals_ui32(_env, x, 1))
            return k;
        }
      }
    }
    return prime_count;
  }

  /*
  __host__ static instance_t *generate_instances(uint32_t count) {
    instance_t *instances=(instance_t *)malloc(sizeof(instance_t)*count);
    int         index;

    mpz_t n;
    mpz_init(n);
    mpz_set_str(n,
        "0x13e0ab623eecd3771fc8f7306d00dc2475262f183746"
        "b69447ac48f39c61f3f097f6a4060fc2708e7a43a0ee2f"
        "29ab3908fd7426fc40b3485eb71ae3bee43aa60fe3d6dc"
        "b79fd846bfd03232de416fad5a576015ec3c75787cc907"
        "35a65feef58500ee3e85d331bde59f1a693cf11c219ae7"
        "82fe06d1f66b2936f139f8c23019a304c9545d6a4e9dfb"
        "6ad477b7c1a5a7e67ac0c4605ba7a7539dd43f09b67c70"
        "4b88d30eb6b6257d65f695", 0);

    for(index=0;index<count;index++) {
      //random_words(instances[index].candidate._limbs, params::BITS/32);
      from_mpz(n, instances[index].candidate._limbs, params::BITS/32);
      mpz_add_ui(n, n, 2);
      instances[index].candidate._limbs[0] |= 1;
      instances[index].passed=0;
    }
    return instances;
  }
  */

  __host__ static int32_t verify_first(instance_t *instances, uint32_t instance_count, uint32_t *primes, uint32_t prime_count) {
    int   index;
    mpz_t candidate;
    mpz_init(candidate);

    for(index=0;index<instance_count;index++) {
      if (instances[index].passed==prime_count) {
        to_mpz(candidate, instances[index].candidate._limbs, params::BITS/32);
        if (mpz_probab_prime_p(candidate, prime_count) == 0) {
          printf("MISMATCH AT INDEX: %d\n", index);
          printf("prime count=%d\n", instances[index].passed);
        }
        return index;
      }
    }
    return -1;
  }

  __host__ static void verify_results(instance_t *instances, uint32_t instance_count, uint32_t *primes, uint32_t prime_count) {
    int   index, total=0;
    mpz_t candidate;
    bool  gmp_prime, xmp_prime, match=true;

    mpz_init(candidate);

    for(index=0;index<instance_count;index++) {
      to_mpz(candidate, instances[index].candidate._limbs, params::BITS/32);
      xmp_prime=(instances[index].passed==prime_count);
      if (xmp_prime) {
        total++;

        gmp_prime=(mpz_probab_prime_p(candidate, prime_count)!=0);
        printf("N + %d\n", index);

        if(gmp_prime!=xmp_prime) {
          printf("MISMATCH AT INDEX: %d\n", index);
          printf("prime count=%d\n", instances[index].passed);
          match=false;
        }
      }
    }
    if(match)
      printf("All results matched\n");
    printf("%d probable primes found in %d random numbers\n", total, instance_count);
    printf("Based on an approximation of the prime gap, we would expect %0.1f primes\n", ((float)instance_count)*2/(0.69315f*params::BITS));
  }
};


template<class params>
__global__ void kernel_miller_rabin(cgbn_error_report_t *report, typename miller_rabin_t<params>::instance_t *instances, uint32_t instance_count, uint32_t *primes, uint32_t prime_count) {
  int32_t instance=(blockIdx.x*blockDim.x + threadIdx.x)/params::TPI;

  if(instance>=instance_count)
    return;

  typedef miller_rabin_t<params> local_mr_t;

  local_mr_t                     mr(cgbn_report_monitor, report, instance);
  typename local_mr_t::bn_t      candidate;
  uint32_t                       passed;

  cgbn_load(mr._env, candidate, &(instances[instance].candidate));

  passed=mr.miller_rabin(candidate, primes, prime_count);

  instances[instance].passed=passed;
}

uint32_t *generate_primes(uint32_t count) {
  uint32_t *list=(uint32_t *)malloc(sizeof(uint32_t)*count);
  int       test, current, index;

  // generate a list of primes
  list[0]=2;
  current=3;
  index=1;
  while(index<count) {
    for(test=1;test<index;test++)
      if(current%list[test]==0)
        break;
    if(test==index)
      list[index++]=current;
    current=current+2;
  }
  return list;
}

template<class params>
int32_t run_test(mpz_t &center, int sign, std::vector<int32_t> offsets, uint32_t prime_count) {
  typedef typename miller_rabin_t<params>::instance_t instance_t;

  instance_t          *gpuInstances;
  cgbn_error_report_t *report;
  uint32_t            *primes, *gpuPrimes;
  uint32_t              TPB=(params::TPB==0) ? 128 : params::TPB;
  uint32_t              TPI=params::TPI, IPB=TPB/TPI;

  primes=generate_primes(prime_count);

  size_t instance_count = offsets.size();
  if (instance_count == 0)
    return -1;

  instance_t *instances = (instance_t *) malloc(sizeof(instance_t)*instance_count);

  mpz_t prime_test;
  mpz_init(prime_test);
  for (size_t i = 0; i < instance_count; i++) {
     auto off = offsets[i];
     assert(off > 0);
     if (sign == 1) {
        mpz_add_ui(prime_test, center, off);
     } else {
        mpz_sub_ui(prime_test, center, off);
     }
     from_mpz(prime_test, instances[i].candidate._limbs, params::BITS/32);
  }
  mpz_clear(prime_test);

  //printf("Copying primes and instances to the GPU ...\n");
  CUDA_CHECK(cudaSetDevice(0));
  CUDA_CHECK(cudaMalloc((void **)&gpuPrimes, sizeof(uint32_t)*prime_count));
  CUDA_CHECK(cudaMemcpy(gpuPrimes, primes, sizeof(uint32_t)*prime_count, cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMalloc((void **)&gpuInstances, sizeof(instance_t)*instance_count));
  CUDA_CHECK(cudaMemcpy(gpuInstances, instances, sizeof(instance_t)*instance_count, cudaMemcpyHostToDevice));

  // create a cgbn_error_report for CGBN to report back errors
  CUDA_CHECK(cgbn_error_report_alloc(&report));

  //printf("Running GPU kernel ...\n");

  size_t BLOCKS_PER_CALL = 16;

  size_t offset = 0;
  int32_t result = 0;
  while (true) {
      //    total_blocks = ceil(instance_count / IPB) and threads = TPB
      if (offset >= instance_count) {
          // Didn't find result
          result = -1;
          break;
      }

      size_t blocks = min(BLOCKS_PER_CALL, (instance_count - offset + IPB - 1) / IPB);
      kernel_miller_rabin<params><<<blocks, TPB>>>(report, gpuInstances + offset, instance_count - offset, gpuPrimes, prime_count);

      // error report uses managed memory, so we sync the device (or stream) and check for cgbn errors
      CUDA_CHECK(cudaDeviceSynchronize());
      CGBN_CHECK(report);

      // copy the instances back from gpuMemory
      size_t num_to_copy = min(blocks * IPB, (instance_count - offset));
      CUDA_CHECK(cudaMemcpy(instances + offset, gpuInstances + offset, sizeof(instance_t)* num_to_copy, cudaMemcpyDeviceToHost));

      result = miller_rabin_t<params>::verify_first(instances + offset, num_to_copy, primes, prime_count);
      if (result >= 0) {
          result += offset;
          //printf("  found %d | %d => %d\n", sign, result, offsets[result]);
          break;
      }

      offset += num_to_copy;
      //printf("  >%ld (%d)\n", offset, offsets[offset]);
  }

  // clean up
  free(primes);
  free(instances);
  CUDA_CHECK(cudaFree(gpuPrimes));
  CUDA_CHECK(cudaFree(gpuInstances));
  CUDA_CHECK(cgbn_error_report_free(report));

  return result;
}
