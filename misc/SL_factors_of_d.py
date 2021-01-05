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

import gmpy2
import sympy
import matplotlib.pyplot as plt
import numpy as np

P = 907

SL = 2 * 10 ** 7
for d in [2, 3, 5, 6, 30]:
        K = gmpy2.primorial(P) // d
        # Ignore 1 which 'messes' up percent factor graph
        test = [i for i in range(2, SL+1) if gmpy2.gcd(i, K) == 1]

        print (d, len(test), test[:20])

        count0 = [ gmpy2.gcd(t, d) == 1 for t in test]
        cum_count = np.cumsum(count0)


        if True:
            # Show percent of X divisible by d
            percent_factor = cum_count / np.arange(1, len(test)+1)
        else:
            # Show local percent of X divisible by d
            group_size = 100
            cum_count = np.array([cum_count[i] - cum_count[max(0, i - group_size)]
                for i in range(len(cum_count))])
            percent_factor = cum_count / group_size


        #scatter = plt.scatter(test, percent_factor, label=str(d), marker='.')
        #c = scatter.get_fc()[0]

        factors_d = sympy.factorint(d)
        assert all(v == 1 for v in factors_d.values())
        factors_d = sorted(factors_d.keys())

        expected = sympy.prod([1 - 1/p for p in factors_d])

        line = plt.plot(test, percent_factor, label=str(d)) #, marker='.')
        c = line[0].get_c()

        plt.axhline(y=expected, color=c)

        if gmpy2.is_prime(d):
            plt.axvline(x=d * P, color=c)

        plt.ylabel("Fraction of X with gcd(X, d) = 1")

        count_after = len([i for i in test if gmpy2.gcd(K + i, d) == 1])
        print (f"{d}\t{count_after}/{len(test)} = {count_after/len(test):.4f} vs {expected:.4f}\t {2*count_after}")


plt.title(f'K={P}#/d')
plt.xlabel("X with (X, K) == 1")

plt.axvline(x=P, color='black')

plt.legend(loc='upper right')
plt.xscale('log')
plt.show()
