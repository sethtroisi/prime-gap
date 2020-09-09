import gmpy2
import re
import time

failed = '''
	BOTH SIDES PRIME: 65309*8009#/1110 2160 130058  => gap = ???
	BOTH SIDES PRIME: 75061*8009#/1110 54212 71922  => gap = 28278
	BOTH SIDES PRIME: 78889*8009#/1110 67638 62614  => gap = 10852
	BOTH SIDES PRIME: 84647*8009#/1110 23398 105474 => gap = 58032
	BOTH SIDES PRIME: 38083*8009#/1110 12800 115314 => gap = 72638
	BOTH SIDES PRIME: 38821*8009#/1110 56418 73366  => gap = 56410
	BOTH SIDES PRIME: 82211*9001#/1110 61062 63962  => gap = 22316
	BOTH SIDES PRIME: 38407*9001#/1110 54452 71754  => gap = 84422
	BOTH SIDES PRIME: 37859*9001#/1110 11988 118292 => gap = 17108
	BOTH SIDES PRIME: 1207*10007#/1110 58508 73734  => gap = 71630
	BOTH SIDES PRIME: 1889*8009#/1110 162 131402    => gap = 73680
	BOTH SIDES PRIME: 1133*9001#/1110 60 126146     => gap = 55302
    BOTH SIDES PRIME: 23083*8009#/1110 25934 104964 => gap = 28598
	BOTH SIDES PRIME: 1027*9001#/1110 30720 93538   => gap = 70558
	BOTH SIDES PRIME: 8497*9001#/1110 7104 124834   => gap = 61390
	BOTH SIDES PRIME: 3497*10007#/1110 120 115574   => gap = 27482
	BOTH SIDES PRIME: 44321*8009#/1110 18750 105596 => gap = 54428
	BOTH SIDES PRIME: 16477*8009#/1110 40844 91362  => gap = 798
    BOTH SIDES PRIME: 7063*13001#/1110 19980 111974 => gap 85406
	BOTH SIDES PRIME: 10163*8009#/1110 8880 123314  => gap 8904
	BOTH SIDES PRIME: 3659*10007#/1110 96 117146    => gap 31394
'''

success = '''
    BOTH SIDES PRIME: 64003*12007#/1110 47954 82944 => gap 130898
'''

checks = '''


'''

for test in checks.strip().split("\n"):
    m, p, d, l, h = map(int, re.search(r"(\d+).(\d+)..(\d+) (\d+) (\d+)", test).groups())
    N = m * gmpy2.primorial(p) // d
    low = N - l
    high = N + h

    print (f"Testing {m}*{p}#/{d} (-{l}, +{h})")
    t0 = time.time()
    assert gmpy2.is_prime(low)
    assert gmpy2.is_prime(high)
    t1 = time.time()

    print ("\tverified endpoints {:.2f} seconds".format(t1-t0))

    z = gmpy2.next_prime(low)
    t2 = time.time()

    print ("\t next_prime {}, {}   {:.1f} seconds".format(
        z == high, z - low, t2 - t1))
