import os
import re

import gmpy2
import subprocess


UNKNOWN_FILENAME_RE = re.compile(
    "^(\d+)_(\d+)_(\d+)_(\d+)_s(\d+)_l(\d+)M(.m2)?(?:.missing)?.txt")


def verify_args(args, fn_extension):
    if args.unknown_filename:
        fn = args.unknown_filename
        if not os.path.exists(fn):
            print ("\"{}\" doesn't exist".format(fn))
            exit(1)
        match = UNKNOWN_FILENAME_RE.match(os.path.basename(fn))
        if not match:
            print ("\"{}\" doesn't match unknown file format".format(fn))
            exit(1)

        ms, p, d, mi, sl, sr, m2 = match.groups()
        args.mstart = int(ms)
        args.minc   = int(mi)
        args.p      = int(p)
        args.d      = int(d)
        args.sieve_length = int(sl)
        args.sieve_range  = int(sr)
        args.method2 = (m2 == ".m2")

    args.sieve_range *= 10 ** 6

    # TODO Should this be moved out?
    #if not os.path.exists(args.prime_db):
    #    print ("prime-db \"{}\" doesn't exist".format(args.prime_db))
    #    exit(1)

    for arg in ('mstart', 'minc', 'p', 'd', 'sieve_length', 'sieve_range'):
        if arg not in args or args.__dict__[arg] in (None, 0):
            print ("Missing required argument", arg)
            exit(1)

    # TODO handle .M2 and .missing
    fn = "{}_{}_{}_{}_s{}_l{}M{}".format(
        args.mstart, args.p, args.d, args.minc,
        args.sieve_length, args.sieve_range // 10 ** 6,
        ".m2" if args.method2 else "")
    fn += fn_extension

    if args.unknown_filename:
        basename = os.path.basename(args.unknown_filename)
        assert basename.startswith(fn), (fn, args.unknown_filename)
    else:
        args.unknown_filename = fn


def openPFGW_is_prime(strn):
    # Overhead of subprocess calls seems to be ~0.03
    s = subprocess.getstatusoutput("./pfgw64 -f0 -q" + strn)
    assert s[1].startswith('PFGW'), s
    return s[0] == 0


def is_prime(num, strnum, dist):
    # TODO print log of which library is being used.
    #if gmpy2.num_digits(num, 2) > 5000:
    #    return openPFGW_is_prime(strnum + str(dist))

    return gmpy2.is_prime(num)

