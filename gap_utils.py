import contextlib
import logging
import os
import re
import subprocess
import sys

import gmpy2


UNKNOWN_FILENAME_RE = re.compile(
    "^(\d+)_(\d+)_(\d+)_(\d+)_s(\d+)_l(\d+)M(.m1)?(?:.missing)?.txt")



class TeeLogger:
    def __init__(self, fn, sysout):
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)-15s %(levelname)s: %(message)s",
            handlers=[
                logging.StreamHandler(sys.stdout),
                logging.FileHandler(fn, mode="w"),
            ]
        )
        self.logger = logging.getLogger()

    def write(self, msg):
        if msg and not msg.isspace():
            self.logger.info(msg)

    def flush(self):
        self.logger.flush()


def logger_context(args):
    # contextlib.nullcontext() requires python 3.7
    if not args.save_logs:
        context = contextlib.suppress()

    assert args.unknown_filename
    log_fn_base = args.unknown_filename + '.log'
    for num in range(0, 5):
        log_fn = log_fn_base
        if num > 0:
            log_fn += "." + str(num)

        if not os.path.exists(log_fn):
            print("Saving logs to '{}'".format(log_fn))
            return contextlib.redirect_stdout(TeeLogger(log_fn, sys.stdout))

    assert False, "log file '{}' already exists x5".format(log_fn_base)
    return context


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
        args.method1 = (m2 == ".m1")

    args.sieve_range *= 10 ** 6

    if 'search_db' in args and args.search_db:
        assert os.path.exists(args.search_db), (
            "Prime Search Database ('{}') doesn't exist".format(args.search_db))

    for arg in ('mstart', 'minc', 'p', 'd', 'sieve_length', 'sieve_range'):
        if arg not in args or args.__dict__[arg] in (None, 0):
            print ("Missing required argument", arg)
            exit(1)

    fn = "{}_{}_{}_{}_s{}_l{}M{}".format(
        args.mstart, args.p, args.d, args.minc,
        args.sieve_length, args.sieve_range // 10 ** 6,
        ".m1" if args.method2 else "")
    fn += fn_extension

    if args.unknown_filename:
        basename = os.path.basename(args.unknown_filename)
        assert basename.startswith(fn), (fn, args.unknown_filename)
    else:
        args.unknown_filename = fn


def K_and_stats(args):
    P = args.p
    D = args.d

    assert P <= 80000
    K = gmpy2.primorial(P)
    K, r = divmod(K, D)
    assert r == 0

    K_digits = gmpy2.num_digits(K, 10)
    K_bits   = gmpy2.num_digits(K, 2)
    K_log    = float(gmpy2.log(K))

    return K, K_digits, K_bits, K_log


def openPFGW_is_prime(strn):
    # Overhead of subprocess calls seems to be ~0.03
    # Process seems to use more than 1 thread, accounting for this gmp is quite competitive.
    s = subprocess.getstatusoutput("./pfgw64 -f0 -q" + strn)
    assert s[1].startswith('PFGW'), s
    return s[0] == 0


def is_prime(num, strnum, dist):
    # TODO print log of which library is being used.
    #if gmpy2.num_digits(num, 2) > 5000:
    #    return openPFGW_is_prime(strnum + str(dist))

    return gmpy2.is_prime(num)

