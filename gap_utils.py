#!/usr/bin/env python3
#
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
        return contextlib.suppress()

    assert args.unknown_filename

    if os.path.isdir("logs"):
        transform_unknown_filename(unknown_fn, "logs", "log")

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


def transform_unknown_filename(unknown_fn, directory, extension):
    '''Return a new path similiar to unknown_fn with corrected extension
    and in direcotry (if exists) otherwise without directory'''
    fn = os.path.splitext(os.path.basename(unknown_fn))[0]
    if not extension.startswith("."):
        extension = "." + extension
    fn += extension

    if os.path.isdir(directory):
        return os.path.join(directory, fn)
    return fn


def parse_unknown_filename(fn):
    '''None or (p,d,ms,mi,sl,mp,m1)'''
    match = UNKNOWN_FILENAME_RE.match(os.path.basename(fn))
    return tuple(map(lambda s: int(s) if s and s.isdigit() else s, match.groups())) if match else match


def generate_unknown_filename(p, d, ms, mi, sl, mp, method1 = False, ext = ".txt"):
    base = f"{p}_{d}_{ms}_{mi}_s{sl}_l{mp}M"
    return base + (".m1" * (method1 == True)) + ext


def verify_args(args):
    if args.unknown_filename:
        fn = args.unknown_filename
        if not (os.path.exists(fn) or os.path.exists("unknowns/" + fn)):
            print ("\"{}\" doesn't exist".format(fn))
            exit(1)

        match = parse_unknown_filename(fn)
        if not match:
            print (f"{fn!r} doesn't match unknown file format")
            exit(1)

        p, d, ms, mi, sl, mp, m1 = match
        args.p              = p
        args.d              = d
        args.mstart         = ms
        args.minc           = mi
        args.sieve_length   = sl
        args.max_prime      = mp
        args.method1        = (m1 == ".m1")

    args.max_prime   *= 10 ** 6

    if 'search_db' in args and args.search_db:
        assert os.path.exists(args.search_db), (
            "Prime Search Database ('{}') doesn't exist".format(args.search_db))

    for arg in ('p', 'd', 'mstart', 'minc', 'sieve_length', 'max_prime'):
        if arg not in args or args.__dict__[arg] in (None, 0):
            print ("Missing required argument", arg)
            exit(1)

    fn = generate_unknown_filename(
        args.p, args.d,
        args.mstart, args.minc,
        args.sieve_length, args.max_prime // 10 ** 6,
        args.method1)

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


def parse_unknown_line(line):
    unknowns = [[], []]

    start, c_l, c_h = line.split(b" | ")

    match = re.match(b"^([0-9]+) : -([0-9]+) \+([0-9]+)", start)
    assert match, start
    mtest, unknown_l, unknown_u = map(int, match.groups())

    # Check if rle or raw
    rle = b" " not in c_l[:20]
    if rle:
        def unrle(sign, digits):
            # Read digits(bytes) in pairs (see save_unknowns_method2)
            values = []
            accum = 0
            for i in range(0, len(digits)//2):
                delta = 128 * digits[2*i] + digits[2*i+1]
                accum += delta
                values.append(sign * accum)
            return values

        unknowns[0] = unrle(-1, c_l)
        unknowns[1] = unrle(+1, c_h[:-1])
    else:
        unknowns[0] = list(map(int,c_l.split(b" ")))
        unknowns[1] = list(map(int,c_h.split(b" ")))

    unknown_l_test = len(unknowns[0])
    unknown_u_test = len(unknowns[1])
    assert unknown_l == unknown_l_test, (unknown_l, unknown_l_test, "\t", start)
    assert unknown_u == unknown_u_test, (unknown_u, unknown_u_test, "\t", start)

    return mtest, unknown_l, unknown_u, unknowns


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
