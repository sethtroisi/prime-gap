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

import array
import contextlib
import logging
import os
import re
import sys


UNKNOWN_FILENAME_RE = re.compile(
    r"^(\d+)_(\d+)_(\d+)_(\d+)_s(\d+)_l(\d+)M(.m1)?(?:.missing)?.txt")


class TeeLogger:
    def __init__(self, fn, out):
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)-15s %(levelname)s: %(message)s",
            handlers=[
                logging.StreamHandler(out),
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
    if not args.save_logs:
        return contextlib.nullcontext()

    assert args.unknown_filename
    unknown_path = transform_unknown_filename(args.unknown_filename, "logs", "log")

    for num in range(0, 5):
        log_fn = unknown_path
        if num > 0:
            log_fn += "." + str(num)

        if not os.path.exists(log_fn):
            print("Saving logs to '{}'".format(log_fn))
            return contextlib.redirect_stdout(TeeLogger(log_fn, sys.stdout))

    assert False, "log file '{}' already exists x5".format(unknown_path)


def transform_unknown_filename(unknown_fn, directory, extension):
    """Return a new path similar to unknown_fn with corrected extension
    and in directory (if exists) otherwise without directory"""
    fn = os.path.splitext(os.path.basename(unknown_fn))[0]
    if not extension.startswith("."):
        extension = "." + extension
    fn += extension

    if os.path.isdir(directory):
        return os.path.join(directory, fn)
    return fn


def parse_unknown_filename(fn):
    """None or (p,d,ms,mi,sl,mp,m1)"""
    match = UNKNOWN_FILENAME_RE.match(os.path.basename(fn))
    return tuple(map(lambda s: int(s) if s and s.isdigit() else s, match.groups())) if match else match


def generate_unknown_filename(p, d, ms, mi, sl, mp, method1=False, ext=".txt"):
    base = f"{p}_{d}_{ms}_{mi}_s{sl}_l{mp}M"
    return base + (".m1" * (method1 is True)) + ext


def verify_args(args):
    if args.unknown_filename:
        fn = args.unknown_filename
        if not (os.path.exists(fn) or os.path.exists("unknowns/" + fn)):
            print("\"{}\" doesn't exist".format(fn))
            exit(1)

        match = parse_unknown_filename(fn)
        if not match:
            print(f"{fn!r} doesn't match unknown file format")
            exit(1)

        p, d, ms, mi, sl, mp, m1 = match
        (args.p, args.d, args.mstart, args.minc,
         args.sieve_length, args.max_prime, _) = match
        args.method1 = (m1 == ".m1")

    args.max_prime *= 10 ** 6

    if 'search_db' in args and args.search_db:
        assert os.path.exists(args.search_db), (
            "Prime Search Database ('{}') doesn't exist".format(args.search_db))

    for arg in ('p', 'd', 'mstart', 'minc', 'sieve_length', 'max_prime'):
        if arg not in args or args.__dict__[arg] in (None, 0):
            print("Missing required argument", arg)
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


def parse_unknown_line(line):
    unknowns = [[], []]

    start, c_l, c_h = line.split(b" | ")

    match = re.match(rb"^([0-9]+) : -([0-9]+) \+([0-9]+)", start)
    assert match, start
    m_test, unknown_l, unknown_h = map(int, match.groups())

    # Check if rle or raw
    rle = b" " not in c_l[:20]
    if rle:
        def accum_rle(sign, digits):
            assert len(digits) % 2 == 0 or digits[-1] == ord("\n")
            # Read digits(bytes) in pairs (see save_unknowns_method2)
            values = array.array('l')
            accum = 0
            for i in range(0, len(digits)//2):
                delta = 128 * (digits[2*i] - 48) + (digits[2*i+1] - 48)
                accum += delta
                values.append(sign * accum)
            return values

        assert len(c_l) == 2 * unknown_l
        assert len(c_h) == 2 * unknown_h + (c_h[-1] == ord("\n"))
        unknowns[0] = accum_rle(-1, c_l)
        unknowns[1] = accum_rle(+1, c_h)
    else:
        unknowns[0] = array.array('l', map(int, c_l.split(b" ")))
        unknowns[1] = array.array('l', map(int, c_h.split(b" ")))

    unknown_l_test = len(unknowns[0])
    unknown_h_test = len(unknowns[1])
    assert unknown_l == unknown_l_test, (unknown_l, unknown_l_test, "\t", start)
    assert unknown_h == unknown_h_test, (unknown_h, unknown_h_test, "\t", start)

    return m_test, unknown_l, unknown_h, unknowns


def convert_to_rle_line(parts):
    """Convert (start, c_l, c_h, (unknowns_l, unknowns_h)) to rle encoded line"""

    assert len(parts) == 4, parts

    line_bytes = array.array('B')
    line_bytes.extend("{} : -{} +{}".format(*parts[:3]).encode())

    for side in parts[3]:
        line_bytes.extend(b" | ")

        last = 0
        for x in side:
            x = abs(x)
            delta = x - last
            last = x
            assert 0 <= delta < (128 * 128)
            upper = 48 + (delta // 128)
            lower = 48 + (delta % 128)
            line_bytes.append(upper)
            line_bytes.append(lower)

    return line_bytes.tobytes()


