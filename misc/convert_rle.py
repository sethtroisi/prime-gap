#!/usr/bin/env python3
#
# Copyright 2021 Seth Troisi
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

"""Convert a non RLE file to RLE encoding"""

import argparse
import math
import os.path
import platform
import re
import sys

import tqdm

import sys
# Hack to allow import of gap_utils
sys.path.append(".")

import gap_utils
import misc_utils

def get_arg_parser():
    parser = argparse.ArgumentParser('convert_rle')

    parser.add_argument(
        '-u', '--unknown-filename', type=str, required=True,
        help="determine p, d, mstart, minc, sieve-length, and max-prime"
             " from unknown-results filename")

    parser.add_argument(
        '--bitcompress', action="store_true",
        help="save to bitcompress encoding")

    return parser

def convert_plaintext_to_encoding(args):
    fn_with_unk = gap_utils.transform_unknown_filename(
        args.unknown_filename, "unknowns", "txt")
    if os.path.exists(fn_with_unk):
        args.unknown_filename = fn_with_unk

    converted_fn = args.unknown_filename + "2"

    print(f"\tLoading unknowns from {args.unknown_filename!r}")
    print(f"\tSaving  unknowns to   {converted_fn!r}")

    P, D, ms, mi, sl, mp, m1 = gap_utils.parse_unknown_filename(args.unknown_filename)
    num_m = misc_utils.count_num_m(ms, mi, D)

    if args.bitcompress:
        K = 1
        for p in range(2, P + 1):
            if all(p % j != 0 for j in range(2, int(p ** 0.5 + 0.01) + 1)):
                K *= p
        # assert K == gmpy2.primorial(P)
        K //= D
        coprime_X = [x for x in range(-sl, sl+1) if math.gcd(K, x) == 1]
        converter = lambda parts: gap_utils.convert_to_bitcompressed_line(K, D, coprime_X, parts)
    else:
        converter = gap_utils.convert_to_rle_line

    with open(args.unknown_filename, "rb") as read_file, open(converted_fn, "wb") as write_file:
        first_line = read_file.readline().strip()
        # go back to start of file
        read_file.seek(0, 0)

        converted_first = converter(gap_utils.parse_unknown_line(first_line))
        if converted_first == first_line:
            print(converted_first)
            print("\nFile is already RLE!\n")
            exit(1)

        print(first_line[:80])
        print(converted_first[:80])

        for line in tqdm.tqdm(read_file, total=num_m):
            parts = gap_utils.parse_unknown_line(line)
            write_file.write(converter(parts))
            write_file.write(b"\n")


if __name__ == "__main__":
    if platform.python_implementation() != 'PyPy':
        print("Much faster with PyPy (alternate python implementation in linux)")

    parser = get_arg_parser()
    args = parser.parse_args()

    convert_plaintext_to_encoding(args)
