#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""

import sys
import os
import argparse
import string
import re
import random
from fasta_utils import format_seq, read_fasta_file_handle

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Replace all the Ns by random nucleotides.')
    parser.add_argument('-i', '--input_fasta', metavar='input',
                        type=argparse.FileType('r'), default='-',
                        help='input fasta file')
    parser.add_argument('-o', '--output_fasta', metavar='output',
                        type=argparse.FileType('w'), default='-',
                        help='ouput fasta file')
    args = parser.parse_args()

    ACGT = ['A', 'C', 'G', 'T']

    for header, sequence in read_fasta_file_handle(args.input_fasta):
        sequence = re.sub(r'[^ACGT]', lambda x: random.choice(ACGT), sequence.upper())
        args.output_fasta.write(">{0}\n{1}\n".format(header, format_seq(sequence)))


