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
    parser.add_argument('-n', '--max_consec_n',
                        metavar = 'MAXN',
                        type = int,
                        default = 5,
                        help = 'Maximum number of consecutive Ns allowed in the sequence. '
                               'Default is %(default)s')
    args = parser.parse_args()

    rejected_seq_num = 0

    for header, sequence in read_fasta_file_handle(args.input_fasta):
        n_matches = re.findall(r'[^ACGT]+', sequence.upper())
        max_length_n = 0
        if n_matches:
            max_length_n = max((len(m) for m in n_matches))
        if max_length_n <= args.max_consec_n:
            sequence = re.sub(r'[^ACGT]', 'A', sequence.upper())
            args.output_fasta.write(">{0}\n{1}\n".format(header, format_seq(sequence)))
        else:
            rejected_seq_num += 1
            #~ print('\nSeq is rejected ({} max consecutive Ns)'.format(max_length_n))
            #~ print(">{0}\n{1}\n".format(header, format_seq(sequence)))
    
    sys.stderr.write('{} sequences were rejected\n'.format(rejected_seq_num))
