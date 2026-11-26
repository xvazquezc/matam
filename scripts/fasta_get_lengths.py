#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""

import argparse
from fasta_utils import read_fasta_file_handle


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Get the length of every sequences from a fasta file.')
    parser.add_argument('-i', '--input_fasta',
                        metavar='INFA',
                        type=argparse.FileType('r'),
                        default='-',
                        help='input fasta file')
    parser.add_argument('-o', '--output_tab',
                        metavar='OUTTAB',
                        type=argparse.FileType('w'),
                        default='-',
                        help='ouput tab file (seqid TAB length)')
    args = parser.parse_args()

    for header, sequence in read_fasta_file_handle(args.input_fasta):
        seq_id = header.split()[0]
        args.output_tab.write("{0}\t{1}\n".format(seq_id, len(sequence)))
