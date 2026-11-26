#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""

import sys
import os
import argparse
import string
import re
from fasta_utils import format_seq, read_fasta_file_handle

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Filter a fasta file based on sequence name.')
    parser.add_argument('-i', '--input_fasta', metavar='input',
                        type=argparse.FileType('r'), required=True,
                        help='input fasta file')
    parser.add_argument('-o', '--output_fasta', metavar='output',
                        type=argparse.FileType('w'), required=True,
                        help='ouput fasta file')
    args = parser.parse_args()

    for header, sequence in read_fasta_file_handle(args.input_fasta):
        sequence_id = header.split()[0]
        args.output_fasta.write(">{0}\n{1}\n".format(sequence_id, format_seq(sequence)))

