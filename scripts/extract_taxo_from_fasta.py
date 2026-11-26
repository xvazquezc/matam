#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import re
from fasta_utils import read_fasta_file_handle

if __name__ == '__main__':

    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input_fasta',
                        metavar='TAB',
                        type=argparse.FileType('r'),
                        default='-',
                        help='input fasta file')
    parser.add_argument('-o', '--output_tab',
                        metavar='TAB',
                        type=argparse.FileType('w'),
                        default='-',
                        help='output tab file')
    args = parser.parse_args()

    for header, sequence in read_fasta_file_handle(args.input_fasta):
        header_tab = header.split()
        ref_id = header_tab[0]
        taxo = 'NULL'
        if len(header_tab) > 1:
            taxo = '_'.join(header_tab[1:])
        args.output_tab.write('{0}\t{1}\n'.format(ref_id, taxo))
