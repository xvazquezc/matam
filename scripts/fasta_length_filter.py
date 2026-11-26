#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fasta_length_filter

Description: Filter a fasta file based on sequence length

  fasta_length_filter.py -i input.fa -o output.fa -m 300
  fasta_length_filter.py -i input.fa -o output.fa -M 1000

-----------------------------------------------------------------------

Author: This software is written and maintained by Pierre Pericard
(pierre.pericard@ed.univ-lille1.fr)
Created: 2016-04-12
Last Modified: 2016-04-12
Licence: GNU GPL 3.0

Copyright 2016 Pierre Pericard

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
from fasta_utils import format_seq, read_fasta_file_handle


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Filter a fasta file based on sequence length.')
    parser.add_argument('-i', '--input_fasta', metavar='input',
                        type=argparse.FileType('r'), default='-',
                        help='input fasta file')
    parser.add_argument('-o', '--output_fasta', metavar='output',
                        type=argparse.FileType('w'), default='-',
                        help='ouput fasta file')
    parser.add_argument('-m', '--min_length', metavar='MIN',
                        type=int, default=0,
                        help='Minimum sequence length')
    parser.add_argument('-M', '--max_length', metavar='MAX',
                        type=int, default=0,
                        help='Maximum sequence length')
    args = parser.parse_args()

    # Code is duplicated here to prevent to have to test args.max_length many times
    if args.max_length:
        for header, sequence in read_fasta_file_handle(args.input_fasta):
            if len(sequence) >= args.min_length and len(sequence) <= args.max_length:
                args.output_fasta.write(">{0}\n{1}\n".format(header, format_seq(sequence)))
    else:
        for header, sequence in read_fasta_file_handle(args.input_fasta):
            if len(sequence) >= args.min_length:
                args.output_fasta.write(">{0}\n{1}\n".format(header, format_seq(sequence)))
