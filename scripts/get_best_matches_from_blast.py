#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import random
from tab_utils import read_tab_file_handle_sorted


if __name__ == '__main__':

    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    # -i / --input_tab
    parser.add_argument('-i', '--input_tab',
                        metavar='INBLAST',
                        type=argparse.FileType('r'),
                        default='-',
                        help='Input blast tab file. '
                             'Assuming sorted by subject, '
                             'eval, score, and qcov')
    # -o / --output_tab
    parser.add_argument('-o', '--output_tab',
                        metavar='OUTBLAST',
                        type=argparse.FileType('w'),
                        default='-',
                        help='Ouput filtered blast tab file')
    # -p / --threshold_percent
    parser.add_argument('-p', '--threshold_percent',
                        metavar='FLOAT',
                        type=float,
                        default=1.0,
                        help='Score threshold percent')
    # --random_best
    parser.add_argument('--random_best',
                        action='store_true',
                        help='Select one random match among all '
                             'best matches for each subject')

    args = parser.parse_args()

    # Initialize random seed
    random.seed()

    # Reading blast tab file
    for tab_list in read_tab_file_handle_sorted(args.input_tab, 0):
        # Get query id and best blast score
        query_id = tab_list[0][0]
        best_score = int(tab_list[0][11])
        # Keep only best matches
        best_lines_list = list()
        for tab in tab_list:
            blast_score = int(tab[11])
            if blast_score >= best_score * args.threshold_percent:
                best_lines_list.append('\t'.join(tab))
        # Print best tabs
        if args.random_best:
            args.output_tab.write('{0}\n'.format(random.choice(best_lines_list)))
        else:
            for line in best_lines_list:
                args.output_tab.write('{0}\n'.format(line))
