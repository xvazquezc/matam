#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import logging
from fasta_utils import read_fasta_file_handle
from sam_utils import parse_cigar
from tab_utils import read_tab_file_handle_sorted


# Create logger
logger = logging.getLogger(__name__)

if __name__ == '__main__':

    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    # -i / --input_sam
    parser.add_argument('-i', '--input_sam',
                        metavar = 'INSAM',
                        type = argparse.FileType('r'),
                        default = '-',
                        help = 'Input sam file, sorted by reference and position')
    # -o / --output_cov
    parser.add_argument('-o', '--output_cov',
                        metavar = 'OUTCOV',
                        type = argparse.FileType('w'),
                        default = '-',
                        help = 'Output ref coverage tab file')
    # -r / --references
    parser.add_argument('-r', '--references',
                        metavar = 'REF',
                        type = argparse.FileType('r'),
                        required = True,
                        help = 'References fasta file')
    # --histo
    parser.add_argument('--histo',
                        metavar = 'FILE',
                        type = str,
                        default = 'ref_coverage_histogram.svg',
                        help = 'References coverage histogram')
    # -v / --verbose
    parser.add_argument('-v', '--verbose',
                        action = 'store_true',
                        help = 'Increase verbosity')
    # --debug
    parser.add_argument('--debug',
                        action = 'store_true',
                        help = 'Debug mode')

    args = parser.parse_args()

    # Set logging
    # create console handler
    ch = logging.StreamHandler()
    #
    if args.debug:
        logger.setLevel(logging.DEBUG)
        # create formatter for debug level
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    else:
        if args.verbose:
            logger.setLevel(logging.INFO)
        else:
            logger.setLevel(logging.WARNING)
        # create default formatter
        formatter = logging.Formatter('%(levelname)s - %(message)s')
    # add the formatter to the console handler
    ch.setFormatter(formatter)
    # add the handler to logger
    logger.addHandler(ch)

    # Get ref seqs and initialize positions count
    logger.info('Reading reference sequences')
    ref_seq_dict = dict()
    for header, seq in read_fasta_file_handle(args.references):
        seqid = header.split()[0]
        ref_seq_dict[seqid] = seq.upper()

    #
    logger.info('Reading through input sam file')
    #
    coverage_list = list()
    frequency_count_dict = defaultdict(int)

    # Reading sam file reference by reference
    for alignment_tabs_list in read_tab_file_handle_sorted(args.input_sam, 2):
        ref_id = alignment_tabs_list[0][2]
        ref_len = len(ref_seq_dict[ref_id])
        total_aligned_nt = 0
        for alignment_tab in alignment_tabs_list:
            cigar = alignment_tab[5]
            for operation, count in parse_cigar(cigar):
                if operation == 'M':
                    total_aligned_nt += count
        # Compute reference coverage
        ref_coverage = total_aligned_nt/ref_len
        #
        frequency_count_dict[ref_coverage] += 1
        coverage_list.append(ref_coverage)

    # Writing ref coverage tab
    logger.info('Writing output ref coverage tab file')
    sorted_ref_coverage_frequency_tuple_tuple = tuple(sorted(frequency_count_dict.items(), key=lambda x: x[0]))
    for frequency_tuple in sorted_ref_coverage_frequency_tuple_tuple:
        args.output_cov.write('{}\t{}\n'.format(*frequency_tuple))

    # Plotting histogram
    logger.info('Plotting histogram')
    plt.hist(coverage_list, range=(0,500), bins=250)
    plt.title("Ref Coverage")
    plt.xlabel("Coverage")
    plt.ylabel("Frequency")
    plt.ylim([0,200])
    plt.savefig(args.histo)
    plt.show()