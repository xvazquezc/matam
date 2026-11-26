#!/usr/bin/env python3

import sys
import re
import random
import argparse
import numpy as np
from fasta_utils import read_fasta_file_handle
from sam_utils import parse_cigar
from tab_utils import group_rows, read_tab_file_handle_sorted


def get_alignment_length_on_ref(cigar):
    """
    """
    return sum((count for operation, count in parse_cigar(cigar) if operation in ('M','D','N','=','X')))


def compute_ref_coverage(ref_sam_tab_list, ref_length):
    """
    """
    # Init ref coverage list
    ref_coverage_list = np.zeros(ref_length, dtype=int)
    # Read alignments on the ref
    for alignment_tab in ref_sam_tab_list:
        starting_pos = int(alignment_tab[3]) - 1 #SAM positions are 1-based
        cigar = alignment_tab[5]
        alignment_length_on_ref = get_alignment_length_on_ref(cigar)
        end_pos = starting_pos + alignment_length_on_ref - 1 #We want the last mapped nucleotide, not the first free one
        ref_coverage_list[starting_pos:end_pos+1] += 1
    #
    return ref_coverage_list


def sample_by_depth(sam_handler, fasta_ref_handler, threshold, sampled_out_sam_handler):
    """
    sample the SAM file based on the coverage at each pos
    """
    # Read reference sequences and store them in a dict
    ref_seq_dict = dict()
    for header, seq in read_fasta_file_handle(fasta_ref_handler):
        seqid = header.split()[0]
        ref_seq_dict[seqid] = seq.upper()
    # Read SAM by ref
    # OPTI: instead of storing all SAM records for each ref in mem, we could read twice the file block corresponding to this ref
    for ref_sam_tab_list in read_tab_file_handle_sorted((l for l in sam_handler if l[0]!='@'), 2):
        ref_id = ref_sam_tab_list[0][2]
        ref_length = len(ref_seq_dict[ref_id])
        # Compute ref coverage list
        ref_coverage_list = compute_ref_coverage(ref_sam_tab_list, ref_length)
        #
        for pos_sam_tab_list in group_rows(ref_sam_tab_list, 3):
            starting_pos = int(pos_sam_tab_list[0][3]) - 1 #SAM positions are 1-based
            ref_coverage = ref_coverage_list[starting_pos]
            if ref_coverage <= threshold:
                for alignment_tab in pos_sam_tab_list:
                    print('{}'.format('\t'.join(alignment_tab)), file=sampled_out_sam_handler)
            else:
                num_alignments_to_remove = ref_coverage - threshold
                random_alignments_order_indices = [i for i in range(len(pos_sam_tab_list))]
                random.shuffle(random_alignments_order_indices)
                #
                for i in random_alignments_order_indices:
                    alignment_tab = pos_sam_tab_list[i]
                    if num_alignments_to_remove > 0:
                        starting_pos = int(alignment_tab[3]) - 1 #SAM positions are 1-based
                        cigar = alignment_tab[5]
                        alignment_length_on_ref = get_alignment_length_on_ref(cigar)
                        end_pos = starting_pos + alignment_length_on_ref - 1 #We want the last mapped nucleotide, not the first free one
                        if min(ref_coverage_list[starting_pos:end_pos+1]) <= threshold:
                            print('{}'.format('\t'.join(alignment_tab)), file=sampled_out_sam_handler)
                        else:
                            ref_coverage_list[starting_pos:end_pos+1] -= 1
                            num_alignments_to_remove -= 1
                    else:
                        print('{}'.format('\t'.join(alignment_tab)), file=sampled_out_sam_handler)


if __name__ == '__main__':
    # Arguments parsing
    parser = argparse.ArgumentParser(description='')
    # -i / --input_sam
    parser.add_argument('-i', '--input_sam',
                        metavar='INSAM',
                        type=argparse.FileType('r'),
                        default='-',
                        help='Input sam file, sorted by reference and position')
    # -o / --output_sam
    parser.add_argument('-o', '--output_sam',
                        metavar='OUTSAM',
                        type=argparse.FileType('w'),
                        default='-',
                        help='Output filtered sam file')
    # -r / --references
    parser.add_argument('-r', '--references',
                        metavar='REF',
                        type=argparse.FileType('r'),
                        required=True,
                        help='References fasta file')
    # -c / --cov_threshold
    parser.add_argument('-c', '--cov_threshold',
                        metavar='COV',
                        type=int,
                        default=50,
                        help='Identity threshold. '
                             'Default is %(default)s')
    # --seed
    parser.add_argument('--seed',
                        metavar='INT',
                        type=int,
                        help='Seed for the random choice of alignments to remove '
                             '(default: system randomness)')
    #
    args = parser.parse_args()
    #
    random.seed(args.seed)
    sample_by_depth(args.input_sam, args.references, args.cov_threshold, args.output_sam)