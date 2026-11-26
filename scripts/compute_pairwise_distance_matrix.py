#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""

import argparse
import sys
from fasta_utils import read_fasta_file_handle


def compute_distance(seq_i, seq_j, semi_global_mode):
    """
    Compute the distance between 2 sequences
    """
    distance = -1
    
    match = 0
    mismatch = 0
    indel = 0
    indel_tmp = 0
    
    for i in range(len(seq_i)):
        a = seq_i[i]
        b = seq_j[i]
        
        if a != '-' or b != '-':
            if a != '-' and b != '-':
                if a == b:
                    match += 1
                else:
                    mismatch += 1
                if indel_tmp or not semi_global_mode: # We just got out of a serie of indels
                    indel += indel_tmp
                    indel_tmp = 0
            else:
                if match or mismatch or not semi_global_mode: # We are in the alignment
                    indel_tmp += 1
    
    if not semi_global_mode and indel_tmp:
        indel += indel_tmp
    
    length = match + mismatch + indel
    if length:
        distance = float(match) / float(length)
    
    return distance
    


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Compute all pairwise distances.')
    parser.add_argument('-i', '--input_fasta', 
                        action = 'store',
                        metavar = 'INFA',
                        type = argparse.FileType('r'), 
                        default = '-',
                        help = 'input fasta file')
    parser.add_argument('-o', '--output_tab', 
                        action = 'store',
                        metavar = 'OUTTAB',
                        type = argparse.FileType('w'), 
                        default = '-',
                        help = 'input tab file')
    parser.add_argument('--semi_global', 
                        action='store_true',
                        help='Compute distance only in the common aligned region (dont consider clipped sequences)')
    args = parser.parse_args()

    seq_list = list()

    for header, sequence in read_fasta_file_handle(args.input_fasta):
        seq_list.append((header, sequence))
    
    for i in range(len(seq_list)-1):
        for j in range(i+1, len(seq_list)):
            seq_i_id = seq_list[i][0]
            seq_j_id = seq_list[j][0]
            seq_i = seq_list[i][1]
            seq_j = seq_list[j][1]
            
            if len(seq_i) != len(seq_j):
                sys.stderr.write('WARNING: {} and {} doesnt have the same length\n'.format(seq_i_id, seq_j_id))
            
            args.output_tab.write('{0}\t{1}\t'.format(seq_i_id, seq_j_id))
            distance = compute_distance(seq_i, seq_j, args.semi_global)
            args.output_tab.write('{:.2f}\n'.format(distance*100.0))












