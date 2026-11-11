#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
import time
import logging
from binary_utils import Binary
import sortmerna_index
import matam_version

# Create logger
logger = logging.getLogger(__name__)

# Get program filename
program_filename = os.path.basename(sys.argv[0])

# Get MATAM root dir absolute path
index_ref_db_bin = os.path.realpath(sys.argv[0])
matam_script_dir = os.path.dirname(index_ref_db_bin)
matam_root_dir = os.path.dirname(matam_script_dir)

# Get all dependencies bin


class DefaultHelpParser(argparse.ArgumentParser):
    """
    This is a slightly modified argparse parser to display the full help
    on parser error instead of only usage
    """
    def error(self, message):
        sys.stderr.write('\nError: %s\n\n' % message)
        self.print_help()
        sys.exit(2)


def parse_arguments():
    """
    Parse the command line, and check if arguments are correct
    """
    # Initiate argument parser
    parser = DefaultHelpParser(description='Index ref db',
                               # to precisely format help display
                               formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=80))

    # --version
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(matam_version.version()))

    # Main parameters
    group_main = parser.add_argument_group('Main parameters')
    # -i / --input_ref_db
    group_main.add_argument('-i', '--input_ref_db',
                            action = 'store',
                            metavar = 'DBPATH',
                            type = str,
                            help = 'MATAM ref db. '
                                   'Default is $MATAM_DIR/db/SILVA_123_SSURef_rdNs_NR95')
    # -v / --verbose
    group_main.add_argument('-v', '--verbose',
                            action = 'store_true',
                            help = 'Increase verbosity')

    # Performance parameters
    group_perf = parser.add_argument_group('Performance')
    # --cpu
    group_perf.add_argument('--cpu',
                            action = 'store',
                            metavar = 'CPU',
                            type = int,
                            default = 1,
                            help = argparse.SUPPRESS)
    # --max_memory
    group_perf.add_argument('--max_memory',
                            action = 'store',
                            metavar = 'MAXMEM',
                            type = int,
                            default = 10000,
                            help = 'Maximum memory to use (in MBi). '
                                   'Default is %(default)s MBi')

    args = parser.parse_args()

    # Set default values
    if not args.input_ref_db:
        args.input_ref_db = os.path.join(matam_root_dir, 'db', 'SILVA_123_SSURef_rdNs_NR95')

    # Get absolute path for all arguments
    args.input_ref_db = os.path.abspath(args.input_ref_db)

    #
    return args


if __name__ == '__main__':

    # Arguments parsing
    args = parse_arguments()

    # Init error code
    error_code = 0

    # Set logging
    # create console handler
    ch = logging.StreamHandler()
    #
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

    # Init list of tmp files to delete at the end
    to_rm_filepath_list = list()

    ##############################################
    # Set all files and directories names + paths

    ref_db_basepath = args.input_ref_db
    ref_db_dir, ref_db_basename = os.path.split(ref_db_basepath)

    complete_ref_db_basename = ref_db_basename + '.complete'
    complete_ref_db_basepath = os.path.join(ref_db_dir, complete_ref_db_basename)
    complete_ref_db_filename = complete_ref_db_basename + '.fasta'
    complete_ref_db_filepath = os.path.join(ref_db_dir, complete_ref_db_filename)

    complete_ref_db_taxo_filename = complete_ref_db_basename + '.taxo.tab'
    complete_ref_db_taxo_filepath = os.path.join(ref_db_dir, complete_ref_db_taxo_filename)

    clustered_ref_db_basename = ref_db_basename + '.clustered'
    clustered_ref_db_basepath = os.path.join(ref_db_dir, clustered_ref_db_basename)
    clustered_ref_db_filename = clustered_ref_db_basename + '.fasta'
    clustered_ref_db_filepath = os.path.join(ref_db_dir, clustered_ref_db_filename)

    ######################################################
    # SortMeRNA indexing of complete and clustered ref db

    # Same indexes, parameters and records as matam_db_preprocessing.py
    # (see sortmerna_index.py), so matam_assembly.py can check and reuse them.
    sortmerna_bin = Binary.assert_which('sortmerna')
    index_params = sortmerna_index.index_params(args.max_memory)

    def call_sortmerna(cmd_line, verbose):
        logger.debug('CMD: {0}'.format(cmd_line))
        returncode = subprocess.call(cmd_line, shell=True)
        if verbose:
            sys.stdout.write('\n')
        return returncode

    for label, ref_db_filepath, ref_db_basepath in (
            ('complete', complete_ref_db_filepath, complete_ref_db_basepath),
            ('clustered', clustered_ref_db_filepath, clustered_ref_db_basepath)):
        logger.info('Indexing {0} ref db'.format(label))
        smr_dir = ref_db_basepath + '.smr'
        try:
            sortmerna_index.write_params(smr_dir, index_params)
            with sortmerna_index.locked(smr_dir):
                sortmerna_index.build(sortmerna_bin, ref_db_filepath, os.path.join(smr_dir, 'idx'),
                                      index_params, verbose=args.verbose, call=call_sortmerna)
        except (sortmerna_index.IndexBuildError, OSError) as error:
            logger.error('Could not index {0} ref db: {1}'.format(label, error))
            error_code += 1

    ###############
    # Exit program

    exit(error_code)
