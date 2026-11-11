#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import collections
import shlex
import re
import subprocess
import time
import logging
from binary_utils import Binary
import sortmerna_index
import discarded_refs_report

# Create logger
logger = logging.getLogger(__name__)

# Get program filename
program_filename = os.path.basename(sys.argv[0])

# Get MATAM root dir absolute path
matam_db_prepro_bin = os.path.realpath(sys.argv[0])
matam_bin_dir = os.path.dirname(matam_db_prepro_bin)
matam_root_dir = os.path.dirname(matam_bin_dir)

# Get all dependencies bin
matam_script_dir = os.path.join(matam_root_dir, 'scripts')
extract_taxo_bin = os.path.join(matam_script_dir, 'extract_taxo_from_fasta.py')
replace_Ns_bin = os.path.join(matam_script_dir, 'replace_Ns_by_As.py')

DEFAULT_CLUSTERING_ID_THRESHOLD = 0.95
DEFAULT_DOMAINS = ('archaea', 'bacteria', 'eukaryota')

# Define a null file handle
FNULL = open(os.devnull, 'w')


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
    parser = DefaultHelpParser(description='MATAM db preprocessing',
                               # to precisely format help display
                               formatter_class=lambda prog: argparse.HelpFormatter(prog, width=120, max_help_position=80))

    # Main parameters
    group_main = parser.add_argument_group('Main parameters')
    # -i / --input_ref
    group_main.add_argument('-i', '--input_ref',
                            action = 'store',
                            metavar = 'INREF',
                            type = str,
                            required = True,
                            help = 'Input reference file (fasta format). '
                                   'Silva-formated taxonomies will be used if available')
    # -d / --db_dir
    group_main.add_argument('-d', '--db_dir',
                            action = 'store',
                            metavar = 'DBDIR',
                            type = str,
                            default = os.getcwd(),
                            help = 'Database output directory. '
                                   'Default is cwd')
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
                            help = 'Max number of CPU to use. '
                                   'Default is %(default)s cpu')
    # --max_memory
    group_perf.add_argument('--max_memory',
                            action = 'store',
                            metavar = 'MAXMEM',
                            type = int,
                            default = 10000,
                            help = 'Maximum memory to use (in MBi). '
                                   'Default is %(default)s MBi')
    # --skip_index
    group_perf.add_argument('--skip_index',
                            action = 'store_true',
                            help = 'Do not build the SortMeRNA indexes of the db. '
                                   'matam_assembly.py then builds them on first use: '
                                   'in the db directory if it is writable (once, '
                                   'shared by later runs), otherwise in each run')

    # Advanced parameters
    group_adv = parser.add_argument_group('Advanced parameters')
    # -m / --min_length
    group_adv.add_argument('-m', '--min_length',
                           action = 'store',
                           metavar = 'MINLGTH',
                           type = int,
                           default = None,
                           help = 'Ref sequences minimum length')
    # -M / --max_length
    group_adv.add_argument('-M', '--max_length',
                           action = 'store',
                           metavar = 'MAXLGTH',
                           type = int,
                           default = None,
                           help = 'Ref sequences maximum length')
    # -n / --max_consecutive_n
    group_adv.add_argument('-n', '--max_consecutive_n',
                           action = 'store',
                           metavar = 'MAXN',
                           type = int,
                           default = 5,
                           help = 'Maximum nb of consecutive Ns a sequence is allowed to have. '
                                  'Default is %(default)s. Setting it to 0 will remove all '
                                  'sequences with Ns. Ns in accepted sequences will be replaced '
                                  'by As')
    # --clustering_id_threshold
    group_adv.add_argument('--clustering_id_threshold',
                           action = 'store',
                           metavar = 'REAL',
                           type = float,
                           default = None,
                           help = 'Identity threshold for clustering. '
                                  'Default is {0}'.format(DEFAULT_CLUSTERING_ID_THRESHOLD))
    # --close_relative_id
    group_adv.add_argument('--close_relative_id',
                           action = 'store',
                           metavar = 'REAL',
                           type = float,
                           default = None,
                           help = 'Identity for a discarded reference sequence (removed by '
                                  'cleaning, or a --clustered_ref centroid not in the cleaned '
                                  'db) to have a close relative left among the remaining '
                                  'references; those without one are reported. 0 skips the '
                                  'search. Default is the clustering identity '
                                  '(--clustering_id_threshold, {0} with --clustered_ref)'.format(
                                      DEFAULT_CLUSTERING_ID_THRESHOLD))
    # --clustered_ref
    group_adv.add_argument('--clustered_ref',
                           action = 'store',
                           metavar = 'FASTA',
                           type = str,
                           help = 'Reference sequences already clustered at the identity '
                                  'of your choice (e.g. the clustered db of an earlier '
                                  'preprocessing) to use as the clustered db instead of '
                                  'clustering INREF. Only their IDs are used: the sequences are '
                                  'taken from the cleaned INREF, so the IDs must be INREF IDs. '
                                  'SILVA NR99 itself is not a clustered set in this sense: it is '
                                  'only made non-redundant at 99%% identity')
    # --by_domain
    group_adv.add_argument('--by_domain',
                           action = 'store_true',
                           help = 'Perform clustering separately for each domain (first rank '
                                  'of the SILVA-formatted taxonomy)')
    # --domains
    group_adv.add_argument('--domains',
                           action = 'store',
                           metavar = 'STR',
                           type = str,
                           nargs = '+',
                           default = None,
                           help = 'Domains to cluster the DB for with --by_domain; sequences of '
                                  'other domains are left out of the clustered db. '
                                  'Default is {0}'.format(' '.join(DEFAULT_DOMAINS)))
    # Former names of --by_domain and --domains, still accepted
    group_adv.add_argument('--by_kingdom', dest = 'by_kingdom_deprecated',
                           action = 'store_true', help = argparse.SUPPRESS)
    group_adv.add_argument('--kingdoms', dest = 'kingdoms_deprecated', metavar = 'STR',
                           type = str, nargs = '+', default = None, help = argparse.SUPPRESS)
    # -o / --out_db_name
    group_adv.add_argument('--out_db_name',
                           action = 'store',
                           metavar = 'OUTNAME',
                           type = str,
                           help = 'Output MATAM db name. '
                                  'Default is composed from parameters')
    # --keep_tmp
    group_adv.add_argument('--keep_tmp',
                            action = 'store_true',
                            help = 'Do not remove tmp files')
    # --debug
    group_adv.add_argument('--debug',
                            action = 'store_true',
                            help = 'Output debug infos')

    #
    args = parser.parse_args()

    # Arguments checking
    if args.by_kingdom_deprecated:
        sys.stderr.write('Warning: --by_kingdom is deprecated, use --by_domain\n')
        args.by_domain = True
    if args.kingdoms_deprecated:
        sys.stderr.write('Warning: --kingdoms is deprecated, use --domains\n')
        if args.domains:
            parser.error('use --domains only (--kingdoms is its former name)')
        args.domains = args.kingdoms_deprecated
    if args.domains is None:
        args.domains = list(DEFAULT_DOMAINS)
    if args.clustered_ref:
        if args.clustering_id_threshold is not None or args.by_domain:
            parser.error('--clustered_ref replaces the clustering: it cannot be used with '
                         '--clustering_id_threshold or --by_domain')
        if not os.path.isfile(args.clustered_ref):
            parser.error('--clustered_ref: no such file: {0}'.format(args.clustered_ref))
        args.clustered_ref = os.path.abspath(args.clustered_ref)
    if args.clustering_id_threshold is None:
        args.clustering_id_threshold = DEFAULT_CLUSTERING_ID_THRESHOLD
    if args.close_relative_id is None:
        args.close_relative_id = args.clustering_id_threshold
    if args.close_relative_id < 0 or args.close_relative_id > 1:
        parser.error('--close_relative_id must be in [0,1]')
    if args.clustering_id_threshold < 0 or args.clustering_id_threshold > 1:
        parser.print_help()
        raise Exception("clustering id threshold not in range [0,1]")

    # Set debug parameters
    if args.debug:
        args.verbose = True
        args.keep_tmp = True

    # Get absolute path for all arguments
    args.input_ref = os.path.abspath(args.input_ref)
    args.db_dir = os.path.abspath(args.db_dir)

    #
    return args


def print_intro(args):
    """
    Print the introduction
    """

    sys.stdout.write("""
#################################
    MATAM db pre-processing
#################################\n\n""")

    # Retrieve complete cmd line
    cmd_line = '{binpath} '.format(binpath=matam_db_prepro_bin)

    # Verbose
    if args.verbose:
        cmd_line += '--verbose '

    # Debug
    if args.debug:
        cmd_line += '--debug '

    # Performance
    cmd_line += """--cpu {cpu} --max_memory {memory} \
""".format(cpu=args.cpu,
           memory=args.max_memory)

    # Advanced parameters
    if args.min_length:
        cmd_line += '--min_length {} '.format(args.min_length)
    if args.max_length:
        cmd_line += '--max_length {} '.format(args.max_length)

    cmd_line += '--max_consecutive_n {0} '.format(args.max_consecutive_n)

    if args.clustered_ref:
        cmd_line += '--clustered_ref {0} '.format(args.clustered_ref)
    else:
        cmd_line += '--clustering_id_threshold {0} '.format(args.clustering_id_threshold)
    cmd_line += '--close_relative_id {0} '.format(args.close_relative_id)

    if args.by_domain:
        cmd_line += '--by_domain --domains '
        for domain in args.domains:
            cmd_line += '{0} '.format(domain)

    if args.skip_index:
        cmd_line += '--skip_index '

    if args.keep_tmp:
        cmd_line += '--keep_tmp '

    if args.out_db_name:
        cmd_line += '--out_db_name {0} '.format(args.out_db_name)

    # Main parameters
    cmd_line += '--db_dir {0} '.format(args.db_dir)

    cmd_line += '--input_ref {0}'.format(args.input_ref)

    # Print cmd line
    sys.stdout.write('CMD: {0}\n\n'.format(cmd_line))

    return 0


def fasta_ids(seqkit_bin, fasta, cpu):
    """The sequence IDs of a fasta file, in file order."""
    output = subprocess.run(
        [seqkit_bin, 'seq', '--name', '--only-id', '-j', str(cpu), fasta],
        stdout=subprocess.PIPE, check=True, universal_newlines=True).stdout
    return output.split()


def fasta_header_taxonomy(seqkit_bin, fasta, cpu):
    """{id: taxonomy} from the headers of a fasta file, as extract_taxo_from_fasta.py makes it."""
    output = subprocess.run([seqkit_bin, 'seq', '--name', '-j', str(cpu), fasta],
                            stdout=subprocess.PIPE, check=True, universal_newlines=True).stdout
    taxonomy = {}
    for header in output.splitlines():
        fields = header.split()
        if fields:
            taxonomy[fields[0]] = '_'.join(fields[1:]) or discarded_refs_report.UNKNOWN
    return taxonomy


def select_clustered_refs(seqkit_bin, clustered_ref, cleaned_fasta, ids_filepath,
                          out_fasta, cpu):
    """Write to out_fasta the cleaned reference sequences whose IDs are in clustered_ref.

    The sequences come from the cleaned complete db, so the clustered db gets
    the same cleaning (U->T, Ns, headers, length filters) and IDs that are
    always found in the complete db and its taxonomies. They keep the order of
    the cleaned db (decreasing length), like vsearch centroids.
    Returns the IDs given (unique, in order) and the set of IDs selected.
    Raises ValueError if none is found.
    """
    given_ids = list(dict.fromkeys(fasta_ids(seqkit_bin, clustered_ref, cpu)))  # unique, in order
    if not given_ids:
        raise ValueError('{0} contains no sequences'.format(clustered_ref))
    with open(ids_filepath, 'w') as handle:
        handle.write('\n'.join(given_ids) + '\n')
    subprocess.run([seqkit_bin, 'grep', '--pattern-file', ids_filepath, '-j', str(cpu),
                    '-o', out_fasta, cleaned_fasta], check=True)

    found_ids = set(fasta_ids(seqkit_bin, out_fasta, cpu))
    missing = [seq_id for seq_id in given_ids if seq_id not in found_ids]
    if not found_ids:
        raise ValueError('none of the {0} IDs of {1} is in the reference db (e.g. {2})'.format(
            len(given_ids), clustered_ref, ', '.join(given_ids[:3])))
    if missing:
        logger.warning('{0} of the {1} clustered sequences are not in the cleaned reference db '
                       '(not in the input, or removed by the length or N filters) and are left '
                       'out, e.g. {2}. Their clusters have no representative in the clustered db '
                       '(their taxa are summarised below); a clustered set made from the cleaned '
                       'db (e.g. the clustered db of an earlier preprocessing) avoids this'.format(
                           len(missing), len(given_ids), ', '.join(missing[:5])))
    logger.info('Clustered db: {0} sequences selected from the reference db'.format(
        len(found_ids)))
    return given_ids, found_ids


def nearest_remaining_refs(vsearch_bin, seqkit_bin, ids, sources, remaining_fasta, tmp_basepath,
                           cpu):
    """{id: (nearest remaining reference or NO_HIT, identity in %)} of the given sequences.

    The sequences are taken from the ``sources`` fasta files (first one that
    has them) and aligned to ``remaining_fasta`` with the semi-global
    alignment of the clustering; a hit must cover QUERY_COVERAGE of the
    sequence (otherwise free end gaps let a few matching columns score 100%).
    The best of the first 16 accepted hits is kept; candidates are tried in
    order of k-mer similarity, giving up after 1024 rejected ones.
    """
    report = discarded_refs_report
    queries_filepath = tmp_basepath + '.fasta'
    hits_filepath = tmp_basepath + '.tsv'
    pattern_filepath = tmp_basepath + '.ids'
    wanted = set(ids)
    try:
        with open(queries_filepath, 'w') as queries:
            for source in sources:
                if not wanted:
                    break
                with open(pattern_filepath, 'w') as handle:
                    handle.write('\n'.join(sorted(wanted)) + '\n')
                found = subprocess.run(
                    [seqkit_bin, 'grep', '--pattern-file', pattern_filepath, '-j', str(cpu), source],
                    stdout=subprocess.PIPE, check=True, universal_newlines=True).stdout
                queries.write(found)
                wanted -= {line[1:].split()[0] for line in found.splitlines() if line.startswith('>')}
        cmd = [vsearch_bin, '--usearch_global', queries_filepath, '--db', remaining_fasta,
               '--id', '0.5', '--query_cov', str(report.QUERY_COVERAGE), '--gapopen', '20I/0E',
               '--maxaccepts', '16', '--maxrejects', '1024', '--maxhits', '1',
               '--userout', hits_filepath, '--userfields', 'query+target+id', '--output_no_hits',
               '--threads', str(cpu), '--quiet']
        logger.debug('CMD: {0}'.format(' '.join(cmd)))
        subprocess.run(cmd, check=True)
        nearest = {}
        with open(hits_filepath) as handle:
            for line in handle:
                query, target, identity = line.rstrip('\n').split('\t')[:3]
                nearest[query] = (report.NO_HIT, 0.0) if target == '*' else (target, float(identity))
    finally:
        for path in (queries_filepath, hits_filepath, pattern_filepath):
            if os.path.exists(path):
                os.remove(path)
    for seq_id in wanted:  # sequence not found in the sources
        nearest[seq_id] = (report.NO_HIT, 0.0)
    return nearest


def report_discarded_refs(seqkit_bin, report_basepath, input_taxo_filepath, cleaned_fasta,
                          clustered_fasta, cpu, clustered_ref=None, given_ids=(),
                          domain_ids=None, input_fasta=None, close_relative_id=0.0):
    """Find the reference sequences that don't reach the clustered db; report and log their taxa.

    See discarded_refs_report.py. ``given_ids`` are the --clustered_ref IDs,
    ``domain_ids`` the IDs clustered with --by_domain (None without it).
    Unless ``close_relative_id`` is 0, the nearest remaining reference of each
    discarded sequence (other domains excepted) is searched in the cleaned db.
    """
    report = discarded_refs_report
    taxonomy = report.read_taxonomy(input_taxo_filepath)    # every input sequence
    cleaned_ids = set(fasta_ids(seqkit_bin, cleaned_fasta, cpu))
    clustered_ids = set(fasta_ids(seqkit_bin, clustered_fasta, cpu))
    discarded = collections.OrderedDict(
        (seq_id, report.REMOVED_BY_CLEANING) for seq_id in taxonomy if seq_id not in cleaned_ids)
    if clustered_ref:
        header_taxonomy = None
        for seq_id in given_ids:
            if seq_id in clustered_ids:
                continue
            if seq_id in taxonomy:
                discarded[seq_id] = report.CENTROID_REMOVED_BY_CLEANING
            else:
                if header_taxonomy is None:
                    header_taxonomy = fasta_header_taxonomy(seqkit_bin, clustered_ref, cpu)
                discarded[seq_id] = report.CENTROID_NOT_IN_INPUT
                taxonomy[seq_id] = header_taxonomy.get(seq_id, report.UNKNOWN)
    if domain_ids is not None:
        for seq_id in fasta_ids(seqkit_bin, cleaned_fasta, cpu):
            if seq_id not in domain_ids:
                discarded[seq_id] = report.OUTSIDE_DOMAINS
    nearest, close_identity = None, None
    searched = [seq_id for seq_id, reason in discarded.items() if reason in report.SEARCHED_REASONS]
    if close_relative_id > 0 and searched:
        logger.info('Searching the nearest remaining reference of {0} discarded sequences'.format(
            len(searched)))
        nearest = nearest_remaining_refs(
            Binary.assert_which('vsearch'), seqkit_bin, searched,
            [path for path in (input_fasta, clustered_ref) if path], cleaned_fasta,
            report_basepath + '.discarded_search.tmp', cpu)
        close_identity = close_relative_id * 100
    summary = report.write_report(report_basepath, discarded, taxonomy, cleaned_ids, clustered_ids,
                                  nearest, close_identity)
    report.summarise(summary, logger)
    return summary


def rm_files(filepath_list):
    """
    Try to delete all files in filepath_list
    """
    for filepath in filepath_list:
        try:
            logger.debug('rm {0}'.format(filepath))
            os.remove(filepath)
        except:
            pass


if __name__ == '__main__':

    # Arguments parsing
    args = parse_arguments()

    # Print intro infos
    print_intro(args)

    # Init error code
    error_code = 0

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

    # Init list of tmp files to delete at the end
    to_rm_filepath_list = list()

    # Resolve multithreaded tool binaries
    seqkit_bin = Binary.assert_which('seqkit')

    ##############################################
    # Set all files and directories names + paths

    complete_ref_db_filepath = args.input_ref
    complete_ref_db_filename = os.path.basename(complete_ref_db_filepath)
    complete_ref_db_basename = os.path.splitext(complete_ref_db_filename)[0]

    try:
        if not os.path.exists(args.db_dir):
            logger.debug('mkdir {0}'.format(args.db_dir))
            os.makedirs(args.db_dir)
    except OSError:
        logger.exception('Could not create output directory {0}'.format(args.db_dir))
        raise

    complete_ref_db_taxo_filename = complete_ref_db_basename + '.taxo.tab'
    complete_ref_db_taxo_filepath = os.path.join(args.db_dir, complete_ref_db_taxo_filename)

    cleaned_complete_ref_db_basename = complete_ref_db_basename
    cleaned_complete_ref_db_filename = cleaned_complete_ref_db_basename + '.cleaned.fasta'
    cleaned_complete_ref_db_filepath = os.path.join(args.db_dir, cleaned_complete_ref_db_filename)

    clustering_id_threshold_int = int(args.clustering_id_threshold * 100)

    vsearch_output_basename = complete_ref_db_basename + '.vsearch_'
    vsearch_output_basename += '{0}'.format(clustering_id_threshold_int)
    if args.by_domain:
        vsearch_output_basename += '_by_domain'
    vsearch_output_filename = vsearch_output_basename + '.fasta'
    vsearch_output_filepath = os.path.join(args.db_dir, vsearch_output_filename)

    vsearch_centroids_basename = vsearch_output_basename + '.centroids'
    vsearch_centroids_basepath = os.path.join(args.db_dir, vsearch_output_basename)
    vsearch_centroids_filename = vsearch_centroids_basename + '.fasta'
    vsearch_centroids_filepath = os.path.join(args.db_dir, vsearch_centroids_filename)

    if args.clustered_ref:
        # Named after the given clustered set, e.g. SILVA_138_NR95.fasta -> SILVA_138_NR95
        clustered_ref_db_basename = os.path.splitext(os.path.basename(args.clustered_ref))[0]
    else:
        clustered_ref_db_basename = cleaned_complete_ref_db_basename + '_NR{0}'.format(clustering_id_threshold_int)
        if args.by_domain:
            clustered_ref_db_basename += '_bd'
    clustered_ref_db_basepath = os.path.join(args.db_dir, clustered_ref_db_basename)
    clustered_ref_db_filename = clustered_ref_db_basename + '.fasta'
    clustered_ref_db_filepath = os.path.join(args.db_dir, clustered_ref_db_filename)

    # This is the output MATAM db basepath to pass to matam_assembly.py
    output_ref_db_basename = clustered_ref_db_basename
    if args.out_db_name:
        output_ref_db_basename = os.path.basename(args.out_db_name)
    output_ref_db_basepath = os.path.join(args.db_dir, output_ref_db_basename)
    # This is the output MATAM db file names
    # For the complete db fasta file
    output_complete_ref_db_basename = output_ref_db_basename + '.complete'
    output_complete_ref_db_basepath = os.path.join(args.db_dir, output_complete_ref_db_basename)
    output_complete_ref_db_filename = output_complete_ref_db_basename + '.fasta'
    output_complete_ref_db_filepath = os.path.join(args.db_dir, output_complete_ref_db_filename)
    # For the complete db taxo file
    output_complete_ref_db_taxo_filename = output_complete_ref_db_basename + '.taxo.tab'
    output_complete_ref_db_taxo_filepath = os.path.join(args.db_dir, output_complete_ref_db_taxo_filename)
    # For the clustered db fasta file
    output_clustered_ref_db_basename = output_ref_db_basename + '.clustered'
    output_clustered_ref_db_basepath = os.path.join(args.db_dir, output_clustered_ref_db_basename)
    output_clustered_ref_db_filename = output_clustered_ref_db_basename + '.fasta'
    output_clustered_ref_db_filepath = os.path.join(args.db_dir, output_clustered_ref_db_filename)

    ##############################################################
    # Ref DB pre-processing (cleaning, extracting taxo, indexing)

    logger.info('Starting ref db pre-processing')

    # Extract taxo from ref DB and sort by ref id
    logger.info('Extracting taxonomies from reference DB')

    cmd_line = extract_taxo_bin + ' -i ' + complete_ref_db_filepath + ' | sort -k1,1 > '
    cmd_line += complete_ref_db_taxo_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # TO DO, maybe, one day:
    # Trim Ns from both sides

    # Convert Us in Ts
    # Option: Either filter out seq with Ns or replace Ns with random nucl
    # Option: Filter too small or too long sequences
    # Sort sequences by decreasing length
    logger.info('Cleaning reference db')

    # Headers: turn tabs (and other whitespace runs) into single spaces.
    # SortMeRNA takes the reference name up to the first space only, so a tab
    # would end up inside the SAM RNAME and shift every following column.
    cmd_line = 'cat ' + complete_ref_db_filepath
    cmd_line += ' | sed -e "/^>/s/[[:space:]]\\+/ /g" -e "/^>/s/ $//"'
    cmd_line += ' | sed "/^>/!s/U/T/g" | sed "/^>/!s/u/t/g" | sed "/^>/!s/ //g"'
    cmd_line += ' | ' + replace_Ns_bin + ' -n {0} '.format(args.max_consecutive_n)
    if args.min_length or args.max_length:
        cmd_line += ' | ' + seqkit_bin + ' seq'
        if args.min_length:
            cmd_line += ' -m ' + str(args.min_length)
        if args.max_length:
            cmd_line += ' -M ' + str(args.max_length)
        cmd_line += ' -j ' + str(args.cpu)
    cmd_line += ' | ' + seqkit_bin + ' sort -l -r -j ' + str(args.cpu) + ' > '
    cmd_line += cleaned_complete_ref_db_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    ####################
    # Ref DB clustering

    if not args.clustered_ref:
        vsearch_bin = Binary.assert_which('vsearch')
        logger.info('Starting ref db clustering')

    # Use the given clustered set, or perform clustering, either by domain or at once
    if args.clustered_ref:
        logger.info('Using the clustered sequences of {0} (no clustering)'.format(
            args.clustered_ref))
        vsearch_output_filepath = os.path.join(
            args.db_dir, complete_ref_db_basename + '.clustered_ref.fasta')
        clustered_ref_ids_filepath = os.path.join(
            args.db_dir, complete_ref_db_basename + '.clustered_ref.ids')
        to_rm_filepath_list.append(clustered_ref_ids_filepath)
        try:
            clustered_ref_given_ids, _ = select_clustered_refs(seqkit_bin, args.clustered_ref,
                                  cleaned_complete_ref_db_filepath,
                                  clustered_ref_ids_filepath, vsearch_output_filepath, args.cpu)
        except (ValueError, subprocess.CalledProcessError) as error:
            logger.error('Cannot use --clustered_ref: {0}'.format(error))
            sys.exit(1)
    elif args.by_domain:

        # Perform clustering by domain
        clustered_domain_ids = set()
        for domain in args.domains:

            # Define by-domain files names + paths
            cleaned_complete_ref_db_domain_basename = cleaned_complete_ref_db_basename + '.' + domain
            cleaned_complete_ref_db_domain_filename = cleaned_complete_ref_db_domain_basename + '.fasta'
            cleaned_complete_ref_db_domain_filepath = os.path.join(args.db_dir, cleaned_complete_ref_db_domain_filename)

            vsearch_output_domain_basename = cleaned_complete_ref_db_domain_basename + '.vsearch_'
            vsearch_output_domain_basename += '{0}'.format(clustering_id_threshold_int)
            vsearch_output_domain_filename = vsearch_output_domain_basename + '.fasta'
            vsearch_output_domain_filepath = os.path.join(args.db_dir, vsearch_output_domain_filename)

            # Extract the sequences of the domain
            logger.info('Extracting sequences from the {0} domain'.format(domain))

            # The domain is the first rank of the taxonomy that follows the ID
            # ('>ID Bacteria;Pseudomonadota;...'), not any occurrence of its name.
            cmd_line = seqkit_bin + ' grep -n -r -i -p ' + shlex.quote(
                r'^\S+\s+' + re.escape(domain) + '(;|$)')
            cmd_line += ' -j ' + str(args.cpu)
            cmd_line += ' ' + cleaned_complete_ref_db_filepath
            cmd_line += ' > ' + cleaned_complete_ref_db_domain_filepath

            logger.debug('CMD: {0}'.format(cmd_line))
            error_code += subprocess.call(cmd_line, shell=True)
            clustered_domain_ids.update(
                fasta_ids(seqkit_bin, cleaned_complete_ref_db_domain_filepath, args.cpu))

            # Clutering with vsearch
            # Aplly a null-penalty value to left and right gaps (--gaopoen "20I/0E") to use vsearch as a semi-global aligner
            logger.info('Clustering {0} sequences @ {1} pct id'.format(domain, clustering_id_threshold_int))

            cmd_line = vsearch_bin + ' --cluster_fast'
            cmd_line += ' ' + cleaned_complete_ref_db_domain_filepath
            cmd_line += ' --centroids ' + vsearch_output_domain_filepath
            cmd_line += ' --id {0:.2f}'.format(args.clustering_id_threshold)
            cmd_line += ' --gapopen "20I/0E"'
            cmd_line += ' --threads ' + str(args.cpu)

            logger.debug('CMD: {0}'.format(cmd_line))
            if args.verbose:
                error_code += subprocess.call(cmd_line, shell=True)
            else:
                # Needed because VSearch doesnt have a verbose option
                # and output everything to stderr
                error_code += subprocess.call(cmd_line, shell=True, stdout=FNULL, stderr=FNULL)

            # Concatenate vsearch outputs
            cmd_line = 'cat ' + vsearch_output_domain_filepath + ' >> '
            cmd_line += vsearch_output_filepath

            logger.debug('CMD: {0}'.format(cmd_line))
            error_code += subprocess.call(cmd_line, shell=True)

            # Tag tmp files for removal
            to_rm_filepath_list.append(cleaned_complete_ref_db_domain_filepath)
            to_rm_filepath_list.append(vsearch_output_domain_filepath)

    else:
        # Clutering with vsearch
        # Aplly a null-penalty value to left and right gaps (--gaopoen "20I/0E") to use vsearch as a semi-global aligner
        logger.info('Clustering sequences @ {0} pct id'.format(clustering_id_threshold_int))
        cmd_line = vsearch_bin + ' --cluster_fast'
        cmd_line += ' ' + cleaned_complete_ref_db_filepath
        cmd_line += ' --centroids ' + vsearch_output_filepath
        cmd_line += ' --id {0:.2f}'.format(args.clustering_id_threshold)
        cmd_line += ' --gapopen "20I/0E"'
        cmd_line += ' --threads ' + str(args.cpu)

        logger.debug('CMD: {0}'.format(cmd_line))
        if args.verbose:
            error_code += subprocess.call(cmd_line, shell=True)
        else:
            # Needed because VSearch doesnt have a verbose option
            # and output everything to stderr
            error_code += subprocess.call(cmd_line, shell=True, stdout=FNULL, stderr=FNULL)
    vsearch_centroids_filepath = vsearch_output_filepath

    # Clean fasta headers (strip description, keep ID only)
    cmd_line = seqkit_bin + ' seq --only-id -j ' + str(args.cpu)
    cmd_line += ' -o ' + clustered_ref_db_filepath + ' ' + vsearch_centroids_filepath

    logger.debug('CMD: {0}'.format(cmd_line))
    error_code += subprocess.call(cmd_line, shell=True)

    # Tag tmp files for removal
    to_rm_filepath_list.append(vsearch_output_filepath)
    to_rm_filepath_list.append(vsearch_centroids_filepath)

    # Safeguard: report the reference sequences that don't reach the clustered
    # db (removed by cleaning, --clustered_ref centroids not in the cleaned db,
    # other domains with --by_domain) and summarise their taxonomy
    try:
        report_discarded_refs(
            seqkit_bin, output_ref_db_basepath, complete_ref_db_taxo_filepath,
            cleaned_complete_ref_db_filepath, clustered_ref_db_filepath, args.cpu,
            clustered_ref=args.clustered_ref,
            given_ids=clustered_ref_given_ids if args.clustered_ref else (),
            domain_ids=clustered_domain_ids if args.by_domain else None,
            input_fasta=complete_ref_db_filepath, close_relative_id=args.close_relative_id)
    except (OSError, subprocess.CalledProcessError) as error:
        logger.error('Could not report the discarded reference sequences: {0}'.format(error))
        error_code += 1

    ##########################################
    # Renaming output files as MATAM db files

    logger.info('Renaming output files as MATAM db files')

    try:
        os.rename(cleaned_complete_ref_db_filepath, output_complete_ref_db_filepath)
        os.rename(complete_ref_db_taxo_filepath, output_complete_ref_db_taxo_filepath)
        os.rename(clustered_ref_db_filepath, output_clustered_ref_db_filepath)
    except OSError:
        logger.exception('Could not rename tmp files into MATAM db files')
        raise

    ######################################################
    # SortMeRNA indexing of complete and clustered ref db

    # matam_assembly.py uses <db>.{complete,clustered}.smr/idx if its record
    # (written by sortmerna_index.build) matches the fasta and the indexing
    # parameters, and otherwise builds the index itself with the parameters
    # stored here (also with --skip_index); see sortmerna_index.py.
    index_params = sortmerna_index.index_params(args.max_memory)
    sortmerna_bin = Binary.assert_which('sortmerna')

    def call_sortmerna(cmd_line, verbose):
        logger.debug('CMD: {0}'.format(cmd_line))
        returncode = subprocess.call(cmd_line, shell=True)
        if verbose:
            sys.stdout.write('\n')
        return returncode

    for label, ref_db_filepath, ref_db_basepath in (
            ('complete', output_complete_ref_db_filepath, output_complete_ref_db_basepath),
            ('clustered', output_clustered_ref_db_filepath, output_clustered_ref_db_basepath)):
        smr_dir = ref_db_basepath + '.smr'
        idx_dir = os.path.join(smr_dir, 'idx')
        try:
            sortmerna_index.write_params(smr_dir, index_params)
        except OSError as error:
            logger.error('Could not store the indexing parameters of {0} ref db: {1}'.format(
                label, error))
            error_code += 1
            continue
        if args.skip_index:
            logger.info('Not indexing {0} ref db (--skip_index); matam_assembly.py will build '
                        'the index on first use'.format(label))
            continue
        logger.info('Indexing {0} ref db'.format(label))
        t0_wall = time.time()
        try:
            with sortmerna_index.locked(smr_dir):
                sortmerna_index.build(sortmerna_bin, ref_db_filepath, idx_dir, index_params,
                                      verbose=args.verbose, call=call_sortmerna)
        except (sortmerna_index.IndexBuildError, OSError) as error:
            logger.error('Could not index {0} ref db: {1}'.format(label, error))
            error_code += 1
            continue
        logger.info('Indexed {0} ref db in {1:.1f} seconds wall time: {2}'.format(
            label, time.time() - t0_wall, idx_dir))

    ###############
    # Exit program

    sys.stdout.write('Output MATAM db: {0}\n'.format(output_ref_db_basepath))

    # Exit if everything went ok
    if not error_code:
        # Try to remove all tmp files
        # won't crash if it cannot
        if not args.keep_tmp:
            sys.stdout.write('\n')
            logger.info('Removing tmp files')
            rm_files(to_rm_filepath_list)
        #
        sys.stdout.write('\n{0} terminated with no error\n'.format(program_filename))
        exit(0)
    # Deal with errors
    else:
        sys.stdout.write('\n{0} terminated with some errors. '.format(program_filename))
        if args.verbose:
            sys.stdout.write('Check the log for additional infos\n')
        else:
            sys.stdout.write('Rerun the program using --verbose or --debug option\n')
        exit(1)
