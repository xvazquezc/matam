#!/usr/bin/env python3
import os
import sys
import gzip
import subprocess
import logging
from collections import defaultdict
import multiprocessing
import argparse

from fasta_utils import read_fasta_file_handle, format_seq
from fastq_utils import read_fastq_file_handle

from assembler_factory import AssemblerFactory

logger = logging.getLogger(__name__)


def _open_fastq(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def is_empty(fpath):
    try:
        return not (os.path.getsize(fpath) > 0)
    except OSError:
        pass
    return True

def extract_interleaved_paired_reads_by_component(interleaved_fastq, read_metanode_component_filepath):
    """Keep only complete pairs whose two reads belong to one component.

    SortMeRNA's ``--paired_in`` output is pair-preserving but interleaved.  It
    is used only as the hand-off from mapping; component files and SGA inputs
    are always emitted as separate R1/R2 FastQ files.
    """
    if not os.path.isfile(interleaved_fastq):
        logger.fatal('The input reads file does not exists:%s' % interleaved_fastq)
        sys.exit("An error occured. Can't extract component's reads")

    if not os.path.isfile(read_metanode_component_filepath):
        logger.fatal('The file storing the correspondance between reads and components does not exists:%s' % read_metanode_component_filepath)
        sys.exit("An error occured. Can't extract component's reads")

    # Reading read --> component file
    logger.debug('Reading read-->component from {}'.format(read_metanode_component_filepath))
    read_component_dict = dict()
    with open(read_metanode_component_filepath, 'r') as read_metanode_component_fh:
        read_component_dict = {t[0]:t[2] for t in (l.split() for l in read_metanode_component_fh) if t[2] != 'NULL'}

    # Storing reads for each component
    logger.debug('Storing complete pairs by component from {}'.format(interleaved_fastq))
    component_reads_dict = defaultdict(lambda: ([], []))
    kept_pairs = 0
    discarded_pairs = 0
    with _open_fastq(interleaved_fastq) as fastq_fh:
        reads = iter(read_fastq_file_handle(fastq_fh))
        while True:
            try:
                forward_read = next(reads)
            except StopIteration:
                break
            try:
                reverse_read = next(reads)
            except StopIteration:
                logger.fatal('SortMeRNA paired output ends with an orphan read: %s',
                             interleaved_fastq)
                sys.exit("An error occured. Can't extract component's reads")

            forward_component = read_component_dict.get(forward_read[0])
            reverse_component = read_component_dict.get(reverse_read[0])
            if (forward_component is not None and
                    forward_component == reverse_component):
                component_reads_dict[forward_component][0].append(forward_read)
                component_reads_dict[forward_component][1].append(reverse_read)
                kept_pairs += 1
            else:
                # A pair with an unselected mate or mates in distinct graph
                # components cannot preserve the original component semantics.
                discarded_pairs += 1
    logger.info('Kept %s complete same-component read pairs; discarded %s',
                kept_pairs, discarded_pairs)
    return component_reads_dict


def extract_separate_paired_reads_by_component(forward_fastq, reverse_fastq,
                                                read_metanode_component_filepath):
    """Keep complete, same-component pairs from separate R1/R2 files."""
    for fastq in (forward_fastq, reverse_fastq):
        if not os.path.isfile(fastq):
            logger.fatal('The input reads file does not exists:%s' % fastq)
            sys.exit("An error occured. Can't extract component's reads")
    if not os.path.isfile(read_metanode_component_filepath):
        logger.fatal('The file storing the correspondance between reads and components does not exists:%s' % read_metanode_component_filepath)
        sys.exit("An error occured. Can't extract component's reads")

    with open(read_metanode_component_filepath, 'r') as component_fh:
        read_component_dict = {
            fields[0]: fields[2]
            for fields in (line.split() for line in component_fh)
            if fields[2] != 'NULL'
        }

    component_reads_dict = defaultdict(lambda: ([], []))
    kept_pairs = 0
    discarded_pairs = 0
    with _open_fastq(forward_fastq) as forward_fh, _open_fastq(reverse_fastq) as reverse_fh:
        forward_reads = iter(read_fastq_file_handle(forward_fh))
        reverse_reads = iter(read_fastq_file_handle(reverse_fh))
        while True:
            try:
                forward_read = next(forward_reads)
            except StopIteration:
                try:
                    next(reverse_reads)
                except StopIteration:
                    break
                logger.fatal('Paired FastQ files have different record counts')
                sys.exit("An error occured. Can't extract component's reads")
            try:
                reverse_read = next(reverse_reads)
            except StopIteration:
                logger.fatal('Paired FastQ files have different record counts')
                sys.exit("An error occured. Can't extract component's reads")

            forward_component = read_component_dict.get(forward_read[0])
            reverse_component = read_component_dict.get(reverse_read[0])
            if (forward_component is not None and
                    forward_component == reverse_component):
                component_reads_dict[forward_component][0].append(forward_read)
                component_reads_dict[forward_component][1].append(reverse_read)
                kept_pairs += 1
            else:
                discarded_pairs += 1
    logger.info('Kept %s complete same-component read pairs; discarded %s',
                kept_pairs, discarded_pairs)
    return component_reads_dict


def extract_lca_by_component(components_lca_filepath):
    if not os.path.isfile(components_lca_filepath):
        logger.fatal('The file storing the correspondance between lca info and components does not exists:%s' % components_lca_filepath)
        sys.exit("An error occured. Can't extract component's lca")

    # Reading components LCA and storing them in a dict
    logger.debug('Reading components LCA assignment from {0}'.format(components_lca_filepath))
    component_lca_dict = dict()
    with open(components_lca_filepath, 'r') as component_lca_fh:
        component_lca_dict = {t[0]:t[1] for t in (l.split() for l in component_lca_fh) if len(t) == 2}
    return component_lca_dict


def save_components(component_reads_dict, directory):
    """
    Take a dict (key=component_id, value=[header,seq,qual])
    and save each component to a file into the given directory:
    directory/
        component%s_R1.fq
        component%s_R2.fq

    Return a dict (key=component_id, value=(forward_fastq, reverse_fastq)).
    """

    try:
        os.mkdir(directory)
    except FileExistsError as fee:
        pass

    components_fq = {}
    for component_id, (forward_reads, reverse_reads) in component_reads_dict.items():
        forward_path = os.path.join(directory, "component%s_R1.fq" % component_id)
        reverse_path = os.path.join(directory, "component%s_R2.fq" % component_id)
        components_fq[component_id] = (forward_path, reverse_path)

        with open(forward_path, 'w') as forward_fh, open(reverse_path, 'w') as reverse_fh:
            for (header, seq, qual) in forward_reads:
                forward_fh.write('@{}\n{}\n+\n{}\n'.format(header, seq, qual))
            for (header, seq, qual) in reverse_reads:
                reverse_fh.write('@{}\n{}\n+\n{}\n'.format(header, seq, qual))
    return components_fq


def isfastq(filepath):
    """
    Determine if filepath is a fastq file based on the extension
    """

    suffix = ('.fq', '.fastq')
    return any( [ filepath.endswith(s) for s in suffix] )


def isfasta(filepath):
    """
    Determine if filepath is a fasta file based on the extension
    """

    suffix = ('.fa', '.fasta', '.fna')
    return any( [filepath.endswith(s) for s in suffix] )


def nucleotidic_number(fastx):
    """
    Compute the total number of nucleotides in the fastx file
    """

    parser = None
    if isfasta(fastx):
        parser = read_fasta_file_handle
    elif isfastq(fastx):
        parser = read_fastq_file_handle

    if parser is None:
        logger.fatal("Can't dertermine whether this file is a fasta or fastq")
        sys.exit('Aborting assembly step')


    count = 0
    with _open_fastq(fastx) as fastx_handle:
        for rec in parser(fastx_handle):
            seq = rec[1]
            count += len(seq)
    return count

def estimate_coverage(forward_fq, reverse_fq, contigs_fa):
    """
    estimated_cov = reads_nt/contigs_nt
    """

    reads_nt = nucleotidic_number(forward_fq) + nucleotidic_number(reverse_fq)
    contigs_nt = nucleotidic_number(contigs_fa)
    estimated_cov = None

    try:
        estimated_cov = reads_nt/contigs_nt
        #logger.debug("Estimated coverage:%s" % estimated_cov)

    except ZeroDivisionError:
        #logger.warning("Can't estimate the coverage for 0 length contigs:%s" % contigs_fa)
        pass

    return estimated_cov

def assemble_component(assembler_name, forward_fastq, reverse_fastq, workdir,
                       read_correction, cpu, coverage_threshold):

    try:
        if read_correction != 'auto' and coverage_threshold is not None:
            logger.warning("Coverage_threshold %s makes no sense when read_correction is not auto. This argument will be ignored" % coverage_threshold)

        logger.debug('Assembling: %s, %s' % (forward_fastq, reverse_fastq))
        assembler_factory = AssemblerFactory()
        assembler = assembler_factory.get(assembler_name)
        assembler.build_command_line(forward_fastq, reverse_fastq, workdir,
                                     read_correction, cpu)

        fasta_file = assembler.run()
        if is_empty(fasta_file): return

        estimated_cov = estimate_coverage(forward_fastq, reverse_fastq, fasta_file)
        logger.debug("Estimated coverage:%s, %s" % (estimated_cov, forward_fastq))

        # Re-run the assembly with error correction activated when read_correction == auto
        if read_correction == 'auto' and estimated_cov is not None and estimated_cov > coverage_threshold:
            assembler.build_command_line(forward_fastq, reverse_fastq, workdir,
                                         'yes', cpu)
            fasta_file = assembler.run()
            if is_empty(fasta_file): return
            estimated_cov2 = estimate_coverage(forward_fastq, reverse_fastq, fasta_file)
            logger.debug("Estimated coverage, before:%s, after:%s, component:%s, cov_threshold:%s" % (estimated_cov, estimated_cov2, forward_fastq, coverage_threshold))

            # suspicious when estimate_cov is not None and estimated_cov2 is None
            if estimated_cov2 is None:
                logger.warning("0 length contigs from reads component: %s" % forward_fastq)
    except SystemExit as se:
        logger.fatal('MATAM failed to assemble component: %s', forward_fastq)
        raise Exception # SystemExit cause the pool to hang, prefer a generic exception instead

    return fasta_file


def concat_components_fasta_with_lca(assembled_components_fasta, contigs_fasta, component_lca_dict):
    if os.path.isfile(contigs_fasta):
        #logger.debug("Remove old contig fasta file:%s" % contigs_fasta)
        os.unlink(contigs_fasta)

    contigs_fh = open(contigs_fasta, 'w')
    contig_count = 0
    for component_id, component_fasta in assembled_components_fasta.items():
        component_lca = component_lca_dict.get(component_id, 'NULL')
        with open(component_fasta, 'r') as assembler_contigs_fh:
            for header, seq in read_fasta_file_handle(assembler_contigs_fh):
                if len(seq):
                    contig_count += 1
                    contigs_fh.write('>{0} component={1} '.format(contig_count, component_id))
                    contigs_fh.write('lca={0}\n{1}\n'.format(component_lca, format_seq(seq)))

    contigs_fh.close()


def _get_workdir(fq):
    """
    Convenient function to build a workdir from the name of the fastqfile
    """
    assembler_wkdir_basename, _ = os.path.splitext(fq)
    return '%s_assembly_wkdir' % assembler_wkdir_basename


def assemble_components(assembler_name, components_dict, components_lca_filepath,
                        out_contigs_fasta, workdir, cpu, read_correction,
                        coverage_threshold):
    logger.info("Save components to fastq files")
    components_reads_fq = save_components(components_dict, workdir).items()
    assembled_components_fasta = {}

    logger.info("Assemble components")

    # Foreach component, build the parameters used by assemble_component and save
    # them into a list to be able to apply a map function on it
    params = []
    component_id_list = []
    for component_id, (forward_fq, reverse_fq) in components_reads_fq:
        params.append((assembler_name, forward_fq, reverse_fq,
                       _get_workdir(forward_fq), read_correction, 1,
                       coverage_threshold))
        component_id_list.append(component_id)

    with multiprocessing.Pool(processes=cpu) as pool:
        fasta_list = pool.starmap(assemble_component, params)


    # Make the correspondance between the component_id and the fasta file
    assembled_components_fasta = dict(zip(component_id_list, fasta_list))

    # Filter out empty components
    total_cpnt_nb = len(assembled_components_fasta)
    assembled_components_fasta = { cpnt:fasta for (cpnt,fasta) in assembled_components_fasta.items() if fasta }
    non_empty_cpnt_nb = len(assembled_components_fasta)
    if non_empty_cpnt_nb != total_cpnt_nb:
        logger.warning('Number of non empty components: %s/%s', non_empty_cpnt_nb, total_cpnt_nb)

    lca_dict = extract_lca_by_component(components_lca_filepath)
    logger.debug("Pool components contigs into: %s" % out_contigs_fasta)
    concat_components_fasta_with_lca(assembled_components_fasta,
                                     out_contigs_fasta, lca_dict)


def assemble_all_components(assembler_name, interleaved_fastq,
                            read_metanode_component_filepath, components_lca_filepath,
                            out_contigs_fasta, workdir,
                            cpu, read_correction, coverage_threshold):
    components_dict = extract_interleaved_paired_reads_by_component(
        interleaved_fastq, read_metanode_component_filepath)
    assemble_components(assembler_name, components_dict, components_lca_filepath,
                        out_contigs_fasta, workdir, cpu, read_correction,
                        coverage_threshold)


if __name__ == '__main__':

    FORMAT = '%(asctime)s - %(processName)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.DEBUG, format=FORMAT)

    parser = argparse.ArgumentParser(description='Assemble components')
    parser.add_argument('-a', '--assembler',
                        choices=[a.name() for a in AssemblerFactory.ASSEMBLER_ENGINES],
                        help="Select the assembler to be used. Default is %(default)s",
                        default="SGA")
    parser.add_argument('--forward',
                        type=argparse.FileType('r'),
                        help='forward reads fastq file',
                        required=True)
    parser.add_argument('--reverse',
                        type=argparse.FileType('r'),
                        help='reverse reads fastq file',
                        required=True)
    parser.add_argument('-m', '--reads_metanode',
                        type=argparse.FileType('r'),
                        help='This file makes the correspondance between reads and the components',
                        required=True)
    parser.add_argument('-l', '--components_lca',
                        type=argparse.FileType('r'),
                        help='This file make the correspondance between lca and the components',
                        required=True)
    parser.add_argument('-w', '--workdir',
                        action = 'store',
                        type = str,
                        help = 'Working  directory')
    parser.add_argument('-o', '--output_fasta',
                        help='output fasta file',
                        required=True)
    parser.add_argument('--cpu',
                        action = 'store',
                        metavar = 'CPU',
                        type = int,
                        default = 1,
                        help = 'Max number of CPU to use. '
                        'Default is %(default)s cpu')
    parser.add_argument('--read_correction',
                        action = 'store',
                        type = str,
                        choices = ['no', 'yes', 'auto'],
                        default = 'no',
                        help = 'Set the assembler read correction step. '
                        'Default is %(default)s')
    parser.add_argument( '--contig_coverage_threshold',
                              action = 'store',
                              type = int,
                              choices = [20, 50],
                              default = 20,
                              help = "When the contig's coverage is sufficient (default %(default)s, then run the assembler with read_correction enabled for this contig.")

    args = parser.parse_args()

    components = extract_separate_paired_reads_by_component(
        args.forward.name, args.reverse.name, args.reads_metanode.name)
    assemble_components(args.assembler, components, args.components_lca.name,
                        args.output_fasta, args.workdir, args.cpu,
                        args.read_correction, args.contig_coverage_threshold)
