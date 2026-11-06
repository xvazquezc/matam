#!/usr/bin/env python3

import os
import sys
import logging
import runner
import collections
import tempfile
import shutil
import re
import shlex

from fasta_clean_name import read_fasta_file_handle, format_seq
from sam_mates import iter_separate_pairs, mate_of, open_text, tag_mates

logger = logging.getLogger(__name__)


def stream_file_command(filepath):
    """Return a shell command that streams either a plain or gzip file."""
    quoted_path = shlex.quote(filepath)
    return 'gzip -cd ' + quoted_path if filepath.endswith('.gz') else 'cat ' + quoted_path

def index_ref(sortmerna_bin_path, input_fasta_ref_path, workdir, verbose=False):
    logger.info('--- Indexing scaffolds ---')
    logger.debug('File to index:  %s' % input_fasta_ref_path)
    cmd_line = '{bin} --ref {input} --workdir {workdir} --task 5'.format(
        bin=sortmerna_bin_path, input=input_fasta_ref_path, workdir=workdir)

    if verbose:
        cmd_line += ' -v'

    runner.logged_check_call(cmd_line, verbose=verbose)


def reads_mapping(sortmerna_bin, fasta_ref_path, workdir, forward_path, reverse_path, output_basepath, best, evalue, cpu, verbose=False):

    logger.info('--- Reads mapping against scaffolds ---')
    logger.debug('Database:%s' % workdir)
    reads = '{} --reads {}'.format(shlex.quote(forward_path), shlex.quote(reverse_path))
    parameters = { 'fasta_ref': fasta_ref_path, 'workdir': workdir, 'reads': reads, 'output': output_basepath,
                   'best': best, 'evalue': evalue, 'cpu': cpu }

    # Only the SAM output is consumed: its SEQ field lets sam_mates.py tell
    # the mates of a pair apart (BLAST tabular output cannot). Do not create
    # unused FASTQ or BLAST reports; keep the SAM compressed.
    cmd_line = '{bin} --ref {fasta_ref} --workdir {workdir} --reads {reads} --aligned {output}' \
               ' --sam --zip-out true' \
               ' --num_alignments {best} -e {evalue:.2e} --threads {cpu}'.format(bin=sortmerna_bin, **parameters)
    if verbose:
        cmd_line += ' -v'

    runner.logged_check_call(cmd_line, verbose=verbose)


def best_alignments_by_mate(tagged_sam_lines, threshold=0.99):
    """Yield ((read, mate), scaffold) for each read's best alignments.

    Records must be grouped per read, as sam_mates.tag_mates writes them.
    An alignment is kept when its AS score is >= threshold * the read's best
    score (0.99 tolerates the few nucleotides SortMeRNA soft-clips at ends).
    """
    def flush(key, group):
        best = max(score for score, _ in group)
        for score, scaffold in group:
            if score >= best * threshold:
                yield key, scaffold

    key, group = None, []
    for line in tagged_sam_lines:
        if not line.strip() or line.startswith('@'):
            continue
        fields = line.rstrip('\n').split('\t')
        current = (fields[0], mate_of(fields[1]))
        score = re.search(r'\tAS:i:(-?\d+)', line)
        if score is None:
            raise ValueError('SAM record has no AS:i score: {!r}'.format(line.strip()))
        if current != key:
            if group:
                yield from flush(key, group)
            key, group = current, []
        group.append((int(score.group(1)), fields[2]))
    if group:
        yield from flush(key, group)


def abundance_from_assignments(assignments):
    """Compute scaffold abundances from (read, scaffold) assignments.

    A read is any hashable key; in MATAM it is (QNAME, mate), so each mate of
    a pair counts as one read. A read hitting n distinct scaffolds adds 1/n to
    each of them; hitting the same scaffold several times counts once.
    """
    scaffolds_by_read = collections.defaultdict(list)
    reads_by_scaffold = collections.defaultdict(list)
    for read, scaffold in assignments:
        scaffolds_by_read[read].append(scaffold)
        reads_by_scaffold[scaffold].append(read)

    #asign a weight for each read
    read_weight = {}
    for read, scaffolds in scaffolds_by_read.items():
        uniq_scaffolds = set(scaffolds)
        weight = 1 / len(uniq_scaffolds)
        read_weight[read] = weight

        if len(scaffolds) != len(uniq_scaffolds):
            # keep scaffolds names where this read map more than once
            more_than_once = {k:v for k,v in collections.Counter(scaffolds).items() if v > 1}
            logger.debug('%s is mapped more than once on the same scaffold (%s) but it will contribute \
to the abundance of this scaffold only as 1 weight where weight=1/uniq_scaffolds_nb=%s' % (read, more_than_once, weight))

    # compute abundance for each scaffold depending on read weight
    abundance_by_scaffold = collections.defaultdict(lambda: 0)
    for scaffold, reads in reads_by_scaffold.items():
        for read in set(reads):
            abundance_by_scaffold[scaffold]+= read_weight[read]

        abundance_by_scaffold[scaffold] = round(abundance_by_scaffold[scaffold], 2)
    return abundance_by_scaffold


def abundance_calculation(blast_path):
    """Abundance from a BLAST tabular file of best matches (read in column 1)."""
    def assignments():
        with open(blast_path, 'r') as f:
            for line in f:
                if line.strip() == '': continue
                read, scaffold = [ s.strip() for s in line.split('\t')[:2] ]
                yield read, scaffold
    return abundance_from_assignments(assignments())


def get_abundance_by_scaffold(sortmerna_bin,
                              input_fasta_ref, forward_fastq, reverse_fastq,
                              best=10, evalue=1e-05,
                              cpu=4,
                              output_dir_basepath="/tmp/",
                              verbose=False,
                              keep_tmp=False):
    """Map the paired reads (R1/R2 files) on the scaffolds and count each mate.

    The same files are used to tell the mates apart in the SAM output.
    """
    output_dir_basepath = os.path.join(output_dir_basepath , '') # add a trailing slash
    outdir = tempfile.mkdtemp(dir=output_dir_basepath, prefix='abundance_')

    #index ref
    smr_workdir = os.path.join(outdir, 'smr_workdir')

    index_ref(sortmerna_bin, input_fasta_ref, smr_workdir, verbose=verbose)

    #reads mapping
    filtered_basepath = os.path.join(outdir, 'filt_prefix')
    reads_mapping(sortmerna_bin, input_fasta_ref, smr_workdir, forward_fastq, reverse_fastq, filtered_basepath, best, evalue, cpu, verbose=verbose)

    #best matches per (read, mate), then abundance
    sam_path = '%s.sam.gz' % filtered_basepath
    with open_text(sam_path) as sam:
        tagged = tag_mates(sam, iter_separate_pairs(forward_fastq, reverse_fastq))
        abundance = abundance_from_assignments(best_alignments_by_mate(tagged, 0.99))

    #delete tempfiles
    if not keep_tmp:
        shutil.rmtree(outdir)

    return abundance


def complete_fasta_with_abundance(input_fasta, output_fasta, abundance):
    in_fasta_handler = open(input_fasta, 'r')
    out_fasta_handler = open(output_fasta, 'w')

    for header, seq in read_fasta_file_handle(in_fasta_handler):
        id = header.split()[0].strip()
        ab = 0
        if id not in abundance:
            logger.warning("Can't find the abundance for:%s. Set to 0 by default." % id)
        else:
            ab = abundance[id]
        header = '{header} count={abundance}'.format(header=header, abundance=ab)
        out_fasta_handler.write( '>{header}\n{seq}\n'.format(header=header, seq=format_seq(seq)) )

    in_fasta_handler.close()
    out_fasta_handler.close()


def get_abundance_from_fasta(fasta, regexp=r'count=(\d+\.\d+|\d+)'):
    abundance = {}
    in_fasta_handler = open(fasta, 'r')
    p = re.compile(regexp)
    for i, [header, _] in enumerate(read_fasta_file_handle(in_fasta_handler)):
        id = header.split()[0].strip()
        m = re.search(p, header)
        if not m:
            logger.fatal("Can't retrieve abundance information:\nfasta_path:%s\nfasta_header:%s" % (fasta, header))
            sys.exit("Cant't retrieve abundance")
        abundance[id] = float(m.group(1))
    in_fasta_handler.close()
    return abundance
