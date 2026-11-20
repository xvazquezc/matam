#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Take the read mapping of a MATAM run from SortMeRNA alignments made earlier.

``matam_assembly.py --read_alignments SAM`` skips the SortMeRNA read mapping
and uses alignments of the input reads made by an earlier run instead, e.g.
the ``marker_alignments.sam.gz`` of a ``--filter_only`` run on the same or a
larger read set. Assembling several subsamples of one read set then needs one
mapping, not one per subsample.

This script writes what the read mapping would have written, under the same
names (``<aligned>_fwd.fq.gz``, ``<aligned>_rev.fq.gz``, ``<aligned>.sam.gz``):

* the input pairs with at least one alignment (as SortMeRNA ``--paired_in``),
  in input order;
* the alignments of those reads, in the order of the SAM file (SortMeRNA
  writes them in input order, which the mate tagging relies on).

Alignments of reads that are not in the input are left out. A read that is in
the input but has no alignment is taken as unaligned, so the SAM file must
come from a mapping that included every input read.

The alignments are not exactly those of a mapping of the input reads alone:
SortMeRNA's minimum alignment score grows with the total length of the reads
mapped together, and ``--num_alignments`` keeps the first alignments that
reach it. A read near the threshold, or with more candidate references than
``--num_alignments``, can therefore get other alignments.
"""

import argparse
import gzip
import shlex
import sys

PROGRAM = 'sortmerna'


class ReadAlignmentsError(ValueError):
    """The SAM file cannot stand in for this run's read mapping."""


def open_text(path, mode='rt'):
    path = str(path)
    if mode.startswith('w') and path.endswith('.gz'):
        return gzip.open(path, mode, compresslevel=1)
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)


def sortmerna_command(sam_path):
    """The SortMeRNA command line from the SAM header (@PG CL), or None."""
    with open_text(sam_path) as handle:
        for line in handle:
            if not line.startswith('@'):
                break
            if not line.startswith('@PG'):
                continue
            fields = dict(field.split(':', 1) for field in line.rstrip('\n').split('\t')[1:]
                          if ':' in field)
            if fields.get('ID') == PROGRAM and 'CL' in fields:
                return fields['CL']
    return None


def _option_values(tokens, option):
    return [tokens[i + 1] for i, token in enumerate(tokens[:-1]) if token == option]


def check_header(sam_path, ref_db_filename, num_alignments, evalue):
    """Check that the SAM file was made against this database with these options.

    The mapping options that change the alignments are compared with the
    @PG line SortMeRNA writes: the reference file name (not its directory, so
    a moved database still matches), --num_alignments and -e.
    """
    command = sortmerna_command(sam_path)
    if command is None:
        raise ReadAlignmentsError(
            '{} has no SortMeRNA @PG header line, so it cannot be checked that it was '
            'made against the same database with the same options'.format(sam_path))
    tokens = shlex.split(command)
    problems = []
    refs = [ref.split(',')[0].rsplit('/', 1)[-1] for ref in _option_values(tokens, '--ref')]
    if refs != [ref_db_filename]:
        problems.append('reference {} instead of {}'.format(
            ', '.join(refs) or 'unknown', ref_db_filename))
    found = _option_values(tokens, '--num_alignments')
    if found != [str(num_alignments)]:
        problems.append('--num_alignments {} instead of {}'.format(
            ', '.join(found) or 'unset', num_alignments))
    found = _option_values(tokens, '-e')
    try:
        same_evalue = len(found) == 1 and float(found[0]) == float('{:.2e}'.format(evalue))
    except ValueError:
        same_evalue = False
    if not same_evalue:
        problems.append('-e {} instead of {:.2e}'.format(', '.join(found) or 'unset', evalue))
    if problems:
        raise ReadAlignmentsError('{} was not made with this run\'s read mapping: {}'.format(
            sam_path, '; '.join(problems)))


def _iter_fastq_records(path):
    with open_text(path) as handle:
        while True:
            record = [handle.readline() for _ in range(4)]
            if not record[0]:
                if any(record[1:]):
                    raise ReadAlignmentsError('Truncated FASTQ record in {}'.format(path))
                return
            if (not all(record) or not record[0].startswith('@')
                    or not record[2].startswith('+')):
                raise ReadAlignmentsError(
                    'Malformed or truncated FASTQ record in {}: {!r}'.format(
                        path, record[0].strip()))
            yield record[0][1:].split()[0], record


def _iter_pairs(forward_path, reverse_path):
    forward = _iter_fastq_records(forward_path)
    reverse = _iter_fastq_records(reverse_path)
    for forward_record in forward:
        reverse_record = next(reverse, None)
        if reverse_record is None:
            raise ReadAlignmentsError('{} has more reads than {}'.format(
                forward_path, reverse_path))
        yield forward_record, reverse_record
    if next(reverse, None) is not None:
        raise ReadAlignmentsError('{} has more reads than {}'.format(
            reverse_path, forward_path))


def _aligned_read_ids(sam_path):
    with open_text(sam_path) as handle:
        return {line.split('\t', 1)[0] for line in handle if not line.startswith('@')}


def select(sam_path, forward_path, reverse_path, aligned_base):
    """Write the aligned input pairs and their alignments; return (pairs, selected, records)."""
    aligned_ids = _aligned_read_ids(sam_path)

    selected_ids = set()
    pairs = selected = 0
    with open_text(aligned_base + '_fwd.fq.gz', 'wt') as forward_out, \
            open_text(aligned_base + '_rev.fq.gz', 'wt') as reverse_out:
        for (forward_id, forward_record), (reverse_id, reverse_record) in _iter_pairs(
                forward_path, reverse_path):
            pairs += 1
            if forward_id in aligned_ids or reverse_id in aligned_ids:
                forward_out.writelines(forward_record)
                reverse_out.writelines(reverse_record)
                selected_ids.update((forward_id, reverse_id))
                selected += 1
    del aligned_ids

    records = 0
    with open_text(sam_path) as sam, open_text(aligned_base + '.sam.gz', 'wt') as sam_out:
        for line in sam:
            if line.startswith('@'):
                sam_out.write(line)
            elif line.split('\t', 1)[0] in selected_ids:
                sam_out.write(line)
                records += 1
    return pairs, selected, records


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--sam', required=True,
                        help='SortMeRNA SAM alignments (optionally gzipped) of a read set '
                             'that includes the input reads')
    parser.add_argument('-1', '--forward', required=True, help='R1 FASTQ (optionally gzipped)')
    parser.add_argument('-2', '--reverse', required=True, help='R2 FASTQ (optionally gzipped)')
    parser.add_argument('--aligned', required=True,
                        help='Output prefix, as SortMeRNA --aligned')
    args = parser.parse_args(argv)
    try:
        pairs, selected, records = select(args.sam, args.forward, args.reverse, args.aligned)
    except ReadAlignmentsError as error:
        sys.exit('ERROR: {}'.format(error))
    print('Selected {} of {} pairs with {} alignments from {}'.format(
        selected, pairs, records, args.sam))


if __name__ == '__main__':
    main()
