#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recover mate identity in SortMeRNA SAM output.

SortMeRNA aligns each mate independently and writes FLAG 0/16 only: both
mates of a CASAVA-named pair share one QNAME, and nothing in the SAM tells
R1 from R2.  MATAM restores the standard SAM pairing flags right after
mapping, so every later step can key reads on (QNAME, mate) without renaming
the user's reads:

* 0x1  read is paired
* 0x40 first mate (R1) / 0x80 second mate (R2)
* 0x8  the other mate has no alignment

Each SAM record is assigned to a mate by comparing its SEQ (reverse
complemented when 0x10 is set) with the R1 and R2 sequences of its pair.  When
the mates have different IDs (legacy ``/1``, ``/2`` naming) the QNAME decides.
SAM records and read pairs are streamed in lockstep: SortMeRNA writes
alignments in input order, so memory use stays constant.  Any inconsistency is
a hard error instead of a silent misassignment.
"""

import argparse
import gzip
import sys

PAIRED = 0x1
MATE_UNMAPPED = 0x8
REVERSE = 0x10
FIRST_MATE = 0x40
SECOND_MATE = 0x80
SECONDARY = 0x100
SUPPLEMENTARY = 0x800
MATE_FLAGS = {1: FIRST_MATE, 2: SECOND_MATE}

_COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


class MateTaggingError(ValueError):
    """The SAM records cannot be matched unambiguously to the read pairs."""


def mate_of(flag):
    """Return 1 for R1, 2 for R2 and 0 for an unpaired record."""
    flag = int(flag)
    first, second = bool(flag & FIRST_MATE), bool(flag & SECOND_MATE)
    if first == second:
        return 0
    return 1 if first else 2


def reverse_complement(sequence):
    return sequence.translate(_COMPLEMENT)[::-1]


def open_text(path, mode='rt'):
    path = str(path)
    if path == '-':
        return sys.stdin if 'r' in mode else sys.stdout
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)


def iter_fastq(handle, source='FASTQ'):
    """Yield (read_id, sequence) from a 4-line FASTQ, validating its structure."""
    while True:
        header = handle.readline()
        if not header:
            return
        sequence, separator, quality = handle.readline(), handle.readline(), handle.readline()
        if not quality or not header.startswith('@') or not separator.startswith('+'):
            raise MateTaggingError('Malformed or truncated FASTQ record in {}: {!r}'.format(
                source, header.strip()))
        yield header[1:].split()[0], sequence.strip()


def iter_separate_pairs(forward_path, reverse_path):
    """Yield ((id1, seq1), (id2, seq2)) from R1/R2 FASTQ files read in step."""
    with open_text(forward_path) as forward, open_text(reverse_path) as reverse:
        forward_records = iter_fastq(forward, forward_path)
        reverse_records = iter_fastq(reverse, reverse_path)
        for first in forward_records:
            second = next(reverse_records, None)
            if second is None:
                raise MateTaggingError('Paired FASTQ files have different record counts')
            yield first, second
        if next(reverse_records, None) is not None:
            raise MateTaggingError('Paired FASTQ files have different record counts')


def _record_groups(sam_lines, header_out):
    """Group contiguous SAM records by QNAME; header lines go to header_out."""
    qname, group = None, []
    for line in sam_lines:
        if not line.strip():
            continue
        if line.startswith('@'):
            header_out.append(line if line.endswith('\n') else line + '\n')
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 11:
            raise MateTaggingError('Truncated SAM record: {!r}'.format(line.strip()))
        flag = int(fields[1])
        if flag & (SECONDARY | SUPPLEMENTARY):
            raise MateTaggingError(
                'Unexpected secondary/supplementary flag {} for {}'.format(flag, fields[0]))
        if fields[0] != qname:
            if group:
                yield qname, group
            qname, group = fields[0], []
        group.append(fields)
    if group:
        yield qname, group


def _sequence_matches(fields, read_sequence):
    sequence = fields[9]
    if sequence == '*':
        return False
    oriented = reverse_complement(sequence) if int(fields[1]) & REVERSE else sequence
    oriented, read_sequence = oriented.upper(), read_sequence.upper()
    if 'H' in fields[5]:  # hard-clipped records only carry part of the read
        return oriented in read_sequence
    return oriented == read_sequence


def _assign(fields, pair):
    (id1, seq1), (id2, seq2) = pair
    qname = fields[0]
    if id1 != id2:
        mate = 1 if qname == id1 else 2
        if not _sequence_matches(fields, seq1 if mate == 1 else seq2):
            raise MateTaggingError('SAM SEQ of {} does not match its read'.format(qname))
        return mate
    first, second = _sequence_matches(fields, seq1), _sequence_matches(fields, seq2)
    if first and second:
        return None  # identical mates: resolved by record order
    if first:
        return 1
    if second:
        return 2
    raise MateTaggingError('SAM record of {} matches neither mate (FLAG {}, {}:{})'.format(
        qname, fields[1], fields[2], fields[3]))


def _resolve_identical_mates(records):
    """Assign records of identical mates by order: R1's alignments come first.

    The second mate starts at the first placement (RNAME, POS, strand) that
    repeats one already given to the first mate.
    """
    seen_first, current = set(), 1
    for record in records:
        if record[1] is not None:
            continue
        fields = record[0]
        placement = (fields[2], fields[3], int(fields[1]) & REVERSE)
        if current == 1 and placement in seen_first:
            current = 2
        if current == 1:
            seen_first.add(placement)
        record[1] = current


def _emit_pair(records):
    _resolve_identical_mates(records)
    present = {mate for _, mate in records}
    for wanted in (1, 2):
        for fields, mate in records:
            if mate != wanted:
                continue
            flag = int(fields[1])
            existing = mate_of(flag)
            if existing and existing != mate:
                raise MateTaggingError('Conflicting mate flag {} for {}'.format(flag, fields[0]))
            flag |= PAIRED | MATE_FLAGS[mate]
            if (3 - mate) not in present:
                flag |= MATE_UNMAPPED
            yield '\t'.join([fields[0], str(flag)] + fields[2:]) + '\n'


def tag_mates(sam_lines, pairs):
    """Yield SAM lines with pairing flags set, records ordered R1 then R2 per pair.

    ``pairs`` yields ((id1, seq1), (id2, seq2)) in the order the reads were
    given to SortMeRNA.  Pairs without any alignment are skipped.
    """
    header = []
    groups = _record_groups(sam_lines, header)
    first_group = next(groups, None)
    for line in header:
        yield line
    yield '@PG\tID:matam_sam_mates\tPN:sam_mates.py\n'
    if first_group is None:
        return

    pairs = iter(pairs)
    pair, records = None, []
    for qname, group in _chain_first(first_group, groups):
        while pair is None or qname not in (pair[0][0], pair[1][0]):
            if records:
                yield from _emit_pair(records)
                records = []
            pair = next(pairs, None)
            if pair is None:
                raise MateTaggingError(
                    'SAM record {} not found in the read pairs (missing or out of order)'.format(qname))
        records.extend([fields, _assign(fields, pair)] for fields in group)
    if records:
        yield from _emit_pair(records)


def _chain_first(first, rest):
    yield first
    yield from rest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('-s', '--sam', required=True,
                        help='SortMeRNA SAM (plain or .gz, "-" for stdin)')
    parser.add_argument('-1', '--forward', required=True,
                        help='R1 FASTQ, in the order given to SortMeRNA')
    parser.add_argument('-2', '--reverse', required=True,
                        help='R2 FASTQ, in the order given to SortMeRNA')
    parser.add_argument('-o', '--output', default='-', help='Output SAM (default: stdout)')
    args = parser.parse_args(argv)

    pairs = iter_separate_pairs(args.forward, args.reverse)

    sam_in = open_text(args.sam)
    sam_out = open_text(args.output, 'wt')
    try:
        sam_out.writelines(tag_mates(sam_in, pairs))
    except MateTaggingError as error:
        sys.stderr.write('ERROR: {}\n'.format(error))
        return 1
    finally:
        if sam_in is not sys.stdin:
            sam_in.close()
        if sam_out is not sys.stdout:
            sam_out.close()
    return 0


if __name__ == '__main__':
    sys.exit(main())
