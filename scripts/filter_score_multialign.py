#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Keep each read's best-scoring alignments.

Records are grouped per read, where a read is (QNAME, mate) as set by
``sam_mates.py``: the two mates of a pair are filtered independently and never
compete with each other.  Records of one read must be contiguous.  Within a
read, alignments are ranked by their ``AS:i`` score, and kept until the first
one that falls below the threshold:

* straight mode: ``threshold * best score``
* geometric mode: ``threshold * previous score``
"""

import argparse
import re
import sys

from sam_mates import mate_of

SCORE_PATTERN = re.compile(r'\tAS:i:(-?\d+)')


def alignment_score(line):
    match = SCORE_PATTERN.search(line)
    if match is None:
        raise ValueError('SAM record has no AS:i score: {!r}'.format(line.strip()))
    return int(match.group(1))


def read_groups(lines):
    """Yield lists of SAM lines sharing the same (QNAME, mate), header skipped."""
    key, group = None, []
    for line in lines:
        if not line.strip() or line.startswith('@'):
            continue
        fields = line.split('\t', 2)
        current = (fields[0], mate_of(fields[1]))
        if current != key:
            if group:
                yield group
            key, group = current, []
        group.append(line if line.endswith('\n') else line + '\n')
    if group:
        yield group


def filter_group(lines, threshold, geometric=False):
    """Return the lines of one read to keep, best score first."""
    ranked = sorted(lines, key=alignment_score, reverse=True)
    kept = ranked[:1]
    reference = alignment_score(ranked[0])
    for line in ranked[1:]:
        score = alignment_score(line)
        if score < reference * threshold:
            break
        kept.append(line)
        if geometric:
            reference = score
    return kept


def filter_alignments(lines, threshold, geometric=False):
    for group in read_groups(lines):
        yield from filter_group(group, threshold, geometric)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('-i', '--input_sam', metavar='input.sam',
                        type=argparse.FileType('r'), default='-',
                        help='input sam file')
    parser.add_argument('-o', '--output_sam', metavar='output.sam',
                        type=argparse.FileType('w'), default='-',
                        help='output sam file')
    parser.add_argument('--geometric', action='store_true',
                        help='Uses an adaptative geometric threshold (straight threshold by default)')
    parser.add_argument('-t', '--threshold', metavar='PCT',
                        type=float, default=0.9,
                        help='Threshold will be either threshold*max score (straight mode) '
                             'or threshold*previous score (geometric mode)')
    args = parser.parse_args(argv)
    args.output_sam.writelines(filter_alignments(args.input_sam, args.threshold, args.geometric))
    return 0


if __name__ == '__main__':
    sys.exit(main())
