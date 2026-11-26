#!/usr/bin/env python3
"""FASTA reading and writing shared by the MATAM scripts."""


def read_fasta_file_handle(fasta_file_handle):
    """
    Parse a fasta file and yield (header, sequence) tuples.

    Each line is stripped of surrounding whitespace and blank lines are
    skipped. The header is the rest of the '>' line; the sequence is its
    following lines joined. An empty file yields nothing. Text before the
    first header is an error (ValueError). The file handle is closed at the end.
    """
    header = None
    seqlines = list()
    for line in (l.strip() for l in fasta_file_handle):
        if not line:
            continue
        if line[0] == '>':
            if header is not None:
                yield (header, ''.join(seqlines))
                seqlines = list()
            header = line[1:]
        elif header is None:
            raise ValueError('FASTA sequence before the first header: {!r}'.format(line[:60]))
        else:
            seqlines.append(line)
    if header is not None:
        yield (header, ''.join(seqlines))
    fasta_file_handle.close()


def format_seq(seq, linereturn=80):
    """
    Format an input sequence in lines of linereturn characters (no final newline)
    """
    buff = list()
    for i in range(0, len(seq), linereturn):
        buff.append("{0}\n".format(seq[i:(i + linereturn)]))
    return ''.join(buff).rstrip()
