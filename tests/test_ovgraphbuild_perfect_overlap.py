"""ovgraphbuild skips the full re-alignment of reads that overlap perfectly.

Each pair of compatible reads is re-aligned de novo (an overlap alignment with
free end gaps) to measure the overlap. When the best alignment score is
reached at a single end cell and equals the length of the diagonal there, the
alignment SeqAn returns can only be that ungapped all-match diagonal, so its
statistics are filled in directly; any other case runs the SeqAn alignment
(see ovgraphbuild/src/perfectOverlap.h). ``--full-realignment`` aligns every
pair with SeqAn as before. Both modes must give byte-identical graphs and the
same statistics.
"""
import os
import random
import re

import pytest

from test_ovgraphbuild_pair_pruning import (BASES, OVGRAPHBUILD, make_inputs, read_bytes,
                                            run, statistics)


pytestmark = pytest.mark.skipif(
    not os.path.isfile(OVGRAPHBUILD) or not os.access(OVGRAPHBUILD, os.X_OK),
    reason='ovgraphbuild binary is not built'
)

REALIGNED = re.compile(r'INFO: Re-aligned read pairs: +(\d+) \((\d+) perfect overlaps')


def realignments(stdout):
    """(re-aligned pairs, of which perfect overlaps that skipped SeqAn)."""
    realigned, perfect = REALIGNED.search(stdout).groups()
    return int(realigned), int(perfect)


def other_statistics(stdout):
    return [line for line in statistics(stdout) if not REALIGNED.match(line)]


def low_complexity_reference(rng, length):
    """Tandem repeats, homopolymers and copied blocks: many equal-score overlaps."""
    parts = []
    while sum(map(len, parts)) < length:
        kind = rng.randrange(4)
        if kind == 0:
            unit = ''.join(rng.choice(BASES) for _ in range(rng.randint(1, 6)))
            parts.append(unit * rng.randint(4, 20))
        elif kind == 1:
            parts.append(rng.choice(BASES) * rng.randint(5, 25))
        elif kind == 2 and parts:
            parts.append(rng.choice(parts))
        else:
            parts.append(''.join(rng.choice(BASES) for _ in range(rng.randint(10, 40))))
    return ''.join(parts)[:length]


def make_tricky_inputs(tmp_path, seed, reads_num=180):
    """Reads of 60-150 bases on low-complexity and random references.

    Besides perfect overlaps they give enclosed reads, duplicates, ties between
    several alignments, substitutions, ambiguous bases and 1-base indels.
    """
    rng = random.Random(seed)
    ref_len = 500
    refs = [''.join(rng.choice(BASES) for _ in range(ref_len))]
    refs += [low_complexity_reference(rng, ref_len) for _ in range(2)]

    reference = tmp_path / 'refs.fasta'
    reference.write_text(''.join('>ref{}\n{}\n'.format(i, s) for i, s in enumerate(refs)))

    lines = []
    previous = None
    for read in range(reads_num):
        if previous and rng.random() < 0.1:
            ref_id, pos, cigar, seq = previous   # duplicate read
        else:
            ref_id = rng.randrange(len(refs))
            read_len = rng.randint(60, 150)
            pos = rng.randrange(ref_len - read_len - 2)
            seq = refs[ref_id][pos:pos + read_len]
            cigar = '{}M'.format(read_len)
            change = rng.random()
            if change < 0.12:
                k = rng.randrange(read_len)
                seq = seq[:k] + rng.choice(BASES.replace(seq[k], '')) + seq[k + 1:]
            elif change < 0.2:
                k = rng.randrange(read_len)
                seq = seq[:k] + 'N' + seq[k + 1:]
            elif change < 0.26:
                seq = seq[:40] + rng.choice(BASES) + seq[40:]
                cigar = '40M1I{}M'.format(read_len - 40)
            elif change < 0.32:
                seq = seq[:40] + seq[41:]
                cigar = '40M1D{}M'.format(read_len - 41)
        previous = (ref_id, pos, cigar, seq)
        flag = 16 if rng.random() < 0.5 else 0
        lines.append('\t'.join(['read{:04d}'.format(read), str(flag), 'ref{}'.format(ref_id),
                                str(pos + 1), '255', cigar, '*', '0', '0', seq, '*']))

    sam = tmp_path / 'reads.sam'
    with sam.open('w') as handle:
        handle.write('@HD\tVN:1.0\n')
        handle.write(''.join('@SQ\tSN:ref{}\tLN:{}\n'.format(i, len(s)) for i, s in enumerate(refs)))
        handle.write('\n'.join(lines) + '\n')
    return reference, sam


def inputs(tmp_path, kind, seed):
    if kind == 'tricky':
        return make_tricky_inputs(tmp_path, seed)
    reference, sam, _ = make_inputs(tmp_path, seed)
    return reference, sam


def assert_same_graph(prefix, other_prefix):
    for suffix in ('.edges.csv', '.nodes.csv', '.asqg'):
        assert read_bytes(str(prefix) + suffix) == read_bytes(str(other_prefix) + suffix), suffix


@pytest.mark.parametrize('kind,seed', [('tricky', 1), ('tricky', 2), ('tricky', 3),
                                       ('families', 1), ('families', 2)])
@pytest.mark.parametrize('id_threshold', ['1', '0.97'])
def test_perfect_overlaps_match_full_realignment(tmp_path, kind, seed, id_threshold):
    reference, sam = inputs(tmp_path, kind, seed)
    full = run(reference, sam, tmp_path / 'full', '--asqg', '--full-realignment',
               '--id_threshold', id_threshold)
    fast = run(reference, sam, tmp_path / 'fast', '--asqg', '--id_threshold', id_threshold)

    assert_same_graph(tmp_path / 'fast', tmp_path / 'full')
    assert other_statistics(fast.stdout) == other_statistics(full.stdout)

    realigned, perfect = realignments(fast.stdout)
    assert realignments(full.stdout) == (realigned, 0)
    assert 0 < perfect <= realigned
    assert read_bytes(str(tmp_path / 'fast') + '.edges.csv').count(b'\n') > 20
    if kind == 'tricky':
        # Low-complexity reads must also take the SeqAn route (ties between
        # alignments), or the comparison proves little.
        assert perfect < realigned
        enclosed = int(re.search(r'enclosed overlaps: +(\d+)', fast.stdout).group(1))
        assert enclosed > 0


@pytest.mark.parametrize('seed', range(1, 9))
@pytest.mark.parametrize('threads', ['1', '4'])
def test_each_perfect_overlap_matches_seqan(tmp_path, seed, threads):
    """--check-perfect-overlaps aligns every shortcut pair with SeqAn too.

    This checks each pair, ties between alignments included, not only the
    ones whose difference would change the graph.
    """
    reference, sam = make_tricky_inputs(tmp_path, seed)
    checked = run(reference, sam, tmp_path / 'checked', '--check-perfect-overlaps',
                  '--threads', threads, '--id_threshold', '0.9')
    assert realignments(checked.stdout)[1] > 0


@pytest.mark.parametrize('threads', [3, 24])
def test_perfect_overlap_threads_and_shards(tmp_path, threads):
    reference, sam = make_tricky_inputs(tmp_path, seed=4)
    full = run(reference, sam, tmp_path / 'full', '--full-realignment', '--threads', '1')
    fast = run(reference, sam, tmp_path / 'fast', '--threads', str(threads))
    assert read_bytes(str(tmp_path / 'fast') + '.edges.csv') == \
        read_bytes(str(tmp_path / 'full') + '.edges.csv')
    realigned, perfect = realignments(fast.stdout)

    merged, metadata = [], []
    for shard in range(3):
        prefix = tmp_path / 'shard-{}'.format(shard)
        run(reference, sam, prefix, '--threads', str(threads),
            '--pair-shard-index', str(shard), '--pair-shard-count', '3')
        with open(str(prefix) + '.edges.csv') as handle:
            merged.extend(handle.readlines()[1:])
        with open(str(prefix) + '.pair-shard.tsv') as handle:
            metadata.append(dict(line.rstrip('\n').split('\t') for line in handle))

    with open(str(tmp_path / 'full') + '.edges.csv') as handle:
        assert merged == handle.readlines()[1:]
    assert sum(int(m['realigned_pair_count']) for m in metadata) == realigned
    assert sum(int(m['perfect_overlap_shortcut_count']) for m in metadata) == perfect
    assert all(m['full_realignment'] == '0' for m in metadata)


def test_dovetail_and_enclosed_overlaps(tmp_path):
    """Exact overlap coordinates of a perfect dovetail, as SeqAn gives them."""
    rng = random.Random(21)
    ref = ''.join(rng.choice(BASES) for _ in range(400))
    reference = tmp_path / 'refs.fasta'
    reference.write_text('>ref0\n{}\n'.format(ref))
    # read0 [1, 120], read1 [61, 180], read2 [71, 150]: read0 overlaps read1 by 60 and read2
    # by 50, and read2 lies inside read1.
    layout = [('read0', 1, 120, 0), ('read1', 61, 120, 0), ('read2', 71, 80, 0)]
    sam = tmp_path / 'reads.sam'
    with sam.open('w') as handle:
        handle.write('@HD\tVN:1.0\n@SQ\tSN:ref0\tLN:400\n')
        for name, pos, length, flag in layout:
            handle.write('\t'.join([name, str(flag), 'ref0', str(pos), '255', '{}M'.format(length),
                                    '*', '0', '0', ref[pos - 1:pos - 1 + length], '*']) + '\n')

    full = run(reference, sam, tmp_path / 'full', '--asqg', '--full-realignment')
    fast = run(reference, sam, tmp_path / 'fast', '--asqg')
    assert_same_graph(tmp_path / 'fast', tmp_path / 'full')
    assert realignments(fast.stdout) == (3, 3)
    assert realignments(full.stdout) == (3, 0)
    with open(str(tmp_path / 'fast') + '.asqg') as handle:
        edges = [line.split('\t')[1] for line in handle if line.startswith('ED')]
    assert edges == ['0 1 60 119 120 0 59 120 0 -1\n', '0 2 70 119 120 0 49 80 0 -1\n']
    # An enclosed overlap is counted but gives no edge (its identity is not measured).
    assert 'enclosed overlaps:  1\n' in fast.stdout


@pytest.mark.parametrize('x_len,y_len,shortcut', [(60, 60, 0), (55, 65, 1), (65, 55, 1)])
def test_tied_overlaps_take_the_seqan_route(tmp_path, x_len, y_len, shortcut):
    """Reads X+Y and Y+X overlap perfectly on Y and, the other way round, on X.

    With |X| == |Y| both overlaps score the same: SeqAn's tie-breaking decides,
    so the pair must not take the shortcut. With different lengths the longer
    overlap is the single best one.
    """
    rng = random.Random(x_len * 100 + y_len)
    x, y = (''.join(rng.choice(BASES) for _ in range(n)) for n in (x_len, y_len))
    ref = x + y + x
    reference = tmp_path / 'refs.fasta'
    reference.write_text('>ref0\n{}\n'.format(ref))
    sam = tmp_path / 'reads.sam'
    read_len = x_len + y_len
    with sam.open('w') as handle:
        handle.write('@HD\tVN:1.0\n@SQ\tSN:ref0\tLN:{}\n'.format(len(ref)))
        for name, pos in (('read0', 0), ('read1', x_len)):
            handle.write('\t'.join([name, '0', 'ref0', str(pos + 1), '255', '{}M'.format(read_len),
                                    '*', '0', '0', ref[pos:pos + read_len], '*']) + '\n')

    full = run(reference, sam, tmp_path / 'full', '--asqg', '--full-realignment')
    fast = run(reference, sam, tmp_path / 'fast', '--asqg', '--check-perfect-overlaps')
    assert_same_graph(tmp_path / 'fast', tmp_path / 'full')
    assert realignments(fast.stdout) == (1, shortcut)
