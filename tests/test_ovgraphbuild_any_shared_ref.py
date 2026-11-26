"""ovgraphbuild --any-shared-ref: a read pair gets an edge when it is compatible
on any reference both reads are aligned to.

By default the first alignment pair on a shared reference (read i's
alignments in SAM order, then read j's) decides the pair: if the reads do not
overlap there, they are incompatible and later shared references are never
looked at. --any-shared-ref keeps looking until one is compatible. Every edge
still has to pass the de novo read-to-read alignment, so the option can only
add edges, and the edges both modes make are written identically.
"""
import os
import random
import subprocess

import pytest

from test_ovgraphbuild_pair_pruning import EDGES_HEADER, make_inputs, read_bytes, run

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OVGRAPHBUILD = os.path.join(ROOT_DIR, 'ovgraphbuild', 'bin', 'ovgraphbuild')

pytestmark = pytest.mark.skipif(
    not os.path.isfile(OVGRAPHBUILD) or not os.access(OVGRAPHBUILD, os.X_OK),
    reason='ovgraphbuild binary is not built'
)

BASES = 'ACGT'


def edges(prefix):
    with open(str(prefix) + '.edges.csv') as handle:
        lines = handle.readlines()
    assert lines[0] == EDGES_HEADER
    return lines[1:]


def write_pair_case(tmp_path, first_ref_for_both):
    """Two reads that overlap by 60 bases on refB. On refA, listed first,
    read2 is aligned 300 bases further on (e.g. a partial or misplaced hit)."""
    rng = random.Random(7)
    ref_a = ''.join(rng.choice(BASES) for _ in range(600))
    ref_b = ''.join(rng.choice(BASES) for _ in range(600))
    reference = tmp_path / 'refs.fasta'
    reference.write_text('>refA x\n{}\n>refB x\n{}\n'.format(ref_a, ref_b))
    read1, read2 = ref_b[0:100], ref_b[40:140]
    order = ['refA', 'refB'] if first_ref_for_both == 'refA' else ['refB', 'refA']
    positions = {'read1': {'refA': 1, 'refB': 1}, 'read2': {'refA': 301, 'refB': 41}}
    sam = tmp_path / 'reads.sam'
    with sam.open('w') as handle:
        handle.write('@HD\tVN:1.0\tSO:queryname\n@SQ\tSN:refA\tLN:600\n@SQ\tSN:refB\tLN:600\n')
        for name, seq in (('read1', read1), ('read2', read2)):
            for ref in order:
                handle.write('\t'.join([name, '0', ref, str(positions[name][ref]), '255', '100M',
                                        '*', '0', '0', seq, '*', 'AS:i:100']) + '\n')
    return reference, sam


def test_an_overlap_on_a_later_shared_reference_needs_the_option(tmp_path):
    reference, sam = write_pair_case(tmp_path, first_ref_for_both='refA')
    run(reference, sam, tmp_path / 'first')
    run(reference, sam, tmp_path / 'any', '--any-shared-ref')
    assert edges(tmp_path / 'first') == []
    assert len(edges(tmp_path / 'any')) == 1


def test_the_default_verdict_depends_on_the_alignment_order(tmp_path):
    reference, sam = write_pair_case(tmp_path, first_ref_for_both='refB')
    run(reference, sam, tmp_path / 'first')
    run(reference, sam, tmp_path / 'any', '--any-shared-ref')
    assert len(edges(tmp_path / 'first')) == 1
    assert edges(tmp_path / 'first') == edges(tmp_path / 'any')


def test_no_compatible_shared_reference_gives_no_edge(tmp_path):
    reference, sam = write_pair_case(tmp_path, first_ref_for_both='refA')
    text = sam.read_text().replace('refB\t41\t', 'refB\t401\t')  # no overlap on refB either
    sam.write_text(text)
    first = run(reference, sam, tmp_path / 'first')
    run(reference, sam, tmp_path / 'any', '--any-shared-ref')
    assert edges(tmp_path / 'first') == edges(tmp_path / 'any') == []
    assert '1 were found incompatible' in first.stdout


def misplace_first_alignments(sam, seed, rate=0.3):
    """Move the first alignment of some reads 200 bases along its reference,
    so that the first shared reference of a pair can disagree with the others."""
    rng = random.Random(seed)
    lines = sam.read_text().splitlines(keepends=True)
    header = [l for l in lines if l.startswith('@')]
    records = [l.split('\t') for l in lines if not l.startswith('@')]
    seen, moved = set(), 0
    for record in records:
        if record[0] in seen:
            continue
        seen.add(record[0])
        if rng.random() < rate:
            record[3] = str((int(record[3]) + 200 - 1) % 380 + 1)
            moved += 1
    sam.write_text(''.join(header) + ''.join('\t'.join(r) for r in records))
    return moved


@pytest.mark.parametrize('seed', [1, 2, 3])
@pytest.mark.parametrize('id_threshold', ['1', '0.97'])
def test_option_only_adds_edges(tmp_path, seed, id_threshold):
    reference, sam, _ = make_inputs(tmp_path, seed)
    assert misplace_first_alignments(sam, seed) > 0
    run(reference, sam, tmp_path / 'first', '--id_threshold', id_threshold)
    run(reference, sam, tmp_path / 'any', '--any-shared-ref', '--id_threshold', id_threshold)

    first, any_ref = edges(tmp_path / 'first'), edges(tmp_path / 'any')
    assert set(first) <= set(any_ref)
    assert len(any_ref) > len(first) > 20
    # Same order: the edges both modes make keep their relative positions.
    assert [e for e in any_ref if e in set(first)] == first
    assert read_bytes(str(tmp_path / 'first') + '.nodes.csv') == \
        read_bytes(str(tmp_path / 'any') + '.nodes.csv')


def test_option_without_disagreeing_references_changes_nothing(tmp_path):
    # References of a family keep the read coordinates: every shared reference agrees.
    reference, sam, _ = make_inputs(tmp_path, seed=6)
    run(reference, sam, tmp_path / 'first')
    run(reference, sam, tmp_path / 'any', '--any-shared-ref')
    assert edges(tmp_path / 'any') == edges(tmp_path / 'first')


@pytest.mark.parametrize('threads', [3, 8])
def test_option_with_threads_and_shards_matches_serial(tmp_path, threads):
    reference, sam, _ = make_inputs(tmp_path, seed=8)
    misplace_first_alignments(sam, 8)
    run(reference, sam, tmp_path / 'serial', '--any-shared-ref', '--threads', '1')
    run(reference, sam, tmp_path / 'threads', '--any-shared-ref', '--threads', str(threads))
    assert read_bytes(str(tmp_path / 'threads') + '.edges.csv') == \
        read_bytes(str(tmp_path / 'serial') + '.edges.csv')

    merged = []
    for shard in range(3):
        prefix = tmp_path / 'shard-{}'.format(shard)
        run(reference, sam, prefix, '--any-shared-ref', '--threads', str(threads),
            '--pair-shard-index', str(shard), '--pair-shard-count', '3')
        merged.extend(edges(prefix))
        with open(str(prefix) + '.pair-shard.tsv') as handle:
            assert 'any_shared_ref\t1\n' in handle.read()
    assert merged == edges(tmp_path / 'serial')
