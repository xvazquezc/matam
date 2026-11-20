"""ovgraphbuild compares only the read pairs that share a reference.

A pair of reads without a common reference can never give an edge: it is
counted as "neither compatible nor incompatible" (and as incompatible) without
any alignment comparison. By default ovgraphbuild skips those pairs and adds
them to the statistics in bulk; ``--all-pairs`` compares every pair as before.
Both modes must give byte-identical graphs and identical statistics.
"""
import itertools
import os
import random
import re
import subprocess

import pytest


ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OVGRAPHBUILD = os.path.join(ROOT_DIR, 'ovgraphbuild', 'bin', 'ovgraphbuild')

pytestmark = pytest.mark.skipif(
    not os.path.isfile(OVGRAPHBUILD) or not os.access(OVGRAPHBUILD, os.X_OK),
    reason='ovgraphbuild binary is not built'
)

EDGES_HEADER = 'Source;Target;Type;Coemitted;Weight;Pid;MultiRef\n'
BASES = 'ACGT'


def mutate(rng, seq, rate):
    return ''.join(rng.choice(BASES.replace(b, '')) if rng.random() < rate else b for b in seq)


def make_inputs(tmp_path, seed, reads_num=260):
    """Random reference families and reads aligned to 1-3 related references.

    References of a family differ by substitutions only, so a read keeps its
    coordinates on each of them. Reads carry sequencing errors, soft clips,
    1-base insertions/deletions and both strands; a few are also aligned to
    another family, so that reference sharing is not only within a family.
    """
    rng = random.Random(seed)
    ref_len, read_len = 600, 100
    families = []
    for _ in range(5):
        base = ''.join(rng.choice(BASES) for _ in range(ref_len))
        families.append([mutate(rng, base, 0.02) for _ in range(3)])
    refs = [(f * 3 + v, seq) for f, fam in enumerate(families) for v, seq in enumerate(fam)]

    reference = tmp_path / 'refs.fasta'
    with reference.open('w') as handle:
        for ref_id, seq in refs:
            handle.write('>ref{} family\n{}\n'.format(ref_id, seq))

    records = []
    for read in range(reads_num):
        family = rng.randrange(len(families))
        source = rng.choice(families[family])
        if rng.random() < 0.2:
            source = mutate(rng, source, 0.01)  # strain absent from the references
        pos = rng.randrange(ref_len - read_len - 2)
        seq = source[pos:pos + read_len]
        if rng.random() < 0.3:
            seq = mutate(rng, seq, 0.01)
        cigar = '{}M'.format(read_len)
        shape = rng.random()
        if shape < 0.08:
            cigar = '5S{}M'.format(read_len - 5)
        elif shape < 0.14:
            seq = seq[:50] + rng.choice(BASES) + seq[50:]
            cigar = '50M1I{}M'.format(read_len - 50)
        elif shape < 0.20:
            seq = seq[:50] + seq[51:]
            cigar = '50M1D{}M'.format(read_len - 51)
        flag = 16 if rng.random() < 0.5 else 0
        ref_ids = rng.sample(range(family * 3, family * 3 + 3), rng.randint(1, 3))
        if rng.random() < 0.1:
            other = rng.choice([f for f in range(len(families)) if f != family])
            ref_ids.append(other * 3 + rng.randrange(3))
        clip = 5 if cigar.startswith('5S') else 0
        for ref_id in ref_ids:
            records.append(('read{:04d}'.format(read), flag, ref_id, pos + 1 + clip, cigar, seq))

    sam = tmp_path / 'reads.sam'
    with sam.open('w') as handle:
        handle.write('@HD\tVN:1.0\tSO:queryname\n')
        for ref_id, seq in refs:
            handle.write('@SQ\tSN:ref{}\tLN:{}\n'.format(ref_id, len(seq)))
        for qname, flag, ref_id, pos, cigar, seq in records:
            handle.write('\t'.join([qname, str(flag), 'ref{}'.format(ref_id), str(pos), '255',
                                    cigar, '*', '0', '0', seq, '*', 'AS:i:100']) + '\n')

    refs_by_read = {}
    for qname, _, ref_id, *_ in records:
        refs_by_read.setdefault(qname, set()).add(ref_id)
    return reference, sam, [refs_by_read[name] for name in sorted(refs_by_read)]


def run(reference, sam, prefix, *extra):
    command = [OVGRAPHBUILD, '--reference', str(reference), '--sam', str(sam), '--csv',
               '--output_basename', str(prefix), '--min_overlap', '50', '--debug', '-v']
    command.extend(extra)
    return subprocess.run(command, check=True, capture_output=True, text=True)


def statistics(stdout):
    """INFO lines of the run, without the ones that describe the work done.

    Progress (single worker only) is rewritten in place with '\\r' (read as a
    newline here) and the last one isn't followed by a newline, so it is cut
    out before splitting the lines.
    """
    stdout = re.sub(r'INFO: Compatibility graph building progress at \S+% *', '', stdout)
    lines = []
    for line in stdout.splitlines():
        if not line.startswith('INFO:'):
            continue
        if line.startswith('INFO: Compared read pairs'):
            continue
        lines.append(line)
    return lines


def compared_pairs(stdout):
    return int(re.search(r'INFO: Compared read pairs: (\d+)', stdout).group(1))


def common_ref_pairs(refs_by_read):
    return sum(1 for a, b in itertools.combinations(refs_by_read, 2) if a & b)


def read_bytes(path):
    with open(path, 'rb') as handle:
        return handle.read()


@pytest.mark.parametrize('seed', [1, 2, 3])
@pytest.mark.parametrize('id_threshold', ['1', '0.97'])
def test_pruned_graph_and_statistics_match_all_pairs(tmp_path, seed, id_threshold):
    reference, sam, refs_by_read = make_inputs(tmp_path, seed)
    full = run(reference, sam, tmp_path / 'full', '--all-pairs', '--id_threshold', id_threshold)
    pruned = run(reference, sam, tmp_path / 'pruned', '--id_threshold', id_threshold)

    edges = read_bytes(str(tmp_path / 'pruned') + '.edges.csv')
    assert edges == read_bytes(str(tmp_path / 'full') + '.edges.csv')
    assert read_bytes(str(tmp_path / 'pruned') + '.nodes.csv') == \
        read_bytes(str(tmp_path / 'full') + '.nodes.csv')
    # The inputs must exercise every outcome, or the comparison proves little.
    assert edges.count(b'\n') > 20
    full_stats = statistics(full.stdout)
    assert statistics(pruned.stdout) == full_stats
    neither = int(re.search(r'and +(\d+) were neither', full.stdout).group(1))
    incompatible = int(re.search(r'(\d+) were found incompatible', full.stdout).group(1))
    assert 0 < neither < incompatible

    reads_num = len(refs_by_read)
    total = reads_num * (reads_num - 1) // 2
    assert compared_pairs(full.stdout) == total
    assert compared_pairs(pruned.stdout) == common_ref_pairs(refs_by_read)
    assert neither == total - common_ref_pairs(refs_by_read)


@pytest.mark.parametrize('threads', [2, 7, 24])
def test_pruned_threads_match_serial_all_pairs(tmp_path, threads):
    reference, sam, _ = make_inputs(tmp_path, seed=4)
    full = run(reference, sam, tmp_path / 'full', '--all-pairs', '--threads', '1')
    pruned = run(reference, sam, tmp_path / 'pruned', '--threads', str(threads))

    assert read_bytes(str(tmp_path / 'pruned') + '.edges.csv') == \
        read_bytes(str(tmp_path / 'full') + '.edges.csv')
    assert statistics(pruned.stdout) == statistics(full.stdout)
    assert not list(tmp_path.glob('*.tmp'))


@pytest.mark.parametrize('shard_count', [3, 8])
def test_pruned_shards_rebuild_serial_graph(tmp_path, shard_count):
    reference, sam, refs_by_read = make_inputs(tmp_path, seed=5)
    run(reference, sam, tmp_path / 'full', '--all-pairs')

    merged, metadata = [], []
    for shard in range(shard_count):
        prefix = tmp_path / 'shard-{}'.format(shard)
        run(reference, sam, prefix, '--threads', '3',
            '--pair-shard-index', str(shard), '--pair-shard-count', str(shard_count))
        with open(str(prefix) + '.edges.csv') as handle:
            lines = handle.readlines()
        assert lines[0] == EDGES_HEADER
        merged.extend(lines[1:])
        with open(str(prefix) + '.pair-shard.tsv') as handle:
            metadata.append(dict(line.rstrip('\n').split('\t') for line in handle))

    with open(str(tmp_path / 'full') + '.edges.csv') as handle:
        assert merged == handle.readlines()[1:]

    reads_num = len(refs_by_read)
    assert sum(int(m['completed_pair_count']) for m in metadata) == reads_num * (reads_num - 1) // 2
    assert sum(int(m['compared_pair_count']) for m in metadata) == common_ref_pairs(refs_by_read)
    assert all(m['all_pairs'] == '0' for m in metadata)


def test_reads_without_common_reference_are_not_compared(tmp_path):
    reference = tmp_path / 'refs.fasta'
    rng = random.Random(9)
    seqs = [''.join(rng.choice(BASES) for _ in range(300)) for _ in range(3)]
    reference.write_text(''.join('>ref{}\n{}\n'.format(i, s) for i, s in enumerate(seqs)))
    # read0, read1 on ref0 (overlapping); read2 on ref1; read3 on ref1 and ref2.
    layout = [('read0', 0, 1), ('read1', 0, 31), ('read2', 1, 1), ('read3', 1, 41), ('read3', 2, 1)]
    sam = tmp_path / 'reads.sam'
    with sam.open('w') as handle:
        handle.write('@HD\tVN:1.0\n')
        for i, s in enumerate(seqs):
            handle.write('@SQ\tSN:ref{}\tLN:{}\n'.format(i, len(s)))
        for name, ref_id, pos in layout:
            seq = seqs[ref_id][pos - 1:pos + 99]
            handle.write('\t'.join([name, '0', 'ref{}'.format(ref_id), str(pos), '255', '100M',
                                    '*', '0', '0', seq, '*']) + '\n')

    full = run(reference, sam, tmp_path / 'full', '--all-pairs')
    pruned = run(reference, sam, tmp_path / 'pruned')

    # Pairs with a common reference: (read0, read1) and (read2, read3).
    assert compared_pairs(pruned.stdout) == 2
    assert compared_pairs(full.stdout) == 6
    assert statistics(pruned.stdout) == statistics(full.stdout)
    assert '4 were neither compatible nor incompatible' in pruned.stdout
    with open(str(tmp_path / 'pruned') + '.edges.csv') as handle:
        edges = [line.split(';')[:2] for line in handle.readlines()[1:]]
    assert edges == [['0', '1'], ['2', '3']]
    assert read_bytes(str(tmp_path / 'pruned') + '.edges.csv') == \
        read_bytes(str(tmp_path / 'full') + '.edges.csv')


def make_lopsided_inputs(tmp_path, lone_reads=200, shared_reads=200):
    """Reads 0..lone_reads-1 each alone on a reference; the others share one.

    Balanced by pair count, the first shards get the lone reads, which have
    the most partners but nothing to compare, and the last shard gets every
    comparison.
    """
    rng = random.Random(11)
    reference = tmp_path / 'refs.fasta'
    sam = tmp_path / 'reads.sam'
    refs = [''.join(rng.choice(BASES) for _ in range(300)) for _ in range(lone_reads + 1)]
    with reference.open('w') as handle:
        handle.write(''.join('>ref{}\n{}\n'.format(i, s) for i, s in enumerate(refs)))
    with sam.open('w') as handle:
        handle.write('@HD\tVN:1.0\n')
        handle.write(''.join('@SQ\tSN:ref{}\tLN:300\n'.format(i) for i in range(len(refs))))
        for read in range(lone_reads + shared_reads):
            ref_id = read if read < lone_reads else lone_reads
            pos = 1 if read < lone_reads else 1 + rng.randrange(200)
            handle.write('\t'.join(['read{:04d}'.format(read), '0', 'ref{}'.format(ref_id), str(pos),
                                    '255', '100M', '*', '0', '0', refs[ref_id][pos - 1:pos + 99],
                                    '*']) + '\n')
    return reference, sam


def shard_metadata(tmp_path, reference, sam, shard_count, *extra):
    merged, metadata = [], []
    for shard in range(shard_count):
        prefix = tmp_path / 'shard-{}'.format(shard)
        run(reference, sam, prefix, '--pair-shard-index', str(shard),
            '--pair-shard-count', str(shard_count), *extra)
        with open(str(prefix) + '.edges.csv') as handle:
            merged.extend(handle.readlines()[1:])
        with open(str(prefix) + '.pair-shard.tsv') as handle:
            metadata.append({k: int(v) for k, v in (line.rstrip('\n').split('\t') for line in handle)})
    return merged, metadata


def test_shards_are_balanced_by_comparisons_not_pairs(tmp_path):
    reference, sam = make_lopsided_inputs(tmp_path)
    run(reference, sam, tmp_path / 'full', '--all-pairs')
    merged, metadata = shard_metadata(tmp_path, reference, sam, 4)

    with open(str(tmp_path / 'full') + '.edges.csv') as handle:
        assert merged == handle.readlines()[1:]
    compared = [m['compared_pair_count'] for m in metadata]
    assert sum(compared) == 200 * 199 // 2
    assert min(compared) > 0.8 * max(compared)    # by pair count: [0, 0, 0, 19900]
    # Work: 1 per read plus one per (later read, shared reference).
    assert metadata[0]['total_work'] == 399 + 200 * 199 // 2
    assert sum(m['assigned_work'] for m in metadata) == metadata[0]['total_work']
    assert sum(m['assigned_pair_count'] for m in metadata) == 400 * 399 // 2


def test_all_pairs_shards_stay_balanced_by_pair_count(tmp_path):
    reference, sam = make_lopsided_inputs(tmp_path)
    _, metadata = shard_metadata(tmp_path, reference, sam, 4, '--all-pairs')

    assert all(m['assigned_work'] == m['assigned_pair_count'] for m in metadata)
    pairs = [m['assigned_pair_count'] for m in metadata]
    assert min(pairs) > 0.95 * max(pairs)
