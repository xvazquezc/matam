import os
import subprocess

import pytest


ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OVGRAPHBUILD = os.path.join(ROOT_DIR, 'ovgraphbuild', 'bin', 'ovgraphbuild')

pytestmark = pytest.mark.skipif(
    not os.path.isfile(OVGRAPHBUILD) or not os.access(OVGRAPHBUILD, os.X_OK),
    reason='ovgraphbuild binary is not built'
)


@pytest.fixture
def graph_inputs(tmp_path):
    reference = tmp_path / 'reference.fasta'
    reference.write_text(
        '>ref1 reference one\n' + 'A' * 200 + '\n'
        '>ref2 reference two\n' + 'C' * 200 + '\n'
    )

    records = [
        ('AAA_read1', 'ref1', 1, 'A' * 100),
        ('AAA_read2', 'ref1', 21, 'A' * 100),
        ('BBB_read3', 'ref1', 41, 'T' * 100),
        ('CCC_read4', 'ref2', 1, 'C' * 100),
        ('CCC_read5', 'ref2', 21, 'C' * 100),
    ]
    sam = tmp_path / 'reads.sam'
    with sam.open('w') as handle:
        handle.write('@HD\tVN:1.0\tSO:queryname\n')
        handle.write('@SQ\tSN:ref1\tLN:200\n')
        handle.write('@SQ\tSN:ref2\tLN:200\n')
        for name, reference_name, position, sequence in records:
            handle.write(
                '{}\t0\t{}\t{}\t255\t100M\t*\t0\t0\t{}\t*\tAS:i:100\n'.format(
                    name, reference_name, position, sequence
                )
            )

    return reference, sam


@pytest.fixture
def many_graph_inputs(tmp_path):
    reference = tmp_path / 'many-reference.fasta'
    reference.write_text('>ref1 reference one\n' + 'A' * 500 + '\n')

    sam = tmp_path / 'many-reads.sam'
    with sam.open('w') as handle:
        handle.write('@HD\tVN:1.0\tSO:queryname\n')
        handle.write('@SQ\tSN:ref1\tLN:500\n')
        for read_index in range(60):
            handle.write(
                'AAA_read{:03d}\t0\tref1\t{}\t255\t100M\t*\t0\t0\t{}\t*\tAS:i:100\n'.format(
                    read_index,
                    read_index + 1,
                    'A' * 100,
                )
            )

    return reference, sam


def run_ovgraphbuild(reference, sam, output_prefix, *extra_args):
    command = [
        OVGRAPHBUILD,
        '--reference', str(reference),
        '--sam', str(sam),
        '--csv',
        '--output_basename', str(output_prefix),
        '--min_overlap', '50',
        '--id_threshold', '1',
    ]
    command.extend(extra_args)
    return subprocess.run(command, check=True, capture_output=True, text=True)


def read_bytes(path):
    with open(path, 'rb') as handle:
        return handle.read()


def read_metadata(path):
    metadata = {}
    with open(path) as handle:
        for line in handle:
            key, value = line.rstrip().split('\t')
            metadata[key] = int(value)
    return metadata


@pytest.mark.parametrize('threads', [2, 4, 8, 24])
def test_threaded_output_is_byte_identical(graph_inputs, tmp_path, threads):
    reference, sam = graph_inputs
    serial = tmp_path / 'serial'
    threaded = tmp_path / 'threaded'

    run_ovgraphbuild(reference, sam, serial, '--threads', '1')
    run_ovgraphbuild(reference, sam, threaded, '--threads', str(threads))

    assert read_bytes(str(threaded) + '.nodes.csv') == read_bytes(str(serial) + '.nodes.csv')
    assert read_bytes(str(threaded) + '.edges.csv') == read_bytes(str(serial) + '.edges.csv')
    assert not list(tmp_path.glob('*.tmp'))


def test_dynamic_chunks_and_thread_cap(many_graph_inputs, tmp_path):
    reference, sam = many_graph_inputs
    output = tmp_path / 'capped'

    result = run_ovgraphbuild(
        reference,
        sam,
        output,
        '--threads', '25',
        '--pair-shard-index', '1',
        '--pair-shard-count', '2',
    )
    metadata = read_metadata(str(output) + '.pair-shard.tsv')

    assert 'Capping graph worker threads at 24' in result.stdout
    assert metadata['requested_threads'] == 25
    assert metadata['thread_limit'] == 24
    assert metadata['threads'] == 24
    assert metadata['chunk_count'] > metadata['threads']
    assert not list(tmp_path.glob('*.tmp'))


@pytest.mark.parametrize('shard_count', [3, 8])
def test_shards_cover_all_pairs_and_reproduce_edge_order(graph_inputs, tmp_path, shard_count):
    reference, sam = graph_inputs
    serial = tmp_path / 'serial'
    run_ovgraphbuild(reference, sam, serial)

    merged_edges = []
    metadata = []
    expected_nodes = read_bytes(str(serial) + '.nodes.csv')

    for shard_index in range(shard_count):
        shard = tmp_path / 'shard-{}'.format(shard_index)
        run_ovgraphbuild(
            reference,
            sam,
            shard,
            '--threads', '2',
            '--pair-shard-index', str(shard_index),
            '--pair-shard-count', str(shard_count),
        )
        assert read_bytes(str(shard) + '.nodes.csv') == expected_nodes
        with open(str(shard) + '.edges.csv') as handle:
            lines = handle.readlines()
        assert lines[0] == 'Source;Target;Type;Coemitted;Weight;Pid;MultiRef\n'
        merged_edges.extend(lines[1:])
        metadata.append(read_metadata(str(shard) + '.pair-shard.tsv'))

    with open(str(serial) + '.edges.csv') as handle:
        serial_edges = handle.readlines()[1:]

    assert merged_edges == serial_edges
    assert sum(item['assigned_pair_count'] for item in metadata) == 10
    assert sum(item['completed_pair_count'] for item in metadata) == 10
    assert all(item['total_pair_count'] == 10 for item in metadata)
    assert metadata[0]['outer_begin'] == 0
    assert metadata[-1]['outer_end'] == 4
    assert all(
        left['outer_end'] == right['outer_begin']
        for left, right in zip(metadata, metadata[1:])
    )
    assert not list(tmp_path.glob('*.tmp'))


def test_rejects_invalid_shard_index(graph_inputs, tmp_path):
    reference, sam = graph_inputs
    output = tmp_path / 'invalid'
    command = [
        OVGRAPHBUILD,
        '--reference', str(reference),
        '--sam', str(sam),
        '--csv',
        '--output_basename', str(output),
        '--pair-shard-index', '2',
        '--pair-shard-count', '2',
    ]

    result = subprocess.run(command, capture_output=True, text=True)

    assert result.returncode != 0
    assert '--pair-shard-index must be smaller' in result.stderr
