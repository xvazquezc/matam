"""ovgraphbuild node identity: a node is one read, keyed by (QNAME, mate)."""

import os
import subprocess

import pytest

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OVGRAPHBUILD = os.path.join(ROOT_DIR, 'ovgraphbuild', 'bin', 'ovgraphbuild')

pytestmark = pytest.mark.skipif(
    not os.path.isfile(OVGRAPHBUILD) or not os.access(OVGRAPHBUILD, os.X_OK),
    reason='ovgraphbuild binary is not built'
)

REF = 'ACGTTGCAAGGCTTACCGATTGACCTGAAGTCCATGGATCCGTAGCTAGGCTAACGTTAGC' * 4
R1, R2 = 0x1 | 0x40, 0x1 | 0x80


def run(tmp_path, records):
    reference = tmp_path / 'ref.fasta'
    reference.write_text('>ref1\n{}\n>ref2\n{}\n'.format(REF, REF))
    sam = tmp_path / 'reads.sam'
    with sam.open('w') as handle:
        handle.write('@HD\tVN:1.0\n@SQ\tSN:ref1\tLN:{0}\n@SQ\tSN:ref2\tLN:{0}\n'.format(len(REF)))
        for qname, flag, ref, pos in records:
            seq = REF[pos - 1:pos - 1 + 60]
            cigar, ref_name, position = '60M', ref, pos
            if flag & 0x4:
                cigar, ref_name, position = '*', '*', 0
            handle.write('\t'.join([qname, str(flag), ref_name, str(position), '255', cigar,
                                    '*', '0', '0', seq, '*', 'AS:i:120']) + '\n')
    prefix = tmp_path / 'graph'
    completed = subprocess.run(
        [OVGRAPHBUILD, '--reference', str(reference), '--sam', str(sam), '--csv',
         '--output_basename', str(prefix), '--min_overlap', '20', '--id_threshold', '1'],
        capture_output=True, text=True)
    nodes = []
    nodes_csv = tmp_path / 'graph.nodes.csv'
    if nodes_csv.exists():
        lines = nodes_csv.read_text().splitlines()
        nodes = [line.split(';') for line in lines[1:]]
        assert lines[0] == 'Id;Label;Mate;Specie'
    return completed, nodes


def test_mates_sharing_a_qname_are_separate_nodes(tmp_path):
    completed, nodes = run(tmp_path, [('pairA', R1, 'ref1', 1), ('pairA', R2 | 0x10, 'ref1', 41)])

    assert completed.returncode == 0, completed.stderr
    assert [(n[1], n[2]) for n in nodes] == [('pairA', '1'), ('pairA', '2')]


def test_all_alignments_of_one_read_form_one_node(tmp_path):
    # One read, two placements: must stay one node (regression check for the
    # moved-from qName in readSamFile).
    completed, nodes = run(tmp_path, [('readA', 0, 'ref1', 1), ('readA', 0, 'ref2', 1),
                                      ('readB', 0, 'ref1', 21)])

    assert completed.returncode == 0, completed.stderr
    assert [(n[1], n[2]) for n in nodes] == [('readA', '0'), ('readB', '0')]


def test_mapped_record_after_an_unmapped_one_starts_its_own_node(tmp_path):
    completed, nodes = run(tmp_path, [('readA', 0, 'ref1', 1), ('readB', 0x4, '*', 0),
                                      ('readB', 0, 'ref1', 21)])

    assert completed.returncode == 0, completed.stderr
    assert [n[1] for n in nodes] == ['readA', 'readB']


def test_non_contiguous_records_of_one_read_are_rejected(tmp_path):
    completed, _ = run(tmp_path, [('pairA', R1, 'ref1', 1), ('pairA', R2, 'ref1', 41),
                                  ('pairA', R1, 'ref2', 1)])

    assert completed.returncode != 0
    assert 'not contiguous' in completed.stderr


def test_flagged_shared_names_build_the_same_graph_as_unique_legacy_names(tmp_path):
    # Overlapping pairs: R1 forward at p, R2 reverse-complemented at p + 30.
    # Mates share a QNAME and differ only by 0x40/0x80, as MATAM tags them;
    # the legacy copy uses unique /1, /2 names without pairing flags.
    shared, legacy = [], []
    for index, start in enumerate((1, 11, 21, 31, 41)):
        name = 'pair{}'.format(index)
        shared += [(name, R1, 'ref1', start), (name, R2 | 0x10, 'ref1', start + 30)]
        legacy += [(name + '/1', 0, 'ref1', start), (name + '/2', 0x10, 'ref1', start + 30)]

    (tmp_path / 'shared').mkdir()
    (tmp_path / 'legacy').mkdir()
    shared_run, shared_nodes = run(tmp_path / 'shared', shared)
    legacy_run, legacy_nodes = run(tmp_path / 'legacy', legacy)

    assert shared_run.returncode == 0 and legacy_run.returncode == 0
    assert len(shared_nodes) == len(legacy_nodes) == 10
    assert [(n[1] + '/' + n[2]) for n in shared_nodes] == [n[1] for n in legacy_nodes]
    shared_edges = (tmp_path / 'shared' / 'graph.edges.csv').read_bytes()
    assert shared_edges.count(b'\n') > 1  # the reads do overlap
    assert shared_edges == (tmp_path / 'legacy' / 'graph.edges.csv').read_bytes()
