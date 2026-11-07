import sys
from pathlib import Path

import pytest


SCRIPTS_DIR = Path(__file__).parent.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS_DIR))
import components_assembly  # noqa: E402
sys.path.pop(0)


def write_fastq(path, read_ids):
    with open(path, 'w') as output:
        for read_id in read_ids:
            output.write('@{}\nACGT\n+\nIIII\n'.format(read_id))


def read_fastq_ids(path):
    with open(path) as input_handle:
        return [line[1:].strip() for line in input_handle if line.startswith('@')]


def extract(tmp_path, forward_ids, reverse_ids, table):
    forward, reverse = tmp_path / 'reads_R1.fastq', tmp_path / 'reads_R2.fastq'
    write_fastq(forward, forward_ids)
    write_fastq(reverse, reverse_ids)
    components = tmp_path / 'read_component.tab'
    components.write_text(table)  # read  mate  metanode  component
    return components_assembly.extract_separate_paired_reads_by_component(
        str(forward), str(reverse), str(components))


def test_a_pair_spanning_two_components_goes_to_both(tmp_path):
    pairs = extract(tmp_path, ['pair-1/1', 'pair-2/1'], ['pair-1/2', 'pair-2/2'],
                    'pair-1/1\t1\tmeta\tcomponent-A\n'
                    'pair-1/2\t2\tmeta\tcomponent-A\n'
                    'pair-2/1\t1\tmeta\tcomponent-A\n'
                    'pair-2/2\t2\tmeta\tcomponent-B\n')
    files = components_assembly.save_components(pairs, str(tmp_path / 'components'))

    assert sorted(files) == ['component-A', 'component-B']
    forward, reverse = files['component-A']
    assert read_fastq_ids(forward) == ['pair-1/1', 'pair-2/1']
    assert read_fastq_ids(reverse) == ['pair-1/2', 'pair-2/2']
    forward, reverse = files['component-B']
    assert read_fastq_ids(forward) == ['pair-2/1']
    assert read_fastq_ids(reverse) == ['pair-2/2']


def test_pair_routing_is_logged(tmp_path, caplog):
    caplog.set_level('INFO')
    extract(tmp_path, ['p1', 'p2', 'p3', 'p4'], ['p1', 'p2', 'p3', 'p4'],
            'p1\t1\tm\tA\np1\t2\tm\tA\n'   # within A
            'p2\t1\tm\tA\np2\t2\tm\tB\n'   # spanning A/B
            'p3\t2\tm\tB\n')                  # one mate in B; p4 in no component
    assert ('Read pairs by component: 1 within one component, 1 spanning two (given to both), '
            '1 with one mate in a component, 1 in no component (not assembled)') in caplog.text


def test_casava_mates_sharing_an_id_are_looked_up_by_mate(tmp_path):
    # CASAVA naming: both mates have the same ID (the comment is not part of it).
    ids = ['M01:1:FC:1:1:1:1', 'M01:1:FC:1:1:1:2']
    pairs = extract(tmp_path,
                    [i + ' 1:N:0:ACGT' for i in ids], [i + ' 2:N:0:ACGT' for i in ids],
                    'M01:1:FC:1:1:1:1\t1\tm0\tA\n'
                    'M01:1:FC:1:1:1:1\t2\tm1\tA\n'   # pair 1: both mates in A -> kept
                    'M01:1:FC:1:1:1:2\t1\tm0\tA\n'
                    'M01:1:FC:1:1:1:2\t2\tm2\tB\n')  # pair 2: split across A/B -> both

    assert sorted(pairs) == ['A', 'B']
    assert [record[0] for record in pairs['A'][0]] == ids
    assert [record[0] for record in pairs['A'][1]] == ids
    assert [record[0] for record in pairs['B'][1]] == ['M01:1:FC:1:1:1:2']


def test_a_pair_with_one_mate_in_a_component_goes_to_that_component(tmp_path):
    pairs = extract(tmp_path, ['pair-1', 'pair-2', 'pair-3'], ['pair-1', 'pair-2', 'pair-3'],
                    'pair-1\t1\tmeta\tcomponent-A\n'
                    'pair-1\t2\tmeta\tcomponent-A\n'
                    'pair-2\t1\tmeta\tcomponent-A\n')   # pair-3: no mate in a component

    assert list(pairs) == ['component-A']
    assert [record[0] for record in pairs['component-A'][0]] == ['pair-1', 'pair-2']
    assert [record[0] for record in pairs['component-A'][1]] == ['pair-1', 'pair-2']


@pytest.mark.parametrize('forward_ids, reverse_ids', [
    (['pair-1', 'pair-2'], ['pair-1', 'pair-3']),   # mates out of sync
    (['pair-1', 'pair-2'], ['pair-1']),             # different record counts
])
def test_inconsistent_pair_files_are_a_fatal_error(tmp_path, forward_ids, reverse_ids):
    with pytest.raises(SystemExit):
        extract(tmp_path, forward_ids, reverse_ids, 'pair-1\t1\tm\tA\npair-1\t2\tm\tA\n')


def test_empty_selection_gives_no_components(tmp_path):
    assert dict(extract(tmp_path, [], [], '')) == {}


def test_component_without_lca_is_labelled_null(tmp_path):
    # Component 2 follows a labelled component, so a stale LCA would leak into it.
    component_1 = tmp_path / 'component1.fa'
    component_2 = tmp_path / 'component2.fa'
    contigs = tmp_path / 'contigs.fa'
    component_1.write_text('>contig-0\nACGT\n>contig-1\n\n')
    component_2.write_text('>contig-0\nGGGG\n')

    components_assembly.concat_components_fasta_with_lca(
        {'1': str(component_1), '2': str(component_2)},
        str(contigs),
        {'1': 'Bacteria;Firmicutes'})

    assert contigs.read_text() == (
        '>1 component=1 lca=Bacteria;Firmicutes\nACGT\n'
        '>2 component=2 lca=NULL\nGGGG\n')


def test_sga_driver_accepts_only_separate_paired_files():
    source = (SCRIPTS_DIR / 'sga_assemble.py').read_text()

    assert "parser.add_argument('-i', '--input_fastq'" not in source
    assert "parser.add_argument('--forward'" in source
    assert "parser.add_argument('--reverse'" in source
    assert 'preprocess -v --pe-mode 1 ' in source


def test_sga_factory_passes_both_component_files_to_driver():
    source = (SCRIPTS_DIR / 'assembler_factory.py').read_text()

    assert "--forward ' + self.forward_fastq" in source
    assert "--reverse ' + self.reverse_fastq" in source
    assert "--paired_end" not in source
