import io
import subprocess
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).parent.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS_DIR))
import compute_lca_from_tab  # noqa: E402
import pe_tables  # noqa: E402
sys.path.pop(0)

R1, R2 = 0x1 | 0x40, 0x1 | 0x80


def write(path, text):
    path.write_text(text)
    return path


def test_read_component_table_keeps_mates_apart_and_drops_excluded_reads(tmp_path):
    components = write(tmp_path / 'components.csv',
                       'Id;Label;Mate;Specie;metanode;component\n'
                       '0;pairA;1;pai;0;0\n'
                       '1;pairA;2;pai;3;1\n'
                       '2;pairB;1;pai;-1;-1\n')
    out = tmp_path / 'read_component.tab'

    stats = pe_tables.write_read_component_table(components, out)

    assert out.read_text() == 'pairA\t1\t0\t0\npairA\t2\t3\t1\n'
    assert stats == pe_tables.GraphStats(nodes=3, reads_in_components=2, components=2)
    assert pe_tables.load_read_components(out) == {('pairA', 1): ('0', '0'), ('pairA', 2): ('3', '1')}


def test_components_without_mate_column_are_unpaired_reads(tmp_path):
    components = write(tmp_path / 'components.csv',
                       'Id;Label;Specie;metanode;component\n0;readA;rea;0;5\n')

    assert list(pe_tables.iter_components_csv(components)) == [('readA', 0, 0, 5)]


def test_complete_taxo_table_joins_on_read_and_mate(tmp_path):
    read_component = write(tmp_path / 'read_component.tab', 'pairA\t1\t0\t0\npairA\t2\t3\t1\n')
    taxo = write(tmp_path / 'taxo.tab', 'refX\tBac;Firm\nrefY\tBac;Prot\n')
    sam = write(tmp_path / 'filt.sam', ''.join(
        '\t'.join([q, str(f), ref, '1', '255', '4M', '*', '0', '0', 'ACGT', '*']) + '\n'
        for q, f, ref in [('pairA', R1, 'refX'), ('pairA', R2, 'refY'),
                          ('pairA', R2, 'refUnknown'), ('pairZ', R1, 'refX')]))
    out = tmp_path / 'complete_taxo.tab'

    rows = pe_tables.write_complete_taxo_table(sam, taxo, read_component, out)

    assert rows == 2
    assert out.read_text() == 'pairA\t1\t0\t0\tBac;Firm\npairA\t2\t3\t1\tBac;Prot\n'


def groups(text, factor, group_by):
    return [[[row[0] + ':' + row[1] for row in group] for group in factor_list]
            for factor_list in compute_lca_from_tab.read_tab_file_handle_sorted(
                io.StringIO(text), factor, group_by, None)]


def test_group_continuing_into_the_next_factor_stays_with_its_factor():
    # Regression for TESTING_PLAN #2: read r2 ends component A and starts B.
    text = 'r1 x A t\nr2 x A t\nr2 y B t\nr3 y B t\n'

    assert groups(text, 2, 0) == [[['r1:x'], ['r2:x']], [['r2:y'], ['r3:y']]]


def test_multi_column_groups_split_mates_of_one_read():
    text = 'pairA 1 C t1\npairA 1 C t2\npairA 2 C t3\n'

    assert groups(text, 2, [0, 1]) == [[['pairA:1', 'pairA:1'], ['pairA:2']]]
    assert groups(text, 2, 0) == [[['pairA:1', 'pairA:1', 'pairA:2']]]


def test_empty_input_yields_no_factor():
    assert groups('', 2, 0) == []


def run_lca(tmp_path, text, *args):
    tab = write(tmp_path / 'in.tab', text)
    out = tmp_path / 'out.tab'
    subprocess.run([sys.executable, str(SCRIPTS_DIR / 'compute_lca_from_tab.py'),
                    '-i', str(tab), '-o', str(out), *args], check=True)
    return out.read_text()


def test_cli_labels_each_component_with_its_own_reads(tmp_path):
    text = ('r1\tm1\tA\tBac;Firm\nr2\tm2\tA\tBac;Firm\n'
            'r2\tm7\tB\tBac;Prot\nr3\tm8\tB\tBac;Prot\n')

    for quorum in ('1.0', '0.51'):
        assert run_lca(tmp_path, text, '-t', '4', '-f', '3', '-g', '1', '-m', quorum) == \
            'A\tBac;Firm\nB\tBac;Prot\n'


def test_cli_counts_each_mate_as_one_vote(tmp_path):
    # pipeline layout: read mate metanode component taxo
    text = ('pairA\t1\t0\tC\tBac;Firm\n'
            'pairA\t2\t0\tC\tBac;Firm\n'
            'pairB\t1\t0\tC\tBac;Prot\n')

    # Two of three mate votes are Firmicutes: above a 0.6 quorum.
    assert run_lca(tmp_path, text, '-t', '5', '-f', '4', '-g', '1,2', '-m', '0.6') == 'C\tBac;Firm\n'
    # Grouping on the QNAME alone merges the pair into one vote: 1 of 2.
    assert run_lca(tmp_path, text, '-t', '5', '-f', '4', '-g', '1', '-m', '0.6') == 'C\tBac\n'


def test_cli_handles_a_component_made_of_one_shared_read(tmp_path):
    text = 'r2\tm2\tA\tBac;Firm\nr2\tm7\tB\tBac;Prot\n'

    assert run_lca(tmp_path, text, '-t', '4', '-f', '3', '-g', '1') == 'A\tBac;Firm\nB\tBac;Prot\n'
