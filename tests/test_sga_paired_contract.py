import sys
from pathlib import Path


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


def test_component_extraction_keeps_only_complete_same_component_pairs(tmp_path):
    selected_reads = tmp_path / 'selected.fastq'
    components = tmp_path / 'read_component.tab'
    output_dir = tmp_path / 'components'
    write_fastq(selected_reads, ['forward-1', 'reverse-1',
                                 'forward-2', 'reverse-2'])
    components.write_text(
        'forward-1 meta component-A\n'
        'reverse-1 meta component-A\n'
        'forward-2 meta component-A\n'
        'reverse-2 meta component-B\n')

    pairs = components_assembly.extract_interleaved_paired_reads_by_component(
        str(selected_reads), str(components))
    files = components_assembly.save_components(pairs, str(output_dir))

    assert list(files) == ['component-A']
    forward, reverse = files['component-A']
    assert read_fastq_ids(forward) == ['forward-1']
    assert read_fastq_ids(reverse) == ['reverse-1']


def test_separate_component_input_uses_pair_positions_and_requires_both_mates(tmp_path):
    forward_reads = tmp_path / 'reads_R1.fastq'
    reverse_reads = tmp_path / 'reads_R2.fastq'
    components = tmp_path / 'read_component.tab'
    write_fastq(forward_reads, ['forward-1', 'forward-2'])
    write_fastq(reverse_reads, ['reverse-1', 'reverse-2'])
    components.write_text(
        'forward-1 meta component-A\n'
        'reverse-1 meta component-A\n'
        'forward-2 meta component-A\n')

    pairs = components_assembly.extract_separate_paired_reads_by_component(
        str(forward_reads), str(reverse_reads), str(components))

    assert list(pairs) == ['component-A']
    assert [record[0] for record in pairs['component-A'][0]] == ['forward-1']
    assert [record[0] for record in pairs['component-A'][1]] == ['reverse-1']


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
