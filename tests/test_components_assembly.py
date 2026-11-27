"""Characterization of components_assembly.assemble_component and the coverage
estimate that decides on read correction (--read_correction auto)."""
import gzip
import logging
import sys
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).parent.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS_DIR))
import components_assembly  # noqa: E402
sys.path.pop(0)


def write_fastq(path, seqs):
    text = ''.join('@r{}\n{}\n+\n{}\n'.format(i, s, 'I' * len(s)) for i, s in enumerate(seqs))
    if str(path).endswith('.gz'):
        with gzip.open(path, 'wt') as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return str(path)


@pytest.fixture
def reads(tmp_path):
    # 2 x 100 nt per mate file: 400 read nucleotides in all
    return (write_fastq(tmp_path / 'c_R1.fq', ['A' * 100, 'C' * 100]),
            write_fastq(tmp_path / 'c_R2.fq', ['G' * 100, 'T' * 100]))


class FakeAssembler:
    """Writes one prepared contig set per run; records each build."""

    def __init__(self, tmp_path, outputs, fail=False):
        self.tmp_path, self.outputs, self.fail = tmp_path, list(outputs), fail
        self.builds = []

    def build_command_line(self, forward, reverse, workdir, read_correction, cpu):
        self.builds.append(read_correction)

    def run(self):
        if self.fail:
            sys.exit('assembler failed')
        fasta = self.tmp_path / 'assembly{}.fasta'.format(len(self.builds))
        fasta.write_text(self.outputs[len(self.builds) - 1])
        return str(fasta)


@pytest.fixture
def fake(monkeypatch, tmp_path):
    def install(*outputs, fail=False):
        assembler = FakeAssembler(tmp_path, outputs, fail)
        monkeypatch.setattr(components_assembly.AssemblerFactory, 'get', lambda self, name: assembler)
        return assembler
    return install


def assemble(reads, tmp_path, read_correction='auto', coverage_threshold=20):
    return components_assembly.assemble_component(
        'SGA', reads[0], reads[1], str(tmp_path / 'wd'), read_correction, 1, coverage_threshold)


def test_high_coverage_reruns_the_assembly_with_read_correction(reads, tmp_path, fake):
    # 400 read nt / 10 contig nt = 40 > 20
    assembler = fake('>1\n' + 'A' * 10 + '\n', '>1\n' + 'A' * 50 + '\n')
    fasta = assemble(reads, tmp_path)
    assert assembler.builds == ['auto', 'yes']
    assert fasta.endswith('assembly2.fasta')


def test_coverage_at_the_threshold_keeps_the_uncorrected_assembly(reads, tmp_path, fake):
    # 400 / 20 = 20, not above 20
    assembler = fake('>1\n' + 'A' * 20 + '\n')
    assert assemble(reads, tmp_path).endswith('assembly1.fasta')
    assert assembler.builds == ['auto']


@pytest.mark.parametrize('mode', ['no', 'yes'])
def test_only_auto_mode_reruns(reads, tmp_path, fake, caplog, mode):
    assembler = fake('>1\nACGT\n')
    with caplog.at_level(logging.WARNING):
        assemble(reads, tmp_path, read_correction=mode, coverage_threshold=20)
    assert assembler.builds == [mode]
    assert 'Coverage_threshold 20 makes no sense' in caplog.text


def test_empty_assembly_gives_none(reads, tmp_path, fake):
    assembler = fake('')
    assert assemble(reads, tmp_path) is None
    assert assembler.builds == ['auto']


def test_empty_corrected_assembly_gives_none(reads, tmp_path, fake):
    fake('>1\nACGT\n', '')
    assert assemble(reads, tmp_path) is None


def test_zero_length_contigs_are_kept_without_rerun(reads, tmp_path, fake):
    """Headers without sequence: the coverage can't be estimated (None), so
    there is no corrected rerun; the file is returned as is."""
    assembler = fake('>1\n\n')
    assert assemble(reads, tmp_path).endswith('assembly1.fasta')
    assert assembler.builds == ['auto']


def test_assembler_exit_becomes_an_exception(reads, tmp_path, fake):
    """SystemExit would hang the multiprocessing pool, so it is re-raised as Exception."""
    fake(fail=True)
    with pytest.raises(Exception) as info:
        assemble(reads, tmp_path)
    assert not isinstance(info.value, SystemExit)


def test_coverage_is_read_nucleotides_over_contig_nucleotides(tmp_path):
    forward = write_fastq(tmp_path / 'r_R1.fq', ['A' * 30])
    reverse = write_fastq(tmp_path / 'r_R2.fq', ['C' * 30])
    contigs = tmp_path / 'contigs.fa'
    contigs.write_text('>1\nACGTACGTAC\n>2\nACGTA\nCGTAC\n')
    assert components_assembly.estimate_coverage(forward, reverse, str(contigs)) == 3.0


def test_coverage_of_zero_length_contigs_is_none(tmp_path):
    forward = write_fastq(tmp_path / 'r_R1.fq', ['A' * 30])
    contigs = tmp_path / 'contigs.fa'
    contigs.write_text('>1\n')
    assert components_assembly.estimate_coverage(forward, forward, str(contigs)) is None


@pytest.mark.parametrize('name', ['r.fq.gz', 'r.fastq.gz'])
def test_nucleotide_count_reads_gzipped_fastq(tmp_path, name):
    assert components_assembly.nucleotidic_number(write_fastq(tmp_path / name, ['ACGT', 'AC'])) == 6


@pytest.mark.parametrize('name', ['c.fa.gz', 'c.fasta.gz', 'c.fna'])
def test_nucleotide_count_reads_fasta_names(tmp_path, name):
    path = tmp_path / name
    text = '>1\nACG\nT\n>2\nAC\n'
    if name.endswith('.gz'):
        with gzip.open(path, 'wt') as fh:
            fh.write(text)
    else:
        path.write_text(text)
    assert components_assembly.nucleotidic_number(str(path)) == 6


def test_nucleotide_count_of_an_unknown_format_stops(tmp_path):
    path = tmp_path / 'reads.txt'
    path.write_text('ACGT\n')
    with pytest.raises(SystemExit):
        components_assembly.nucleotidic_number(str(path))


@pytest.mark.parametrize('id1, id2', [
    ('r1', 'r1'), ('r1/1', 'r1/2'), ('Heavy_1691.1', 'Heavy_1691.2'), ('r1_1', 'r1_2'), ('r1:1', 'r1:2'),
])
def test_pairs_are_split_by_component_for_every_accepted_naming(id1, id2):
    """The input check (validate_paired_fastq) accepts these namings, so the
    component split must too. MarkerMAG's .1/.2 names used to stop every
    assembly with 'Mates are out of sync'."""
    pairs = [((id1, 'ACGT', 'IIII'), (id2, 'TTTT', 'IIII'))]
    components = components_assembly._split_pairs_by_component(pairs, {(id1, 1): '7', (id2, 2): '7'})
    assert dict(components) == {'7': ([(id1, 'ACGT', 'IIII')], [(id2, 'TTTT', 'IIII')])}


@pytest.mark.parametrize('id1, id2', [('r1.1', 'r2.2'), ('r1/1', 'r1.2'), ('r1', 'r2')])
def test_mates_out_of_step_stop_the_split(id1, id2):
    pairs = [((id1, 'ACGT', 'IIII'), (id2, 'TTTT', 'IIII'))]
    with pytest.raises(SystemExit):
        components_assembly._split_pairs_by_component(pairs, {})
