"""Input parsers: FASTQ quality encoding (#3), CIGAR (#4) and mpileup (#5)."""

import gzip
import importlib
import io
import sys
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).parent.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS_DIR))
import fastq_utils  # noqa: E402
import scaffold_contigs  # noqa: E402
sys.path.pop(0)


# --- FASTQ reader and Phred+33 detection (#3) ---------------------------------

def write_fastq(path, quals, gz=False):
    opener = gzip.open if gz else open
    with opener(path, 'wt') as handle:
        for i, qual in enumerate(quals):
            handle.write('@r{0}\n{1}\n+\n{2}\n'.format(i, 'A' * len(qual), qual))
    return str(path)


def test_reader_yields_nothing_for_an_empty_file():
    assert list(fastq_utils.read_fastq_file_handle(io.StringIO(''))) == []


def test_reader_yields_every_record():
    records = list(fastq_utils.read_fastq_file_handle(io.StringIO('@a x\nAC\n+\nII\n@b\nG\n+\n#\n')))
    assert records == [('a', 'AC', 'II'), ('b', 'G', '#')]


@pytest.mark.parametrize('quals', [
    ['#,:FFFJ'],                      # binned NovaSeq Phred+33
    ['IIIIIIII', 'IIIIIIII'],         # uniform high quality (ART, simulated)
    ['!KLMN~', 'IIII'],               # Phred+33 above Q41 ('K' = Q42 ... '~' = Q93)
    ['@@@@', 'JJJJ'],                 # ambiguous but within the Illumina Phred+33 range
], ids=['binned', 'uniform-high', 'above-Q41', 'ambiguous-range'])
def test_phred33_input_is_accepted(tmp_path, quals):
    assert fastq_utils.is_phred33(write_fastq(tmp_path / 'r.fq', quals))


@pytest.mark.parametrize('quals', [
    ['hhhhBBB', 'hhgfe'],             # Illumina 1.5 Phred+64
    ['@ABCDEFGHIJ', 'ghhh'],          # Illumina 1.3 Phred+64
    ['II\x7fI'],                      # outside the FASTQ range
], ids=['illumina-1.5', 'illumina-1.3', 'out-of-range'])
def test_phred64_or_invalid_input_is_rejected(tmp_path, quals):
    assert not fastq_utils.is_phred33(write_fastq(tmp_path / 'r.fq', quals))


def test_empty_fastq_is_not_rejected_and_does_not_crash(tmp_path):
    (tmp_path / 'empty.fq').write_text('')
    assert fastq_utils.is_phred33(str(tmp_path / 'empty.fq'))


def test_gzipped_input_and_read_limit(tmp_path):
    # Only the first read is tested: the Phred+64 tail is never reached.
    path = write_fastq(tmp_path / 'r.fq.gz', ['#FFF'] + ['hhhh'] * 3, gz=True)
    assert fastq_utils.is_phred33(path, number_of_reads_to_test=1)
    assert fastq_utils.is_phred33(path)  # '#' already settles it


# --- CIGAR parsing (#4) -------------------------------------------------------

def _parse_cigar_functions():
    sys.path.insert(0, str(SCRIPTS_DIR))
    try:
        modules = {name: importlib.import_module(name) for name in (
            'compute_contigs_compatibility', 'filter_sam_by_pid', 'sample_sam_by_coverage',
            'compute_ref_coverage_histogram', 'compute_assembly_stats')}
    finally:
        sys.path.pop(0)
    functions = {name: module.parse_cigar for name, module in modules.items()
                 if hasattr(module, 'parse_cigar')}
    functions['compute_contigs_compatibility'] = (
        modules['compute_contigs_compatibility'].SamAlignment.parse_cigar)
    return functions


PARSE_CIGAR = _parse_cigar_functions()
MERGING = [name for name in PARSE_CIGAR if name != 'compute_assembly_stats']


@pytest.mark.parametrize('name', sorted(PARSE_CIGAR))
def test_parse_cigar_reads_operations_and_counts(name):
    assert PARSE_CIGAR[name]('5S10M2I3M1D4M') == [
        ('S', 5), ('M', 10), ('I', 2), ('M', 3), ('D', 1), ('M', 4)]


@pytest.mark.parametrize('name', sorted(MERGING))
def test_parse_cigar_merges_adjacent_repeated_operations(name):
    # Used to raise TypeError ('tuple' object does not support item assignment)
    assert PARSE_CIGAR[name]('5M5M2D1D3M') == [('M', 10), ('D', 3), ('M', 3)]


# --- mpileup read bases (#5) --------------------------------------------------

def test_mpileup_read_bases_are_split_into_events():
    events = list(scaffold_contigs.iter_read_bases('^F.,+2AG-1c*$A'))
    assert events == ['^F', '.', ',', '+AG', '-c', '*', '$', 'A']


def test_unknown_mpileup_character_raises_a_clear_error():
    # Used to raise NameError: ParsingError is not defined
    with pytest.raises(ValueError, match=r"'#' at index 2"):
        list(scaffold_contigs.iter_read_bases('..#'))
