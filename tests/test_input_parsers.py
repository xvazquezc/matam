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



# --- Paired FASTQ input validation (validate_paired_fastq) -------------------

def pair_files(tmp_path, records1, records2, suffix='.fq'):
    """R1/R2 files, each from (id, seq, qual) records or from raw text."""
    paths = []
    for mate, records in ((1, records1), (2, records2)):
        path = tmp_path / 'r{}{}'.format(mate, suffix)
        text = records if isinstance(records, str) else ''.join(
            '@{}\n{}\n+\n{}\n'.format(*r) for r in records)
        opener = gzip.open if suffix.lower().endswith('.gz') else open
        with opener(str(path), 'wt', newline='') as handle:
            handle.write(text)
        paths.append(str(path))
    return paths


def reads(n, name='r{}', qual='IIII'):
    return [(name.format(i), 'ACGT', qual) for i in range(n)]


@pytest.fixture(params=[fastq_utils._BATCH, 3], ids=['batch', 'small-batch'])
def batch(request, monkeypatch):
    """Run each test with the real batch size and with 3-record batches."""
    monkeypatch.setattr(fastq_utils, '_BATCH', request.param)
    return request.param


@pytest.mark.parametrize('name1, name2', [
    ('M01:1:FC:1:1101:{} 1:N:0:ACGT', 'M01:1:FC:1:1101:{} 2:N:0:ACGT'),   # CASAVA 1.8+
    ('r{}/1', 'r{}/2'),                                                  # legacy
    ('SRR1.{}', 'SRR1.{}'),                                              # SRA
    ('SRR1.{}.1', 'SRR1.{}.2'),                                          # fastq-dump --readids
    ('r{}_1', 'r{}_2'),
], ids=['casava', 'legacy', 'sra', 'readids', 'underscore'])
@pytest.mark.parametrize('suffix', ['.fq', '.fastq.gz', '.FQ.GZ'])
def test_valid_pairs_are_counted(tmp_path, batch, name1, name2, suffix):
    forward, reverse = pair_files(tmp_path, reads(7, name1), reads(7, name2), suffix)
    assert fastq_utils.validate_paired_fastq(forward, reverse) == 7


def test_line_ends_and_trailing_blank_lines_are_accepted(tmp_path, batch):
    text = '@a\nAC\n+\nII\r\n@b\r\nGT\r\n+\r\nJJ'               # CRLF, no final newline
    forward, reverse = pair_files(tmp_path, text, text + '\n\n\n')
    assert fastq_utils.validate_paired_fastq(forward, reverse) == 2


def raises(forward, reverse, *fragments):
    with pytest.raises(fastq_utils.FastqFormatError) as error:
        fastq_utils.validate_paired_fastq(forward, reverse)
    for fragment in fragments:
        assert fragment in str(error.value)
    return str(error.value)


def test_record_problems_name_the_file_and_record(tmp_path, batch):
    good = reads(5)
    for bad_record, message in [
            ('@r3\nACGT\n+\nIII\n', 'record 4 has 4 bases but 3 quality values'),
            ('@r3\nACGT\nIIII\n+\n', 'record 4 is not FASTQ'),
            ('r3\nACGT\n+\nIIII\n', 'record 4 is not FASTQ'),
            ('\n@r3\nACGT\n+\nIIII\n', 'record 4'),
    ]:
        text = ''.join('@{}\n{}\n+\n{}\n'.format(*r) for r in good[:3]) + bad_record + \
            '@r4\nACGT\n+\nIIII\n'
        forward, reverse = pair_files(tmp_path, reads(5), text)
        raises(forward, reverse, reverse, message)


def test_a_truncated_last_record_is_an_error(tmp_path, batch):
    text = ''.join('@{}\n{}\n+\n{}\n'.format(*r) for r in reads(4)) + '@r4\nACGT\n'
    forward, reverse = pair_files(tmp_path, reads(5), text)
    raises(forward, reverse, reverse, 'the last record is truncated')


@pytest.mark.parametrize('n1, n2', [(5, 4), (4, 5), (10, 3), (3, 10)])
def test_different_read_numbers_are_an_error(tmp_path, batch, n1, n2):
    # (5, 4): exactly one extra read in R1, which a plain zip() would miss.
    forward, reverse = pair_files(tmp_path, reads(n1), reads(n2))
    longer = forward if n1 > n2 else reverse
    raises(forward, reverse, 'different numbers of reads', longer,
           'the {} reads of the other file'.format(min(n1, n2)))


@pytest.mark.parametrize('name1, name2', [('r{}', 'r{}'), ('r{}/1', 'r{}/2')])
def test_pairs_out_of_step_are_an_error(tmp_path, batch, name1, name2):
    second = reads(8, name2)
    second[5], second[6] = second[6], second[5]          # two R2 records swapped
    forward, reverse = pair_files(tmp_path, reads(8, name1), second)
    raises(forward, reverse, 'not in step: read 6', name1.format(5), name2.format(6))


@pytest.mark.parametrize('qual, message', [
    ('hhgf', 'looks like Phred+64'),
    ('II I', 'quality character outside the FASTQ range'),
    ('II\x7fI', 'quality character outside the FASTQ range'),
])
def test_quality_encoding_is_checked_in_the_first_reads(tmp_path, batch, qual, message):
    first = reads(6, qual='hhhh' if 'Phred' in message else 'IIII')
    first[5] = ('r5', 'ACGT', qual)                     # after the first batch of 3
    forward, reverse = pair_files(tmp_path, first, reads(6, qual=first[0][2]))
    raises(forward, reverse, message)


def test_qualities_after_the_sample_are_not_checked(tmp_path, batch, monkeypatch):
    # The encoding is the same for a whole file: only the first pairs are read for it.
    monkeypatch.setattr(fastq_utils, '_QUALITY_SAMPLE', 4)
    first = reads(6)
    first[5] = ('r5', 'ACGT', 'hhgf')
    forward, reverse = pair_files(tmp_path, first, reads(6))
    assert fastq_utils.validate_paired_fastq(forward, reverse) == 6
    first[3] = ('r3', 'ACGT', 'II\x7fI')                # within the sample
    forward, reverse = pair_files(tmp_path, first, reads(6))
    raises(forward, reverse, 'record 4 has a quality character outside the FASTQ range')


def test_phred33_with_high_qualities_is_accepted(tmp_path, batch):
    # Above 'J' but with a value below ';' somewhere: only Phred+33 has those.
    first = reads(6, qual='KKKK') + [('r6', 'ACGT', '#KKK')]
    forward, reverse = pair_files(tmp_path, first, reads(7, qual='KKKK'))
    assert fastq_utils.validate_paired_fastq(forward, reverse) == 7


def test_empty_input_and_other_file_types_are_errors(tmp_path):
    forward, reverse = pair_files(tmp_path, [], [])
    raises(forward, reverse, 'contain no reads')
    raises(str(tmp_path / 'r1.fasta'), str(tmp_path / 'r2.fq'), 'not a FASTQ file name')


def test_a_corrupt_gzip_file_is_an_error(tmp_path, batch):
    forward, reverse = pair_files(tmp_path, reads(50), reads(50), '.fq.gz')
    data = Path(reverse).read_bytes()
    Path(reverse).write_bytes(data[:len(data) // 2])
    raises(forward, reverse, reverse, 'cannot be read')


@pytest.mark.parametrize('id1, id2, same', [
    ('r', 'r', True), ('r/1', 'r/2', True), ('S.7.1', 'S.7.2', True), ('r_1', 'r_2', True),
    ('r/2', 'r/1', False), ('r5/1', 'r6/2', False), ('r/1', 'r.2', False), ('a1', 'a2', False),
])
def test_is_same_pair(id1, id2, same):
    assert fastq_utils.is_same_pair(id1, id2) is same


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
