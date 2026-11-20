"""SortMeRNA behaviour that MATAM's mate tagging relies on.

sam_mates.py assumes SortMeRNA (7.x):
* writes no SAM pairing flags (FLAG 0/16), so MATAM must add them;
* writes SAM records in input order, each read's records contiguous;
* writes full read sequences in SEQ, so each record matches one mate;
* with --paired_in --out2, writes the aligned pairs as R1/R2 files in step,
  in input order.
If a SortMeRNA release changes this, this test names the broken assumption.
"""

import collections
import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from pe_fixtures import iter_example_pairs, write_pe_reads, write_reference_fasta

SCRIPTS_DIR = Path(__file__).parent.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS_DIR))
import sam_mates  # noqa: E402
sys.path.pop(0)

SORTMERNA = shutil.which('sortmerna')
pytestmark = pytest.mark.skipif(SORTMERNA is None, reason='requires SortMeRNA 7 on PATH')


@pytest.fixture(scope='module')
def casava_mapping(tmp_path_factory):
    base = tmp_path_factory.mktemp('smr_pe')
    reference = write_reference_fasta(base / 'refs.fasta')
    forward, reverse = write_pe_reads(base, 'casava', gz=True)
    prefix = base / 'aligned'
    # Same options as the read-mapping step of matam_assembly.py
    subprocess.run(
        [SORTMERNA, '--ref', str(reference), '--workdir', str(base / 'wd'),
         '--reads', str(forward), '--reads', str(reverse), '--paired_in', '--out2',
         '--aligned', str(prefix), '--fastx', '--sam', '--blast', '1', '--zip-out', 'true',
         '--num_alignments', '10', '-e', '1e-5', '--threads', '4'],
        check=True, capture_output=True, text=True)
    return (Path(str(prefix) + '.sam.gz'),
            Path(str(prefix) + '_fwd.fq.gz'), Path(str(prefix) + '_rev.fq.gz'))


def sam_records(lines):
    return [line.split('\t') for line in lines if line.strip() and not line.startswith('@')]


def test_sortmerna_writes_no_pairing_flags(casava_mapping):
    sam, _, _ = casava_mapping
    with gzip.open(sam, 'rt') as handle:
        flags = {int(record[1]) for record in sam_records(handle)}
    assert flags <= {0, 16}


def test_every_record_is_tagged_with_its_mate(casava_mapping):
    sam, forward, reverse = casava_mapping
    with gzip.open(sam, 'rt') as handle:
        raw = sam_records(handle)
    with gzip.open(sam, 'rt') as handle:
        tagged = sam_records(sam_mates.tag_mates(handle, sam_mates.iter_separate_pairs(forward, reverse)))

    assert len(tagged) == len(raw)
    mates = collections.Counter(sam_mates.mate_of(record[1]) for record in tagged)
    assert set(mates) == {1, 2}
    # Mates of most pairs align, so neither mate dominates the records.
    assert min(mates.values()) > 0.4 * len(tagged)


def test_paired_fastx_output_is_two_files_in_step_and_in_input_order(casava_mapping):
    _, forward, reverse = casava_mapping
    input_order = [core_id for core_id, _, _ in iter_example_pairs()]
    output = [(r1[0], r2[0]) for r1, r2 in sam_mates.iter_separate_pairs(forward, reverse)]

    assert all(r1 == r2 for r1, r2 in output)  # CASAVA: mates share the ID
    positions = [input_order.index(r1) for r1, _ in output]
    assert positions == sorted(positions)


AMBIGUOUS = 'KRYMSWBDHV'


def with_ambiguous_bases(pairs):
    """Every other pair gets an IUPAC code in R1 and a lower-case base in R2."""
    for index, (core_id, (seq1, qual1), (seq2, qual2)) in enumerate(pairs):
        if index % 2 == 0:
            code = AMBIGUOUS[(index // 2) % len(AMBIGUOUS)]
            seq1 = seq1[:40] + code + seq1[41:]
            seq2 = seq2[:60] + seq2[60].lower() + seq2[61:]
        yield core_id, (seq1, qual1), (seq2, qual2)


@pytest.mark.parametrize('style', ['casava', 'legacy'])
def test_ambiguous_read_bases_are_written_as_n_and_still_tagged(tmp_path, style):
    """SortMeRNA writes IUPAC codes as N in SEQ and keeps them in its FASTQ.

    MATAM compares SEQ with the FASTQ reads to tell the mates apart, so the
    comparison must read any base other than A, C, G, T as N.
    """
    reference = write_reference_fasta(tmp_path / 'refs.fasta')
    forward, reverse = write_pe_reads(
        tmp_path, style, pairs=with_ambiguous_bases(iter_example_pairs(every=10)))
    prefix = tmp_path / 'aligned'
    subprocess.run(
        [SORTMERNA, '--ref', str(reference), '--workdir', str(tmp_path / 'wd'),
         '--reads', str(forward), '--reads', str(reverse), '--paired_in', '--out2',
         '--aligned', str(prefix), '--fastx', '--sam', '--blast', '1',
         '--num_alignments', '10', '-e', '1e-5', '--threads', '2'],
        check=True, capture_output=True, text=True)
    sam = Path(str(prefix) + '.sam')
    aligned_reads = (Path(str(prefix) + '_fwd.fq'), Path(str(prefix) + '_rev.fq'))

    with open(sam) as handle:
        raw = sam_records(handle)
    sequences = {seq for r1, r2 in sam_mates.iter_separate_pairs(*aligned_reads)
                 for _, seq in (r1, r2)}
    assert any(set(seq) & set(AMBIGUOUS) for seq in sequences)  # kept in the FASTQ
    assert any(seq != seq.upper() for seq in sequences)
    assert {base for record in raw for base in record[9]} <= set('ACGTN')  # N in SEQ

    with open(sam) as handle:
        tagged = sam_records(sam_mates.tag_mates(handle, sam_mates.iter_separate_pairs(*aligned_reads)))
    assert len(tagged) == len(raw)
    assert {sam_mates.mate_of(record[1]) for record in tagged} == {1, 2}
