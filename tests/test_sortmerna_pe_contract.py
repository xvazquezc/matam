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
