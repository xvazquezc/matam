"""Integration test for reusing one full SortMeRNA mapping across subsamples."""

from collections import Counter
import gzip
from pathlib import Path
import random
import re
import shutil
import subprocess

import pytest


SORTMERNA = shutil.which("sortmerna")
pytestmark = pytest.mark.skipif(
    SORTMERNA is None,
    reason="requires SortMeRNA 7 on PATH",
)


def canonical_pair_id(read_id):
    """Remove the conventional mate suffix used by this test fixture."""
    return re.sub(r"/[12]$", "", read_id.split()[0])


def write_fastq(path, records):
    with path.open("w") as handle:
        for read_id, sequence in records:
            handle.write("@{}\n{}\n+\n{}\n".format(read_id, sequence, "I" * len(sequence)))


def open_text(path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def fastq_records(path):
    records = []
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()
            assert header.startswith("@") and sequence and separator == "+\n" and quality
            records.append((header[1:].strip(), sequence.strip(), quality.strip()))
    return records


def mapping_records(path):
    with open_text(path) as handle:
        return Counter(
            line.rstrip()
            for line in handle
            if line.strip() and not line.startswith("@")
        )


def map_pairs(reference, index_dir, forward, reverse, output_dir, label,
              compress_output):
    """Run the MATAM read-mapping SortMeRNA contract on one paired input."""
    output_prefix = output_dir / label
    workdir = output_dir / "{}_workdir".format(label)
    command = [
        SORTMERNA,
        "--ref", str(reference),
        "--idx-dir", str(index_dir),
        "--workdir", str(workdir),
        "--reads", str(forward), "--reads", str(reverse),
        "--paired_in",
        "--aligned", str(output_prefix),
        "--fastx", "--sam", "--blast", "1",
        "--num_alignments", "1", "-e", "1e-03", "--threads", "1",
    ]
    if compress_output:
        command.extend(["--zip-out", "true"])
    subprocess.run(command, check=True, capture_output=True, text=True)

    # SortMeRNA 7 normalises its paired FASTX output suffix to ``.fq``.
    suffix = ".gz" if compress_output else ""
    sam_path = output_prefix.with_suffix(".sam" + suffix)
    fastq_path = output_prefix.with_suffix(".fq" + suffix)
    assert sam_path.is_file(), "SortMeRNA did not produce {}".format(sam_path)
    assert fastq_path.is_file(), "SortMeRNA did not produce {}".format(fastq_path)
    return sam_path, fastq_path


@pytest.mark.parametrize("compress_output", [False, True])
def test_subsetting_a_full_mapping_matches_mapping_the_selected_pairs(
        tmp_path, compress_output):
    """Selected SAM/FASTQ records must not depend on when the subset is made.

    This is the required invariant before replacing per-percentage SortMeRNA
    runs with one full mapping followed by pair-level subsetting of its outputs.
    """
    random_source = random.Random(17)
    reference_sequence = "".join(random_source.choice("ACGT") for _ in range(1200))
    reference = tmp_path / "reference.fasta"
    reference.write_text(">reference\n{}\n".format(reference_sequence))

    forward_records = []
    reverse_records = []
    for index, start in enumerate((20, 180, 340, 500, 660, 820)):
        pair_id = "pair{:02d}".format(index)
        forward_records.append((pair_id + "/1", reference_sequence[start:start + 100]))
        reverse_records.append((pair_id + "/2", reference_sequence[start + 25:start + 125]))

    full_forward = tmp_path / "full_R1.fastq"
    full_reverse = tmp_path / "full_R2.fastq"
    write_fastq(full_forward, forward_records)
    write_fastq(full_reverse, reverse_records)

    selected_pair_ids = {"pair00", "pair02", "pair05"}
    subset_forward = tmp_path / "subset_R1.fastq"
    subset_reverse = tmp_path / "subset_R2.fastq"
    write_fastq(
        subset_forward,
        [record for record in forward_records
         if canonical_pair_id(record[0]) in selected_pair_ids],
    )
    write_fastq(
        subset_reverse,
        [record for record in reverse_records
         if canonical_pair_id(record[0]) in selected_pair_ids],
    )

    index_dir = tmp_path / "index"
    full_sam, full_fastq = map_pairs(
        reference, index_dir, full_forward, full_reverse, tmp_path, "full",
        compress_output)
    subset_sam, subset_fastq = map_pairs(
        reference, index_dir, subset_forward, subset_reverse, tmp_path, "subset",
        compress_output)

    full_selected_sam = Counter(
        record for record in mapping_records(full_sam)
        if canonical_pair_id(record.split("\t", 1)[0]) in selected_pair_ids
    )
    assert full_selected_sam
    assert full_selected_sam == mapping_records(subset_sam)

    full_selected_fastq = Counter(
        record for record in fastq_records(full_fastq)
        if canonical_pair_id(record[0]) in selected_pair_ids
    )
    assert {canonical_pair_id(record[0]) for record in full_selected_fastq} == selected_pair_ids
    assert full_selected_fastq == Counter(fastq_records(subset_fastq))
