"""SortMeRNA contract behind reusing one mapping for several read subsets.

Slicing a mapping of a larger read set gives the alignments that mapping the
subset alone gives, except near SortMeRNA's minimum alignment score: that score
grows with the total length of the reads mapped together (E-value), and
``--num_alignments`` keeps the first alignments that reach it. MATAM's
``--read_alignments`` relies on the first part and documents the second.
"""

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
              compress_output, num_alignments=1, evalue="1e-03"):
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
        "--num_alignments", str(num_alignments), "-e", evalue, "--threads", "1",
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
    """Away from the score threshold, slicing the full mapping gives the subset's own.

    The reads match the reference exactly, so their scores are far above
    SortMeRNA's minimum score for either read set.
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


def test_alignments_near_the_minimum_score_depend_on_the_reads_mapped_together(tmp_path):
    """The minimum score rises with the reads mapped together, so a weak alignment
    found when a few reads are mapped alone is lost in a larger read set."""
    random_source = random.Random(5)

    def random_sequence(length):
        return "".join(random_source.choice("ACGT") for _ in range(length))

    reference_sequence = random_sequence(1500)
    reference = tmp_path / "reference.fasta"
    reference.write_text(">reference\n{}\n".format(reference_sequence))
    # R1 of probe k starts with k reference bases (score about 2k), the rest is random.
    probes = [("probe{:02d}".format(k),
               reference_sequence[100 + 20 * k:100 + 21 * k] + random_sequence(100 - k),
               random_sequence(100))
              for k in range(20, 40)]
    unrelated = [("other{:05d}".format(i), random_sequence(100), random_sequence(100))
                 for i in range(3000)]

    def aligned_reads_and_minimum_score(label, pairs):
        forward, reverse = tmp_path / (label + "_R1.fq"), tmp_path / (label + "_R2.fq")
        write_fastq(forward, [(name + "/1", r1) for name, r1, _ in pairs])
        write_fastq(reverse, [(name + "/2", r2) for name, _, r2 in pairs])
        sam_path, _ = map_pairs(reference, tmp_path / "index", forward, reverse, tmp_path,
                                label, compress_output=False, num_alignments=10, evalue="1e-05")
        log = (tmp_path / (label + ".log")).read_text()
        minimum = int(re.search(r"Minimal SW score based on E-value = (\d+)", log).group(1))
        return set(mapping_records_by_read(sam_path)), minimum

    alone, alone_minimum = aligned_reads_and_minimum_score("alone", probes)
    together, together_minimum = aligned_reads_and_minimum_score("together", probes + unrelated)

    assert alone_minimum < together_minimum
    assert together < alone, "a probe aligned alone should be lost among more reads"


def mapping_records_by_read(path):
    return Counter(record.split("\t", 1)[0] for record in mapping_records(path))
