"""Characterization of compute_contigs_compatibility.py (--contigs_binning)."""
import subprocess
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).parent.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS_DIR))
import compute_contigs_compatibility as ccc  # noqa: E402
sys.path.pop(0)

# A reference stretch; contigs are pieces of it placed at their 1-based position.
REF = 'ACGTTGCAAGCTTCGATCGGATCCATGCAGTCAGTTACGG'


def contig(name, start, end, mismatch_at=None, ref='ref1'):
    """SAM fields of a contig matching REF[start:end+1] (0-based, inclusive end)."""
    seq = list(REF[start:end + 1])
    if mismatch_at is not None:
        i = mismatch_at - start
        seq[i] = 'A' if seq[i] != 'A' else 'C'
    return [name, '0', ref, str(start + 1), '60', '{}M'.format(len(seq)),
            '*', '0', '0', ''.join(seq), '*']


def alignment(*args, **kwargs):
    return ccc.SamAlignment(contig(*args, **kwargs))


def test_sam_alignment_coordinates_are_0_based_and_inclusive():
    a = ccc.SamAlignment(['c', '0', 'r', '11', '60', '2S8M1I3M2D4M1S',
                          '*', '0', '0', 'A' * 19, '*'])
    assert (a.subject_start, a.subject_end) == (10, 26)   # 8 + 3 + 2 + 4 reference bases
    assert (a.query_start, a.query_end) == (2, 17)        # 2 soft-clipped, 16 aligned, 1 clipped
    assert (a.left_softclip_num, a.right_softclip_num, a.indel_num) == (2, 1, 3)


def test_identical_overlap_is_compatible():
    assert ccc.return_compatibility_status(alignment('a', 0, 19), alignment('b', 10, 29)) == 1


def test_a_mismatch_in_the_overlap_is_incompatible():
    assert ccc.return_compatibility_status(
        alignment('a', 0, 19), alignment('b', 10, 29, mismatch_at=15)) == 2


def test_a_mismatch_outside_the_overlap_is_ignored():
    assert ccc.return_compatibility_status(
        alignment('a', 0, 19), alignment('b', 10, 29, mismatch_at=25)) == 1


def test_contained_contig_is_compared_on_its_whole_length():
    assert ccc.return_compatibility_status(alignment('a', 0, 19), alignment('b', 5, 9)) == 1
    assert ccc.return_compatibility_status(
        alignment('a', 0, 19), alignment('b', 5, 9, mismatch_at=9)) == 2


def test_two_base_overlap_is_compared():
    assert ccc.return_compatibility_status(alignment('a', 0, 9), alignment('b', 8, 17)) == 1
    assert ccc.return_compatibility_status(
        alignment('a', 0, 9), alignment('b', 8, 17, mismatch_at=8)) == 2


def test_one_base_overlap_counts_as_no_overlap():
    """Current behaviour: overlap size is end - start with an inclusive end,
    so contigs sharing exactly one reference base are not compared."""
    assert ccc.return_compatibility_status(
        alignment('a', 0, 9), alignment('b', 9, 18, mismatch_at=9)) == 0


def test_disjoint_contigs_do_not_overlap():
    assert ccc.return_compatibility_status(alignment('a', 0, 9), alignment('b', 20, 29)) == 0


def test_matrix_scan_stops_at_the_first_contig_that_does_not_overlap():
    contigs = [alignment('a', 0, 9), alignment('b', 20, 29), alignment('c', 25, 34)]
    assert ccc.compute_contigs_compatibility_matrix(contigs) == [
        [0, 0, 0],
        [0, 0, 1],
        [0, 0, 0]]


def test_early_exit_skips_later_overlaps_after_a_one_base_contig():
    """Current behaviour: a one-base contig inside a is reported as no
    overlap, which ends the scan for a, so a's conflict with c is never
    recorded. (MATAM contigs are far longer than one base.)"""
    a, b, c = alignment('a', 0, 19), alignment('b', 10, 10), alignment('c', 12, 30, mismatch_at=15)
    assert ccc.return_compatibility_status(a, c) == 2
    assert ccc.compute_contigs_compatibility_matrix([a, b, c]) == [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]]


def test_empty_list_gives_an_empty_matrix():
    assert ccc.compute_contigs_compatibility_matrix([]) == []


def test_compatible_contigs_share_a_bin():
    matrix = [[0, 1, 1],
              [0, 0, 1],
              [0, 0, 0]]
    assert ccc.compute_bin_list(matrix) == [[0, 1, 2]]


def test_incompatible_contig_starts_a_new_bin():
    matrix = [[0, 2],
              [0, 0]]
    assert ccc.compute_bin_list(matrix) == [[0], [1]]


def test_contig_overlapping_no_bin_member_starts_a_new_bin():
    matrix = [[0, 0],
              [0, 0]]
    assert ccc.compute_bin_list(matrix) == [[0], [1]]


def test_secondary_incompatibility_keeps_compatible_contigs_apart():
    """b agrees with a, but a later contig c agrees with a and not with b:
    b goes to its own bin and c joins a."""
    contigs = [alignment('a', 0, 19), alignment('b', 10, 29), alignment('c', 15, 34, mismatch_at=25)]
    matrix = ccc.compute_contigs_compatibility_matrix(contigs)
    assert matrix == [[0, 1, 1],
                      [0, 0, 2],
                      [0, 0, 0]]
    assert ccc.compute_bin_list(matrix) == [[0, 2], [1]]


def run_script(sam_text):
    return subprocess.run([sys.executable, str(SCRIPTS_DIR / 'compute_contigs_compatibility.py')],
                          input=sam_text, capture_output=True, text=True, check=True).stdout


def test_script_keeps_a_reference_whole_when_its_contigs_agree():
    lines = [contig('1', 0, 19), contig('2', 10, 29)]
    out = run_script(''.join('\t'.join(t) + '\n' for t in lines))
    assert out.splitlines() == ['@HD\tVN:1.0', '@SQ\tSN:ref1\tLN:2147483647',
                                '\t'.join(lines[0]), '\t'.join(lines[1])]


def test_script_splits_a_reference_into_bins_of_compatible_contigs():
    lines = [contig('1', 0, 19), contig('2', 10, 29, mismatch_at=15)]
    out = run_script(''.join('\t'.join(t) + '\n' for t in lines)).splitlines()
    assert out[1] == '@SQ\tSN:ref1_0\tLN:2147483647'
    assert out[2].split('\t')[:3] == ['1', '0', 'ref1_0']
    assert out[3] == '@SQ\tSN:ref1_1\tLN:2147483647'
    assert out[4].split('\t')[:3] == ['2', '0', 'ref1_1']


def test_script_accepts_an_empty_sam():
    """No contig aligned: only the header (this used to raise IndexError)."""
    assert run_script('') == '@HD\tVN:1.0\n'
