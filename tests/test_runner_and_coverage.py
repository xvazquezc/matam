import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).parent.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS_DIR))
import runner  # noqa: E402
import sample_sam_by_coverage  # noqa: E402
sys.path.pop(0)


def test_pipeline_failure_is_masked_without_pipefail():
    # Documents the default /bin/sh behaviour that pipefail exists to avoid.
    assert runner.logged_call('false | true') == 0


def test_pipefail_reports_a_failure_anywhere_in_the_pipeline():
    assert runner.logged_call('false | true', pipefail=True) != 0
    assert runner.logged_call('gzip -cd /nonexistent.gz | wc -l', pipefail=True) != 0
    assert runner.logged_call('true | true', pipefail=True) == 0


def test_reference_coverage_counts_each_aligned_reference_position():
    # Two reads on a 10 bp reference: 1-based POS 1 with 4M, POS 3 with 2M1D2M.
    alignments = [['r1', '0', 'ref', '1', '255', '4M'],
                  ['r2', '0', 'ref', '3', '255', '2M1D2M']]

    coverage = sample_sam_by_coverage.compute_ref_coverage(alignments, 10)

    assert coverage.tolist() == [1, 1, 2, 2, 1, 1, 1, 0, 0, 0]
