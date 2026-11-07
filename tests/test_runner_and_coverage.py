import logging
import sys
import time
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).parent.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS_DIR))
import runner  # noqa: E402
import sample_sam_by_coverage  # noqa: E402
sys.path.pop(0)


def test_pipeline_failure_is_masked_only_when_pipefail_is_disabled():
    # Documents the plain /bin/sh behaviour that pipefail exists to avoid.
    assert runner.logged_call('false | true', pipefail=False) == 0


def test_pipefail_is_the_default_and_reports_a_failure_anywhere_in_the_pipeline():
    assert runner.logged_call('false | true') != 0
    assert runner.logged_call('gzip -cd /nonexistent.gz | wc -l') != 0
    assert runner.logged_call('true | true') == 0


def test_check_call_exits_on_a_failure_inside_a_pipeline():
    with pytest.raises(SystemExit):
        runner.logged_check_call('false | cat')


class SlowFirstHandler(logging.Handler):
    """Collects messages; the first one is slow, as a busy log handler can be."""

    def __init__(self):
        super().__init__(logging.INFO)
        self.messages = []

    def emit(self, record):
        if not self.messages:
            time.sleep(0.5)
        self.messages.append(record.getMessage())


def test_verbose_output_is_logged_to_the_end():
    # The output used to be read only while the process was running: output
    # written while the first chunk was being logged, by a process that then
    # exited, was lost.
    handler = SlowFirstHandler()
    root = logging.getLogger()
    old_level = root.level
    root.addHandler(handler)
    root.setLevel(logging.INFO)
    try:
        command = 'printf start; sleep 0.1; printf " END-é"'
        assert runner.logged_call(command, verbose=True) == 0
    finally:
        root.removeHandler(handler)
        root.setLevel(old_level)
    assert ''.join(handler.messages) == 'start END-é'


def test_reference_coverage_counts_each_aligned_reference_position():
    # Two reads on a 10 bp reference: 1-based POS 1 with 4M, POS 3 with 2M1D2M.
    alignments = [['r1', '0', 'ref', '1', '255', '4M'],
                  ['r2', '0', 'ref', '3', '255', '2M1D2M']]

    coverage = sample_sam_by_coverage.compute_ref_coverage(alignments, 10)

    assert coverage.tolist() == [1, 1, 2, 2, 1, 1, 1, 0, 0, 0]
