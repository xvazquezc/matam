"""SortMeRNA step completeness: reuse complete output, clear incomplete or outdated output (#13)."""

import json
import os
import shutil
import sys
import textwrap
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).parent.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS_DIR))
import sortmerna_step  # noqa: E402
sys.path.pop(0)

# Behaves like SortMeRNA 7 for what matters here: refuses to start when run
# state from an earlier run is present in the work directory, keeps run state
# in <workdir>/kvdb, writes <aligned>.sam / .blast / .log, and can fail midway.
FAKE_SORTMERNA = textwrap.dedent('''\
    import os, sys
    args = sys.argv[1:]
    value = lambda name: args[args.index(name) + 1]
    workdir, aligned, reads = value('--workdir'), value('--aligned'), value('--reads')
    with open(os.environ['FAKE_SMR_CALLS'], 'a') as calls:
        calls.write(' '.join(args) + '\\n')
    kvdb = os.path.join(workdir, 'kvdb')
    if os.path.isdir(kvdb) and os.listdir(kvdb):
        sys.exit('kvdb resume aborted: inputs mismatch')
    os.makedirs(kvdb, exist_ok=True)
    open(os.path.join(kvdb, 'state'), 'w').write(reads)
    sam = open(aligned + '.sam', 'w')
    sam.write('@HD\\tVN:1.0\\n')
    if os.environ.get('FAKE_SMR_FAIL'):
        sys.exit(1)  # interrupted: partial SAM and run state left behind
    sam.write('read\\t0\\tref\\t1\\n' * len(open(reads).read()))
    open(aligned + '.blast', 'w').write('read\\tref\\n')
    open(aligned + '.log', 'w').write('done\\n')
    ''')

OUTPUTS = [('', ('.sam',)), ('', ('.blast',))]


@pytest.fixture
def step(tmp_path, monkeypatch):
    fake = tmp_path / 'fake_sortmerna.py'
    fake.write_text(FAKE_SORTMERNA)
    calls = tmp_path / 'calls.txt'
    calls.write_text('')
    monkeypatch.setenv('FAKE_SMR_CALLS', str(calls))
    run_dir = tmp_path / 'run'
    (run_dir / 'workdir').mkdir(parents=True)
    workdir = run_dir / 'sortmerna_scaffolding'
    reads = run_dir / 'workdir' / 'contigs.fasta'
    reads.write_text('>c1\nACGT\n')
    ref = tmp_path / 'refs.fasta'
    ref.write_text('>r1\nACGTACGT\n')

    class Step:
        def __init__(self):
            self.tmp_path, self.calls_file = tmp_path, calls
            self.set_paths(run_dir)

        def set_paths(self, base):
            self.run_dir = base
            self.workdir = str(base / 'sortmerna_scaffolding')
            self.aligned = str(base / 'workdir' / 'contigs.sortmerna_vs_refs')
            self.reads = str(base / 'workdir' / 'contigs.fasta')

        def command(self, evalue='1e-05'):
            return '{} {} --ref {} --workdir {} --reads {} --aligned {} -e {}'.format(
                sys.executable, fake, ref, self.workdir, self.reads, self.aligned, evalue)

        def run(self, run_options='', prepare=None, **kwargs):
            return sortmerna_step.run(self.command(**kwargs), self.workdir, self.aligned,
                                      [self.reads, str(ref)], OUTPUTS, run_options,
                                      prepare=prepare)

        def completed_outputs(self):
            return sortmerna_step.completed_outputs(self.command(), self.workdir, self.aligned,
                                                    [self.reads, str(ref)], OUTPUTS)

        def calls(self):
            return len(self.calls_file.read_text().splitlines())

    return Step()


def test_first_run_maps_records_completion_and_drops_run_state(step):
    outputs = step.run()

    assert outputs == [step.aligned + '.sam', step.aligned + '.blast']
    assert step.calls() == 1
    record = json.loads(Path(step.aligned + sortmerna_step.RECORD_SUFFIX).read_text())
    assert [o['path'] for o in record['outputs']] == ['contigs.sortmerna_vs_refs.sam',
                                                       'contigs.sortmerna_vs_refs.blast']
    assert not os.path.exists(os.path.join(step.workdir, 'kvdb'))


def test_complete_output_for_the_same_inputs_is_reused(step):
    step.run()
    outputs = step.run(run_options=' --threads 8 -v')  # threads/verbosity don't matter

    assert step.calls() == 1
    assert outputs == [step.aligned + '.sam', step.aligned + '.blast']


def test_prepare_runs_only_when_sortmerna_runs_and_its_options_are_not_recorded(step):
    # matam_assembly finds or builds the reference index in prepare: not
    # needed for a reused output, and its location doesn't change the output.
    prepared = []

    def prepare():
        prepared.append(len(prepared))
        return ' --idx-dir {}/idx{} --index 0'.format(step.tmp_path, len(prepared))

    step.run(prepare=prepare)
    assert prepared == [0]
    assert '--idx-dir {}/idx1 --index 0'.format(step.tmp_path) in step.calls_file.read_text()
    assert '--idx-dir' not in Path(step.aligned + sortmerna_step.RECORD_SUFFIX).read_text()

    step.run(prepare=prepare)              # complete output reused
    assert prepared == [0] and step.calls() == 1


@pytest.mark.parametrize('change', ['input', 'command', 'truncated output', 'missing output'])
def test_outdated_or_damaged_output_is_cleared_and_remapped(step, change):
    step.run()
    if change == 'input':
        Path(step.reads).write_text('>c1\nACGTACGTAAAA\n')  # rebuilt contigs
        kwargs = {}
    elif change == 'command':
        kwargs = {'evalue': '1e-10'}
    elif change == 'truncated output':
        Path(step.aligned + '.sam').write_text('@HD\tVN:1.0\n')
        kwargs = {}
    else:
        os.remove(step.aligned + '.blast')
        kwargs = {}

    step.run(**kwargs)

    assert step.calls() == 2


def test_interrupted_run_is_cleared_before_running_again(step, monkeypatch):
    monkeypatch.setenv('FAKE_SMR_FAIL', '1')
    with pytest.raises(SystemExit):
        step.run()
    assert os.listdir(os.path.join(step.workdir, 'kvdb'))      # stale run state left
    assert not os.path.exists(step.aligned + sortmerna_step.RECORD_SUFFIX)

    monkeypatch.delenv('FAKE_SMR_FAIL')
    step.run()  # the fake aborts like SortMeRNA if the old kvdb is still there

    assert step.calls() == 2
    assert sortmerna_step.is_complete(step.aligned, step.command(), [step.reads, str(step.tmp_path / 'refs.fasta')],
                                      step.workdir)


def test_clearing_keeps_other_files_sharing_the_prefix(step):
    step.run()
    downstream = Path(step.aligned + '.blast.best_only.tab')
    downstream.write_text('kept')
    Path(step.reads).write_text('>c1\nTTTT\n')

    step.run()

    assert downstream.read_text() == 'kept'


def test_a_moved_run_directory_is_still_complete(step):
    step.run()
    moved = step.tmp_path / 'moved'
    shutil.copytree(step.run_dir, moved, copy_function=shutil.copy2)  # keeps file times
    step.set_paths(moved)

    step.run()

    assert step.calls() == 1


def test_skipped_step_returns_complete_outputs(step):
    step.run()
    assert step.completed_outputs() == [step.aligned + '.sam', step.aligned + '.blast']


def test_skipped_step_without_a_record_uses_existing_outputs_with_a_warning(step, caplog):
    step.run()
    os.remove(step.aligned + sortmerna_step.RECORD_SUFFIX)  # output of an older MATAM

    assert step.completed_outputs() == [step.aligned + '.sam', step.aligned + '.blast']
    assert 'No completion record' in caplog.text


@pytest.mark.parametrize('damage', ['outdated', 'no output'])
def test_skipped_step_with_incomplete_output_is_an_error(step, damage):
    if damage == 'outdated':
        step.run()
        Path(step.reads).write_text('>c1\nGGGG\n')
    with pytest.raises(sortmerna_step.IncompleteStepError):
        step.completed_outputs()
