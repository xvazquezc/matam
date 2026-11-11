"""SortMeRNA reference indexes: reuse only an index that matches its fasta and parameters.

SortMeRNA 7 reuses any non-empty index files named after the reference file
name: a stale index (fasta rewritten) or one built with other -L/--max_pos/-m
is used silently. sortmerna_index records what an index was built from and
checks it before use; the index is built in the database directory when it is
writable (once, under a lock), otherwise privately in the run.
"""

import json
import re
import logging
import os
import shutil
import stat
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
from pe_fixtures import EXAMPLE_REFS, SCRIPTS_DIR, tool_env  # noqa: E402

sys.path.insert(0, str(SCRIPTS_DIR))
import sortmerna_index  # noqa: E402
sys.path.pop(0)

SORTMERNA = shutil.which('sortmerna', path=tool_env()['PATH'])
pytestmark = pytest.mark.skipif(SORTMERNA is None, reason='missing: sortmerna')
PARAMS = sortmerna_index.index_params(2000)
# Parameters of a database without stored ones (indexed by an older MATAM).
LEGACY = sortmerna_index.index_params(sortmerna_index.LEGACY_MEMORY_MB)


def write_refs(path, records=None):
    """The 16 example references (or the given ones) as a plain fasta."""
    if records is None:
        records = []
        for line in open(str(EXAMPLE_REFS)):
            if line.startswith('>'):
                records.append((line[1:].split()[0], ''))
            else:
                records[-1] = (records[-1][0], records[-1][1] + line.strip())
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(''.join('>{}\n{}\n'.format(name, seq) for name, seq in records))
    return path


@pytest.fixture
def db(tmp_path):
    ref = write_refs(tmp_path / 'db' / 'mini.clustered.fasta')
    return ref, tmp_path / 'db' / 'mini.clustered.smr' / 'idx'


def built(ref, idx_dir, params=PARAMS):
    with sortmerna_index.locked(str(idx_dir.parent)):
        sortmerna_index.build(SORTMERNA, str(ref), str(idx_dir), params)


def check(ref, idx_dir, params=PARAMS):
    return sortmerna_index.check(str(idx_dir), str(ref), SORTMERNA, params)


def test_build_records_the_fasta_and_parameters_and_leaves_no_temporary_files(db):
    ref, idx_dir = db
    built(ref, idx_dir)

    record = json.loads((idx_dir / 'matam_index.json').read_text())
    assert record['reference']['name'] == 'mini.clustered.fasta'
    assert record['reference']['size'] == ref.stat().st_size
    assert record['params'] == {'L': 18, 'max_pos': sortmerna_index.MAX_POS, 'm': 2000}
    assert record['sortmerna_version'] == sortmerna_index.sortmerna_version(SORTMERNA)
    assert {entry['name'] for entry in record['files']} == (
        set(os.listdir(str(idx_dir))) - {'matam_index.json'})
    # Only the index: no kvdb/out from '--task 5', no temporary build directory.
    assert sorted(os.listdir(str(idx_dir.parent))) == ['.matam_index.lock', 'idx']
    assert check(ref, idx_dir) == (True, '')


def test_sortmerna_maps_with_the_built_index_without_reindexing(db, tmp_path):
    ref, idx_dir = db
    built(ref, idx_dir)
    reads = tmp_path / 'reads.fa'
    reads.write_text('>r\n{}\n'.format(ref.read_text().splitlines()[1][100:200]))
    completed = subprocess.run(
        [SORTMERNA, '--ref', str(ref), '--idx-dir', str(idx_dir), '--index', '0',
         '--workdir', str(tmp_path / 'wd'), '--reads', str(reads),
         '--aligned', str(tmp_path / 'al'), '--blast', '1', '-v'],
        env=tool_env(), capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout[-2000:]
    assert 'Skipping indexing' in completed.stdout
    assert (tmp_path / 'al.blast').read_text().startswith('r\t')


def test_fasta_rewritten_with_the_same_size_makes_the_index_stale(db):
    # SortMeRNA itself uses such an index silently.
    ref, idx_dir = db
    built(ref, idx_dir)
    text = ref.read_text()
    ref.write_text(text.replace('ACGT', 'ACGA', 1))
    assert ref.stat().st_size == len(text)

    assert check(ref, idx_dir) == (
        False, 'the fasta changed since the index was built (checksum differs)')


def test_fasta_with_another_size_makes_the_index_stale(db):
    # A shorter fasta with the old index made SortMeRNA segfault.
    ref, idx_dir = db
    built(ref, idx_dir)
    ref.write_text(''.join(ref.read_text().splitlines(True)[:4]))
    assert check(ref, idx_dir) == (
        False, 'the fasta changed since the index was built (size differs)')


def test_a_copied_fasta_with_the_same_content_keeps_the_index_and_refreshes_its_record(db):
    ref, idx_dir = db
    built(ref, idx_dir)
    mtime = ref.stat().st_mtime_ns
    os.utime(str(ref), ns=(mtime + 10**9, mtime + 10**9))

    assert check(ref, idx_dir) == (True, '')
    record = json.loads((idx_dir / 'matam_index.json').read_text())
    assert record['reference']['mtime_ns'] == mtime + 10**9


@pytest.mark.parametrize('other', [{'m': 3000}, {'L': 20}, {'max_pos': 1000}])
def test_other_indexing_parameters_make_the_index_stale(db, other):
    ref, idx_dir = db
    built(ref, idx_dir)
    params = dict(PARAMS, **other)
    usable, reason = check(ref, idx_dir, params)
    assert not usable and reason.startswith('the index was built with other parameters')


def test_missing_or_truncated_index_files_make_the_index_stale(db):
    ref, idx_dir = db
    built(ref, idx_dir)
    kmer = next(idx_dir.glob('*.kmer_0.dat'))
    kmer.write_bytes(kmer.read_bytes()[:10])
    usable, reason = check(ref, idx_dir)
    assert not usable and reason.startswith('the index is incomplete')
    kmer.unlink()
    assert not check(ref, idx_dir)[0]


def legacy_index(ref, idx_dir, *options):
    """An index as MATAM built it before index records: plain 'sortmerna --task 5'."""
    subprocess.run([SORTMERNA, '--ref', str(ref), '--workdir', str(idx_dir.parent), '--task', '5',
                    *options], check=True, capture_output=True, env=tool_env())


def test_an_index_without_a_matam_record_is_checked_against_sortmerna_statistics(db, tmp_path, caplog):
    # Existing databases: rebuilding a full SILVA index would take hours.
    ref, idx_dir = db
    legacy_index(ref, idx_dir)
    defaults = LEGACY
    assert check(ref, idx_dir, defaults) == (False, sortmerna_index.NO_RECORD)
    caplog.set_level(logging.INFO, logger='sortmerna_index')
    local = tmp_path / 'run' / 'idx'

    assert sortmerna_index.ensure(SORTMERNA, str(ref), str(idx_dir), str(local), require=False) == str(idx_dir)

    assert 'no MATAM record but matches the fasta and parameters' in caplog.text
    assert 'Building' not in caplog.text and not local.exists()
    record = json.loads((idx_dir / 'matam_index.json').read_text())
    assert record['params'] == defaults and 'adopted' in record
    assert check(ref, idx_dir, defaults) == (True, '')     # adopted once


def test_a_split_index_from_an_older_matam_is_adopted_with_its_parameters(db, tmp_path):
    # Full SILVA with SortMeRNA's default -m 3072 gives several parts.
    ref, idx_dir = db
    legacy_index(ref, idx_dir, '-m', '0.1')   # 25 kB fasta: 3 parts
    assert 'part_num: 3' in next(idx_dir.glob('*.idx_stats.yaml')).read_text()
    assert sortmerna_index.ensure(SORTMERNA, str(ref), str(idx_dir), str(tmp_path / 'idx'), require=False) == str(idx_dir)
    assert json.loads((idx_dir / 'matam_index.json').read_text())['params'] == LEGACY


@pytest.mark.parametrize('change, reason', [
    ('fasta size', 'size differs'),
    ('sequences', 'sequences differ'),
    ('seed length', 'another seed length'),
    ('truncated', 'incomplete'),
    ('several parts', 'partitioning'),
    ('capped at 1000', 'another --max_pos'),
])
def test_an_index_without_a_record_that_does_not_match_is_rebuilt(db, tmp_path, caplog, change, reason):
    ref, idx_dir = db
    legacy_index(ref, idx_dir, *(['-L', '20'] if change == 'seed length' else
                                 ['-m', '0.1'] if change == 'several parts' else []))
    if change == 'fasta size':
        ref.write_text(''.join(ref.read_text().splitlines(True)[:8]))
    elif change == 'sequences':
        # Same size, one sequence fewer: the second header becomes sequence.
        lines = ref.read_text().splitlines(True)
        merged = lines[1].rstrip('\n') + 'A' * (len(lines[2]) + 1) + lines[3]
        ref.write_text(lines[0] + merged + ''.join(lines[4:]))
        assert ref.stat().st_size == sum(len(line) for line in lines)
    elif change == 'several parts':
        # -m changes the alignments through the number of parts (alignments are
        # selected per part). A database whose parameters say -m 10000 can't
        # use a split index from an older MATAM (-m 3072).
        sortmerna_index.write_params(str(idx_dir.parent), sortmerna_index.index_params(10000))
    elif change == 'capped at 1000':
        # What 'sortmerna --max_pos 1000' (the value in its help) gives on SILVA.
        stats = next(idx_dir.glob('*.idx_stats.yaml'))
        stats.write_text(re.sub(r'max_positions: \d+', 'max_positions: 1000', stats.read_text()))
    elif change == 'truncated':
        pos = next(idx_dir.glob('*.pos_0.dat'))
        pos.write_bytes(pos.read_bytes()[:100])
    caplog.set_level(logging.INFO, logger='sortmerna_index')

    sortmerna_index.ensure(SORTMERNA, str(ref), str(idx_dir), str(tmp_path / 'run' / 'idx'), require=False)

    assert reason in caplog.text and 'Building it there now' in caplog.text
    assert 'adopted' not in json.loads((idx_dir / 'matam_index.json').read_text())


def test_a_failed_build_leaves_no_index_and_no_temporary_files(db, tmp_path):
    ref, idx_dir = db
    fake = tmp_path / 'sortmerna'
    fake.write_text('#!/bin/sh\nexit 3\n')
    fake.chmod(0o755)
    with pytest.raises(sortmerna_index.IndexBuildError, match='return code 3'):
        sortmerna_index.build(str(fake), str(ref), str(idx_dir), PARAMS)
    assert os.listdir(str(idx_dir.parent)) == []


def test_db_params_are_the_stored_ones_or_the_defaults(tmp_path):
    smr_dir = tmp_path / 'db.clustered.smr'
    assert sortmerna_index.db_params(str(smr_dir)) == LEGACY
    sortmerna_index.write_params(str(smr_dir), PARAMS)
    assert sortmerna_index.db_params(str(smr_dir)) == PARAMS


def test_ensure_builds_a_missing_index_in_a_writable_db_once(db, tmp_path, caplog):
    ref, idx_dir = db
    sortmerna_index.write_params(str(idx_dir.parent), PARAMS)
    local = tmp_path / 'run' / 'sortmerna_reads_mapping' / 'idx'
    caplog.set_level(logging.INFO, logger='sortmerna_index')

    assert sortmerna_index.ensure(SORTMERNA, str(ref), str(idx_dir), str(local), require=False) == str(idx_dir)
    assert 'no usable prebuilt index' in caplog.text and 'Building it there now' in caplog.text
    # Built with the database's parameters.
    assert json.loads((idx_dir / 'matam_index.json').read_text())['params'] == PARAMS
    caplog.clear()

    assert sortmerna_index.ensure(SORTMERNA, str(ref), str(idx_dir), str(local), require=False) == str(idx_dir)
    assert 'using the prebuilt index' in caplog.text and 'Building' not in caplog.text
    assert not local.exists()


def test_ensure_rebuilds_a_stale_index_in_place(db, tmp_path, caplog):
    ref, idx_dir = db
    built(ref, idx_dir, LEGACY)
    ref.write_text(''.join(ref.read_text().splitlines(True)[:8]))
    caplog.set_level(logging.INFO, logger='sortmerna_index')

    sortmerna_index.ensure(SORTMERNA, str(ref), str(idx_dir), str(tmp_path / 'local' / 'idx'), require=False)

    assert 'size differs' in caplog.text
    assert check(ref, idx_dir, LEGACY) == (True, '')


@pytest.mark.skipif(hasattr(os, 'geteuid') and os.geteuid() == 0, reason='root ignores permissions')
def test_ensure_builds_a_private_index_in_the_run_when_the_db_is_read_only(db, tmp_path, caplog):
    ref, idx_dir = db
    db_dir = ref.parent
    sortmerna_index.write_params(str(idx_dir.parent), PARAMS)
    before = sorted(str(p) for p in db_dir.rglob('*'))
    local = tmp_path / 'run' / 'sortmerna_reads_mapping' / 'idx'
    caplog.set_level(logging.INFO, logger='sortmerna_index')
    os.chmod(str(idx_dir.parent), stat.S_IRUSR | stat.S_IXUSR)
    os.chmod(str(db_dir), stat.S_IRUSR | stat.S_IXUSR)
    try:
        assert sortmerna_index.ensure(SORTMERNA, str(ref), str(idx_dir), str(local), require=False) == str(local)
        assert 'database directory is not writable' in caplog.text
        assert sorted(str(p) for p in db_dir.rglob('*')) == before
        assert json.loads((local / 'matam_index.json').read_text())['params'] == PARAMS
        caplog.clear()
        # A rerun of the same run reuses its private index.
        assert sortmerna_index.ensure(SORTMERNA, str(ref), str(idx_dir), str(local), require=False) == str(local)
        assert 'using the index built earlier for this run' in caplog.text
    finally:
        os.chmod(str(db_dir), 0o755)
        os.chmod(str(idx_dir.parent), 0o755)


CONCURRENT_ENSURE = textwrap.dedent('''\
    import sys
    sys.path.insert(0, sys.argv[1])
    import sortmerna_index
    print(sortmerna_index.ensure(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], require=False))
    ''')


def test_concurrent_runs_build_a_shared_index_once(db, tmp_path):
    ref, idx_dir = db
    # A SortMeRNA that is slow to index, and logs its index builds.
    calls = tmp_path / 'calls.txt'
    slow = tmp_path / 'sortmerna'
    slow.write_text('#!/bin/sh\ncase "$*" in *"--task 5"*) echo build >> {calls}; sleep 2;; esac\n'
                    'exec {real} "$@"\n'.format(calls=calls, real=SORTMERNA))
    slow.chmod(0o755)
    script = tmp_path / 'ensure.py'
    script.write_text(CONCURRENT_ENSURE)
    runs = [subprocess.Popen([sys.executable, str(script), str(SCRIPTS_DIR), str(slow), str(ref),
                              str(idx_dir), str(tmp_path / 'run{}'.format(n) / 'idx')],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            for n in range(3)]
    results = [run.communicate() + (run.returncode,) for run in runs]

    assert [code for _, _, code in results] == [0, 0, 0], results
    assert [out.strip() for out, _, _ in results] == [str(idx_dir)] * 3
    assert calls.read_text() == 'build\n'
    assert check(ref, idx_dir, LEGACY) == (True, '')


# require=True (the default; matam_assembly.py without --no_require_index):
# the prebuilt index or an error.

def test_by_default_the_prebuilt_index_is_used_without_taking_a_lock(db, tmp_path, monkeypatch):
    ref, idx_dir = db
    built(ref, idx_dir)
    sortmerna_index.write_params(str(idx_dir.parent), PARAMS)
    (idx_dir.parent / sortmerna_index.LOCK_NAME).unlink()

    def no_lock(directory):
        raise AssertionError('requiring the index must not lock ' + directory)
    monkeypatch.setattr(sortmerna_index, 'locked', no_lock)
    local = tmp_path / 'run' / 'idx'

    assert sortmerna_index.ensure(SORTMERNA, str(ref), str(idx_dir), str(local)) == str(idx_dir)
    assert not local.exists()
    assert not (idx_dir.parent / sortmerna_index.LOCK_NAME).exists()


@pytest.mark.parametrize('state, reason', [
    ('missing', 'no index'),
    ('stale', 'the fasta changed since the index was built (size differs)'),
    ('other params', 'the index was built with other parameters'),
])
def test_by_default_a_missing_or_unusable_index_is_never_built(db, tmp_path, state, reason):
    ref, idx_dir = db
    sortmerna_index.write_params(str(idx_dir.parent), PARAMS)
    if state != 'missing':
        built(ref, idx_dir, PARAMS if state == 'stale' else dict(PARAMS, m=3000))
    if state == 'stale':
        ref.write_text(''.join(ref.read_text().splitlines(True)[:4]))
    before = sorted(str(p) for p in ref.parent.rglob('*'))
    local = tmp_path / 'run' / 'idx'

    with pytest.raises(sortmerna_index.IndexMissingError) as error:
        sortmerna_index.ensure(SORTMERNA, str(ref), str(idx_dir), str(local),
                               label='clustered ref db')

    assert reason in str(error.value) and 'clustered ref db' in str(error.value)
    assert sorted(str(p) for p in ref.parent.rglob('*')) == before   # nothing built or changed
    assert not local.exists()


def test_require_adopts_an_index_without_a_record(db):
    ref, idx_dir = db
    legacy_index(ref, idx_dir)
    assert sortmerna_index.prebuilt(SORTMERNA, str(ref), str(idx_dir)) == str(idx_dir)
    assert 'adopted' in json.loads((idx_dir / 'matam_index.json').read_text())


def test_require_reports_a_missing_fasta(tmp_path):
    with pytest.raises(sortmerna_index.IndexMissingError, match='does not exist'):
        sortmerna_index.prebuilt(SORTMERNA, str(tmp_path / 'db.complete.fasta'),
                                 str(tmp_path / 'db.complete.smr' / 'idx'), 'complete ref db')
