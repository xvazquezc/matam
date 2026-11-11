#!/usr/bin/env python3
"""SortMeRNA reference indexes that are known to match their fasta.

SortMeRNA 7 names an index after a hash of the reference file *name* only and,
when it finds non-empty index files of that name, uses them without checking
them against the fasta: a stale index (fasta rewritten under the same name) is
used silently, and can give wrong alignments or crash SortMeRNA. The indexing
options are not checked either, although they change the alignments: ``-L``
and ``--max_pos`` obviously, but also ``-m``, which sets how many parts the
index is split into (alignments are selected per part).

Indexing parameters
    :func:`index_params` is the one definition of the indexing parameters.
    matam_db_preprocessing.py stores the parameters of a database in
    ``<db>.<kind>.smr/matam_index_params.json`` (also with ``--skip_index``),
    and every index of that database, prebuilt or built at mapping time, is
    built with them (:func:`db_params`). A database without stored parameters
    (made by an older MATAM) keeps the SortMeRNA defaults it was indexed with.

Index record
    :func:`build` writes a record (``matam_index.json``) next to the index
    files: the fasta name, size, modification time and SHA-256, the indexing
    parameters, the SortMeRNA version and the size of each index file.
    :func:`check` accepts an index only when the record matches the fasta (the
    checksum is recomputed only when the size matches but the modification
    time does not, e.g. after a copy), the parameters, the SortMeRNA version,
    and the index files are intact. An index without a record (built by an
    older MATAM) is adopted only if SortMeRNA's own index statistics match the
    fasta and the parameters (:func:`adopt`); the record is then added.

Safe builds
    An index is built in a temporary directory beside its final location and
    moved into place, with its record, in one rename, so an index that has a
    record is complete. Builds are serialised with a lock file, so concurrent
    MATAM runs build a shared index once and never read a half-written one.

:func:`ensure` finds the index used by a mapping step. By default
(``require=True``) it must be the prebuilt index in the database directory
(``<db>.<kind>.smr/idx``): :func:`prebuilt` accepts it when it matches the
fasta and the parameters (or adopts one without a record), and otherwise raises
:class:`IndexMissingError`, without building anything or taking a lock. An
unindexed database thus stops a run at once, instead of every run (e.g. every
job of a pipeline) spending hours indexing it.

With ``require=False`` (matam_assembly.py ``--no_require_index``) a missing or
unusable index is built instead:

1. the prebuilt index, when usable, as above;
2. otherwise, when the database directory is writable, it builds the index
   there once, for later runs too;
3. otherwise (read-only database), it builds a private index inside the run
   (the step's SortMeRNA work directory), reused by reruns of the same run.
"""

import contextlib
import fcntl
import hashlib
import json
import logging
import os
import re
import shlex
import shutil
import subprocess
import tempfile

import runner

logger = logging.getLogger(__name__)

# Indexing parameters (SortMeRNA -L, --max_pos, -m), passed explicitly so a
# change of SortMeRNA defaults cannot change the index silently.
# L=18 and max_pos=10000 are the SortMeRNA 2.1b defaults used by upstream
# MATAM; they are also SortMeRNA 7.1's effective defaults (its help says
# max_pos 1000, but an index built without --max_pos is byte-identical to one
# built with --max_pos 10000, and differs from --max_pos 1000).
# -m is the --max_memory of matam_db_preprocessing.py (default DEFAULT_MEMORY_MB).
SEED_LENGTH = 18
MAX_POS = 10000
DEFAULT_MEMORY_MB = 10000

RECORD_NAME = 'matam_index.json'
PARAMS_NAME = 'matam_index_params.json'
LOCK_NAME = '.matam_index.lock'
RECORD_FORMAT = 1
NO_RECORD = 'the index has no MATAM index record ({})'.format(RECORD_NAME)
# -m of the indexes built by MATAM before index records (SortMeRNA's default).
LEGACY_MEMORY_MB = 3072
_version_cache = {}


class IndexBuildError(Exception):
    """SortMeRNA failed to build an index."""


class IndexMissingError(Exception):
    """The database has no usable prebuilt index, and building one is not allowed."""


def index_params(memory_mb=DEFAULT_MEMORY_MB):
    """The indexing parameters, as SortMeRNA option name -> value."""
    return {'L': SEED_LENGTH, 'max_pos': MAX_POS, 'm': int(memory_mb)}


def _params_options(params):
    return ' -L {L} --max_pos {max_pos} -m {m}'.format(**params)


def params_path(smr_dir):
    return os.path.join(smr_dir, PARAMS_NAME)


def write_params(smr_dir, params):
    """Store the indexing parameters of a database (matam_db_preprocessing.py)."""
    os.makedirs(smr_dir, exist_ok=True)
    tmp = params_path(smr_dir) + '.tmp.{}'.format(os.getpid())
    with open(tmp, 'w') as handle:
        json.dump(params, handle, indent=1)
    os.replace(tmp, params_path(smr_dir))


def db_params(smr_dir):
    """The indexing parameters stored for a database.

    A database without stored parameters was indexed by an older MATAM, with
    SortMeRNA's defaults (-m 3072); its indexes keep these parameters, so an
    existing index stays usable and a rebuilt one gives the same alignments.
    """
    try:
        with open(params_path(smr_dir)) as handle:
            params = json.load(handle)
    except (OSError, ValueError):
        return index_params(LEGACY_MEMORY_MB)
    if not isinstance(params, dict) or set(params) != set(index_params()):
        return index_params(LEGACY_MEMORY_MB)
    return params


def sortmerna_version(sortmerna_bin):
    """The SortMeRNA version string (e.g. '7.1.0'), or None if unknown."""
    if sortmerna_bin not in _version_cache:
        try:
            output = subprocess.run([sortmerna_bin, '--version'], stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, universal_newlines=True).stdout
        except OSError:
            output = ''
        match = re.search(r'SortMeRNA version (\S+)', output)
        _version_cache[sortmerna_bin] = match.group(1) if match else None
    return _version_cache[sortmerna_bin]


def sha256sum(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def record_path(idx_dir):
    return os.path.join(idx_dir, RECORD_NAME)


def read_record(idx_dir):
    try:
        with open(record_path(idx_dir)) as handle:
            record = json.load(handle)
    except (OSError, ValueError):
        return None
    return record if isinstance(record, dict) else None


def _write_record(idx_dir, record):
    tmp = record_path(idx_dir) + '.tmp.{}'.format(os.getpid())
    with open(tmp, 'w') as handle:
        json.dump(record, handle, indent=1)
    os.replace(tmp, record_path(idx_dir))


def _index_files(idx_dir):
    return sorted(name for name in os.listdir(idx_dir)
                  if name != RECORD_NAME and os.path.isfile(os.path.join(idx_dir, name)))


def check(idx_dir, ref_fasta, sortmerna_bin, params):
    """Return (True, '') if idx_dir holds a complete index of ref_fasta built
    with params, else (False, reason)."""
    if not os.path.isdir(idx_dir) or not os.listdir(idx_dir):
        return False, 'no index'
    record = read_record(idx_dir)
    if record is None or record.get('format') != RECORD_FORMAT:
        return False, NO_RECORD
    reference = record.get('reference', {})
    if reference.get('name') != os.path.basename(ref_fasta):
        return False, 'the index was built for {}'.format(reference.get('name'))
    stat = os.stat(ref_fasta)
    if stat.st_size != reference.get('size'):
        return False, 'the fasta changed since the index was built (size differs)'
    if record.get('params') != params:
        return False, 'the index was built with other parameters ({} instead of {})'.format(
            record.get('params'), params)
    if not record.get('files'):
        return False, 'the index record lists no index file'
    for entry in record['files']:
        path = os.path.join(idx_dir, entry['name'])
        if not os.path.isfile(path) or os.path.getsize(path) != entry['size']:
            return False, 'the index is incomplete ({} missing or truncated)'.format(entry['name'])
    version = sortmerna_version(sortmerna_bin)
    if record.get('sortmerna_version') != version:
        return False, 'the index was built with SortMeRNA {}, this is SortMeRNA {}'.format(
            record.get('sortmerna_version'), version)
    if stat.st_mtime_ns != reference.get('mtime_ns'):
        # Same size, other timestamp: a copy, or a rewrite. Compare the content.
        if sha256sum(ref_fasta) != reference.get('sha256'):
            return False, 'the fasta changed since the index was built (checksum differs)'
        reference['mtime_ns'] = stat.st_mtime_ns
        try:  # skip the checksum next time; not needed for correctness
            _write_record(idx_dir, record)
        except OSError:
            pass
    return True, ''


def _fasta_counts(path):
    """(number of sequences, total sequence length) of a fasta file."""
    numseq = full_len = 0
    with open(path, 'rb') as handle:
        for line in handle:
            if line.startswith(b'>'):
                numseq += 1
            else:
                full_len += len(line.rstrip(b'\r\n'))
    return numseq, full_len


def _read_idx_stats(path):
    """The fields of a SortMeRNA ``<hash>.idx_stats.yaml`` needed to check the index."""
    text = open(path).read()

    def ints(key):
        return [int(value) for value in re.findall(r'^\s*{}:\s*(\d+)\s*$'.format(key), text, re.M)]

    reference = re.search(r'^reference:\s*(.+?)\s*$', text, re.M)
    prefix = os.path.basename(path)[:-len('.idx_stats.yaml')]
    files = {prefix + '.' + name: int(size) for name, size in re.findall(
        r'^\s*((?:bursttrie|pos|kmer)_\d+\.dat):\s*(\d+)\s*$', text, re.M)}
    files.update({prefix + '.stats': size for size in ints('stats')})
    return {
        'reference': reference.group(1).strip('\'"') if reference else None,
        'numseq': ints('numseq'), 'full_len': ints('full_len'), 'part_num': ints('part_num'),
        'seed_win_len': ints('seed_win_len'), 'seq_part_size': ints('seq_part_size'),
        'max_positions': ints('max_positions'), 'files': files,
    }


def adopt(idx_dir, ref_fasta, sortmerna_bin, params):
    """Accept an index without a MATAM record if SortMeRNA's own statistics match.

    Indexes built before MATAM wrote records (``sortmerna --task 5`` without
    -L/--max_pos/-m, i.e. L=18, max_pos=10000, -m 3072) are checked against
    ``<hash>.idx_stats.yaml``: same fasta name, size, number of sequences and
    total length, seed length, positions within max_pos, a partitioning that
    -m gives (any with -m 3072, the parameters of a database without stored
    parameters; otherwise a single part), and intact index files. The record is then written
    (if the directory is writable), so this is done once.
    Returns (True, '') or (False, reason).
    """
    name = os.path.basename(ref_fasta)
    stats = None
    for entry in sorted(os.listdir(idx_dir)):
        if entry.endswith('.idx_stats.yaml'):
            candidate = _read_idx_stats(os.path.join(idx_dir, entry))
            if candidate['reference'] and os.path.basename(candidate['reference']) == name:
                stats, stats_name = candidate, entry
                break
    if stats is None:
        return False, NO_RECORD + ' and no SortMeRNA index statistics for {}'.format(name)
    size = os.path.getsize(ref_fasta)
    if sum(stats['seq_part_size']) != size:
        return False, NO_RECORD + ' and the fasta changed since it was built (size differs)'
    if stats['seed_win_len'] != [params['L']]:
        return False, NO_RECORD + ' and it was built with another seed length ({})'.format(
            stats['seed_win_len'])
    # Below the cap it was built with, max_positions is the true maximum and the
    # index is the same as with any higher cap. 1000 is SortMeRNA's documented
    # (not effective) default, so an index stopping at 1000 was likely capped
    # there; other caps below max_pos can't be told from a true maximum.
    max_positions = max(stats['max_positions'] or [0])
    if (not stats['max_positions'] or max_positions > params['max_pos']
            or (max_positions == 1000 and params['max_pos'] != 1000)):
        return False, NO_RECORD + ' and it was built with another --max_pos ({})'.format(
            max_positions)
    # The number of parts (set by -m) changes the alignments. An older MATAM
    # used -m 3072; with another -m, only a single-part index (which a larger
    # -m would also give) is the same.
    parts = len(stats['seq_part_size'])
    if stats['part_num'] != [parts] or (params['m'] != LEGACY_MEMORY_MB and (
            parts != 1 or params['m'] < LEGACY_MEMORY_MB)):
        return False, NO_RECORD + ' and its partitioning may differ from -m {}'.format(params['m'])
    if not stats['files']:
        return False, NO_RECORD + ' and its statistics list no index file'
    for file_name, file_size in stats['files'].items():
        path = os.path.join(idx_dir, file_name)
        if not os.path.isfile(path) or os.path.getsize(path) != file_size:
            return False, NO_RECORD + ' and it is incomplete ({} missing or truncated)'.format(
                file_name)
    if [list(_fasta_counts(ref_fasta))] != [stats['numseq'] + stats['full_len']]:
        return False, NO_RECORD + ' and the fasta changed since it was built (sequences differ)'

    stat = os.stat(ref_fasta)
    record = {
        'format': RECORD_FORMAT,
        'reference': {'name': name, 'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns,
                      'sha256': sha256sum(ref_fasta)},
        'params': params,
        'sortmerna_version': sortmerna_version(sortmerna_bin),
        'files': [{'name': file_name, 'size': os.path.getsize(os.path.join(idx_dir, file_name))}
                  for file_name in sorted(list(stats['files']) + [stats_name])],
        'adopted': 'built without a MATAM record; checked against SortMeRNA index statistics',
    }
    try:
        _write_record(idx_dir, record)
    except OSError as error:
        logger.info('Cannot write the MATAM record of the index %s (%s); it will be checked '
                    'again next time', idx_dir, error)
    return True, ''


def build(sortmerna_bin, ref_fasta, idx_dir, params, verbose=False, call=None):
    """Build the index of ref_fasta into idx_dir (replacing any index there) and record it.

    The caller must hold the lock (:func:`locked`) of the directory containing
    idx_dir when other processes may use it. ``call(command, verbose)`` runs
    the SortMeRNA command and returns its exit code (default:
    runner.logged_call).
    """
    parent = os.path.dirname(os.path.abspath(idx_dir))
    os.makedirs(parent, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix='.matam_index_tmp.', dir=parent)
    try:
        new_idx = os.path.join(tmp, 'idx')
        stat = os.stat(ref_fasta)
        reference = {'name': os.path.basename(ref_fasta), 'size': stat.st_size,
                     'mtime_ns': stat.st_mtime_ns, 'sha256': sha256sum(ref_fasta)}
        # --workdir in tmp too: '--task 5' also creates kvdb/ and out/ there.
        cmd_line = '{} --ref {} --workdir {} --idx-dir {} --task 5'.format(
            shlex.quote(sortmerna_bin), shlex.quote(ref_fasta),
            shlex.quote(os.path.join(tmp, 'work')), shlex.quote(new_idx))
        cmd_line += _params_options(params)
        if verbose:
            cmd_line += ' -v'
        returncode = (call or runner.logged_call)(cmd_line, verbose)
        if returncode != 0 or not os.path.isdir(new_idx) or not _index_files(new_idx):
            raise IndexBuildError('SortMeRNA could not index {} (return code {}): {}'.format(
                ref_fasta, returncode, cmd_line))
        if os.stat(ref_fasta).st_size != reference['size'] or sha256sum(ref_fasta) != reference['sha256']:
            raise IndexBuildError('{} changed while it was being indexed'.format(ref_fasta))
        _write_record(new_idx, {
            'format': RECORD_FORMAT,
            'reference': reference,
            'params': params,
            'sortmerna_version': sortmerna_version(sortmerna_bin),
            'files': [{'name': name, 'size': os.path.getsize(os.path.join(new_idx, name))}
                      for name in _index_files(new_idx)],
        })
        if os.path.lexists(idx_dir):
            os.rename(idx_dir, os.path.join(tmp, 'old'))
        os.rename(new_idx, idx_dir)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


@contextlib.contextmanager
def locked(directory):
    """Hold an exclusive lock on directory (created if needed) for index builds."""
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, LOCK_NAME), 'a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def _writable(directory):
    """True if an index can be built in directory (created if needed)."""
    try:
        os.makedirs(directory, exist_ok=True)
    except OSError:
        return False
    return os.access(directory, os.W_OK | os.X_OK)


def prebuilt(sortmerna_bin, ref_fasta, shared_idx_dir, label='reference'):
    """Return shared_idx_dir if it holds a usable index of ref_fasta; never builds one.

    The index must match the fasta and the database's parameters
    (:func:`check`), or be an index without a MATAM record that
    :func:`adopt` accepts. Takes no lock, so it also works where file locks
    are not supported. Raises IndexMissingError (with the reason) otherwise.
    """
    if not os.path.isfile(ref_fasta):
        raise IndexMissingError('the {} fasta {} does not exist'.format(label, ref_fasta))
    params = db_params(os.path.dirname(os.path.abspath(shared_idx_dir)))
    usable, reason = check(shared_idx_dir, ref_fasta, sortmerna_bin, params)
    if usable:
        logger.info('SortMeRNA index of the %s: using the prebuilt index %s', label, shared_idx_dir)
        return shared_idx_dir
    if reason == NO_RECORD:
        usable, reason = adopt(shared_idx_dir, ref_fasta, sortmerna_bin, params)
        if usable:
            logger.info('SortMeRNA index of the %s: using the prebuilt index %s, which has no '
                        'MATAM record but matches the fasta and parameters', label, shared_idx_dir)
            return shared_idx_dir
    raise IndexMissingError('no usable prebuilt SortMeRNA index of the {} {} in {}: {}'.format(
        label, ref_fasta, shared_idx_dir, reason))


def ensure(sortmerna_bin, ref_fasta, shared_idx_dir, local_idx_dir, verbose=False,
           label='reference', require=True):
    """Return the directory of an index of ref_fasta.

    shared_idx_dir is the database's index (``<db>.<kind>.smr/idx``). With
    require=True (the default), only it is used (:func:`prebuilt`): nothing is
    built, and IndexMissingError is raised when it is missing or unusable.
    With require=False, a missing or unusable index is built: in
    shared_idx_dir when possible, else in local_idx_dir (the run's private
    fallback), with the database's parameters (:func:`db_params`); raises
    IndexBuildError if SortMeRNA fails.
    """
    if require:
        return prebuilt(sortmerna_bin, ref_fasta, shared_idx_dir, label)
    shared_dir = os.path.dirname(os.path.abspath(shared_idx_dir))
    params = db_params(shared_dir)
    usable, reason = check(shared_idx_dir, ref_fasta, sortmerna_bin, params)
    if usable:
        logger.info('SortMeRNA index of the %s: using the prebuilt index %s', label, shared_idx_dir)
        return shared_idx_dir

    def adopted():
        """Adopt an index from before MATAM records; updates reason otherwise."""
        nonlocal reason
        if reason != NO_RECORD:
            return False
        usable, reason = adopt(shared_idx_dir, ref_fasta, sortmerna_bin, params)
        if usable:
            logger.info('SortMeRNA index of the %s: using the prebuilt index %s, which has no '
                        'MATAM record but matches the fasta and parameters (record added)',
                        label, shared_idx_dir)
        return usable

    if _writable(shared_dir):
        try:
            with locked(shared_dir):
                # Another run may have built it while this one waited for the lock.
                usable, reason = check(shared_idx_dir, ref_fasta, sortmerna_bin, params)
                if usable:
                    logger.info('SortMeRNA index of the %s: using the index %s, built by a '
                                'concurrent run', label, shared_idx_dir)
                    return shared_idx_dir
                if adopted():
                    return shared_idx_dir
                logger.warning('SortMeRNA index of the %s: no usable prebuilt index in %s (%s). '
                               'Building it there now; later runs will reuse it.',
                               label, shared_idx_dir, reason)
                build(sortmerna_bin, ref_fasta, shared_idx_dir, params, verbose)
                logger.info('SortMeRNA index of the %s: built %s', label, shared_idx_dir)
                return shared_idx_dir
        except OSError as error:
            logger.warning('SortMeRNA index of the %s: cannot build the index in the database '
                           'directory (%s)', label, error)
    elif adopted():
        return shared_idx_dir

    with locked(os.path.dirname(os.path.abspath(local_idx_dir))):
        usable, _ = check(local_idx_dir, ref_fasta, sortmerna_bin, params)
        if usable:
            logger.info('SortMeRNA index of the %s: using the index built earlier for this run, '
                        '%s', label, local_idx_dir)
            return local_idx_dir
        logger.warning('SortMeRNA index of the %s: no usable prebuilt index in %s (%s) and the '
                       'database directory is not writable. Building a private index for this '
                       'run in %s; every run on this database will rebuild it. Index the '
                       'database with matam_db_preprocessing.py (without --skip_index) in a '
                       'writable directory to avoid this.',
                       label, shared_idx_dir, reason, local_idx_dir)
        build(sortmerna_bin, ref_fasta, local_idx_dir, params, verbose)
    return local_idx_dir
