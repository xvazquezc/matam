#!/usr/bin/env python3
"""Run SortMeRNA steps so that reruns and ``--resume_from`` are safe.

SortMeRNA 7 keeps run state in ``<workdir>/kvdb`` (and ``<workdir>/readb``).
When it finds state from a run with different inputs it aborts
(``kvdb resume aborted: inputs mismatch``), so a rerun or a resumed MATAM run
whose inputs changed (e.g. rebuilt contigs) failed. SortMeRNA has no option to
force an overwrite.

A step is therefore treated as follows:

* after SortMeRNA succeeds, a completion record (``<aligned>.matam_complete.json``)
  stores the command, the input files (size and modification time) and the
  required outputs (size); the now useless run state is removed;
* when the step comes up again and the record matches the current command,
  inputs and outputs, the output is complete and is reused;
* otherwise (no record, missing or truncated output, changed inputs or
  command) the step's SortMeRNA outputs and run state are removed and
  SortMeRNA runs from scratch.

Options that don't change the result (``--threads``, ``-v``, the location of
the reference index) are passed as
``run_options`` and are not part of the record. Paths inside the run (the
SortMeRNA work directory and the directory of the outputs) are stored relative
to it, so a run directory that was moved or copied can still be resumed.
"""

import json
import logging
import os
import shutil

import runner

logger = logging.getLogger(__name__)

RECORD_SUFFIX = '.matam_complete.json'
RUN_STATE_DIRS = ('kvdb', 'readb')
# What SortMeRNA 7 can write for an --aligned prefix (``--zip-out`` adds .gz).
_OUTPUT_STEMS = ('', '_fwd', '_rev')
_OUTPUT_EXTENSIONS = ('.sam', '.blast', '.fq', '.fastq', '.fa', '.fasta', '.log')


class IncompleteStepError(Exception):
    """A SortMeRNA step that is skipped has no complete output to continue from."""


def find_sortmerna_output(basepath, extensions):
    """Return the one SortMeRNA output matching a prefix and file type.

    SortMeRNA 7 writes ``.fq`` for FASTQ output regardless of whether the
    input was named ``.fq`` or ``.fastq``.  ``--zip-out`` additionally appends
    ``.gz`` to all report files.  Do not infer the output name from the input
    name: discover the generated file instead.
    """
    candidates = [
        basepath + extension + compression
        for extension in extensions
        for compression in ('', '.gz')
    ]
    found = [candidate for candidate in candidates if os.path.isfile(candidate)]
    if len(found) == 1:
        return found[0]
    if not found:
        raise FileNotFoundError(
            'SortMeRNA did not produce an expected output: {}'.format(
                ', '.join(candidates)))
    raise RuntimeError(
        'SortMeRNA produced ambiguous outputs for {}: {}'.format(
            basepath, ', '.join(found)))


def record_path(aligned_base):
    return aligned_base + RECORD_SUFFIX


def _relative(path, base_dir):
    """Path relative to base_dir when inside it, absolute otherwise."""
    path = os.path.abspath(path)
    if os.path.commonpath([path, base_dir]) == base_dir:
        return os.path.relpath(path, base_dir)
    return path


def _base_dir(aligned_base):
    return os.path.dirname(os.path.abspath(aligned_base))


def _canonical_command(command, workdir, aligned_base):
    """The command with the run's own directories replaced by placeholders."""
    replacements = {}
    for path, placeholder in ((workdir, '<sortmerna_workdir>'), (_base_dir(aligned_base), '<output_dir>')):
        replacements[os.path.abspath(path)] = placeholder
        replacements[path] = placeholder
    for path in sorted(replacements, key=len, reverse=True):
        command = command.replace(path, replacements[path])
    return command


def _signature(path, base_dir):
    stat = os.stat(path)
    return {'path': _relative(path, base_dir), 'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns}


def _output_candidates(aligned_base):
    return [aligned_base + stem + extension + compression
            for stem in _OUTPUT_STEMS
            for extension in _OUTPUT_EXTENSIONS
            for compression in ('', '.gz')]


def _discover(aligned_base, required_outputs):
    """Paths of the required outputs, given as (prefix suffix, extensions) pairs."""
    return [find_sortmerna_output(aligned_base + suffix, extensions)
            for suffix, extensions in required_outputs]


def _read_record(aligned_base):
    try:
        with open(record_path(aligned_base)) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def is_complete(aligned_base, command, inputs, workdir):
    """True if the recorded run used this command and these inputs, and its outputs are intact."""
    record = _read_record(aligned_base)
    base_dir = _base_dir(aligned_base)
    if (not isinstance(record, dict)
            or record.get('command') != _canonical_command(command, workdir, aligned_base)):
        return False
    try:
        if record.get('inputs') != [_signature(path, base_dir) for path in inputs]:
            return False
    except OSError:  # an input is missing
        return False
    outputs = record.get('outputs')
    if not outputs:
        return False
    for output in outputs:
        path = os.path.join(base_dir, output['path'])
        if not os.path.isfile(path) or os.path.getsize(path) != output['size']:
            return False
    return True


def clear(workdir, aligned_base):
    """Remove the step's completion record, SortMeRNA outputs and run state."""
    removed = []
    for path in [record_path(aligned_base)] + _output_candidates(aligned_base):
        if os.path.lexists(path):
            os.remove(path)
            removed.append(path)
    for name in RUN_STATE_DIRS:
        path = os.path.join(workdir, name)
        if os.path.isdir(path):
            shutil.rmtree(path)
            removed.append(path)
    return removed


def run(command, workdir, aligned_base, inputs, required_outputs, run_options='', verbose=False,
        prepare=None):
    """Run a SortMeRNA step unless a complete output for these inputs exists.

    ``command`` must contain everything that determines the output, and
    ``run_options`` the rest (threads, verbosity). Returns the paths of the
    required outputs, in the order of ``required_outputs``.

    ``prepare``, if given, is called only when SortMeRNA has to run, and
    returns more run options (e.g. finds or builds the reference index and
    returns ``--idx-dir``); a reused output doesn't need it.
    """
    if is_complete(aligned_base, command, inputs, workdir):
        logger.info('Reusing complete SortMeRNA output (same command and inputs): %s', aligned_base)
        return _discover(aligned_base, required_outputs)

    removed = clear(workdir, aligned_base)
    if removed:
        logger.info('Removing incomplete or outdated SortMeRNA output before running: %s',
                    ', '.join(removed))
    if prepare is not None:
        run_options = prepare() + run_options
    runner.logged_check_call(command + run_options, verbose=verbose)

    outputs = _discover(aligned_base, required_outputs)
    base_dir = _base_dir(aligned_base)
    record = {
        'command': _canonical_command(command, workdir, aligned_base),
        'inputs': [_signature(path, base_dir) for path in inputs],
        'outputs': [{'path': _relative(path, base_dir), 'size': os.path.getsize(path)}
                    for path in outputs],
    }
    with open(record_path(aligned_base), 'w') as handle:
        json.dump(record, handle, indent=1)
    # The run state only serves to resume an interrupted SortMeRNA run; the
    # record now covers this step, and stale state would block the next run.
    for name in RUN_STATE_DIRS:
        shutil.rmtree(os.path.join(workdir, name), ignore_errors=True)
    return outputs


def completed_outputs(command, workdir, aligned_base, inputs, required_outputs):
    """Outputs of a step that is skipped (e.g. by ``--resume_from``).

    Raises IncompleteStepError when the step's output is missing, incomplete
    or was made from other inputs. Output from a MATAM version that wrote no
    completion record is accepted, with a warning, if the required files exist.
    """
    if is_complete(aligned_base, command, inputs, workdir):
        return _discover(aligned_base, required_outputs)
    if _read_record(aligned_base) is None:
        try:
            outputs = _discover(aligned_base, required_outputs)
        except (FileNotFoundError, RuntimeError) as error:
            raise IncompleteStepError(str(error))
        logger.warning('No completion record for the SortMeRNA output %s; '
                       'its completeness cannot be checked, using it as is', aligned_base)
        return outputs
    raise IncompleteStepError(
        'The SortMeRNA output {} is incomplete or was made with other inputs or '
        'options'.format(aligned_base))
