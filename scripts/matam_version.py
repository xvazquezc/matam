#!/usr/bin/env python3
"""MATAM version, for --version and the logs.

The release number is in the VERSION file at the root of the repository
(empty until the first release). On the commit tagged ``v<VERSION>`` the
version is that number; on any other commit the short commit hash is added
(``<VERSION>+<hash>``), or is the whole version while VERSION is empty. Outside
a git checkout of MATAM (no git, or no ``.git`` of its own), it is VERSION
alone, or 'unknown'.
"""

import os
import subprocess

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
VERSION_FILE = os.path.join(ROOT_DIR, 'VERSION')


def _git(*args):
    """Output of a git command run in the MATAM repository, or None if it fails."""
    try:
        completed = subprocess.run(['git', '-C', ROOT_DIR, *args], stdout=subprocess.PIPE,
                                   stderr=subprocess.DEVNULL, universal_newlines=True,
                                   timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    output = completed.stdout.strip()
    return output if completed.returncode == 0 and output else None


def release():
    """The release number in VERSION ('' if empty or missing)."""
    try:
        with open(VERSION_FILE) as handle:
            return handle.read().strip()
    except OSError:
        return ''


def commit():
    """Short hash of the checked-out commit, or None outside a git checkout of MATAM."""
    toplevel = _git('rev-parse', '--show-toplevel')
    # A MATAM copy inside another repository must not report that one's commit.
    if toplevel is None or os.path.realpath(toplevel) != ROOT_DIR:
        return None
    return _git('rev-parse', '--short', 'HEAD')


def version():
    number = release()
    short_hash = commit()
    if short_hash is None:
        return number or 'unknown'
    if number:
        tagged = _git('rev-parse', '--verify', '--quiet', 'refs/tags/v{}^{{commit}}'.format(number))
        if tagged is not None and tagged == _git('rev-parse', 'HEAD'):
            return number
        return '{}+{}'.format(number, short_hash)
    return short_hash


if __name__ == '__main__':
    print(version())
