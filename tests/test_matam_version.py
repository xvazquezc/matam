"""MATAM --version: the VERSION number, plus the short commit hash off the tagged commit."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
from pe_fixtures import SCRIPTS_DIR  # noqa: E402

pytestmark = pytest.mark.skipif(shutil.which('git') is None, reason='missing: git')
GIT_ENV = dict(os.environ, GIT_AUTHOR_NAME='t', GIT_AUTHOR_EMAIL='t@t', GIT_COMMITTER_NAME='t',
               GIT_COMMITTER_EMAIL='t@t', GIT_CONFIG_GLOBAL=os.devnull, GIT_CONFIG_NOSYSTEM='1')


def git(repo, *args):
    return subprocess.run(['git', '-C', str(repo), *args], check=True, env=GIT_ENV,
                          stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()


def matam_copy(root, number, init_git=True):
    """A MATAM-like tree: VERSION and scripts/matam_version.py, optionally a git repo."""
    (root / 'scripts').mkdir(parents=True)
    shutil.copy(str(SCRIPTS_DIR / 'matam_version.py'), str(root / 'scripts'))
    (root / 'VERSION').write_text(number)
    if init_git:
        git(root, 'init', '-q')
        commit(root)
    return root


def commit(repo, message='c'):
    git(repo, 'add', '-A')
    git(repo, 'commit', '-q', '--allow-empty', '-m', message)
    return git(repo, 'rev-parse', '--short', 'HEAD')


def version(root, cwd=None):
    return subprocess.run([sys.executable, str(root / 'scripts' / 'matam_version.py')],
                          check=True, stdout=subprocess.PIPE, universal_newlines=True,
                          cwd=str(cwd or root), env=GIT_ENV).stdout.strip()


def test_empty_version_file_gives_the_short_commit_hash(tmp_path):
    root = matam_copy(tmp_path / 'matam', '')
    assert version(root) == git(root, 'rev-parse', '--short', 'HEAD')


def test_the_tagged_commit_gives_the_version_number(tmp_path):
    root = matam_copy(tmp_path / 'matam', '2.0.0\n')
    git(root, 'tag', 'v2.0.0')
    assert version(root) == '2.0.0'
    git(root, 'tag', '-d', 'v2.0.0')
    git(root, 'tag', '-a', 'v2.0.0', '-m', 'annotated')      # annotated tags too
    assert version(root) == '2.0.0'


def test_other_commits_add_their_short_hash(tmp_path):
    root = matam_copy(tmp_path / 'matam', '2.0.0')
    short = git(root, 'rev-parse', '--short', 'HEAD')
    assert version(root) == '2.0.0+' + short              # no tag yet
    git(root, 'tag', 'v2.0.0')
    short = commit(root, 'after the release')
    assert version(root) == '2.0.0+' + short              # a commit after the tag
    git(root, 'tag', 'v1.9.0')
    assert version(root) == '2.0.0+' + short              # tagged, but not as VERSION


def test_runs_from_any_directory_and_through_a_symlink(tmp_path):
    root = matam_copy(tmp_path / 'matam', '')
    (tmp_path / 'bin').mkdir()
    link = tmp_path / 'bin' / 'matam_version.py'
    link.symlink_to(root / 'scripts' / 'matam_version.py')
    out = subprocess.run([sys.executable, str(link)], check=True, stdout=subprocess.PIPE,
                         universal_newlines=True, cwd=str(tmp_path), env=GIT_ENV).stdout.strip()
    assert out == git(root, 'rev-parse', '--short', 'HEAD')


def test_outside_a_git_checkout_only_the_version_file_counts(tmp_path):
    assert version(matam_copy(tmp_path / 'a', '2.0.0', init_git=False)) == '2.0.0'
    assert version(matam_copy(tmp_path / 'b', '', init_git=False)) == 'unknown'


def test_a_copy_inside_another_repository_does_not_report_its_commit(tmp_path):
    outer = tmp_path / 'project'
    outer.mkdir()
    git(outer, 'init', '-q')
    root = matam_copy(outer / 'vendor' / 'matam', '2.0.0', init_git=False)
    commit(outer)
    assert version(root) == '2.0.0'


@pytest.mark.parametrize('script', ['matam_assembly.py', 'matam_db_preprocessing.py',
                                    'index_ref_db.py'])
def test_entry_points_print_the_version(script):
    out = subprocess.run([sys.executable, str(SCRIPTS_DIR / script), '--version'], check=True,
                         stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()
    expected = subprocess.run([sys.executable, str(SCRIPTS_DIR / 'matam_version.py')], check=True,
                              stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()
    assert out == '{} {}'.format(script, expected) and expected != 'unknown'
