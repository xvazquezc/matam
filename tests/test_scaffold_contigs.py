"""Characterization of scaffold_contigs.find_called_base (scaffold consensus).

MATAM runs ``samtools mpileup`` without ``-f``, so read bases are letters
(never ``.``/``,``) and the reference column is N.
"""
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).parent.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS_DIR))
from scaffold_contigs import find_called_base, iter_read_bases  # noqa: E402
sys.path.pop(0)


def test_read_bases_are_tokenised_with_markers_and_indels():
    assert list(iter_read_bases('^]A+2CGa-1T*$')) == ['^]', 'A', '+CG', 'a', '-T', '*', '$']


def test_unknown_read_base_is_an_error():
    with pytest.raises(ValueError, match='Unrecognised mpileup read base'):
        list(iter_read_bases('A#'))


@pytest.mark.parametrize('read_bases, expected', [
    ('AAC', 'A'),
    ('aaC', 'A'),          # case is ignored
    ('^]AA$C', 'A'),       # read start/end markers are not bases
    ('NNA', 'N'),          # N is a base like any other
    ('A-2CGA-2CGC', 'A'),  # deletions after this position are not counted here
])
def test_majority_base_is_called(read_bases, expected):
    assert find_called_base('N', read_bases, 3) == expected


def test_majority_deletion_calls_no_base():
    assert find_called_base('N', '**A', 3) == ''


@pytest.mark.parametrize('read_bases, expected', [
    ('AT', 'T'),
    ('CG', 'G'),
    ('A*', 'A'),   # '*' sorts before letters, so a base beats a tied deletion
])
def test_ties_are_settled_alphabetically_last_wins(read_bases, expected):
    assert find_called_base('N', read_bases, 2) == expected


def test_insertion_carried_by_half_the_coverage_is_added():
    assert find_called_base('N', 'A+2CGA+2CGAA', 4) == 'ACG'


def test_insertion_below_half_the_coverage_is_dropped():
    assert find_called_base('N', 'A+2CGA+2CGAAA', 5) == 'A'


def test_most_common_insertion_wins_and_ties_go_to_the_last_alphabetically():
    assert find_called_base('N', 'A+1CA+1CA+1G', 3) == 'AC'
    assert find_called_base('N', 'A+1CA+1G', 2) == 'AG'


def test_insertion_after_a_deleted_position_is_kept():
    assert find_called_base('N', '*+2TT*+2TT', 2) == 'TT'


def test_coverage_argument_not_the_token_count_sets_the_insertion_threshold():
    # Same bases; only the reported coverage changes.
    assert find_called_base('N', 'A+1GA', 2) == 'AG'
    assert find_called_base('N', 'A+1GA', 3) == 'A'


def test_scaffolds_break_at_any_uncovered_position(tmp_path):
    """max_N_string_length is 1, so a gap is never filled with Ns: the
    N-padding branch in the main loop is unreachable."""
    mpileup = tmp_path / 'x.mpileup'
    mpileup.write_text(''.join(
        'ref\t{}\tN\t1\t{}\t~\n'.format(pos, base)
        for pos, base in [(1, 'A'), (2, 'C'), (4, 'G'), (5, 'T')]) +
        'other\t6\tN\t1\tA\t~\n')
    out = subprocess.run([sys.executable, str(SCRIPTS_DIR / 'scaffold_contigs.py'), '-i', str(mpileup)],
                         capture_output=True, text=True, check=True).stdout
    assert out == '>1\nAC\n>2\nGT\n>3\nA\n'
