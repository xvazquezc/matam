import sys
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).parent.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS_DIR))
import filter_score_multialign as fsm  # noqa: E402
sys.path.pop(0)

R1 = 0x1 | 0x40
R2 = 0x1 | 0x80


def sam(qname, flag, ref, score):
    return '\t'.join([qname, str(flag), ref, '1', '255', '100M', '*', '0', '0',
                      'A' * 100, '*', 'AS:i:{}'.format(score), 'NM:i:0']) + '\n'


def kept(lines, threshold=0.9, geometric=False):
    return [(l.split('\t')[0], int(l.split('\t')[1]), l.split('\t')[2])
            for l in fsm.filter_alignments(lines, threshold, geometric)]


@pytest.mark.parametrize('geometric', [False, True])
def test_weaker_mate_is_not_dropped_by_the_stronger_mate(geometric):
    # Verified SortMeRNA case: R1 AS 300, R2 AS 240 (< 0.9 * 300), same QNAME.
    lines = [sam('pair', R1, 'refA', 300), sam('pair', R2 | 0x10, 'refA', 240)]

    assert kept(lines, geometric=geometric) == [('pair', R1, 'refA'), ('pair', R2 | 0x10, 'refA')]


def test_straight_threshold_is_relative_to_the_best_alignment_of_each_mate():
    lines = [sam('pair', R1, 'r1', 100), sam('pair', R1, 'r2', 92), sam('pair', R1, 'r3', 85),
             sam('pair', R2, 'r4', 50), sam('pair', R2, 'r5', 46), sam('pair', R2, 'r6', 44)]

    assert [r for _, _, r in kept(lines)] == ['r1', 'r2', 'r4', 'r5']


def test_geometric_threshold_follows_the_previous_kept_score():
    # 100 -> 92 -> 85 is kept geometrically (85 >= 0.9 * 92) but not straight.
    lines = [sam('read', 0, 'r3', 85), sam('read', 0, 'r1', 100), sam('read', 0, 'r2', 92),
             sam('read', 0, 'r4', 70)]

    assert [r for _, _, r in kept(lines, geometric=True)] == ['r1', 'r2', 'r3']
    assert [r for _, _, r in kept(lines, geometric=False)] == ['r1', 'r2']


def test_unpaired_reads_and_headers_are_handled():
    lines = ['@HD\tVN:1.0\n', sam('a', 0, 'r1', 10), sam('b', 16, 'r1', 20)]

    assert kept(lines) == [('a', 0, 'r1'), ('b', 16, 'r1')]


def test_missing_score_is_an_error():
    line = '\t'.join(['a', '0', 'r1', '1', '255', '4M', '*', '0', '0', 'ACGT', '*']) + '\n'
    with pytest.raises(ValueError, match='no AS:i score'):
        list(fsm.filter_alignments([line], 0.9))
