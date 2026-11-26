"""Characterization of compute_lca_from_tab.compute_lca (component LCA, --quorum)."""
import sys
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).parent.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS_DIR))
from compute_lca_from_tab import compute_lca  # noqa: E402
sys.path.pop(0)

# One inner list per read: the taxonomies of its alignments (each read weighs 1).


def test_strict_lca_is_the_common_prefix_of_every_taxonomy():
    assert compute_lca([['A;B;C'], ['A;B;D'], ['A;B;C']], 1.0) == 'A;B'


def test_strict_lca_counts_every_alignment_of_a_read():
    assert compute_lca([['A;B;C', 'A;E;F'], ['A;B;C']], 1.0) == 'A'


def test_strict_lca_without_a_common_rank_is_root():
    assert compute_lca([['Bacteria;X'], ['Archaea;Y']], 1.0) == 'Root'


def test_quorum_keeps_the_deepest_taxon_reached_by_enough_reads():
    # 2 of 3 reads reach A;B;C, 2/3 >= 0.51
    assert compute_lca([['A;B;C'], ['A;B;C'], ['A;E;F']], 0.51) == 'A;B;C'


def test_quorum_falls_back_to_a_shallower_rank():
    # A;B;C 1/3, A;B;D 1/3, A;B 2/3
    assert compute_lca([['A;B;C'], ['A;B;D'], ['A;E;F']], 0.51) == 'A;B'


def test_a_read_with_several_alignments_shares_its_weight():
    # A;B;C: 0.5 + 1 = 1.5 of 2 reads; A;B;D: 0.5
    reads = [['A;B;C', 'A;B;D'], ['A;B;C']]
    assert compute_lca(reads, 0.51) == 'A;B;C'
    assert compute_lca(reads, 0.8) == 'A;B'


def test_quorum_mode_without_any_taxon_reaching_the_quorum_is_root():
    assert compute_lca([['A;B'], ['C;D'], ['E;F']], 0.51) == 'Root'


@pytest.mark.parametrize('reads, expected', [
    ([['A;B'], ['A;C']], 'A;C'),
    ([['A;C'], ['A;B']], 'A;B'),
])
def test_tied_siblings_below_0_51_are_decided_by_input_order(reads, expected):
    """Why --quorum must be at least 0.51: at 0.5 two sibling taxa can both
    reach the quorum, and the one seen last wins."""
    assert compute_lca(reads, 0.5) == expected


def test_taxon_reaching_full_proportion_beats_its_ancestors():
    # Equal weights: the deeper taxon is listed after its ancestors and wins.
    assert compute_lca([['A;B;C'], ['A;B;C']], 0.51) == 'A;B;C'
