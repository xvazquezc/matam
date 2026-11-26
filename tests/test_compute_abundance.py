import os
import sys
import tempfile

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.join(CURRENT_DIR, '..', 'scripts')
sys.path.append(SCRIPTS_DIR)


SAMPLE_DIR = os.path.join(CURRENT_DIR, 'sample')

from compute_abundance import abundance_calculation
import pytest

@pytest.mark.parametrize('blast,expected_abundance',

    [    # Basic test
        ['scaffolds.blast',
         {'scaff1': 1.75,
          'scaff2': 3.75,
          'scaff3': 1.25,
          'scaff4': 1.25
         }
        ],
        # Test the case where a read can be found several times on the same scaffolds.
        # this behavior is tolerated but is not intented to occur often
        ['scaffolds_multiple_reads.blast',
         { '159': 0.5,
           '161': 1,
           '175': 0.5,
           '240': 3
         },
        ]
    ]
)

def test_abundance_calculation(blast, expected_abundance):
    blast_file = os.path.join(SAMPLE_DIR, blast)
    abundance = abundance_calculation(blast_file)

    assert set(abundance.keys()) == set(expected_abundance.keys())
    assert sorted(abundance.values()) == pytest.approx(sorted(expected_abundance.values()))


from compute_abundance import abundance_from_assignments, best_alignments_by_mate

R1, R2 = 0x1 | 0x40, 0x1 | 0x80


def tagged_sam(qname, flag, scaffold, score):
    return '\t'.join([qname, str(flag), scaffold, '1', '255', '100M', '*', '0', '0',
                      'A' * 100, '*', 'AS:i:{}'.format(score)]) + '\n'


def test_best_alignments_are_selected_per_mate():
    lines = ['@HD\tVN:1.0\n',
             tagged_sam('pair', R1, 's1', 200), tagged_sam('pair', R1, 's2', 199),
             tagged_sam('pair', R1, 's3', 150),
             # R2 scores far below R1 but is judged against its own best score
             tagged_sam('pair', R2, 's3', 120)]

    assert list(best_alignments_by_mate(lines, 0.99)) == [
        (('pair', 1), 's1'), (('pair', 1), 's2'), (('pair', 2), 's3')]


def test_each_mate_counts_as_one_read():
    # pairA: R1 splits between s1 and s2, R2 is on s1; pairB: both mates on s2.
    assignments = [(('pairA', 1), 's1'), (('pairA', 1), 's2'), (('pairA', 2), 's1'),
                   (('pairB', 1), 's2'), (('pairB', 2), 's2')]

    abundance = abundance_from_assignments(assignments)

    assert dict(abundance) == {'s1': 1.5, 's2': 2.5}
    assert sum(abundance.values()) == 4  # four mates


from compute_abundance import complete_fasta_with_abundance, get_abundance_from_fasta


def test_abundance_round_trip_through_the_fasta_headers(tmp_path):
    scaffolds = tmp_path / 'scaffolds.fa'
    scaffolds.write_text('>1 lca=Bacteria;X\n' + 'A' * 100 + '\n>2\nCCGG\n>3\nTT\n')
    out = tmp_path / 'scaffolds_abundance.fa'

    complete_fasta_with_abundance(str(scaffolds), str(out), {'1': 2.5, '2': 7, '3': 1 / 3})

    assert out.read_text().splitlines()[:3] == ['>1 lca=Bacteria;X count=2.5', 'A' * 80, 'A' * 20]
    assert get_abundance_from_fasta(str(out)) == {'1': 2.5, '2': 7.0, '3': 1 / 3}


def test_scaffold_without_an_abundance_gets_0(tmp_path, caplog):
    scaffolds = tmp_path / 'scaffolds.fa'
    scaffolds.write_text('>1\nACGT\n>2\nACGT\n')
    out = tmp_path / 'scaffolds_abundance.fa'

    complete_fasta_with_abundance(str(scaffolds), str(out), {'1': 4})

    assert get_abundance_from_fasta(str(out)) == {'1': 4.0, '2': 0.0}
    assert "Can't find the abundance for:2" in caplog.text


def test_reading_a_fasta_without_abundances_stops(tmp_path):
    fasta = tmp_path / 'scaffolds.fa'
    fasta.write_text('>1\nACGT\n')
    with pytest.raises(SystemExit):
        get_abundance_from_fasta(str(fasta))
