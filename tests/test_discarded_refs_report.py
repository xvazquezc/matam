"""Report of the reference sequences that don't reach the clustered db (discarded_refs_report.py)."""

import csv
import logging
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).parent.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS_DIR))
import discarded_refs_report as report  # noqa: E402
sys.path.pop(0)

TAXONOMY = {
    'A1': 'Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_subtilis',
    'A2': 'Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_cereus',
    'B1': 'Bacteria;Firmicutes;Bacilli;Lactobacillales;Listeriaceae;Listeria;Listeria_sp.',
    'C1': 'Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;Methanobacteriaceae;'
          'Methanobrevibacter;M._smithii',
    'E1': 'Eukaryota;Amorphea;Obazoa;Opisthokonta;Holozoa;Choanozoa;Metazoa;'
          'Animalia;Arthropoda;Insecta;Drosophila;Drosophila_melanogaster',
    'N1': 'NULL',
}


def rows(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def run(tmp_path, discarded, complete_ids=(), clustered_ids=(), nearest=None, close=95.0, caplog=None):
    summary = report.write_report(str(tmp_path / 'db'), discarded, dict(TAXONOMY), set(complete_ids),
                                  set(clustered_ids), nearest, close if nearest is not None else None)
    if caplog is not None:
        with caplog.at_level(logging.INFO):
            report.summarise(summary)
    return summary


def test_each_discarded_sequence_gets_its_nearest_remaining_reference(tmp_path):
    run(tmp_path, {'A1': report.CENTROID_REMOVED_BY_CLEANING, 'B1': report.REMOVED_BY_CLEANING,
                   'X9': report.CENTROID_NOT_IN_INPUT},
        nearest={'A1': ('A2', 97.0), 'B1': ('C1', 80.5), 'X9': (report.NO_HIT, 0.0)})

    assert rows(tmp_path / 'db.discarded.tsv') == [
        {'id': 'A1', 'reason': report.CENTROID_REMOVED_BY_CLEANING, 'taxonomy': TAXONOMY['A1'],
         'nearest_ref': 'A2', 'nearest_ref_identity': '97.0',
         'nearest_ref_lineage': 'Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus',
         'close_relative_left': 'yes'},
        {'id': 'B1', 'reason': report.REMOVED_BY_CLEANING, 'taxonomy': TAXONOMY['B1'],
         'nearest_ref': 'C1', 'nearest_ref_identity': '80.5',
         'nearest_ref_lineage': TAXONOMY['C1'].rsplit(';', 1)[0], 'close_relative_left': 'no'},
        {'id': 'X9', 'reason': report.CENTROID_NOT_IN_INPUT, 'taxonomy': 'NULL',
         'nearest_ref': '*', 'nearest_ref_identity': 'NA', 'nearest_ref_lineage': 'NULL',
         'close_relative_left': 'no'},
    ]


def test_closeness_is_by_identity_not_by_taxon_name(tmp_path):
    # A1 keeps a namesake (A2) but it is distant; B1 loses its genus name but a
    # close relative of another genus remains.
    run(tmp_path, {'A1': report.REMOVED_BY_CLEANING, 'B1': report.REMOVED_BY_CLEANING},
        complete_ids=['A2', 'C1'], nearest={'A1': ('A2', 88.0), 'B1': ('A2', 96.0)})

    taxa = {(row['level'], row['taxon']): row for row in rows(tmp_path / 'db.discarded_taxa.tsv')}
    bacillus = taxa[('6', 'Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus')]
    assert (bacillus['without_close_relative'], bacillus['left_in_complete_db']) == ('1', '1')
    listeria = taxa[('6', 'Bacteria;Firmicutes;Bacilli;Lactobacillales;Listeriaceae;Listeria')]
    assert (listeria['without_close_relative'], listeria['left_in_complete_db']) == ('0', '0')
    assert taxa[('2', 'Bacteria;Firmicutes')]['without_close_relative'] == '1'


def test_levels_follow_the_taxonomy_and_eukaryotes_stop_at_the_domain(tmp_path):
    taxonomy = {'S144': 'Bacteria;Pseudomonadati;Spirochaetota;Spirochaetia;Spirochaetales;'
                        'Borreliaceae;Borrelia;Borrelia burgdorferi',         # SILVA 144: 8 fields
                'E1': TAXONOMY['E1']}
    report.write_report(str(tmp_path / 'db'), {'S144': report.REMOVED_BY_CLEANING,
                                               'E1': report.REMOVED_BY_CLEANING}, taxonomy, set(), set())

    levels = {row['taxon'].split(';')[-1]: (row['level'], row['lineage'])
              for row in rows(tmp_path / 'db.discarded_taxa.tsv')}
    assert levels['Borreliaceae'] == ('6', 'no') and levels['Borrelia'] == ('7', 'yes')
    assert 'Borrelia burgdorferi' not in levels                    # the organism name is no taxon
    assert levels['Eukaryota'] == ('1', 'yes') and 'Drosophila' not in levels


def test_summary_warns_about_sequences_without_a_close_relative(tmp_path, caplog):
    run(tmp_path, {'A1': report.CENTROID_REMOVED_BY_CLEANING, 'B1': report.CENTROID_REMOVED_BY_CLEANING},
        nearest={'A1': ('A2', 99.0), 'B1': ('C1', 80.2)}, caplog=caplog)

    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert any('(centroid removed by cleaning): 2 (Bacteria 2); without a close relative left in '
               'the references (>= 95% identity over >= 90% of the sequence): 1, '
               'e.g. Listeria B1 (nearest: 80.2% Methanobrevibacter)' in w for w in warnings), warnings
    assert any('db.discarded_taxa.tsv' in w for w in warnings)


def test_cleaning_losses_with_close_relatives_are_info_only(tmp_path, caplog):
    run(tmp_path, {'A2': report.REMOVED_BY_CLEANING}, nearest={'A2': ('A1', 99.1)}, caplog=caplog)

    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]
    assert 'without a close relative left in the references (>= 95% identity' in caplog.text


def test_other_domains_are_reported_without_a_search(tmp_path, caplog):
    run(tmp_path, {'C1': report.OUTSIDE_DOMAINS}, nearest={}, caplog=caplog)

    assert rows(tmp_path / 'db.discarded.tsv')[0]['close_relative_left'] == 'NA'
    assert 'WARNING' in caplog.text and 'not clustered, by choice of --domains' in caplog.text


def test_search_disabled_reports_na(tmp_path, caplog):
    run(tmp_path, {'A1': report.CENTROID_REMOVED_BY_CLEANING}, caplog=caplog)

    assert rows(tmp_path / 'db.discarded.tsv')[0]['close_relative_left'] == 'NA'
    assert rows(tmp_path / 'db.discarded_taxa.tsv')[0]['without_close_relative'] == 'NA'
    assert 'close relatives not searched' in caplog.text


def test_nothing_discarded_still_writes_the_files(tmp_path, caplog):
    run(tmp_path, {}, complete_ids=TAXONOMY, clustered_ids=['A1'], caplog=caplog)

    assert rows(tmp_path / 'db.discarded.tsv') == []
    assert rows(tmp_path / 'db.discarded_taxa.tsv') == []
    assert 'Discarded reference sequences: none' in caplog.text


def test_sequences_without_taxonomy_are_counted_but_have_no_taxa(tmp_path):
    summary = run(tmp_path, {'N1': report.REMOVED_BY_CLEANING})

    assert summary['domains'][report.REMOVED_BY_CLEANING] == {'NULL': 1}
    assert rows(tmp_path / 'db.discarded_taxa.tsv') == []
