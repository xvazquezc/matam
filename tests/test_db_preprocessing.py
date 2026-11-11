"""matam_db_preprocessing: reference headers with tabs (#9), SortMeRNA indexes."""

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
from pe_fixtures import (EXAMPLE_REFS, SCRIPTS_DIR, build_mini_db, tool_env, write_pe_reads,  # noqa: E402
                         write_reference_fasta)

TOOLS = ('sortmerna', 'vsearch', 'seqkit')
MISSING = [tool for tool in TOOLS if shutil.which(tool, path=tool_env()['PATH']) is None]
pytestmark = pytest.mark.skipif(bool(MISSING), reason='missing: {}'.format(', '.join(MISSING)))


@pytest.fixture(scope='module')
def tab_header_db(tmp_path_factory):
    # The bundled references use '>ID<TAB>1459bp<TAB>Species' headers.
    base = tmp_path_factory.mktemp('db_tabs')
    assert '\t' in EXAMPLE_REFS.read_text().splitlines()[0]
    subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / 'matam_db_preprocessing.py'),
         '-i', str(EXAMPLE_REFS), '-d', str(base), '--cpu', '2',
         '--max_memory', '2000', '--out_db_name', 'tabs'],
        check=True, env=tool_env(), capture_output=True, text=True)
    return base


def headers(path):
    return [line[1:].rstrip('\n') for line in open(path) if line.startswith('>')]


def test_complete_db_headers_have_no_tabs_and_keep_the_description(tab_header_db):
    complete = headers(tab_header_db / 'tabs.complete.fasta')
    assert len(complete) == 16
    assert not any('\t' in header or header != header.strip() for header in complete)
    assert 'ABA_X81660 1459bp Acinetobacter baumannii' in complete


def test_reference_ids_agree_across_db_files(tab_header_db):
    complete_ids = {header.split(' ')[0] for header in headers(tab_header_db / 'tabs.complete.fasta')}
    taxo_ids = {line.split('\t')[0] for line in open(tab_header_db / 'tabs.complete.taxo.tab')}
    clustered_ids = set(headers(tab_header_db / 'tabs.clustered.fasta'))
    assert complete_ids == taxo_ids
    assert clustered_ids <= complete_ids


def test_sortmerna_rname_is_the_reference_id(tab_header_db, tmp_path):
    # Scaffolding maps against the complete DB; a tab in the header used to
    # end up inside RNAME and shift every following SAM column.
    pairs = write_pe_reads(tmp_path, 'legacy', every=200)
    prefix = tmp_path / 'aligned'
    subprocess.run(
        ['sortmerna', '--ref', str(tab_header_db / 'tabs.complete.fasta'),
         '--workdir', str(tmp_path / 'wd'), '--reads', str(pairs[0]),
         '--aligned', str(prefix), '--sam', '--num_alignments', '1', '--threads', '2'],
        check=True, env=tool_env(), capture_output=True, text=True)
    taxo_ids = {line.split('\t')[0] for line in open(tab_header_db / 'tabs.complete.taxo.tab')}
    records = [line.rstrip('\n').split('\t') for line in open(str(prefix) + '.sam')
               if not line.startswith('@')]
    assert records
    for fields in records:
        assert fields[2] in taxo_ids
        assert fields[3].isdigit(), fields[:6]                       # POS
        assert re.fullmatch(r'(\d+[MIDNSHP=X])+', fields[5]), fields[:6]  # CIGAR


def _sortmerna_index():
    sys.path.insert(0, str(SCRIPTS_DIR))
    try:
        import sortmerna_index
    finally:
        sys.path.pop(0)
    return sortmerna_index


@pytest.mark.parametrize('kind', ['complete', 'clustered'])
def test_db_indexes_are_prebuilt_and_recorded_with_the_max_memory(tab_header_db, kind):
    sortmerna_index = _sortmerna_index()
    smr_dir = tab_header_db / 'tabs.{}.smr'.format(kind)
    params = sortmerna_index.db_params(str(smr_dir))
    assert params == sortmerna_index.index_params(2000)   # -m is --max_memory
    assert sortmerna_index.check(str(smr_dir / 'idx'), str(tab_header_db / 'tabs.{}.fasta'.format(kind)),
                                 shutil.which('sortmerna', path=tool_env()['PATH']), params) == (True, '')
    # No SortMeRNA run state (kvdb/, out/) left in the database.
    assert sorted(os.listdir(str(smr_dir))) == ['.matam_index.lock', 'idx', 'matam_index_params.json']


def test_skip_index_builds_no_index_but_stores_the_indexing_parameters(tmp_path):
    sortmerna_index = _sortmerna_index()
    db = build_mini_db(tmp_path / 'db', extra_args=['--skip_index'])
    for kind in ('complete', 'clustered'):
        smr_dir = Path(str(db) + '.{}.smr'.format(kind))
        assert Path(str(db) + '.{}.fasta'.format(kind)).is_file()
        assert os.listdir(str(smr_dir)) == ['matam_index_params.json']
        assert sortmerna_index.db_params(str(smr_dir)) == sortmerna_index.index_params(2000)


def _fasta(path):
    """{id: sequence} of a fasta file, in file order."""
    records, seq_id = {}, None
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('>'):
            seq_id = line[1:].split()[0]
            records[seq_id] = ''
        else:
            records[seq_id] += line
    return records


def _preprocess(tmp_path, *args):
    return subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / 'matam_db_preprocessing.py'),
         '-i', str(EXAMPLE_REFS), '-d', str(tmp_path / 'db'), '--cpu', '2',
         '--max_memory', '2000', '-v', '--skip_index', *args],
        env=tool_env(), capture_output=True, text=True)


def test_clustered_ref_replaces_clustering_with_the_cleaned_sequences_of_its_ids(tmp_path):
    refs = _fasta(EXAMPLE_REFS)
    chosen = list(refs)[3:8]
    # A clustered set made elsewhere: other headers, RNA alphabet, lower case,
    # its own order, and an ID the reference db doesn't have.
    clustered = tmp_path / 'SILVA_X_NR95.fasta'
    clustered.write_text(''.join('>{} centroid\n{}\n'.format(seq_id, refs[seq_id].lower().replace('t', 'u'))
                                 for seq_id in reversed(chosen)) + '>NOT_IN_DB\nACGT\n')

    result = _preprocess(tmp_path, '--clustered_ref', str(clustered))

    assert result.returncode == 0, result.stderr
    db = tmp_path / 'db' / 'SILVA_X_NR95'                     # named after the clustered set
    complete = _fasta(str(db) + '.complete.fasta')
    selected = _fasta(str(db) + '.clustered.fasta')
    assert list(selected) == [seq_id for seq_id in complete if seq_id in chosen]  # db order
    assert all(selected[seq_id] == complete[seq_id] for seq_id in selected)      # cleaned sequences
    assert headers(str(db) + '.clustered.fasta') == list(selected)                # IDs only
    assert 'Clustering sequences' not in result.stderr                            # no vsearch
    assert '1 of the 6 clustered sequences are not in the cleaned reference db' in result.stderr
    assert 'NOT_IN_DB' in result.stderr
    assert not list((tmp_path / 'db').glob('*.clustered_ref.*'))                 # tmp files removed


def test_clustered_ref_without_any_reference_id_is_an_error(tmp_path):
    clustered = tmp_path / 'other.fasta'
    clustered.write_text('>X1\nACGT\n>X2\nACGT\n')

    result = _preprocess(tmp_path, '--clustered_ref', str(clustered))

    assert result.returncode != 0
    assert 'none of the 2 IDs' in result.stderr


@pytest.mark.parametrize('other', [['--clustering_id_threshold', '0.97'], ['--by_domain']])
def test_clustered_ref_cannot_be_combined_with_clustering_options(tmp_path, other):
    clustered = tmp_path / 'c.fasta'
    clustered.write_text('>X1\nACGT\n')

    result = _preprocess(tmp_path, '--clustered_ref', str(clustered), *other)

    assert result.returncode != 0
    assert '--clustered_ref replaces the clustering' in result.stderr


@pytest.mark.parametrize('options', [['--by_domain'], ['--by_kingdom']])  # former name still accepted
def test_by_domain_clusters_each_domain_by_the_first_taxonomy_rank(tmp_path, options):
    refs = write_reference_fasta(tmp_path / 'refs.fasta')   # 16 Bacteria
    # A eukaryote whose name contains another domain must only be clustered as a eukaryote.
    random_seq = ''.join('ACGT'[(i * 7 + i // 3) % 4] for i in range(1400))
    with open(refs, 'a') as handle:
        handle.write('>EUK1 Eukaryota;Opisthokonta;Fungi;Ascomycota;Some bacteria endosymbiont\n')
        handle.write(random_seq + '\n')

    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / 'matam_db_preprocessing.py'),
         '-i', str(refs), '-d', str(tmp_path / 'db'), '--cpu', '2', '--max_memory', '2000',
         '-v', '--skip_index', *options],
        env=tool_env(), capture_output=True, text=True)

    assert result.returncode == 0, result.stderr
    clustered = headers(tmp_path / 'db' / 'refs_NR95_bd.clustered.fasta')
    assert clustered.count('EUK1') == 1
    assert len(clustered) == len(set(clustered))
    assert 'Clustering eukaryota sequences' in result.stderr
    assert ('deprecated' in result.stderr) == (options == ['--by_kingdom'])


def _discarded(db_basepath):
    import csv
    with open(str(db_basepath) + '.discarded.tsv') as handle:
        return {row['id']: row for row in csv.DictReader(handle, delimiter='\t')}


@pytest.fixture
def refs_with_losses(tmp_path):
    """16 Bacteria + 2 sequences the N filter removes + an archaeon.

    NRUN1 is unrelated to every other reference (own taxa, random sequence);
    CLOSE1 is ABA_X81660 with an N run and 3 substitutions, so a close relative remains.
    """
    import random
    refs = write_reference_fasta(tmp_path / 'refs.fasta')
    rng = random.Random(7)
    unrelated = ''.join(rng.choice('ACGT') for _ in range(1400))
    archaeon = ''.join(rng.choice('ACGT') for _ in range(1400))
    aba = ''.join(line.strip() for line in open(EXAMPLE_REFS).read().split('>ABA_X81660', 1)[1]
                  .split('>', 1)[0].splitlines()[1:])
    with open(refs, 'a') as handle:
        handle.write('>NRUN1 Bacteria;Testphylum;Testclass;Testorder;Testfamily;Testgenus;Test sp.\n')
        handle.write(unrelated[:600] + 'N' * 20 + unrelated[600:] + '\n')
        handle.write('>CLOSE1 Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;'
                     'Moraxellaceae;Acinetobacter;Acinetobacter sp.\n')
        close = list(aba[:500] + 'N' * 12 + aba[512:])
        for position in (100, 900, 1200):                # ~99.8% identity
            close[position] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[close[position]]
        handle.write(''.join(close) + '\n')
        handle.write('>ARC1 Archaea;Euryarchaeota;Methanobacteria;Methanobacteriales;'
                     'Methanobacteriaceae;Methanobrevibacter;M. smithii\n')
        handle.write(archaeon + '\n')
    return refs


def _run_preprocessing(refs, db_dir, *args):
    return subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / 'matam_db_preprocessing.py'), '-i', str(refs),
         '-d', str(db_dir), '--cpu', '2', '--max_memory', '2000', '-v', '--skip_index', *args],
        env=tool_env(), capture_output=True, text=True)


def test_clustering_reports_the_removed_sequences_and_whether_a_close_relative_is_left(
        refs_with_losses, tmp_path):
    result = _run_preprocessing(refs_with_losses, tmp_path / 'db')

    assert result.returncode == 0, result.stderr
    discarded = _discarded(tmp_path / 'db' / 'refs_NR95')
    assert {seq_id: row['reason'] for seq_id, row in discarded.items()} == {
        'NRUN1': 'removed by cleaning', 'CLOSE1': 'removed by cleaning'}
    assert discarded['NRUN1']['taxonomy'] == \
        'Bacteria;Testphylum;Testclass;Testorder;Testfamily;Testgenus;Test_sp.'
    assert discarded['NRUN1']['close_relative_left'] == 'no'
    assert (discarded['CLOSE1']['nearest_ref'], discarded['CLOSE1']['close_relative_left']) == \
        ('ABA_X81660', 'yes')
    assert float(discarded['CLOSE1']['nearest_ref_identity']) >= 99
    assert ('WARNING - Discarded reference sequences (removed by cleaning): 2 (Bacteria 2); without '
            'a close relative left in the references (>= 95% identity over >= 90% of the sequence): '
            '1, e.g. Testgenus NRUN1') in result.stderr


def test_close_relative_identity_is_configurable_and_the_search_can_be_skipped(refs_with_losses, tmp_path):
    strict = _run_preprocessing(refs_with_losses, tmp_path / 'strict', '--close_relative_id', '1.0',
                                '--out_db_name', 'db')
    skipped = _run_preprocessing(refs_with_losses, tmp_path / 'skipped', '--close_relative_id', '0',
                                 '--out_db_name', 'db')

    assert strict.returncode == skipped.returncode == 0, strict.stderr + skipped.stderr
    assert 'over >= 90% of the sequence): 2' in strict.stderr       # 100% needed
    assert 'close relatives not searched' in skipped.stderr
    assert _discarded(tmp_path / 'skipped' / 'db')['CLOSE1']['close_relative_left'] == 'NA'


def test_clustered_ref_reports_the_discarded_centroids_with_their_taxa(refs_with_losses, tmp_path):
    clustered = tmp_path / 'given_centroids.fasta'
    clustered.write_text('>NRUN1\nACGT\n>ARC1\nACGT\n>ABA_X81660\nACGT\n'
                         '>ELSEWHERE1 Bacteria;Otherphylum;C;O;F;G;S\nACGT\n')

    result = _run_preprocessing(refs_with_losses, tmp_path / 'db', '--clustered_ref', str(clustered))

    assert result.returncode == 0, result.stderr
    discarded = _discarded(tmp_path / 'db' / 'given_centroids')
    assert {seq_id: row['reason'] for seq_id, row in discarded.items()} == {
        'NRUN1': 'centroid removed by cleaning', 'CLOSE1': 'removed by cleaning',
        'ELSEWHERE1': 'centroid not in the input'}
    assert discarded['ELSEWHERE1']['taxonomy'] == 'Bacteria;Otherphylum;C;O;F;G;S'
    assert discarded['NRUN1']['close_relative_left'] == 'no'
    assert discarded['ELSEWHERE1']['close_relative_left'] == 'no'   # 4 bp: no hit
    assert 'WARNING - Discarded reference sequences (centroid removed by cleaning): 1' in result.stderr


def test_by_domain_reports_the_sequences_of_other_domains(refs_with_losses, tmp_path):
    result = _run_preprocessing(refs_with_losses, tmp_path / 'db', '--by_domain', '--domains', 'bacteria')

    assert result.returncode == 0, result.stderr
    discarded = _discarded(tmp_path / 'db' / 'refs_NR95_bd')
    assert {seq_id: row['reason'] for seq_id, row in discarded.items()} == {
        'NRUN1': 'removed by cleaning', 'CLOSE1': 'removed by cleaning', 'ARC1': 'outside --domains'}
    assert 'Discarded reference sequences (outside --domains): 1 (Archaea 1)' in result.stderr
