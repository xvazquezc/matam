"""End-to-end paired-end test: MATAM results must not depend on read naming.

The same read pairs are run three times, named in the CASAVA style (gzipped),
the SRA style and the legacy ``/1``, ``/2`` style. Mates are the counting
unit, so every variant must give the same overlap graph, the same pairs
handed to SGA, the same scaffolds and the same abundances.
"""

import glob
import gzip
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from pe_fixtures import NAMING_STYLES, ROOT_DIR, SCRIPTS_DIR, build_mini_db, tool_env, write_pe_reads

OVGRAPHBUILD = ROOT_DIR / 'ovgraphbuild' / 'bin' / 'ovgraphbuild'
COMPONENTSEARCH = ROOT_DIR / 'componentsearch' / 'componentsearch'
# Keep every 4th pair: ~1,450 pairs (~2,900 mates). The overlap graph is
# quadratic in mates, so the full example takes minutes per run.
PAIR_SUBSAMPLING = 4


def _missing_requirements():
    missing = [str(p) for p in (OVGRAPHBUILD, COMPONENTSEARCH) if not p.is_file()]
    env_path = tool_env()['PATH']
    missing += [b for b in ('sortmerna', 'sga', 'samtools', 'vsearch', 'seqkit', 'ktImportText')
                if shutil.which(b, path=env_path) is None]
    return missing


pytestmark = pytest.mark.skipif(
    bool(_missing_requirements()),
    reason='missing MATAM binaries: {}'.format(', '.join(_missing_requirements())))


def single(pattern):
    matches = glob.glob(pattern)
    assert len(matches) == 1, 'expected one match for {}: {}'.format(pattern, matches)
    return Path(matches[0])


def read_fasta(path):
    records, header, seq = [], None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq)))
                header, seq = line[1:], []
            elif line:
                seq.append(line)
    if header is not None:
        records.append((header, ''.join(seq)))
    return records


def summarise_run(out_dir, log_text):
    workdir = out_dir / 'workdir'
    nodes = single(str(workdir / '*.ovgb_*.nodes.csv'))
    edges = single(str(workdir / '*.ovgb_*.edges.csv'))
    components = single(str(workdir / '*.components.csv'))
    with open(nodes) as handle:
        node_header = handle.readline().rstrip('\n').split(';')
        node_rows = [line.rstrip('\n').split(';') for line in handle if line.strip()]
    with open(edges) as handle:
        edge_count = sum(1 for line in handle if line.strip()) - 1
    with open(components) as handle:
        component_header = handle.readline().rstrip('\n').split(';')
        component_ids = {line.rstrip('\n').split(';')[-1] for line in handle if line.strip()}
    kept = re.search(r'Read pairs by component: (\d+) within one component, (\d+) spanning two \(given to both\), (\d+) with one mate in a component, (\d+) in no component', log_text)
    scaffolds = read_fasta(out_dir / 'final_assembly.fa')
    abundance = sorted(
        (seq, float(re.search(r'count=([0-9.]+)', header).group(1)))
        for header, seq in scaffolds)
    return {
        'node_header': node_header,
        'node_rows': node_rows,
        'node_count': len(node_rows),
        'edge_count': edge_count,
        'component_header': component_header,
        'component_count': len({c for c in component_ids if not c.startswith('-')}),
        'pair_routing': kept.groups() if kept else None,
        'marker_reads': tuple(int(n) for n in re.search(
            r'Identified as marker: (\d+) / (\d+) reads', log_text).groups()),
        'scaffold_sequences': sorted(seq for _, seq in scaffolds),
        'abundance': abundance,
    }


@pytest.fixture(scope='module')
def pe_workspace(tmp_path_factory):
    base = tmp_path_factory.mktemp('matam_pe')
    return base, build_mini_db(base / 'db')


def run_matam(base, db, style, label, *extra_args, return_log=False):
    forward, reverse = write_pe_reads(base / 'reads', style, every=PAIR_SUBSAMPLING,
                                      gz=(style == 'casava'))
    out_dir = base / 'out_{}_{}'.format(label, style)
    completed = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / 'matam_assembly.py'),
         '--forward', str(forward), '--reverse', str(reverse),
         '-d', str(db), '-o', str(out_dir),
         '--cpu', str(min(8, multiprocessing.cpu_count())), '--max_memory', '2000', '-v',
         # ComponentSearch seeds from the clock by default.
         '--seed', '42', *extra_args],
        env=tool_env(), capture_output=True, text=True)
    log_text = completed.stdout + completed.stderr
    assert completed.returncode == 0, '{} run failed:\n{}'.format(style, log_text[-4000:])
    if return_log:
        return summarise_run(out_dir, log_text), log_text
    return summarise_run(out_dir, log_text)


@pytest.fixture(scope='module')
def pe_runs(pe_workspace):
    base, db = pe_workspace
    return {style: run_matam(base, db, style, 'default') for style in NAMING_STYLES}


@pytest.fixture(scope='module')
def pe_runs_coverage(pe_workspace):
    """Runs with reference-coverage subsampling, which re-sorts the tagged SAM."""
    base, db = pe_workspace
    return {style: run_matam(base, db, style, 'coverage', '--coverage_threshold', '20')
            for style in ('casava', 'legacy')}


def test_read_counts_are_in_mates(pe_runs):
    for style, run in pe_runs.items():
        selected, total = run['marker_reads']
        assert total == 2 * len(range(0, 5825, PAIR_SUBSAMPLING)), style
        assert 0 < selected <= total, style


def test_every_naming_style_produces_scaffolds(pe_runs):
    for style, run in pe_runs.items():
        assert run['scaffold_sequences'], '{} produced no scaffolds'.format(style)


def test_graph_nodes_are_mates_with_mate_column(pe_runs):
    for style, run in pe_runs.items():
        assert run['node_header'] == ['Id', 'Label', 'Mate', 'Specie'], style
        assert run['component_header'][:4] == ['Id', 'Label', 'Mate', 'Specie'], style
        node_keys = [(row[1], row[2]) for row in run['node_rows']]
        assert all(mate in ('1', '2') for _, mate in node_keys), style
        assert len(node_keys) == len(set(node_keys)), '{}: duplicate (read, mate) nodes'.format(style)


@pytest.mark.parametrize('style', [s for s in NAMING_STYLES if s != 'legacy'])
def test_graph_is_independent_of_read_naming(pe_runs, style):
    reference, run = pe_runs['legacy'], pe_runs[style]
    assert run['node_count'] == reference['node_count']
    assert run['edge_count'] == reference['edge_count']
    assert run['component_count'] == reference['component_count']
    assert run['pair_routing'] == reference['pair_routing']


@pytest.mark.parametrize('style', [s for s in NAMING_STYLES if s != 'legacy'])
def test_scaffolds_and_abundance_are_independent_of_read_naming(pe_runs, style):
    reference, run = pe_runs['legacy'], pe_runs[style]
    assert run['scaffold_sequences'] == reference['scaffold_sequences']
    assert run['abundance'] == reference['abundance']


def test_coverage_subsampling_keeps_reads_whole_and_naming_independent(pe_runs_coverage, pe_runs):
    casava, legacy = pe_runs_coverage['casava'], pe_runs_coverage['legacy']
    node_keys = [(row[1], row[2]) for row in casava['node_rows']]
    assert len(node_keys) == len(set(node_keys))
    assert casava['node_count'] < pe_runs['casava']['node_count']  # subsampling removed reads
    for key in ('node_count', 'edge_count', 'component_count', 'pair_routing',
                'scaffold_sequences', 'abundance'):
        assert casava[key] == legacy[key], key


def fastq_records(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        lines = [line.rstrip('\n') for line in handle]
    return [(lines[i][1:].split()[0], lines[i + 1]) for i in range(0, len(lines), 4)]


def test_filter_only_publishes_selected_pairs_as_r1_r2_files(pe_workspace):
    """--filter_only output contract used by MarkerMAG-nf."""
    base, db = pe_workspace
    forward, reverse = write_pe_reads(base / 'reads_filter', 'casava', every=PAIR_SUBSAMPLING, gz=True)
    out_dir = base / 'out_filter_only'
    completed = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / 'matam_assembly.py'),
         '--forward', str(forward), '--reverse', str(reverse),
         '-d', str(db), '-o', str(out_dir), '--cpu', '2', '--max_memory', '2000',
         '--filter_only'],
        env=tool_env(), capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr[-4000:]

    selected_r1 = fastq_records(out_dir / 'marker_reads_R1.fq.gz')
    selected_r2 = fastq_records(out_dir / 'marker_reads_R2.fq.gz')
    assert selected_r1 and len(selected_r1) == len(selected_r2)
    assert [r for r, _ in selected_r1] == [r for r, _ in selected_r2]   # in step
    input_r1 = set(fastq_records(forward))
    input_r2 = set(fastq_records(reverse))
    assert set(selected_r1) <= input_r1 and set(selected_r2) <= input_r2  # mates not swapped
    assert not (out_dir / 'final_assembly.fa').exists()


def _matam(forward, reverse, db, out_dir, *extra_args, keep_tmp=True):
    completed = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / 'matam_assembly.py'),
         '--forward', str(forward), '--reverse', str(reverse),
         '-d', str(db), '-o', str(out_dir), '--cpu', str(min(8, multiprocessing.cpu_count())),
         '--max_memory', '2000', '-v', '--seed', '42', *(['--keep_tmp'] if keep_tmp else []),
         *extra_args],
        env=tool_env(), capture_output=True, text=True)
    return completed.returncode, completed.stdout + completed.stderr


REUSED = 'Reusing complete SortMeRNA output'
CLEARED = 'Removing incomplete or outdated SortMeRNA output'


@pytest.fixture(scope='module')
def resumable_run(pe_workspace):
    """A finished run with intermediate files kept, for resuming in place (#13)."""
    base, db = pe_workspace
    forward, reverse = write_pe_reads(base / 'reads_resume', 'casava', every=PAIR_SUBSAMPLING, gz=True)
    out_dir = base / 'out_resume'
    code, log = _matam(forward, reverse, db, out_dir)
    assert code == 0, log[-4000:]
    return forward, reverse, db, out_dir, (out_dir / 'final_assembly.fa').read_bytes()


@pytest.mark.parametrize('step, expected_sortmerna_logs', [
    # Same reads: the read mapping is reused; contigs are rebuilt, so the
    # scaffolding mapping is cleared and redone.
    ('reads_mapping', [REUSED, CLEARED]),
    ('contigs_assembly', [CLEARED]),
    # Same contigs: the scaffolding mapping is reused.
    ('scaffolding', [REUSED]),
])
def test_resume_in_place_reuses_or_redoes_sortmerna_steps(resumable_run, step, expected_sortmerna_logs):
    # Before #13 every resume up to scaffolding failed with SortMeRNA's
    # 'kvdb resume aborted: inputs mismatch'.
    forward, reverse, db, out_dir, final_assembly = resumable_run
    (out_dir / 'final_assembly.fa').unlink()

    code, log = _matam(forward, reverse, db, out_dir, '--resume_from', step)

    assert code == 0, log[-4000:]
    assert re.findall('|'.join((REUSED, CLEARED)), log) == expected_sortmerna_logs
    assert (out_dir / 'final_assembly.fa').read_bytes() == final_assembly


def test_rerun_into_the_same_directory_with_other_reads_redoes_the_mapping(resumable_run, tmp_path):
    forward, reverse, db, out_dir, _ = resumable_run
    # Same file names (so the same SortMeRNA output names), other read pairs
    other_forward, other_reverse = write_pe_reads(tmp_path / 'reads_resume', 'casava',
                                                  every=2 * PAIR_SUBSAMPLING, gz=True)
    assert other_forward.name == forward.name
    rerun_dir = tmp_path / 'out'
    shutil.copytree(out_dir, rerun_dir, symlinks=True, copy_function=shutil.copy2)

    code, log = _matam(other_forward, other_reverse, db, rerun_dir)

    assert code == 0, log[-4000:]
    assert re.findall('|'.join((REUSED, CLEARED)), log) == [CLEARED, CLEARED]
    selected = fastq_records(rerun_dir / 'marker_reads_R1.fq.gz')
    assert {read for read, _ in selected} <= {read for read, _ in fastq_records(other_forward)}


def test_resume_past_an_incomplete_read_mapping_is_a_clear_error(resumable_run, tmp_path):
    forward, reverse, db, out_dir, _ = resumable_run
    broken_dir = tmp_path / 'broken'
    shutil.copytree(out_dir, broken_dir, symlinks=True, copy_function=shutil.copy2)
    sam = next((broken_dir / 'workdir').glob('*.sortmerna_vs_*_b10_m10.sam.gz'))
    sam.write_bytes(sam.read_bytes()[:100])   # truncated, as if interrupted

    code, log = _matam(forward, reverse, db, broken_dir, '--resume_from', 'alignments_filtering')

    assert code != 0
    assert 'is incomplete or was made with other inputs or options' in log
    assert '--resume_from reads_mapping' in log


# SortMeRNA reference index: prebuilt by matam_db_preprocessing.py, or, with
# --no_require_index, built by matam_assembly.py when the database was made with --skip_index.
PREBUILT = 'using the prebuilt index'
BUILT_IN_DB = 'Building it there now; later runs will reuse it'
BUILT_IN_RUN = 'the database directory is not writable. Building a private index for this run'


def test_mapping_uses_the_prebuilt_index(pe_workspace):
    base, db = pe_workspace
    forward, reverse = write_pe_reads(base / 'reads_prebuilt', 'casava', every=PAIR_SUBSAMPLING, gz=True)
    code, log = _matam(forward, reverse, db, base / 'out_prebuilt', '--filter_only')
    assert code == 0, log[-4000:]
    assert 'SortMeRNA index of the clustered ref db: ' + PREBUILT in log
    assert BUILT_IN_DB not in log and BUILT_IN_RUN not in log
    assert not (base / 'out_prebuilt' / 'sortmerna_reads_mapping' / 'idx').exists()


def test_unindexed_writable_db_is_indexed_in_place_once(tmp_path):
    db = build_mini_db(tmp_path / 'db', extra_args=['--skip_index'])
    idx_dir = Path(str(db) + '.clustered.smr') / 'idx'
    assert not idx_dir.exists()
    forward, reverse = write_pe_reads(tmp_path / 'reads', 'casava', every=PAIR_SUBSAMPLING, gz=True)

    code, log = _matam(forward, reverse, db, tmp_path / 'out1', '--filter_only', '--no_require_index')
    assert code == 0, log[-4000:]
    assert BUILT_IN_DB in log
    assert (idx_dir / 'matam_index.json').is_file()

    code, log = _matam(forward, reverse, db, tmp_path / 'out2', '--filter_only', '--no_require_index')
    assert code == 0, log[-4000:]
    assert PREBUILT in log and BUILT_IN_DB not in log
    assert (fastq_records(tmp_path / 'out1' / 'marker_reads_R1.fq.gz')
            == fastq_records(tmp_path / 'out2' / 'marker_reads_R1.fq.gz'))


def _tree(directory):
    return sorted((str(p.relative_to(directory)), p.stat().st_size if p.is_file() else None)
                  for p in directory.rglob('*'))


@pytest.mark.skipif(hasattr(os, 'geteuid') and os.geteuid() == 0, reason='root ignores permissions')
def test_unindexed_read_only_db_gives_the_same_assembly_with_private_indexes(pe_runs, tmp_path_factory):
    base = tmp_path_factory.mktemp('matam_pe_readonly_db')
    db = build_mini_db(base / 'db', extra_args=['--skip_index'])
    db_dir = db.parent
    # Not about SortMeRNA: scaffolding runs 'samtools view -T <db>.complete.fasta',
    # which writes <db>.complete.fasta.fai next to the fasta if it is missing.
    subprocess.run(['samtools', 'faidx', str(db) + '.complete.fasta'], check=True, env=tool_env())
    before = _tree(db_dir)
    for path in [db_dir] + [p for p in db_dir.iterdir() if p.is_dir()]:
        path.chmod(0o555)
    try:
        summary, log = run_matam(base, db, 'casava', 'readonly_db', '--no_require_index', return_log=True)
    finally:
        for path in [db_dir] + [p for p in db_dir.iterdir() if p.is_dir()]:
            path.chmod(0o755)

    assert log.count(BUILT_IN_RUN) == 2   # clustered (reads) and complete (scaffolding)
    out_dir = base / 'out_readonly_db_casava'
    for step in ('sortmerna_reads_mapping', 'sortmerna_scaffolding'):
        assert (out_dir / step / 'idx' / 'matam_index.json').is_file()
    assert _tree(db_dir) == before
    # Same assembly as with the prebuilt indexes (same indexing parameters).
    assert summary['scaffold_sequences'] == pe_runs['casava']['scaffold_sequences']
    assert summary['abundance'] == pe_runs['casava']['abundance']


# Default (no --no_require_index): use only the prebuilt indexes, else stop before any step.
REQUIRE_ERROR = 'No usable prebuilt SortMeRNA index (see --no_require_index)'


def test_by_default_only_the_prebuilt_indexes_are_used(pe_workspace):
    base, db = pe_workspace
    forward, reverse = write_pe_reads(base / 'reads_require', 'casava', every=PAIR_SUBSAMPLING, gz=True)
    code, log = _matam(forward, reverse, db, base / 'out_require')
    assert code == 0, log[-4000:]
    assert re.search(r'MATAM version: \S+', log)          # provenance in every run log
    assert 'SortMeRNA index of the clustered ref db: ' + PREBUILT in log
    assert 'SortMeRNA index of the complete ref db: ' + PREBUILT in log
    assert BUILT_IN_DB not in log and BUILT_IN_RUN not in log
    assert (base / 'out_require' / 'final_assembly.fa').exists()


@pytest.mark.parametrize('extra', [[], ['--filter_only']])
def test_by_default_an_unindexed_db_stops_the_run_before_any_step(tmp_path, extra):
    db = build_mini_db(tmp_path / 'db', extra_args=['--skip_index'])
    forward, reverse = write_pe_reads(tmp_path / 'reads', 'casava', every=PAIR_SUBSAMPLING, gz=True)
    before = _tree(db.parent)

    code, log = _matam(forward, reverse, db, tmp_path / 'out', *extra)

    assert code != 0
    assert REQUIRE_ERROR in log
    assert 'index_ref_db.py -i {}'.format(db) in log
    assert '=== Reads mapping' not in log and '=== Input' not in log
    assert not (tmp_path / 'out' / 'sortmerna_reads_mapping').exists()
    assert _tree(db.parent) == before            # nothing indexed in the db


def test_by_default_the_scaffolding_index_is_checked_before_the_read_mapping(tmp_path):
    db = build_mini_db(tmp_path / 'db')
    shutil.rmtree(str(db) + '.complete.smr')     # only the complete db (scaffolding) unindexed
    forward, reverse = write_pe_reads(tmp_path / 'reads', 'casava', every=PAIR_SUBSAMPLING, gz=True)

    code, log = _matam(forward, reverse, db, tmp_path / 'out')
    assert code != 0 and REQUIRE_ERROR in log
    assert 'complete ref db' in log and '=== Reads mapping' not in log

    # --filter_only maps against the clustered db only.
    code, log = _matam(forward, reverse, db, tmp_path / 'out_filter', '--filter_only')
    assert code == 0, log[-4000:]
    assert (tmp_path / 'out_filter' / 'marker_reads_R1.fq.gz').exists()


def test_index_ref_db_makes_an_unindexed_db_usable(tmp_path):
    db = build_mini_db(tmp_path / 'db', extra_args=['--skip_index'])
    subprocess.run([sys.executable, str(SCRIPTS_DIR / 'index_ref_db.py'), '-i', str(db),
                    '--max_memory', '2000'], check=True, env=tool_env(), capture_output=True)
    forward, reverse = write_pe_reads(tmp_path / 'reads', 'casava', every=PAIR_SUBSAMPLING, gz=True)

    code, log = _matam(forward, reverse, db, tmp_path / 'out', '--filter_only')
    assert code == 0, log[-4000:]
    assert PREBUILT in log


# --stop_after (MarkerMAG-nf stops after scaffolding: it doesn't use abundances).

def _scaffolds(path):
    return sorted(seq for _, seq in read_fasta(path))


def _abundance(path):
    return sorted((seq, float(re.search(r'count=([0-9.]+)', header).group(1)))
                  for header, seq in read_fasta(path))


def test_stop_after_scaffolding_gives_the_scaffolds_and_resumes_to_the_abundances(pe_workspace, pe_runs):
    base, db = pe_workspace
    # The reads of pe_runs['casava'].
    forward, reverse = write_pe_reads(base / 'reads_stop', 'casava', every=PAIR_SUBSAMPLING, gz=True)
    out_dir = base / 'out_stop_scaffolding'

    code, log = _matam(forward, reverse, db, out_dir, '--stop_after', 'scaffolding', keep_tmp=False)

    assert code == 0, log[-4000:]
    assert '=== Abundance calculation' not in log
    assert 'Stopping after scaffolding (--stop_after)' in log
    final_assembly = out_dir / 'final_assembly.fa'
    assert os.readlink(str(final_assembly)) == 'workdir/scaffolds.NR.min_500bp.fa'
    assert _scaffolds(final_assembly) == pe_runs['casava']['scaffold_sequences']
    assert not any('count=' in header for header, _ in read_fasta(final_assembly))
    workdir = out_dir / 'workdir'
    assert not list(workdir.glob('*.abd.fa'))
    assert not list(workdir.glob('*.mpileup'))           # tmp files removed as usual

    # The run continues where it stopped.
    code, log = _matam(forward, reverse, db, out_dir, '--resume_from', 'abundance_calculation',
                       keep_tmp=False)
    assert code == 0, log[-4000:]
    assert os.readlink(str(final_assembly)).endswith('.abd.fa')
    assert _abundance(final_assembly) == pe_runs['casava']['abundance']


def test_stop_after_an_early_step_keeps_what_the_next_steps_read(pe_workspace, pe_runs_coverage):
    base, db = pe_workspace
    forward, reverse = write_pe_reads(base / 'reads_stop_early', 'casava', every=PAIR_SUBSAMPLING,
                                      gz=True)
    out_dir = base / 'out_stop_early'
    # With --coverage_threshold the filtered SAM is a tmp file, also read by
    # the LCA labelling (graph_compaction).
    coverage = ('--coverage_threshold', '20')

    code, log = _matam(forward, reverse, db, out_dir, *coverage, '--stop_after',
                       'alignments_filtering', keep_tmp=False)

    assert code == 0, log[-4000:]
    assert 'continue with --resume_from overlap_graph_building' in log
    assert '=== Overlap-graph building' not in log
    workdir = out_dir / 'workdir'
    assert list(workdir.glob('*.scr_filt_geo_90pct.sam'))
    assert not list(workdir.glob('*.nodes.csv')) and not (out_dir / 'final_assembly.fa').exists()

    code, log = _matam(forward, reverse, db, out_dir, *coverage, '--resume_from',
                       'overlap_graph_building', keep_tmp=False)
    assert code == 0, log[-4000:]
    assert _scaffolds(out_dir / 'final_assembly.fa') == pe_runs_coverage['casava']['scaffold_sequences']
    assert not list(workdir.glob('*.scr_filt_geo_90pct.sam'))   # removed at the end


def test_stop_before_scaffolding_does_not_need_the_complete_index(tmp_path):
    db = build_mini_db(tmp_path / 'db')
    shutil.rmtree(str(db) + '.complete.smr')
    forward, reverse = write_pe_reads(tmp_path / 'reads', 'casava', every=PAIR_SUBSAMPLING, gz=True)

    code, log = _matam(forward, reverse, db, tmp_path / 'out', '--stop_after', 'contigs_assembly')

    assert code == 0, log[-4000:]
    assert 'continue with --resume_from scaffolding' in log
    assert '=== Scaffolding' not in log


@pytest.mark.parametrize('options, message', [
    (['--resume_from', 'scaffolding', '--stop_after', 'contigs_assembly'],
     '--stop_after contigs_assembly comes before --resume_from scaffolding'),
    (['--filter_only', '--stop_after', 'scaffolding'], 'do not combine it with --stop_after'),
])
def test_stop_after_arguments_are_checked(tmp_path, options, message):
    code, log = _matam(tmp_path / 'r1.fq', tmp_path / 'r2.fq', tmp_path / 'db', tmp_path / 'out',
                       *options)
    assert code == 2 and message in log


def test_seed_zero_seeds_the_graph_compaction_and_is_logged(pe_workspace, tmp_path):
    # 'if args.seed:' used to treat --seed 0 as no seed (clock-seeded ComponentSearch).
    base, db = pe_workspace
    forward, reverse = write_pe_reads(base / 'reads_seed0', 'casava', every=PAIR_SUBSAMPLING, gz=True)
    code, log = _matam(forward, reverse, db, tmp_path / 'out', '--seed', '0', '--debug',
                       '--stop_after', 'graph_compaction')
    assert code == 0, log[-4000:]
    assert re.search(r'CMD: \S*componentsearch .*-r 0 ', log)
    assert re.search(r'matam_assembly\.py .*--seed 0 ', log)     # the logged command line
