"""End-to-end paired-end test: MATAM results must not depend on read naming.

The same read pairs are run three times, named in the CASAVA style (gzipped),
the SRA style and the legacy ``/1``, ``/2`` style. Mates are the counting
unit, so every variant must give the same overlap graph, the same pairs
handed to SGA, the same scaffolds and the same abundances.
"""

import glob
import gzip
import multiprocessing
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


def run_matam(base, db, style, label, *extra_args):
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
