"""--read_alignments: take the read mapping from alignments made by an earlier run.

MarkerMAG-nf assembles several subsamples of one read set. The read set is
mapped once (``--filter_only``, which publishes ``marker_alignments.sam.gz``)
and every subsample run takes its mapping from those alignments.
"""

import gzip
import multiprocessing
import shutil
import subprocess
import sys

import pytest

from pe_fixtures import SCRIPTS_DIR, build_mini_db, iter_example_pairs, tool_env, write_pe_reads

sys.path.insert(0, str(SCRIPTS_DIR))
import read_alignments  # noqa: E402

PG = ('@PG\tID:sortmerna\tVN:1.0\tCL:/usr/bin/sortmerna --ref /db/mini16.clustered.fasta '
      '--reads r1.fq --reads r2.fq --paired_in --out2 --aligned /x/out --fastx --sam --blast 1 '
      '--zip-out true --num_alignments 10 -e 1.00e-05 --idx-dir /db/idx --threads 8\n')


def sam_record(qname, flag=0, ref='ref1', pos=1):
    return '{}\t{}\t{}\t{}\t255\t4M\t*\t0\t0\tACGT\tIIII\tAS:i:8\n'.format(qname, flag, ref, pos)


def write_text(path, text):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt') as handle:
        handle.write(text)
    return path


def read_text(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        return handle.read()


def fastq(names):
    return ''.join('@{}\nACGT\n+\nIIII\n'.format(name) for name in names)


def fastq_ids(path):
    return [line[1:].split()[0] for i, line in enumerate(read_text(path).splitlines()) if i % 4 == 0]


def sam_body(path):
    return [line for line in read_text(path).splitlines(True) if not line.startswith('@')]


# --- selection ------------------------------------------------------------------

@pytest.mark.parametrize('gz', [False, True])
def test_selects_the_aligned_pairs_and_their_alignments(tmp_path, gz):
    """CASAVA names: both mates share the QNAME. The SAM covers a larger read set."""
    suffix = '.gz' if gz else ''
    forward = write_text(tmp_path / ('r1.fq' + suffix), fastq(
        ['a 1:N:0:1', 'b 1:N:0:1', 'c 1:N:0:1', 'd 1:N:0:1']))
    reverse = write_text(tmp_path / ('r2.fq' + suffix), fastq(
        ['a 2:N:0:1', 'b 2:N:0:1', 'c 2:N:0:1', 'd 2:N:0:1']))
    records = [sam_record('x'), sam_record('a'), sam_record('a', 16, 'ref2'),
               sam_record('y'), sam_record('c', 16), sam_record('z')]
    sam = write_text(tmp_path / ('full.sam' + suffix), '@HD\tVN:1.0\tSO:unsorted\n' + PG
                     + ''.join(records))

    counts = read_alignments.select(str(sam), str(forward), str(reverse), str(tmp_path / 'out'))

    assert counts == (4, 2, 3)
    assert fastq_ids(tmp_path / 'out_fwd.fq.gz') == ['a', 'c']
    assert fastq_ids(tmp_path / 'out_rev.fq.gz') == ['a', 'c']
    assert read_text(tmp_path / 'out.sam.gz').startswith('@HD\tVN:1.0\tSO:unsorted\n' + PG)
    assert sam_body(tmp_path / 'out.sam.gz') == [records[1], records[2], records[4]]


def test_a_pair_with_only_its_second_mate_aligned_is_selected(tmp_path):
    """Legacy /1 /2 names: each mate has its own QNAME; the whole pair is kept."""
    forward = write_text(tmp_path / 'r1.fq', fastq(['a/1', 'b/1', 'c/1']))
    reverse = write_text(tmp_path / 'r2.fq', fastq(['a/2', 'b/2', 'c/2']))
    records = [sam_record('a/1'), sam_record('a/2', 16), sam_record('b/2'), sam_record('q/1')]
    sam = write_text(tmp_path / 'full.sam', PG + ''.join(records))

    counts = read_alignments.select(str(sam), str(forward), str(reverse), str(tmp_path / 'out'))

    assert counts == (3, 2, 3)
    assert fastq_ids(tmp_path / 'out_fwd.fq.gz') == ['a/1', 'b/1']
    assert fastq_ids(tmp_path / 'out_rev.fq.gz') == ['a/2', 'b/2']
    assert sam_body(tmp_path / 'out.sam.gz') == records[:3]


def test_mate_files_of_different_lengths_are_rejected(tmp_path):
    forward = write_text(tmp_path / 'r1.fq', fastq(['a/1', 'b/1']))
    reverse = write_text(tmp_path / 'r2.fq', fastq(['a/2']))
    sam = write_text(tmp_path / 'full.sam', PG + sam_record('a/1'))
    with pytest.raises(read_alignments.ReadAlignmentsError, match='more reads than'):
        read_alignments.select(str(sam), str(forward), str(reverse), str(tmp_path / 'out'))


# --- header check ---------------------------------------------------------------------

def test_header_check_accepts_the_same_database_and_options(tmp_path):
    sam = write_text(tmp_path / 'full.sam.gz', PG + sam_record('a'))
    read_alignments.check_header(str(sam), 'mini16.clustered.fasta', 10, 1e-5)


@pytest.mark.parametrize('ref_db_filename, best, evalue, problem', [
    ('other.clustered.fasta', 10, 1e-5, 'reference mini16.clustered.fasta instead of other'),
    ('mini16.clustered.fasta', 5, 1e-5, '--num_alignments 10 instead of 5'),
    ('mini16.clustered.fasta', 10, 1e-3, '-e 1.00e-05 instead of 1.00e-03'),
])
def test_header_check_rejects_another_mapping(tmp_path, ref_db_filename, best, evalue, problem):
    sam = write_text(tmp_path / 'full.sam', PG + sam_record('a'))
    with pytest.raises(read_alignments.ReadAlignmentsError, match=problem):
        read_alignments.check_header(str(sam), ref_db_filename, best, evalue)


def test_header_check_needs_the_sortmerna_command(tmp_path):
    sam = write_text(tmp_path / 'full.sam', '@HD\tVN:1.0\n' + sam_record('a'))
    with pytest.raises(read_alignments.ReadAlignmentsError, match='no SortMeRNA @PG'):
        read_alignments.check_header(str(sam), 'mini16.clustered.fasta', 10, 1e-5)


# --- MATAM runs -------------------------------------------------------------------------

def _missing_tools():
    path = tool_env()['PATH']
    return [tool for tool in ('sortmerna', 'sga', 'samtools', 'vsearch', 'seqkit')
            if shutil.which(tool, path=path) is None]


needs_tools = pytest.mark.skipif(bool(_missing_tools()),
                                 reason='missing tools: {}'.format(', '.join(_missing_tools())))

# Every 4th example pair (~1,450 pairs), as in test_functional_pe.py.
EVERY = 4


def matam(forward, reverse, db, out_dir, *extra_args):
    completed = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / 'matam_assembly.py'),
         '--forward', str(forward), '--reverse', str(reverse), '-d', str(db), '-o', str(out_dir),
         '--cpu', str(min(8, multiprocessing.cpu_count())), '--max_memory', '2000', '-v',
         '--seed', '42', '--keep_tmp', *extra_args],
        env=tool_env(), capture_output=True, text=True)
    return completed.returncode, completed.stdout + completed.stderr


def workdir_file(out_dir, pattern):
    matches = sorted((out_dir / 'workdir').glob(pattern))
    assert len(matches) == 1, matches
    return matches[0]


def scaffolds(out_dir):
    lines = (out_dir / 'final_assembly.fa').read_text().split()
    return sorted(line for line in lines if not line.startswith('>'))


@pytest.fixture(scope='module')
def mapped_set(tmp_path_factory):
    """A read set mapped once with --filter_only, as MarkerMAG-nf's MATAM_FILTER does."""
    base = tmp_path_factory.mktemp('read_alignments')
    db = build_mini_db(base / 'db')
    forward, reverse = write_pe_reads(base / 'reads', 'casava', every=EVERY, gz=True)
    code, log = matam(forward, reverse, db, base / 'filter', '--filter_only')
    assert code == 0, log[-4000:]
    sam = base / 'filter' / 'marker_alignments.sam.gz'
    assert sam.is_symlink() and sam.resolve().is_file()
    return base, db, forward, reverse, sam


@needs_tools
def test_the_same_reads_give_the_same_assembly_as_mapping_them(mapped_set):
    base, db, forward, reverse, sam = mapped_set
    code, log = matam(forward, reverse, db, base / 'mapped', '--stop_after', 'scaffolding')
    assert code == 0, log[-4000:]
    code, log = matam(forward, reverse, db, base / 'reused', '--stop_after', 'scaffolding',
                      '--read_alignments', str(sam))
    assert code == 0, log[-4000:]
    assert 'Using the alignments in {} (--read_alignments)'.format(sam) in log

    filtered = '*.scr_filt_geo_90pct.sam'
    assert (workdir_file(base / 'reused', filtered).read_bytes()
            == workdir_file(base / 'mapped', filtered).read_bytes())
    assert scaffolds(base / 'reused') and scaffolds(base / 'reused') == scaffolds(base / 'mapped')


@needs_tools
def test_a_subsample_takes_its_reads_alignments_from_the_full_mapping(mapped_set):
    base, db, _, _, sam = mapped_set
    # Every 2nd pair of the mapped set.
    sub_forward, sub_reverse = write_pe_reads(base / 'reads_sub', 'casava', gz=True, pairs=(
        pair for i, pair in enumerate(iter_example_pairs(EVERY)) if i % 2 == 0))
    out_dir = base / 'sub'
    code, log = matam(sub_forward, sub_reverse, db, out_dir, '--read_alignments', str(sam),
                      '--stop_after', 'alignments_filtering')
    assert code == 0, log[-4000:]

    sub_ids = set(fastq_ids(sub_forward))
    selected = fastq_ids(out_dir / 'marker_reads_R1.fq.gz')
    assert selected and set(selected) <= sub_ids
    assert fastq_ids(out_dir / 'marker_reads_R2.fq.gz') == selected
    # The subset's records of the full mapping, in the same order.
    assert sam_body(out_dir / 'marker_alignments.sam.gz') == [
        line for line in sam_body(sam) if line.split('\t', 1)[0] in sub_ids]
    # Every read pair of the subset that has an alignment is selected.
    aligned = {line.split('\t', 1)[0] for line in sam_body(sam)}
    assert set(selected) == sub_ids & aligned


@needs_tools
def test_a_run_from_read_alignments_resumes(mapped_set):
    base, db, forward, reverse, sam = mapped_set
    out_dir = base / 'resume'
    code, log = matam(forward, reverse, db, out_dir, '--read_alignments', str(sam),
                      '--stop_after', 'alignments_filtering')
    assert code == 0, log[-4000:]

    code, log = matam(forward, reverse, db, out_dir, '--read_alignments', str(sam),
                      '--resume_from', 'alignments_filtering', '--stop_after', 'overlap_graph_building')
    assert code == 0, log[-4000:]
    assert '=== Reads mapping' not in log

    # Resuming without the option: the mapping in the run directory is not the
    # one this command would make, so it is not used.
    code, log = matam(forward, reverse, db, out_dir,
                      '--resume_from', 'alignments_filtering', '--stop_after', 'alignments_filtering')
    assert code != 0
    assert 'Resume from the read mapping step' in log


@needs_tools
def test_alignments_made_with_other_mapping_options_are_rejected(mapped_set):
    base, db, forward, reverse, sam = mapped_set
    code, log = matam(forward, reverse, db, base / 'other_best', '--read_alignments', str(sam),
                      '--best', '5', '--stop_after', 'alignments_filtering')
    assert code != 0
    assert "was not made with this run's read mapping: --num_alignments 10 instead of 5" in log


@needs_tools
def test_the_clustered_database_index_is_not_needed(mapped_set, tmp_path):
    base, db, forward, reverse, sam = mapped_set
    db_copy = tmp_path / 'db'
    shutil.copytree(str(db.parent), str(db_copy))
    shutil.rmtree(str(db_copy / 'mini16.clustered.smr'))
    code, log = matam(forward, reverse, db_copy / 'mini16', tmp_path / 'out',
                      '--read_alignments', str(sam), '--stop_after', 'alignments_filtering')
    assert code == 0, log[-4000:]
