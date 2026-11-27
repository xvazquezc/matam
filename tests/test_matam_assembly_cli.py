"""Characterization of matam_assembly.parse_arguments, the small matam_assembly
helpers, and the assembler factory / SGA command line."""
import gzip
import importlib
import os
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).parent.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS_DIR))
matam_assembly = importlib.import_module('matam_assembly')
sortmerna_step = importlib.import_module('sortmerna_step')
assembler_factory = importlib.import_module('assembler_factory')
fasta_utils = importlib.import_module('fasta_utils')
sys.path.pop(0)

REQUIRED = ['--forward', 'r1.fq', '--reverse', 'r2.fq', '-d', 'db/silva']


def parse(*extra):
    return matam_assembly.parse_arguments(REQUIRED + list(extra))


# --- parse_arguments ---------------------------------------------------------

def test_paths_are_made_absolute_and_out_dir_defaults(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    args = parse()
    assert args.forward == str(tmp_path / 'r1.fq')
    assert args.reverse == str(tmp_path / 'r2.fq')
    assert args.ref_db == str(tmp_path / 'db' / 'silva')
    assert args.out_dir == str(tmp_path / 'matam_assembly')


def test_workdir_defaults_to_out_dir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert parse('-o', 'out').workdir == str(tmp_path / 'out')
    assert parse('-o', 'out', '--workdir', 'wd').workdir == str(tmp_path / 'wd')


def test_debug_turns_on_verbose_and_keep_tmp():
    args = parse('--debug')
    assert args.verbose and args.keep_tmp
    args = parse()
    assert not args.verbose and not args.keep_tmp


def test_coverage_threshold_applies_only_to_auto_read_correction():
    assert parse().contig_coverage_threshold == 20
    assert parse('--read_correction', 'auto', '--contig_coverage_threshold', '50').contig_coverage_threshold == 50
    assert parse('--read_correction', 'no').contig_coverage_threshold is None
    assert parse('--read_correction', 'yes', '--contig_coverage_threshold', '50').contig_coverage_threshold is None


@pytest.mark.parametrize('option, value, message', [
    ('--score_threshold', '-0.01', 'score threshold not in range [0,1]'),
    ('--score_threshold', '1.01', 'score threshold not in range [0,1]'),
    ('--min_identity', '1.5', 'min identity not in range [0,1]'),
    ('--quorum', '0.5', 'quorum not in range [0.51,1]'),
    ('--quorum', '0', 'quorum not in range [0.51,1]'),
    ('--quorum', '1.01', 'quorum not in range [0.51,1]'),
    ('--rdp_cutoff', '-1', 'cutoff not in range [0,1]'),
])
def test_out_of_range_values_are_rejected(capsys, option, value, message):
    with pytest.raises(SystemExit) as error:
        parse(option, value)
    assert error.value.code == 2
    assert 'error: ' + message in capsys.readouterr().err


@pytest.mark.parametrize('script, argv, message', [
    ('matam_assembly.py', ['--forward', 'r1.fq', '--reverse', 'r2.fq', '-d', 'db', '--quorum', '0.5'],
     'error: quorum not in range [0.51,1]'),
    ('matam_assembly.py', ['--forward', 'r1.fq', '-d', 'db'],
     'error: the following arguments are required: --reverse'),
    ('matam_db_preprocessing.py', ['-i', 'refs.fasta', '--clustering_id_threshold', '1.5'],
     'error: clustering id threshold not in range [0,1]'),
    ('index_ref_db.py', [], 'error: the following arguments are required'),
])
def test_argument_errors_are_short_and_end_with_the_message(tmp_path, script, argv, message):
    """Usage line and the error as the last line, exit code 2: no traceback and
    no full help scrolling the error out of sight (--help still gives it)."""
    result = subprocess.run([sys.executable, str(SCRIPTS_DIR / script), *argv],
                            cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 2
    assert result.stdout == ''
    lines = result.stderr.strip().splitlines()
    assert lines[0].startswith('usage: ')
    assert message in lines[-1]
    assert 'Traceback' not in result.stderr
    assert len(lines) < 25


@pytest.mark.parametrize('option, value', [
    ('--score_threshold', '0'), ('--score_threshold', '1'),
    ('--min_identity', '0'), ('--min_identity', '1'),
    ('--quorum', '0.51'), ('--quorum', '1'),
    ('--rdp_cutoff', '0'), ('--rdp_cutoff', '1'),
])
def test_range_limits_are_accepted(option, value):
    assert getattr(parse(option, value), option[2:]) == float(value)


def test_any_shared_ref_is_off_by_default():
    assert parse().any_shared_ref is False
    assert parse('--any_shared_ref').any_shared_ref is True


def test_quorum_default_is_the_lowest_accepted_value():
    assert parse().quorum == 0.51


def test_unknown_contig_coverage_threshold_is_rejected(capsys):
    with pytest.raises(SystemExit):
        parse('--contig_coverage_threshold', '30')


# --- helpers -------------------------------------------------------------------

def test_fasta_stats_of_an_empty_file_are_minus_one(tmp_path):
    fasta = tmp_path / 'empty.fa'
    fasta.write_text('')
    stats = matam_assembly.compute_fasta_stats(str(fasta))
    assert (stats.seq_num, stats.total_nt) == (0, 0)
    assert stats.get_avg_length() == stats.get_min_length() == stats.get_max_length() == -1


def test_fasta_stats(tmp_path):
    fasta = tmp_path / 'contigs.fa'
    fasta.write_text('>1\nACGT\nAC\n>2\nACGTACGT\n')
    stats = matam_assembly.compute_fasta_stats(str(fasta))
    assert (stats.seq_num, stats.total_nt) == (2, 14)
    assert (stats.get_min_length(), stats.get_avg_length(), stats.get_max_length()) == (6, 7, 8)


@pytest.mark.parametrize('name', ['r.fq', 'r.fq.gz'])
def test_fastq_line_count(tmp_path, name):
    text = '@a\nACGT\n+\nIIII\n@b\nAC\n+\nII\n'
    path = tmp_path / name
    if name.endswith('.gz'):
        with gzip.open(path, 'wt') as fh:
            fh.write(text)
    else:
        path.write_text(text)
    assert matam_assembly.fastq_line_count(str(path)) == 8


@pytest.mark.parametrize('seq, expected', [
    ('', ''),
    ('ACGT', 'ACGT'),
    ('A' * 80, 'A' * 80),
    ('A' * 81, 'A' * 80 + '\nA'),
])
def test_format_seq_wraps_at_80(seq, expected):
    assert fasta_utils.format_seq(seq) == expected


def test_both_plain_and_gzipped_sortmerna_outputs_are_ambiguous(tmp_path):
    base = str(tmp_path / 'aligned')
    for name in (base + '.sam', base + '.sam.gz'):
        Path(name).touch()
    with pytest.raises(RuntimeError, match='ambiguous outputs'):
        sortmerna_step.find_sortmerna_output(base, ('.sam',))


def test_fq_and_fastq_outputs_are_ambiguous(tmp_path):
    base = str(tmp_path / 'aligned_fwd')
    for name in (base + '.fq', base + '.fastq.gz'):
        Path(name).touch()
    with pytest.raises(RuntimeError, match='ambiguous outputs'):
        sortmerna_step.find_sortmerna_output(base, ('.fq', '.fastq'))


def test_missing_sortmerna_output_lists_the_names_looked_for(tmp_path):
    base = str(tmp_path / 'aligned')
    with pytest.raises(FileNotFoundError, match=r'aligned\.sam\.gz'):
        sortmerna_step.find_sortmerna_output(base, ('.sam',))


# --- assembler factory and SGA command line --------------------------------------

@pytest.fixture
def fake_binaries(monkeypatch):
    monkeypatch.setattr(assembler_factory.binary_utils.Binary, 'assert_which',
                        staticmethod(lambda name: '/bin/' + name))


def test_unknown_assembler_is_a_key_error(fake_binaries):
    with pytest.raises(KeyError, match="Valid keys: \\['SGA'\\]"):
        assembler_factory.AssemblerFactory().get('MEGAHIT')


def test_sga_is_the_only_assembler(fake_binaries):
    assert isinstance(assembler_factory.AssemblerFactory().get('SGA'), assembler_factory.SGA)


@pytest.mark.parametrize('read_correction, corrected', [('no', False), ('auto', False), ('yes', True)])
def test_sga_command_disables_correction_except_in_yes_mode(fake_binaries, read_correction, corrected):
    sga = assembler_factory.AssemblerFactory().get('SGA')
    sga.build_command_line('c_R1.fq', 'c_R2.fq', 'wd', read_correction, 3)
    assert ('--no_correction' in sga.cmd_line) is not corrected
    assert '/bin/sga_assemble.py --forward c_R1.fq --reverse c_R2.fq' in sga.cmd_line
    assert '-o wd/assembly.fasta --sga_bin /bin/sga' in sga.cmd_line
    assert ' --cpu 3 --tmp_dir wd/tmp >> wd/assembly.log 2>&1' in sga.cmd_line
    assert sga.fasta_file == os.path.join('wd', 'assembly.fasta')


def test_run_before_build_stops(fake_binaries):
    with pytest.raises(SystemExit):
        assembler_factory.AssemblerFactory().get('SGA').run()
