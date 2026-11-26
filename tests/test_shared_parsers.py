"""The shared FASTA, table and CIGAR parsers (fasta_utils, tab_utils,
sam_utils), which replaced 19 local copies with diverging edge cases."""
import ast
import io
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).parent.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS_DIR))
from fasta_utils import format_seq, read_fasta_file_handle  # noqa: E402
from sam_utils import parse_cigar  # noqa: E402
from tab_utils import group_rows, read_tab_file_handle_sorted  # noqa: E402
sys.path.pop(0)


def fasta(text):
    return list(read_fasta_file_handle(io.StringIO(text)))


# --- read_fasta_file_handle -------------------------------------------------------

def test_records_with_multiline_sequences():
    assert fasta('>a desc\nACGT\nAC\n>b\nGG\n') == [('a desc', 'ACGTAC'), ('b', 'GG')]


@pytest.mark.parametrize('text', ['', '\n', '\n\n  \n'])
def test_empty_input_yields_nothing(text):
    """Most of the old copies yielded ('', '') here, which crashed every
    caller that does header.split()[0]."""
    assert fasta(text) == []


def test_blank_lines_whitespace_and_crlf_are_ignored():
    assert fasta('\n>a  \r\n AC GT \r\n\r\n  GG\n\n>b\nT') == [('a', 'AC GTGG'), ('b', 'T')]


def test_header_without_sequence_is_kept():
    assert fasta('>a\n>b\nAC\n>c\n') == [('a', ''), ('b', 'AC'), ('c', '')]


def test_empty_header_is_kept():
    assert fasta('>\nAC\n') == [('', 'AC')]


@pytest.mark.parametrize('text', ['ACGT\n>a\nGG\n', 'ACGT\n'])
def test_sequence_before_the_first_header_is_an_error(text):
    """It used to be silently prepended to the first record (or returned
    with an empty header when the file had no header at all)."""
    with pytest.raises(ValueError, match='before the first header'):
        fasta(text)


def test_the_handle_is_closed():
    handle = io.StringIO('>a\nAC\n')
    list(read_fasta_file_handle(handle))
    assert handle.closed


# --- format_seq -------------------------------------------------------------------

def test_format_seq_line_length_is_configurable():
    assert format_seq('ACGTACG', 3) == 'ACG\nTAC\nG'


# --- read_tab_file_handle_sorted --------------------------------------------------

def tab(text, index=0):
    return list(read_tab_file_handle_sorted(io.StringIO(text), index))


def test_rows_are_grouped_by_runs_of_the_column():
    assert tab('q1 r1\nq1 r2\nq2 r1\nq1 r3\n') == [
        [['q1', 'r1'], ['q1', 'r2']], [['q2', 'r1']], [['q1', 'r3']]]


def test_grouping_column_is_selectable():
    assert tab('a x 1\nb x 2\nc y 3\n', 1) == [
        [['a', 'x', '1'], ['b', 'x', '2']], [['c', 'y', '3']]]


@pytest.mark.parametrize('text', ['', '\n \n'])
def test_empty_table_yields_nothing(text):
    """The old copies yielded one empty group, and every caller then failed
    on tab_list[0]."""
    assert tab(text) == []


def test_tabs_spaces_and_blank_lines():
    assert tab('q1\tr1\t5\n\n q1  r2 6 \nq2\tr1\t7\n') == [
        [['q1', 'r1', '5'], ['q1', 'r2', '6']], [['q2', 'r1', '7']]]


def test_any_iterable_of_lines_is_accepted():
    lines = (l for l in ['@HD\tVN:1.0\n', 'r1 0 ref1\n', 'r2 0 ref1\n'] if l[0] != '@')
    assert list(read_tab_file_handle_sorted(lines, 2)) == [[['r1', '0', 'ref1'], ['r2', '0', 'ref1']]]


def test_the_table_handle_is_closed():
    handle = io.StringIO('a 1\n')
    list(read_tab_file_handle_sorted(handle))
    assert handle.closed


def test_group_rows_groups_split_rows_by_column():
    rows = [['r1', '0', 'ref', '5'], ['r2', '0', 'ref', '5'], ['r3', '0', 'ref', '9']]
    assert list(group_rows(rows, 3)) == [rows[:2], rows[2:]]
    assert list(group_rows([], 3)) == []


# --- parse_cigar ------------------------------------------------------------------

@pytest.mark.parametrize('cigar, expected', [
    ('5S10M2I3M1D4M', [('S', 5), ('M', 10), ('I', 2), ('M', 3), ('D', 1), ('M', 4)]),
    ('5M5M2D1D3M', [('M', 10), ('D', 3), ('M', 3)]),
    ('M', [('M', 1)]),
    ('100=1X2N', [('=', 100), ('X', 1), ('N', 2)]),
    ('', []),
])
def test_parse_cigar(cigar, expected):
    assert parse_cigar(cigar) == expected


# --- no copies left ------------------------------------------------------------------

SHARED = {'read_fasta_file_handle', 'format_seq', 'parse_cigar', 'read_tab_file_handle_sorted'}
OWN = {'fasta_utils.py', 'tab_utils.py', 'sam_utils.py'}
# A different function under the same name: groups by factor, then by read.
DIFFERENT = {('compute_lca_from_tab.py', 'read_tab_file_handle_sorted')}


def test_no_script_defines_its_own_copy():
    copies = []
    for path in sorted(SCRIPTS_DIR.glob('*.py')):
        if path.name in OWN:
            continue
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.FunctionDef) and node.name in SHARED \
                    and (path.name, node.name) not in DIFFERENT:
                copies.append('{}:{}'.format(path.name, node.lineno))
    assert copies == []


@pytest.mark.parametrize('script', ['fasta_length_filter.py', 'replace_Ns_by_As.py',
                                    'remove_redundant_sequences.py', 'fasta_clean_name.py'])
def test_pipeline_fasta_scripts_accept_an_empty_fasta(tmp_path, script):
    """These wrote a '>' record with no name (or crashed) on an empty input."""
    empty = tmp_path / 'empty.fa'
    empty.write_text('')
    out = tmp_path / 'out.fa'
    args = {'fasta_length_filter.py': ['-i', str(empty), '-o', str(out), '-m', '0'],
            'replace_Ns_by_As.py': ['-i', str(empty), '-o', str(out)],
            'remove_redundant_sequences.py': ['-i', str(empty), '-o', str(out)],
            'fasta_clean_name.py': ['-i', str(empty), '-o', str(out)]}[script]
    subprocess.run([sys.executable, str(SCRIPTS_DIR / script)] + args, check=True,
                   capture_output=True, text=True)
    assert out.read_text() == ''


@pytest.mark.parametrize('script', ['get_best_matches_from_blast.py', 'generate_scaffolding_blast.py',
                                    'compute_contigs_compatibility.py'])
def test_pipeline_table_scripts_accept_an_empty_table(script):
    """These failed on tab_list[0] when the old reader yielded one empty group
    (e.g. no contig with a hit in the scaffolding step)."""
    result = subprocess.run([sys.executable, str(SCRIPTS_DIR / script)], input='',
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
