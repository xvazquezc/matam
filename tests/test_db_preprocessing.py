"""matam_db_preprocessing: reference headers with tabs (#9)."""

import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
from pe_fixtures import EXAMPLE_REFS, SCRIPTS_DIR, tool_env, write_pe_reads  # noqa: E402

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
