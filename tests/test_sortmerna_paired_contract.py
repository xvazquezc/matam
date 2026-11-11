from pathlib import Path
import importlib
import sys


SCRIPTS_DIR = Path(__file__).parent.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS_DIR))
matam_assembly = importlib.import_module('matam_assembly')
compute_abundance = importlib.import_module('compute_abundance')
sys.path.pop(0)


def test_paired_mapping_keeps_complete_pairs_in_separate_fastx_files():
    """The MarkerMAG-nf filter wrapper relies on this MATAM/SortMeRNA contract."""
    source = (Path(__file__).parent.parent / 'scripts' / 'matam_assembly.py').read_text()
    mapping_stage = source.split("# Everything that determines the mapping output", 1)[1].split(
        "logger.info('=== Reads mapping against ref db ===')", 1)[0]

    assert "--reads ' + args.forward + ' --reads ' + args.reverse" in mapping_stage
    assert "cmd_line += ' --paired_in --out2'" in mapping_stage


def test_selected_pairs_are_published_under_stable_names(tmp_path):
    """MarkerMAG-nf reads <out_dir>/marker_reads_R1/R2 after --filter_only."""
    assert matam_assembly.marker_reads_filepath(
        str(tmp_path), 'R1', '/x/aligned_fwd.fq.gz') == str(tmp_path / 'marker_reads_R1.fq.gz')
    assert matam_assembly.marker_reads_filepath(
        str(tmp_path), 'R2', '/x/aligned_rev.fq') == str(tmp_path / 'marker_reads_R2.fq')


def test_public_cli_requires_separate_paired_input_files():
    source = (Path(__file__).parent.parent / 'scripts' / 'matam_assembly.py').read_text()
    parser_source = source.split('def parse_arguments():', 1)[1].split(
        'def print_intro(args):', 1)[0]

    assert "group_main.add_argument('-i', '--input_fastq'" not in parser_source
    assert "'--forward'" in parser_source
    assert "'--reverse'" in parser_source
    assert parser_source.count('required = True') >= 2


def test_sortmerna_output_discovery_accepts_fq_fastq_and_gzip(tmp_path):
    prefix = tmp_path / 'aligned'

    fastq = Path(str(prefix) + '.fq.gz')
    sam = Path(str(prefix) + '.sam.gz')
    fastq.touch()
    sam.touch()

    assert matam_assembly.find_sortmerna_output(
        str(prefix), ('.fq', '.fastq')) == str(fastq)
    assert matam_assembly.find_sortmerna_output(
        str(prefix), ('.sam',)) == str(sam)
    assert matam_assembly.stream_file_command(str(sam)).startswith('gzip -cd ')

    fastq.unlink()
    fastq = Path(str(prefix) + '.fastq')
    fastq.touch()
    assert matam_assembly.find_sortmerna_output(
        str(prefix), ('.fq', '.fastq')) == str(fastq)


def test_abundance_streams_compressed_blast_reports():
    command = compute_abundance.stream_file_command('/tmp/reads.blast.gz')

    assert command.startswith('gzip -cd ')
