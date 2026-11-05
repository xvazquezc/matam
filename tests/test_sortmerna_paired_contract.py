from pathlib import Path


def test_paired_mapping_keeps_complete_pairs_in_fastx_output():
    """The MarkerMAG-nf filter wrapper relies on this MATAM/SortMeRNA contract."""
    source = (Path(__file__).parent.parent / 'scripts' / 'matam_assembly.py').read_text()
    mapping_stage = source.split("logger.info('=== Reads mapping against ref db ===')", 1)[1]

    assert "--reads ' + args.forward + ' --reads ' + args.reverse" in mapping_stage
    assert "cmd_line += ' --paired_in'" in mapping_stage


def test_public_cli_requires_separate_paired_input_files():
    source = (Path(__file__).parent.parent / 'scripts' / 'matam_assembly.py').read_text()
    parser_source = source.split('def parse_arguments():', 1)[1].split(
        'def print_intro(args):', 1)[0]

    assert "group_main.add_argument('-i', '--input_fastq'" not in parser_source
    assert "'--forward'" in parser_source
    assert "'--reverse'" in parser_source
    assert parser_source.count('required = True') >= 2
