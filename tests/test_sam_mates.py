import gzip
import sys
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).parent.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS_DIR))
import sam_mates  # noqa: E402
sys.path.pop(0)

REF = 'ACGTTGCAAGGCTTACCGATTGACCTGAAGTCCATGGATCCGTAGCTAGGCTAACGTTAGC'
R1_SEQ = REF[0:20]
R2_SEQ = sam_mates.reverse_complement(REF[30:50])   # R2 reads the opposite strand


def sam_line(qname, flag, ref, pos, seq, score=40):
    return '\t'.join([qname, str(flag), ref, str(pos), '255', '{}M'.format(len(seq)),
                      '*', '0', '0', seq, '*', 'AS:i:{}'.format(score)]) + '\n'


def records(lines):
    return [line.rstrip('\n').split('\t') for line in lines if not line.startswith('@')]


def tag(sam, pairs):
    return records(sam_mates.tag_mates(sam, pairs))


def test_shared_qname_mates_are_told_apart_by_sequence():
    # SortMeRNA output: both mates share the QNAME; R2 aligns reverse-complemented.
    sam = [
        '@HD\tVN:1.0\tSO:unsorted\n',
        sam_line('pair1', 0, 'ref1', 1, R1_SEQ),
        sam_line('pair1', 16, 'ref1', 31, REF[30:50]),
    ]
    pairs = [(('pair1', R1_SEQ), ('pair1', R2_SEQ))]

    tagged = list(sam_mates.tag_mates(sam, pairs))

    assert tagged[0] == '@HD\tVN:1.0\tSO:unsorted\n'
    assert tagged[1].startswith('@PG\tID:matam_sam_mates')
    out = records(tagged)
    assert [(r[0], int(r[1])) for r in out] == [('pair1', 0x1 | 0x40), ('pair1', 0x1 | 0x10 | 0x80)]
    assert [sam_mates.mate_of(r[1]) for r in out] == [1, 2]


def test_records_are_regrouped_r1_first_and_unaligned_mate_is_flagged():
    sam = [
        sam_line('pair1', 16, 'ref1', 31, REF[30:50]),   # R2 listed first
        sam_line('pair1', 0, 'ref1', 1, R1_SEQ),
        sam_line('pair1', 0, 'ref2', 7, R1_SEQ),
        sam_line('pair2', 0, 'ref1', 1, R1_SEQ),         # pair2: only R1 aligned
    ]
    pairs = [(('pair1', R1_SEQ), ('pair1', R2_SEQ)),
             (('pair2', R1_SEQ), ('pair2', 'GGGGGGGGGGGGGGGGGGGG'))]

    out = tag(sam, pairs)

    assert [(r[0], r[2], sam_mates.mate_of(r[1])) for r in out] == [
        ('pair1', 'ref1', 1), ('pair1', 'ref2', 1), ('pair1', 'ref1', 2), ('pair2', 'ref1', 1)]
    assert not int(out[0][1]) & sam_mates.MATE_UNMAPPED
    assert int(out[3][1]) & sam_mates.MATE_UNMAPPED


def test_pairs_without_alignments_are_skipped():
    sam = [sam_line('pair3', 0, 'ref1', 1, R1_SEQ)]
    pairs = [(('pair1', 'AAAA'), ('pair1', 'CCCC')),
             (('pair2', 'AAAA'), ('pair2', 'CCCC')),
             (('pair3', R1_SEQ), ('pair3', R2_SEQ))]

    assert [r[0] for r in tag(sam, pairs)] == ['pair3']


def test_legacy_suffixed_names_are_assigned_by_qname():
    sam = [sam_line('pair1/1', 0, 'ref1', 1, R1_SEQ),
           sam_line('pair1/2', 16, 'ref1', 31, REF[30:50])]
    pairs = [(('pair1/1', R1_SEQ), ('pair1/2', R2_SEQ))]

    out = tag(sam, pairs)

    assert [(r[0], sam_mates.mate_of(r[1])) for r in out] == [('pair1/1', 1), ('pair1/2', 2)]


def test_identical_mates_fall_back_to_record_order():
    sam = [sam_line('pair1', 0, 'ref1', 1, R1_SEQ),
           sam_line('pair1', 0, 'ref2', 9, R1_SEQ),
           sam_line('pair1', 0, 'ref1', 1, R1_SEQ),     # same placement again -> R2
           sam_line('pair1', 0, 'ref2', 9, R1_SEQ)]
    pairs = [(('pair1', R1_SEQ), ('pair1', R1_SEQ))]

    assert [sam_mates.mate_of(r[1]) for r in tag(sam, pairs)] == [1, 1, 2, 2]


@pytest.mark.parametrize('sam, pairs, message', [
    ([sam_line('pair1', 0, 'ref1', 1, 'TTTTTTTTTTTTTTTTTTTT')],
     [(('pair1', R1_SEQ), ('pair1', R2_SEQ))], 'matches neither mate'),
    ([sam_line('pair2', 0, 'ref1', 1, R1_SEQ), sam_line('pair1', 0, 'ref1', 1, R1_SEQ)],
     [(('pair1', R1_SEQ), ('pair1', R2_SEQ)), (('pair2', R1_SEQ), ('pair2', R2_SEQ))],
     'not found in the read pairs'),
    ([sam_line('pair1', 256, 'ref1', 1, R1_SEQ)],
     [(('pair1', R1_SEQ), ('pair1', R2_SEQ))], 'secondary/supplementary'),
])
def test_inconsistent_input_fails_loudly(sam, pairs, message):
    with pytest.raises(sam_mates.MateTaggingError, match=message):
        tag(sam, pairs)


def write_fastq(path, records_, gz=False):
    opener = gzip.open if gz else open
    with opener(path, 'wt') as handle:
        for header, seq in records_:
            handle.write('@{}\n{}\n+\n{}\n'.format(header, seq, 'I' * len(seq)))


def test_cli_reads_gzipped_paired_files(tmp_path):
    sam = tmp_path / 'aligned.sam.gz'
    with gzip.open(sam, 'wt') as handle:
        handle.write(sam_line('A00123:8:FC:1:1101:1000:1', 0, 'ref1', 1, R1_SEQ))
        handle.write(sam_line('A00123:8:FC:1:1101:1000:1', 16, 'ref1', 31, REF[30:50]))
    forward, reverse = tmp_path / 'R1.fq.gz', tmp_path / 'R2.fq.gz'
    write_fastq(forward, [('A00123:8:FC:1:1101:1000:1 1:N:0:ACGT', R1_SEQ)], gz=True)
    write_fastq(reverse, [('A00123:8:FC:1:1101:1000:1 2:N:0:ACGT', R2_SEQ)], gz=True)
    out = tmp_path / 'tagged.sam'

    assert sam_mates.main(['-s', str(sam), '-1', str(forward), '-2', str(reverse), '-o', str(out)]) == 0
    with open(out) as handle:
        assert [sam_mates.mate_of(r[1]) for r in records(handle)] == [1, 2]


def test_cli_reports_pair_files_out_of_step(tmp_path, capsys):
    sam = tmp_path / 'aligned.sam'
    sam.write_text(sam_line('pair2', 0, 'ref1', 1, R1_SEQ))
    forward, reverse = tmp_path / 'R1.fq', tmp_path / 'R2.fq'
    write_fastq(forward, [('pair1 1:N:0:A', R1_SEQ), ('pair2 1:N:0:A', R1_SEQ)])
    write_fastq(reverse, [('pair1 2:N:0:A', R2_SEQ)])

    assert sam_mates.main(['-s', str(sam), '-1', str(forward), '-2', str(reverse),
                           '-o', str(tmp_path / 'o.sam')]) == 1
    assert 'different record counts' in capsys.readouterr().err
