"""Paired-end test fixtures built from the bundled 16-species example.

The bundled reads are ART-simulated pairs, interleaved, with legacy
``/1`` and ``/2`` suffixes.  Real Illumina data names both mates
identically (CASAVA 1.8+: ``@ID 1:N:0:IDX`` / ``@ID 2:N:0:IDX``), and SRA
dumps often have identical IDs with no comment at all.  These helpers
rewrite the bundled pairs into each convention while keeping the same core
ID, so that every variant has the same read order after MATAM sorts by name.
"""

import gzip
import os
import subprocess
import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = ROOT_DIR / 'scripts'
EXAMPLE_DIR = ROOT_DIR / 'examples' / '16sp_simulated_dataset'
EXAMPLE_READS = EXAMPLE_DIR / '16sp.art_HS25_pe_100bp_50x.fq'
EXAMPLE_REFS = EXAMPLE_DIR / '16sp.fasta'
EXAMPLE_TAXO = EXAMPLE_DIR / '16sp.taxo.tab'

NAMING_STYLES = ('casava', 'sra', 'legacy')


def iter_example_pairs(every=1):
    """Yield (core_id, (seq1, qual1), (seq2, qual2)) for every n-th pair."""
    with open(EXAMPLE_READS) as handle:
        index = 0
        while True:
            first = [handle.readline() for _ in range(4)]
            if not first[0]:
                return
            second = [handle.readline() for _ in range(4)]
            if not all(second):
                raise ValueError('Example reads end with an orphan record')
            id1, id2 = first[0][1:].split()[0], second[0][1:].split()[0]
            if not (id1.endswith('/1') and id2.endswith('/2') and id1[:-2] == id2[:-2]):
                raise ValueError('Unexpected example pair: {} {}'.format(id1, id2))
            if index % every == 0:
                yield (id1[:-2],
                       (first[1].strip(), first[3].strip()),
                       (second[1].strip(), second[3].strip()))
            index += 1


def mate_header(core_id, mate, style):
    if style == 'casava':
        return '{} {}:N:0:ACGTACGT'.format(core_id, mate)
    if style == 'sra':
        return core_id
    if style == 'legacy':
        return '{}/{}'.format(core_id, mate)
    raise ValueError('Unknown naming style: {}'.format(style))


def write_pe_reads(out_dir, style, every=1, gz=False, pairs=None):
    """Write R1/R2 FASTQ files in the requested naming style.

    Returns (forward_path, reverse_path).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    suffix = '.fq.gz' if gz else '.fq'
    forward = out_dir / '{}_R1{}'.format(style, suffix)
    reverse = out_dir / '{}_R2{}'.format(style, suffix)
    opener = gzip.open if gz else open
    if pairs is None:
        pairs = iter_example_pairs(every)
    with opener(forward, 'wt') as fwd, opener(reverse, 'wt') as rev:
        for core_id, (seq1, qual1), (seq2, qual2) in pairs:
            fwd.write('@{}\n{}\n+\n{}\n'.format(mate_header(core_id, 1, style), seq1, qual1))
            rev.write('@{}\n{}\n+\n{}\n'.format(mate_header(core_id, 2, style), seq2, qual2))
    return forward, reverse


def write_reference_fasta(path):
    """Write the 16 references with SILVA-style ``>ID taxonomy`` headers."""
    taxonomy = {}
    with open(EXAMPLE_TAXO) as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) == 2:
                taxonomy[fields[0]] = fields[1]
    with open(EXAMPLE_REFS) as handle, open(path, 'w') as out:
        for line in handle:
            if line.startswith('>'):
                ref_id = line[1:].split()[0]
                out.write('>{} {}\n'.format(ref_id, taxonomy[ref_id[:3]]))
            else:
                out.write(line)
    return Path(path)


def tool_env():
    """Environment exposing the MATAM scripts and this interpreter's tools."""
    env = dict(os.environ)
    path = [str(SCRIPTS_DIR), str(Path(sys.executable).parent), env.get('PATH', '')]
    env['PATH'] = os.pathsep.join(path)
    return env


def build_mini_db(out_dir, cpu=2):
    """Build a MATAM database from the 16 example references.

    Returns the database prefix to pass to ``matam_assembly.py -d``.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    refs = write_reference_fasta(out_dir / 'refs.fasta')
    subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / 'matam_db_preprocessing.py'),
         '-i', str(refs), '-d', str(out_dir), '--cpu', str(cpu),
         '--max_memory', '2000', '--out_db_name', 'mini16'],
        check=True, env=tool_env(), capture_output=True, text=True)
    return out_dir / 'mini16'
