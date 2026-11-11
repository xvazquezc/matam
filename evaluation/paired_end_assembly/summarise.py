#!/usr/bin/env python3
"""Summarise MATAM runs made by run_compare.sh.

Usage: summarise.py OUT_DIR [PREFIX ...]   (default prefixes: up fork)

For each run OUT_DIR/<prefix>_<seed> it reports graph size, read usage,
contigs and final assembly, and evaluates the contigs and the final assembly
against the true references (OUT_DIR/refs.fasta) with this repository's
evaluate_assembly.py (needs exonerate on PATH). Prints one column per run,
then the mean per prefix. Evaluation files go to OUT_DIR/eval/.
"""
import glob
import os
import re
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def fasta_lengths(path):
    lengths, current = [], None
    for line in open(path):
        if line.startswith('>'):
            if current is not None:
                lengths.append(current)
            current = 0
        else:
            current += len(line.strip())
    if current is not None:
        lengths.append(current)
    return lengths


def n50(lengths):
    total, acc = sum(lengths), 0
    for length in sorted(lengths, reverse=True):
        acc += length
        if acc * 2 >= total:
            return length
    return 0


def evaluate(out_dir, fasta, tag):
    eval_dir = os.path.join(out_dir, 'eval', tag)
    shutil.rmtree(eval_dir, ignore_errors=True)
    os.makedirs(eval_dir)
    local = os.path.join(eval_dir, os.path.basename(fasta))
    shutil.copy(fasta, local)
    env = dict(os.environ, PATH=os.path.join(REPO, 'scripts') + os.pathsep + os.environ['PATH'])
    subprocess.run([sys.executable, os.path.join(REPO, 'scripts', 'evaluate_assembly.py'),
                    '-r', os.path.join(out_dir, 'refs.fasta'), '-i', local],
                   check=True, env=env, capture_output=True)
    text = open(glob.glob(os.path.join(eval_dir, '*.assembly.stats'))[0]).read()
    get = lambda pattern: re.search(pattern, text).group(1)
    per_ref = re.findall(r'\n\t\t(\S+)\t([\d.]+)%\t', text.split('Per-reference Stats:')[1])
    return {
        'error': float(get(r'Assembly error rate\s+=\s+([\d.]+)%')),
        'ref_cov': float(get(r'Total ref coverage\s+=\s+([\d.]+)%')),
        'refs_90': sum(1 for _, cov in per_ref if float(cov) >= 90.0),
        'refs_hit': sum(1 for _, cov in per_ref if float(cov) > 0.0),
        'per_ref': {ref: float(cov) for ref, cov in per_ref},
    }


def run_stats(out_dir, run):
    log = open(os.path.join(out_dir, run + '.log')).read()
    stats = {}
    m = re.search(r'Overlap graph stats: (\d+) nodes, (\d+) edges', log)
    stats['nodes'], stats['edges'] = int(m.group(1)), int(m.group(2))
    stats['components'] = int(re.search(r'Compressed graph: (\d+) components', log).group(1))

    components_csv = glob.glob(os.path.join(out_dir, run, 'workdir', '*.components.csv'))[0]
    with open(components_csv) as handle:
        col = handle.readline().rstrip('\n').split(';').index('component')
        stats['mates_in_components'] = sum(1 for line in handle
                                           if int(line.rstrip('\n').split(';')[col]) >= 0)

    # Reads given to SGA, by the pair rule of the version that ran
    routed = re.search(r'Read pairs by component: (\d+) within one component, (\d+) spanning two '
                       r'\(given to both\), (\d+) with one mate', log)
    same_component_only = re.search(r'Kept (\d+) complete same-component read pairs', log)
    if routed:                       # current fork: pairs go to every component of a mate
        within, spanning, one_mate = map(int, routed.groups())
        stats['reads_to_sga'] = 2 * (within + 2 * spanning + one_mate)
    elif same_component_only:        # fork 82c990d..264cf57: same-component pairs only
        stats['reads_to_sga'] = 2 * int(same_component_only.group(1))
    else:                            # upstream: every read of a component, unpaired
        stats['reads_to_sga'] = stats['mates_in_components']

    m = re.search(r'Number of non empty components: (\d+)/(\d+)', log)
    if m:
        stats['assembled_components'], stats['components_with_reads'] = map(int, m.groups())
    else:  # MATAM logs this only when some component assembled nothing
        stats['assembled_components'] = stats['components_with_reads'] = stats['components']

    contigs_fasta = os.path.join(out_dir, run, 'workdir', 'contigs.fasta')
    contigs = fasta_lengths(contigs_fasta)
    stats['contigs'], stats['contigs_bp'], stats['contigs_n50'] = len(contigs), sum(contigs), n50(contigs)
    final_fasta = os.path.join(out_dir, run, 'final_assembly.fa')
    scaffolds = fasta_lengths(final_fasta)
    stats['scaffolds'], stats['scaffolds_bp'] = len(scaffolds), sum(scaffolds)
    stats['scaffolds_n50'] = n50(scaffolds)
    stats['scaffolds_max'] = max(scaffolds) if scaffolds else 0
    for key, value in evaluate(out_dir, contigs_fasta, run + '_contigs').items():
        stats['contigs_' + key] = value
    for key, value in evaluate(out_dir, final_fasta, run + '_final').items():
        stats['final_' + key] = value
    stats['wall_s'] = float(open(os.path.join(out_dir, run + '.time')).read().split()[0])
    return stats


ROWS = [
    ('nodes', 'Graph nodes (mates)'), ('edges', 'Graph edges'), ('components', 'Components'),
    ('mates_in_components', 'Mates in components'), ('reads_to_sga', 'Reads given to SGA'),
    ('components_with_reads', 'Components given to SGA'),
    ('assembled_components', 'Components with contigs'),
    ('contigs', 'Contigs'), ('contigs_bp', 'Contigs total bp'), ('contigs_n50', 'Contig N50'),
    ('contigs_error', 'Contig error rate %'), ('contigs_ref_cov', 'Ref coverage by contigs %'),
    ('scaffolds', 'Final scaffolds'), ('scaffolds_bp', 'Final total bp'),
    ('scaffolds_n50', 'Final N50'), ('scaffolds_max', 'Longest scaffold'),
    ('final_error', 'Final error rate %'), ('final_ref_cov', 'Ref coverage by final %'),
    ('final_refs_90', 'Refs >=90% covered'), ('final_refs_hit', 'Refs hit'),
    ('wall_s', 'Wall time (s)'),
]


def main(argv):
    out_dir = os.path.abspath(argv[1])
    prefixes = argv[2:] or ['up', 'fork']
    runs = sorted((os.path.basename(path)[:-len('.log')]
                   for prefix in prefixes
                   for path in glob.glob(os.path.join(out_dir, prefix + '_*.log'))),
                  key=lambda run: (prefixes.index(run.rsplit('_', 1)[0]), run))
    results = {run: run_stats(out_dir, run) for run in runs}

    print('metric\t' + '\t'.join(runs))
    for key, label in ROWS:
        print(label + '\t' + '\t'.join(str(results[run][key]) for run in runs))

    by_prefix = defaultdict(list)
    for run in runs:
        by_prefix[run.rsplit('_', 1)[0]].append(run)
    print('\nmetric\t' + '\t'.join('{} (mean of {})'.format(p, len(by_prefix[p])) for p in prefixes))
    for key, label in ROWS:
        means = [statistics.mean(results[run][key] for run in by_prefix[p]) for p in prefixes]
        print(label + '\t' + '\t'.join('{:.2f}'.format(mean) for mean in means))

    print('\nPer-reference coverage by the final assembly (%)')
    for run in runs:
        per_ref = results[run]['final_per_ref']
        print(run + '\t' + ' '.join('{}={:.0f}'.format(ref.split('_')[0], cov)
                                    for ref, cov in sorted(per_ref.items())))


if __name__ == '__main__':
    main(sys.argv)
