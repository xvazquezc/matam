"""componentsearch: graph compaction and components identification.

componentsearch was rewritten to hold the read graph once, in compressed form
(8 bytes per edge instead of about 50). Without an edge filter (-E 0) it must
give the outputs of upstream ComponentSearch2 (f705924) for the same seed,
byte for byte (glibc's rand).

The edge filter is new. Upstream was given -N instead of -E, counted links per
outside read instead of per meta node, and removed the wrong neighbor on one
side only; on production graphs it removed none of the weight-1 edges. Now -E T
removes the meta edges of at most T read links, from both sides.
"""
import collections
import hashlib
import os
import random
import re
import subprocess

import pytest


ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COMPONENTSEARCH = os.path.join(ROOT_DIR, 'componentsearch', 'componentsearch')

pytestmark = pytest.mark.skipif(
    not os.path.isfile(COMPONENTSEARCH) or not os.access(COMPONENTSEARCH, os.X_OK),
    reason='componentsearch binary is not built'
)

EDGES_HEADER = 'Source;Target;Type;Coemitted;Weight;Pid;MultiRef\n'
OUTPUTS = ('components', 'metaNodes', 'metaEdges')


def write_graph(tmp_path, reads_num, edges, edges_header=EDGES_HEADER, edge_format='{};{};Undirected;1;1;1;0\n'):
    nodes = tmp_path / 'nodes.csv'
    nodes.write_text('Id;Label;Mate;Specie\n'
                     + ''.join('{};read{};0;rea\n'.format(i, i) for i in range(reads_num)))
    edges_csv = tmp_path / 'edges.csv'
    edges_csv.write_text(edges_header + ''.join(edge_format.format(a, b) for a, b in edges))
    return nodes, edges_csv


def run(tmp_path, nodes, edges, *args, check=True):
    base = str(tmp_path / 'cpts')
    proc = subprocess.run([COMPONENTSEARCH, '-n', str(nodes), '-e', str(edges), '-b', base] + list(args),
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if check:
        assert proc.returncode == 0, proc.stderr
    return proc


def output(tmp_path, name):
    return (tmp_path / 'cpts.{}.csv'.format(name)).read_text()


def interval_graph(seed, reads_num=1500, genomes=6, link=0.9):
    """Reads along a few genomes, linked when they overlap by >= 50 bases
    (with probability `link`), plus a few chimeric links between genomes.
    Edges are listed source < target, in order, as ovgraphbuild writes them."""
    rng = random.Random(seed)
    reads = []
    for _ in range(reads_num):
        start = int(rng.random() * 1400)
        reads.append((int(rng.random() * genomes), start, start + 80 + int(rng.random() * 70)))
    edges = set()
    for i in range(reads_num):
        for j in range(i + 1, reads_num):
            a, b = reads[i], reads[j]
            if a[0] == b[0] and min(a[2], b[2]) - max(a[1], b[1]) >= 50 and rng.random() < link:
                edges.add((i, j))
    for _ in range(10):
        a, b = int(rng.random() * reads_num), int(rng.random() * reads_num)
        if a != b:
            edges.add((min(a, b), max(a, b)))
    return sorted(edges)


# A hub (2) with a branch, a finger, a chain and a triangle, and an isolated read.
SMALL_EDGES = [(0, 1), (1, 2), (2, 3), (3, 4), (2, 5), (5, 6), (2, 7),
               (8, 9), (9, 10), (8, 10)]

SMALL_EXPECTED = {
    'stdout': '\n--- Loading ---\n12 nodes loaded \n\n--- Run algorithms ---\n-> BFS...\n-> Contraction...\n9 nodes in the contracted graph\n6 edges in the contracted graph\n-> Absorb fingers...\n7 nodes in the contracted graph\n4 edges in the contracted graph\n\n--- Filtering ---\n-> Filtering nodes...\n7 nodes in the filtered graph\n-> Filtering edges...\n4 edges in the filtered graph\n\n--- Splicing ---\n\n--- Saving results ---\n-> Saving components\n-> Saving meta graph\n\n--- Program ended ---\n\n',
    'components': 'Id;Label;Mate;Specie;metanode;component\n0;read0;0;rea;0;0\n1;read1;0;rea;1;0\n2;read2;0;rea;2;0\n3;read3;0;rea;3;0\n4;read4;0;rea;4;0\n5;read5;0;rea;2;0\n6;read6;0;rea;2;0\n7;read7;0;rea;2;0\n8;read8;0;rea;7;1\n9;read9;0;rea;7;1\n10;read10;0;rea;7;1\n11;read11;0;rea;8;2\n',
    'metaNodes': 'node;size\n0;1\n1;1\n2;4\n3;1\n4;1\n7;3\n8;1\n',
    'metaEdges': 'from;to\n0;1\n1;2\n2;3\n3;4\n',
}


def test_small_graph_matches_upstream(tmp_path):
    nodes, edges = write_graph(tmp_path, 12, SMALL_EDGES)
    proc = run(tmp_path, nodes, edges, '-r', '42', '-N', '0', '-E', '0')
    assert proc.stdout == SMALL_EXPECTED['stdout']
    for name in OUTPUTS:
        assert output(tmp_path, name) == SMALL_EXPECTED[name], name


# sha256 of (components, metaNodes, metaEdges) per (graph seed, -r, -N, -E).
# -N 0 -E 0: from the upstream binary. Other -E 0 entries: upstream with its
# edge filter switched off (equal to upstream where that filter removed
# nothing, e.g. (2, 42, 1, 0) = upstream -N 1 -E 1). -E 1 and 2: the new edge
# filter, from a build checked against weights counted independently.
INTERVAL_EXPECTED = {
    (1, 42, 0, 0): ('0949ec71eed760f5119b389156a0e56308f919f2bc7a35d1f437de985175c205',
                    '640cb92e60b3c4a4f781a3936822a3c01c69d3c340d01956e64411816e7d21b7',
                    'fc8ec118536eb5b494afc9847fe39d16b2adf3b6e7c25078f992e3a83274aab4'),  # lines [1501, 116, 80]
    (1, 42, 1, 0): ('e06d984c42d3073179b207dda0d1cbb176feddbd03e54ac0a628a6d426ac3289',
                    '15bcac9de5108f4e534b9e6ae1d7adea10dc0675c32671283aeaadd28320f6c3',
                    '0e4b4480a7187e21cf223799d1099f1a9a3e0e79a9bf8971374c16f880c5c97a'),  # lines [1501, 104, 95]
    (1, 3, 2, 0): ('b921dbc629da7ee1451056cbf4ff212bb1f70c2d039aaac0b87530286a4a34e0',
                   '2e8019a8f1b876a8f9892cad3e156d0d1b65275b1e4c2b87b211ab0b4ba83e6b',
                   '67eee1f4a27bd46ea01018972ab973edfb464933a2c5c7eda3289e1328be08ce'),  # lines [1501, 101, 89]
    (1, 3, 2, 2): ('b921dbc629da7ee1451056cbf4ff212bb1f70c2d039aaac0b87530286a4a34e0',
                   '2e8019a8f1b876a8f9892cad3e156d0d1b65275b1e4c2b87b211ab0b4ba83e6b',
                   '67eee1f4a27bd46ea01018972ab973edfb464933a2c5c7eda3289e1328be08ce'),  # lines [1501, 101, 89]
    (2, 42, 0, 0): ('c7d41babc6f4e1a55b1f6bb6e1591eb22ecf33cfd60d06a5b9650b0c25e2dfe4',
                    '82895789704e2ee793622c207849ef984dc665415f0463584de335fb9f41ea25',
                    '0a1d551e2c067c477c9f8145da390df5676c5622ec3ac26f79c6471c2efa638a'),  # lines [1501, 112, 78]
    (2, 42, 1, 0): ('67931526af489c26416e71bcf7ff36c92d3e5afb22122e224865c48bc91014ed',
                    '3e7b33f656a5b039237f3387e60469b05aec98438db277d530b44ea96312a92b',
                    'a0ccc3e81bcba213fb5d6425f531a4484d025c2de559dd57a47acd17aeb65048'),  # lines [1501, 102, 78]
    (2, 42, 1, 1): ('d82413ad1c84f4248d71cd0df5df342a63d5fa0977b3b8b273c1c3f004049e4b',
                    '3e7b33f656a5b039237f3387e60469b05aec98438db277d530b44ea96312a92b',
                    '2660befe6a3e6018eec43e5a8d0a8934f844c3898baaa8975091c532a70e546f'),  # lines [1501, 102, 81]
    (2, 3, 2, 0): ('533f86b2d964f4c01d1ee2a38506e40a32292c37dd7a121edf0c3ae1dc54cc29',
                   '20dee6bf17d62a1eff1d2e39861c8c70893161d4a13c539ac585dec3919b822a',
                   '7f9279fcfa7d37360f7800c4a356568829c1edfba7e576cdc11ed58979fc9b16'),  # lines [1501, 101, 83]
}


@pytest.mark.parametrize('graph_seed,seed,min_reads,min_edges', sorted(INTERVAL_EXPECTED))
def test_interval_graph_outputs(tmp_path, graph_seed, seed, min_reads, min_edges):
    nodes, edges = write_graph(tmp_path, 1500, interval_graph(graph_seed))
    run(tmp_path, nodes, edges, '-r', str(seed), '-N', str(min_reads), '-E', str(min_edges))
    digests = tuple(hashlib.sha256(output(tmp_path, name).encode()).hexdigest() for name in OUTPUTS)
    assert digests == INTERVAL_EXPECTED[(graph_seed, seed, min_reads, min_edges)]


def meta_weights(tmp_path, edges):
    """Read links between two meta nodes, from the read -> meta node column."""
    meta = [int(line.split(';')[-2]) for line in output(tmp_path, 'components').splitlines()[1:]]
    weights = collections.Counter()
    for a, b in edges:
        x, y = meta[a], meta[b]
        if x >= 0 and y >= 0 and x != y:
            weights[(min(x, y), max(x, y))] += 1
    return weights


def filtered_edges_nb(stdout):
    return int(re.search(r'(\d+) edges in the filtered graph', stdout).group(1))


@pytest.mark.parametrize('min_reads,min_edges', [(0, 1), (1, 1), (1, 2), (2, 0), (2, 3), (3, 1)])
def test_edge_filter_removes_weak_meta_edges(tmp_path, min_reads, min_edges):
    """-E T removes the meta edges of at most T read links, from both meta
    nodes; -N only filters meta nodes. Weights are counted here from the
    meta node of each read (the edge filter doesn't change meta nodes)."""
    edges = interval_graph(4, link=0.15)   # many meta edges of weight 1-3
    nodes, edges_csv = write_graph(tmp_path, 1500, edges)
    args = ['-r', '5', '-N', str(min_reads)]
    unfiltered = run(tmp_path, nodes, edges_csv, *(args + ['-E', '0']))
    weights = meta_weights(tmp_path, edges)
    assert filtered_edges_nb(unfiltered.stdout) == len(weights)
    meta_nodes = output(tmp_path, 'metaNodes')

    filtered = run(tmp_path, nodes, edges_csv, *(args + ['-E', str(min_edges)]))
    assert output(tmp_path, 'metaNodes') == meta_nodes
    kept = sum(1 for weight in weights.values() if weight > min_edges)
    if min_edges:
        assert 0 < kept < len(weights)
    # Half the neighbor entries: a one-sided removal would show here.
    assert filtered_edges_nb(filtered.stdout) == kept
    for line in output(tmp_path, 'metaEdges').splitlines()[1:]:
        a, b = map(int, line.split(';'))
        assert weights[(a, b)] > min_edges


def test_edge_columns_found_by_name(tmp_path):
    """Source and Target are found by name, in any order and case."""
    edges = interval_graph(3, reads_num=400)
    nodes, default = write_graph(tmp_path, 400, edges)
    run(tmp_path, nodes, default, '-r', '5', '-N', '1')
    expected = [output(tmp_path, name) for name in OUTPUTS]

    swapped = tmp_path / 'swapped.csv'
    swapped.write_text('weight;TARGET;source\n' + ''.join('1;{};{}\n'.format(b, a) for a, b in edges))
    run(tmp_path, nodes, swapped, '-r', '5', '-N', '1')
    assert [output(tmp_path, name) for name in OUTPUTS] == expected


def test_same_seed_same_components(tmp_path):
    nodes, edges = write_graph(tmp_path, 1500, interval_graph(2))
    results = []
    for _ in range(2):
        run(tmp_path, nodes, edges, '-r', '11', '-N', '1')
        results.append([output(tmp_path, name) for name in OUTPUTS])
    assert results[0] == results[1]


def test_edge_to_unknown_read_is_an_error(tmp_path):
    nodes, edges = write_graph(tmp_path, 3, [(0, 1), (1, 3)])
    proc = run(tmp_path, nodes, edges, '-r', '1', check=False)
    assert proc.returncode == 1
    assert 'unknown nodes' in proc.stderr


def test_missing_source_or_target_column(tmp_path):
    """As upstream: reported, and the graph is built without edges."""
    nodes, edges = write_graph(tmp_path, 3, [(0, 1)], edges_header='From;To\n', edge_format='{};{}\n')
    proc = run(tmp_path, nodes, edges, '-r', '1', '-N', '0')
    assert 'Source or Target label absent' in proc.stderr
    assert output(tmp_path, 'metaNodes') == 'node;size\n0;1\n1;1\n2;1\n'
    assert output(tmp_path, 'metaEdges') == 'from;to\n'


def test_missing_input_file(tmp_path):
    nodes, _ = write_graph(tmp_path, 1, [])
    proc = run(tmp_path, nodes, tmp_path / 'absent.csv', check=False)
    assert proc.returncode == 1
    assert 'seems to be missing' in proc.stderr
