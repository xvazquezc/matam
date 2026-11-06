#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Read-level tables keyed by (read, mate).

A read is identified by its QNAME and its mate number (1 or 2 for paired
reads, 0 for unpaired reads; see ``sam_mates.py``).  These helpers replace the
shell ``cut``/``join`` pipelines that keyed reads on the QNAME alone, which
merged the two mates of CASAVA-named pairs.

Table formats (tab-separated, no header):

* read-component table: ``read  mate  metanode  component``
* complete taxonomy table: ``read  mate  metanode  component  taxonomy``
  (one row per filtered alignment)
"""

import collections

from sam_mates import mate_of

GraphStats = collections.namedtuple('GraphStats', 'nodes reads_in_components components')


def iter_components_csv(path):
    """Yield (read, mate, metanode, component) from ComponentSearch output.

    The file holds the ovgraphbuild node columns (``Id;Label;Mate;Specie``)
    followed by ``metanode;component``.  Files without a ``Mate`` column
    (older ovgraphbuild) give mate 0.
    """
    with open(path) as handle:
        header = handle.readline().rstrip('\n').split(';')
        columns = {name.lower(): index for index, name in enumerate(header)}
        for required in ('label', 'metanode', 'component'):
            if required not in columns:
                raise ValueError('{} has no {!r} column: {}'.format(path, required, header))
        mate_index = columns.get('mate')
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split(';')
            if len(fields) != len(header):
                raise ValueError('Unexpected number of columns in {}: {!r}'.format(path, line))
            mate = int(fields[mate_index]) if mate_index is not None else 0
            yield (fields[columns['label']], mate,
                   int(fields[columns['metanode']]), int(fields[columns['component']]))


def write_read_component_table(components_csv, out_path):
    """Write reads that belong to a component (component >= 0).

    Returns GraphStats(nodes, reads_in_components, components).
    """
    nodes, reads, components = 0, 0, set()
    with open(out_path, 'w') as out:
        for read, mate, metanode, component in iter_components_csv(components_csv):
            nodes += 1
            if component < 0:
                continue
            reads += 1
            components.add(component)
            out.write('{}\t{}\t{}\t{}\n'.format(read, mate, metanode, component))
    return GraphStats(nodes, reads, len(components))


def load_read_components(path):
    """Return {(read, mate): (metanode, component)} from a read-component table."""
    table = {}
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 4:
                raise ValueError('Unexpected read-component row in {}: {!r}'.format(path, line))
            table[(fields[0], int(fields[1]))] = (fields[2], fields[3])
    return table


def load_ref_taxonomy(path):
    """Return {reference_id: taxonomy} from a two-column taxonomy table."""
    taxonomy = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 2:
                taxonomy[fields[0]] = fields[1]
    return taxonomy


def write_complete_taxo_table(sam_path, ref_taxo_path, read_component_path, out_path):
    """Join filtered alignments with reference taxonomies and read components.

    Alignments of reads outside components, or on references without a
    taxonomy, are skipped. Returns the number of rows written.
    """
    taxonomy = load_ref_taxonomy(ref_taxo_path)
    read_components = load_read_components(read_component_path)
    rows = 0
    with open(sam_path) as sam, open(out_path, 'w') as out:
        for line in sam:
            if not line.strip() or line.startswith('@'):
                continue
            fields = line.split('\t', 3)
            key = (fields[0], mate_of(fields[1]))
            located = read_components.get(key)
            taxo = taxonomy.get(fields[2])
            if located is None or taxo is None:
                continue
            out.write('{}\t{}\t{}\t{}\t{}\n'.format(key[0], key[1], located[0], located[1], taxo))
            rows += 1
    return rows
