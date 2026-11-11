#!/usr/bin/env python3
"""Report the reference sequences that do not reach the clustered db, with their taxonomy.

matam_db_preprocessing.py loses reference sequences on the way to the
clustered db (used for read recruitment and the overlap graph):

* ``removed by cleaning``: input sequences removed by the N or length filters.
  vsearch clusters the cleaned db, so these never compete to be centroids;
* ``centroid removed by cleaning`` / ``centroid not in the input``: with
  ``--clustered_ref``, given centroids that are not in the cleaned db. The
  cluster each one represented has no representative left;
* ``outside --domains``: with ``--by_domain``, cleaned sequences of other
  domains, which are deliberately not clustered.

What matters is whether a discarded sequence leaves a **close relative** among
the remaining references, not whether its taxon *name* survives: another
sequence named alike can be distant, and a close relative can carry another
name. So each discarded sequence (except other domains) is aligned to the
remaining references (the cleaned complete db) and its nearest one is
reported. It has a close relative left if that one has at least the
close-relative identity (default: the clustering identity) over at least
QUERY_COVERAGE of the discarded sequence. vsearch counts N and other IUPAC
codes as matching compatible bases, so masked stretches don't lower identity.
(The clustered db is no reference for this: a centroid is, by construction,
less similar than the clustering identity to every other centroid.)

Taxonomy: taxa are the prefixes of the taxonomy path, by level (1 = domain);
rank names are not used because rank counts differ (SILVA 144 prokaryotes:
domain;kingdom;phylum;class;order;family;genus;organism name). The organism
name (last field) is not a taxon; the remaining full path is the sequence's
lineage (e.g. the genus of a prokaryote). Eukaryote taxonomies are too uneven
to compare level by level, so eukaryotes are reported at the domain level.

Files:

* ``<db>.discarded.tsv``: one row per discarded sequence, with its reason,
  taxonomy, nearest remaining reference, identity, lineage, and whether a
  close relative is left;
* ``<db>.discarded_taxa.tsv``: one row per reason and taxon: discarded
  sequences, how many of them have no close relative left, and (by name) the
  sequences of that taxon left in the complete and clustered db.
"""

import collections
import logging

logger = logging.getLogger(__name__)

REMOVED_BY_CLEANING = 'removed by cleaning'
CENTROID_REMOVED_BY_CLEANING = 'centroid removed by cleaning'
CENTROID_NOT_IN_INPUT = 'centroid not in the input'
OUTSIDE_DOMAINS = 'outside --domains'
REASONS = (REMOVED_BY_CLEANING, CENTROID_REMOVED_BY_CLEANING, CENTROID_NOT_IN_INPUT,
           OUTSIDE_DOMAINS)
# Reasons whose sequences are searched for a close relative (other domains are
# left out of the clustered db on purpose).
SEARCHED_REASONS = (REMOVED_BY_CLEANING, CENTROID_REMOVED_BY_CLEANING, CENTROID_NOT_IN_INPUT)
QUERY_COVERAGE = 0.9
NO_HIT = '*'
UNKNOWN = 'NULL'
EUKARYOTA = 'eukaryota'
LOGGED_EXAMPLES = 5


def read_taxonomy(path):
    """{id: taxonomy string} from a MATAM ``.taxo.tab`` (``id<TAB>taxonomy``)."""
    taxonomy = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if fields[0]:
                taxonomy[fields[0]] = fields[1] if len(fields) > 1 and fields[1] else UNKNOWN
    return taxonomy


def _prefixes(taxo):
    """[(level, taxon)] of a taxonomy string, from the domain (level 1) to the lineage.

    The last field (the organism name of SILVA taxonomies) is left out.
    Eukaryotes only get their domain (see the module docstring).
    """
    if not taxo or taxo == UNKNOWN:
        return []
    ranks = [rank for rank in taxo.split(';') if rank]
    if ranks[0].lower() == EUKARYOTA:
        return [(1, ranks[0])]
    if len(ranks) > 1:
        ranks = ranks[:-1]
    return [(depth + 1, ';'.join(ranks[:depth + 1])) for depth in range(len(ranks))]


def lineage(taxo):
    """The most specific taxon of a taxonomy string (see _prefixes), or UNKNOWN."""
    prefixes = _prefixes(taxo)
    return prefixes[-1][1] if prefixes else UNKNOWN


def _count(ids, taxonomy):
    counts = collections.Counter()
    for seq_id in ids:
        for _, prefix in _prefixes(taxonomy.get(seq_id)):
            counts[prefix] += 1
    return counts


def write_report(basepath, discarded, taxonomy, complete_ids, clustered_ids, nearest=None,
                 close_identity=None):
    """Write the report files and return a summary (see summarise()).

    ``discarded`` maps sequence IDs to one of REASONS; ``taxonomy`` maps IDs to
    taxonomy strings (IDs it lacks are reported with an unknown taxonomy).
    ``nearest`` maps the searched IDs to (nearest remaining reference or NO_HIT,
    identity in %); ``close_identity`` (in %) is the close-relative identity.
    Without them (search disabled) closeness is reported as NA.
    """
    nearest = nearest or {}
    in_complete = _count(complete_ids, taxonomy)
    in_clustered = _count(clustered_ids, taxonomy)
    taxa = collections.Counter()               # (reason, level, taxon) -> discarded
    without = collections.Counter()            # (reason, level, taxon) -> no close relative
    lineages = collections.defaultdict(set)    # reason -> lineages of its discarded sequences
    by_reason = collections.Counter()
    domains = collections.defaultdict(collections.Counter)
    no_relative = collections.defaultdict(list)  # reason -> [(identity, id, lineage, nearest lineage)]
    searched = collections.Counter()
    with open(basepath + '.discarded.tsv', 'w') as handle:
        handle.write('id\treason\ttaxonomy\tnearest_ref\tnearest_ref_identity\t'
                     'nearest_ref_lineage\tclose_relative_left\n')
        for seq_id, reason in discarded.items():
            taxo = taxonomy.get(seq_id, UNKNOWN)
            by_reason[reason] += 1
            prefixes = _prefixes(taxo)
            domains[reason][prefixes[0][1] if prefixes else UNKNOWN] += 1
            if seq_id in nearest and close_identity is not None:
                ref, identity = nearest[seq_id]
                close = ref != NO_HIT and identity >= close_identity
                ref_lineage = lineage(taxonomy.get(ref)) if ref != NO_HIT else UNKNOWN
                searched[reason] += 1
                if not close:
                    no_relative[reason].append((identity, seq_id, lineage(taxo), ref_lineage))
                fields = [ref, '{0:.1f}'.format(identity) if ref != NO_HIT else 'NA',
                          ref_lineage, 'yes' if close else 'no']
            else:
                close, fields = None, ['NA', 'NA', 'NA', 'NA']
            handle.write('\t'.join([seq_id, reason, taxo] + fields) + '\n')
            for level, prefix in prefixes:
                taxa[(reason, level, prefix)] += 1
                if close is False:
                    without[(reason, level, prefix)] += 1
            if prefixes:
                lineages[reason].add(prefixes[-1][1])

    with open(basepath + '.discarded_taxa.tsv', 'w') as handle:
        handle.write('reason\tlevel\ttaxon\tlineage\tdiscarded\twithout_close_relative\t'
                     'left_in_complete_db\tleft_in_clustered_db\n')
        for key, count in sorted(taxa.items(), key=lambda item: (
                REASONS.index(item[0][0]), item[0][1], item[0][2])):
            reason, level, prefix = key
            handle.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(
                reason, level, prefix, 'yes' if prefix in lineages[reason] else 'no', count,
                without[key] if reason in searched else 'NA', in_complete[prefix],
                in_clustered[prefix]))
    return {'by_reason': by_reason, 'domains': domains, 'no_relative': no_relative,
            'searched': searched, 'close_identity': close_identity,
            'files': (basepath + '.discarded.tsv', basepath + '.discarded_taxa.tsv')}


def _short(taxon):
    """The last name of a taxon, for the log (the files have the full paths)."""
    return taxon.split(';')[-1]


def summarise(summary, log=logger):
    """Log the summary of write_report(), one line per reason."""
    if not summary['by_reason']:
        log.info('Discarded reference sequences: none; every input sequence is in the '
                 'cleaned db and could be clustered')
        return
    warned = False
    for reason in REASONS:
        count = summary['by_reason'].get(reason)
        if not count:
            continue
        domains = ', '.join('{0} {1}'.format(domain, number)
                            for domain, number in summary['domains'][reason].most_common())
        message = 'Discarded reference sequences ({0}): {1} ({2})'.format(reason, count, domains)
        if reason == OUTSIDE_DOMAINS:
            message += '; not clustered, by choice of --domains'
            warn = True
        elif summary['searched'].get(reason):
            missing = sorted(summary['no_relative'][reason])
            message += ('; without a close relative left in the references (>= {0:g}% identity '
                        'over >= {1:g}% of the sequence): {2}'.format(
                            summary['close_identity'], QUERY_COVERAGE * 100, len(missing)))
            if missing:
                examples = ['{0} {1} (nearest: {2})'.format(
                    _short(seq_lineage), seq_id,
                    '{0:.1f}% {1}'.format(identity, _short(ref_lineage))
                    if ref_lineage != UNKNOWN or identity else 'none')
                    for identity, seq_id, seq_lineage, ref_lineage in missing[:LOGGED_EXAMPLES]]
                message += ', e.g. {0}'.format('; '.join(examples))
            warn = bool(missing)
        else:
            message += '; close relatives not searched'
            warn = reason != REMOVED_BY_CLEANING
        (log.warning if warn else log.info)(message)
        warned = warned or warn
    (log.warning if warned else log.info)(
        'Discarded reference sequences and their taxa: {0}, {1}'.format(*summary['files']))
