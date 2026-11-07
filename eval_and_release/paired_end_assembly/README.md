# Paired-end assembly in this fork: behaviour, reasoning and comparison with upstream MATAM

This note covers how this fork handles paired-end reads from mapping to
assembly, the problems found on the way, why each fix was chosen, and how the
result compares with upstream MATAM. All numbers come from the bundled
16-species example and can be reproduced with the scripts in this directory
(see [Reproducing](#reproducing)).

**Baselines used throughout**

| Name | What it is |
|---|---|
| **Upstream** | bonsai-team MATAM at `e7cd9b2` (2025-03-28), the fork point. Run from bioconda `matam` 1.6.1, whose scripts are identical to `e7cd9b2`. |
| **Old fork rule** | This fork from `82c990d` to `264cf57`: only read pairs with both mates in the same component reach the assembler. |
| **Current fork** | This fork with the pair routing described in [section 4](#4-assembly-input-pair-routing-across-components). |

## Summary

- The fork assembles paired reads as pairs and keeps R1/R2 in separate files
  throughout. Mate identity is recovered after mapping, so read names never
  have to be changed: CASAVA (`@ID 1:N:0`), SRA (`@ID`) and legacy (`@ID/1`)
  naming give identical results.
- The switch to paired assembly in `82c990d` added a rule that dropped any
  pair whose mates landed in different graph components. On the example this
  threw away **46%** of the reads before assembly. Coverage of the reference
  by the final assembly fell from 97.9% (upstream) to 80.9%, and one species
  (*S. flexneri*) disappeared.
- Components are short (median ~200 bp along a gene), barely longer than a
  sequenced fragment, so mates often straddle a component boundary.
- Merging components linked by mates was tested and rejected: the links chain
  different species together.
- The current fork sends each pair to every component that holds one of its
  mates. Final reference coverage is **98.3%**, with all 16 genes recovered
  at 96–100% and an error rate of ~0%. That matches or slightly exceeds
  upstream.

## 1. Pipeline background

The steps that matter here:

1. **Mapping.** SortMeRNA maps all read pairs against the clustered reference
   DB. With `--paired_in --out2`, a pair is selected if either mate aligns,
   and the selected pairs are written in step to
   `marker_reads_R1.fq.gz` / `marker_reads_R2.fq.gz`.
2. **Alignment filtering** by score, and optional coverage subsampling.
3. **Overlap graph** (ovgraphbuild). One **node** per read (per mate). Two
   nodes share an edge when their alignments to the same reference overlap by
   at least `--min_overlap` bases with matching sequence.
4. **Components** (componentsearch). The graph is walked breadth-first from
   a random start (hence `--seed`). Connected reads at the same depth are
   merged into **metanodes**: stacks of reads covering about the same stretch
   of a gene. Short dead-end branches are absorbed, and metanodes or links
   below `-N` / `-E` are dropped. The connected groups that remain are the
   **components**. Each read gets a component number, or `-1` if it was
   filtered out.
5. **Assembly.** SGA runs once per component on that component's reads.
6. **Scaffolding.** Contigs from all components are mapped to the complete
   reference DB and joined into scaffolds. Redundant and short (< 500 bp)
   scaffolds are removed.

Scale on the example (5,825 pairs, seed 1, current fork): 11,650 nodes →
1,021 metanodes (266 after filtering) → 98 components. A component holds
about 116 reads and spans a median of 204 bp of a gene (interquartile range
180–280 bp). The 16 genes (~24 kb) are therefore cut into ~6 pieces each.
Cutting them is what lets MATAM separate similar species; scaffolding joins
the pieces back together.

## 2. Mate identity

SortMeRNA 7.1 writes FLAG `0`/`16` only, with no pairing bits, and gives
both mates the same QNAME when they share an ID (CASAVA naming). Before this
work every later step identified a read by QNAME alone. As a result:

- the two mates competed in the score filter, and the weaker one was dropped;
- both mates collapsed into one graph node;
- the component lookup could not tell the mates apart;
- abundance depended on how the reads were named.

On the end-to-end test subsample, CASAVA-named input gave 1,457 graph nodes
instead of 2,914.

**Fix.** `scripts/sam_mates.py` restores the standard SAM pairing flags
(0x1, 0x40/0x80, and 0x8 when the mate has no alignment) right after
mapping. It matches each record's SEQ, reverse-complemented when 0x10 is set,
against the R1/R2 output of `--out2`, streaming the SAM and the pairs in
lockstep. Any inconsistency is a hard error. Every later step keys reads on
(QNAME, mate). Verified on all 41,986 records of the full example: every
record assigned, 0 ambiguous.

User input is never renamed. The tagged reads keep their names, and the mate
is carried in the flags and in explicit `mate` table columns (nodes CSV
`Id;Label;Mate;Specie`).

## 3. The assembly-input problem

### 3.1 How each version chooses the reads to assemble

| Version | Rule |
|---|---|
| Upstream (`e7cd9b2`) | `extract_reads_by_component` puts every read of a component into that component's file. SGA always assembles single-end; the wrapper's `--paired_end` branch was a stub that MATAM never passed. A read is dropped only if it is in no component. |
| Old fork rule (`82c990d`) | Reads go to SGA as pairs, so a pair is kept only if **both mates are in the same component**. With CASAVA names the mates were merged into one node at the time, so the check never fired; with distinct mates it does. |
| Current fork | Each pair goes to every component holding one of its mates (section 4). |

### 3.2 Why so many pairs were lost

What happens to the 5,825 pairs under the old fork rule (seed 1):

| Pair outcome | Pairs | Mates sent to SGA |
|---|---|---|
| both mates in the same component | 3,173 (54%) | 6,346 |
| mates in different components | 2,387 (41%) | 0 |
| one mate filtered out of the graph | 259 (4%) | 0 |
| both mates filtered out | 6 | 0 |

The cause is geometry. In this simulated library, the outer span of a pair
on the reference is a median of **149 bp** (interquartile range 143–156); the
100 bp mates overlap each other by about 50 bp. Components span a median of
**204 bp**. A 150 bp fragment fits inside a ~200 bp component only if it
starts in roughly the first 50 bp; otherwise it crosses a boundary and the
pair is lost. The split is not evidence that the mates come from different
places: they overlap and are usually directly linked in the overlap graph.
componentsearch places a boundary between them anyway.

Real libraries (150–250 bp reads, 300–550 bp inserts) have fragments longer
than most components here, so the old rule would probably discard even more
of their data. This hasn't been measured; longer reads may also give longer
components.

## 4. Assembly input: pair routing across components

### 4.1 Rejected: merging components linked by mates

A pair spanning two components suggests the components are continuous, so
the first idea was to merge mate-linked components into one assembly unit.
On seed 1 it fails:

- 38 of the 98 components are multi-species (< 90% of their reads from one
  species): conserved regions shared by up to 8 species.
- There are 194 component links (1–65 pairs each, median 9). Only 525 of the
  2,387 linking pairs join components with the same dominant species.
- Merging transitively through the shared components chains species
  together:

| Merge components linked by at least | Groups left | Largest group | Multi-species groups |
|---|---|---|---|
| 1 pair | 3 | 6,407 reads, 9 species | 2 |
| 5 pairs | 7 | 6,405 reads, 9 species | 3 |
| 10 pairs | 24 | 5,946 reads, 9 species | 4 |
| 20 pairs | 59 | 4,895 reads, 8 species | 10 |

SGA assembles each group as a unit, so merging would invite chimeric contigs
at every threshold.

### 4.2 Chosen: send each pair to every component of its mates

Mate links are valid evidence at the level of a single pair, since a pair is
one molecule from one species. `components_assembly._split_pairs_by_component`
now:

- sends a pair with both mates in one component to that component;
- sends a pair spanning two components, complete, **to both**;
- sends a pair with one mate in a component (the other filtered out or
  unaligned) to that component;
- leaves out only pairs with no mate in any component;
- never merges components.

Each component's assembly thus receives all of its own reads plus the
partner mates that reach across its boundaries. That is extra sequence from
the same molecule, so contigs can extend past the boundary without mixing
species. SGA uses the pairing only in `sga preprocess --pe-mode 1`, which
drops both mates if one fails its filters; every later SGA step treats them
as plain reads.

The routing is logged, for example:

```
Read pairs by component: 3173 within one component, 2387 spanning two (given to both), 259 with one mate in a component, 6 in no component (not assembled)
```

Not tested: keeping same-component pairs as pairs and passing split-pair
mates as unpaired reads. With SGA using pairs only at preprocessing, it is
not expected to differ much from the chosen routing.

## 5. Behaviour comparison

### 5.1 Method

- **Data:** the full bundled example, `examples/16sp_simulated_dataset`
  (5,825 ART pairs, 100 bp reads, 16 species).
- **References:** the 16 references with SILVA-style `>ID taxonomy` headers,
  the same FASTA for every version. Each version builds its own DB with its
  own `matam_db_preprocessing.py`.
- **Input:** upstream gets its documented input, one interleaved file with
  `/1`, `/2` names; the fork gets the same reads as R1/R2 files.
- **Settings:** defaults, `--cpu 8 --max_memory 4000 -v`, seeds 1, 2 and 3
  (componentsearch is random).
- **Evaluation:** contigs and final assemblies of every run go through this
  repository's `scripts/evaluate_assembly.py`, which aligns them to the true
  references with exonerate (`affine:overlap`, exhaustive, best hit). Error
  rate and reference coverage come from its output.
- **Tools:** upstream uses SortMeRNA 2.1b and the bioconda 1.6.1 binaries;
  the fork uses SortMeRNA 7.1 and this repository's ovgraphbuild /
  componentsearch builds.

### 5.2 Results (mean of 3 seeds)

| Metric | Upstream | Old fork rule | Current fork |
|---|---|---|---|
| Graph nodes (mates) | 11,650 | 11,650 | 11,650 |
| Graph edges | 424,569 | 421,351 | 421,351 |
| Components | 107.3 | 102.7 | 102.7 |
| Mates in components | 11,377 | 11,391 | 11,391 |
| **Reads given to SGA** | 11,377 | **6,193** | 16,589 |
| Components given to SGA | 106.3 | 91.3 | 102.7 |
| Components with contigs | 103.7 | 79.3 | 102.3 |
| Contigs | 468 | 271 | 655 |
| Contigs total length | 63.2 kb | 39.7 kb | 97.3 kb |
| Contig N50 | 128 | 138 | 149 |
| Contig error rate | 0.02% | 0.02% | 0.03% |
| Reference covered by contigs | 97.9% | 94.6% | 98.2% |
| Final scaffolds (≥ 500 bp, non-redundant) | 16 | 16.7 | 16 |
| Final total length | 23.6 kb | 19.5 kb | 23.6 kb |
| Final N50 | 1,492 | 1,402 | 1,499 |
| Longest scaffold | 1,536 | 1,506 | 1,538 |
| Final error rate | ~0% | ~0% | ~0% |
| **Reference covered by final assembly** | 97.9% | **80.9%** | **98.3%** |
| **Genes ≥ 90% covered (of 16)** | 16 | **9.7** | **16** |
| Genes hit at all (of 16) | 16 | 15 | 16 |
| Wall time | 20.3 s | 8.5 s | 11.6 s |

Per-seed values are in `summarise.py`'s output (one column per run).

### 5.3 Observations

- **Graph:** all versions build essentially the same graph and place the
  same number of reads in components. The ~0.8% edge difference between
  upstream and the fork comes from the different SortMeRNA versions and
  binaries.
- **Old fork rule:** the 46% read loss shows up as fewer, shorter contigs,
  and 11–13 components that do receive reads produce no contigs (upstream:
  2–4). Most of the coverage is lost at scaffolding: 94.6% by contigs but
  80.9% by the final assembly.
- **Lost species:** under the old rule *S. flexneri* is missing from the
  final assembly in every seed, even though its contigs cover it at 94–95%;
  no *S. flexneri* scaffold survives scaffolding. It is plausibly merged into
  or deduplicated against the near-identical *E. coli* scaffold; this was not
  confirmed. The related Enterobacteriaceae (*E. coli*, *S. flexneri*,
  *K. pneumoniae*, *S. enterica*) keep only 34–46% of their component reads
  as same-component pairs.
- **Current fork:** every component produces contigs. All 16 genes are
  covered at 96–100% in every seed (*S. flexneri* at 97–99%), and the final
  assembly matches upstream on length and N50.
- **Read count:** SGA gets more reads than upstream (16,589 vs 11,377)
  because spanning pairs are counted in both components. That duplication
  is what lets contigs reach across boundaries.
- **Accuracy:** the error rate stays at ~0%, so there is no sign of
  chimeric contigs from the extra reads.
- **Speed:** the current fork is ~3 s slower than the old rule because of the
  extra reads, and still faster than upstream.

### 5.4 Caveats

- One simulated dataset with short fragments (~150 bp). A real dataset with
  longer inserts is needed before generalising, in particular for the
  fragment-versus-component geometry of section 3.2.
- Upstream and the fork differ in tool versions (SortMeRNA 2.1b vs 7.1,
  ovgraphbuild, componentsearch and SGA builds), not only in the pair rule.
  The old-rule versus current-fork comparison uses the same tools, so it
  isolates the rule; the upstream column is context.
- The evaluation is reference-based (exonerate best hits). It measures
  completeness and base accuracy, not taxonomic assignment.

## 6. Tests

Run in the `testmatam` environment: `python -m pytest tests/ -q`.

| Area | Tests |
|---|---|
| Mate tagging | `test_sam_mates.py`: SEQ-based assignment, reverse-strand records, mates without alignment (0x8), identical-sequence ties, order divergence and secondary/supplementary records as errors, gzip input. `test_sortmerna_pe_contract.py`: SortMeRNA 7 assumptions (no pairing flags, `--out2` in step, record order). |
| Filtering | `test_filter_score_multialign.py`: mates no longer compete; straight and geometric thresholds within one mate. |
| Graph | `test_ovgraphbuild_mates.py`: mates sharing a QNAME are two nodes; all alignments of one read form one node; flagged shared names build the same graph as legacy `/1`, `/2` names; non-contiguous records are rejected. |
| Tables and LCA | `test_pe_tables_and_lca.py`: (read, mate) tables, multi-column `-g`, LCA grouping across factors. |
| Assembly input | `test_sga_paired_contract.py`: a spanning pair goes to both components, a pair with one mate in a component goes there, the routing log line, CASAVA names looked up by mate, out-of-sync R1/R2 files are fatal. |
| Abundance | `test_compute_abundance.py`: one count per mate; mates hitting different scaffolds. |
| End to end | `test_functional_pe.py`: builds a 16-reference DB and runs MATAM on CASAVA (gzipped), SRA-style and legacy copies of the same pairs. It requires identical nodes, edges, components, pair routing, scaffolds and abundances across the three namings, also with `--coverage_threshold`, and checks the `--filter_only` R1/R2 output. |

The comparison in section 5 isn't part of the test suite: it needs upstream
MATAM and takes several minutes.

## Reproducing

Prerequisites:
- this fork built (`python build.py`) and its dependencies in an
  environment, here `testmatam`; `exonerate` is needed for the evaluation;
- upstream MATAM in a separate prefix. SortMeRNA 2.1b can't share an
  environment with SortMeRNA 7. Flexible channel priority is needed because
  SortMeRNA 2.1b pins an old zlib:

```bash
mamba create -p /path/to/matam_upstream --channel-priority flexible \
    -c conda-forge -c bioconda matam=1.6.1
```

Upstream versus the current fork (3 seeds each):

```bash
conda activate testmatam
UPSTREAM_PREFIX=/path/to/matam_upstream \
    eval_and_release/paired_end_assembly/run_compare.sh /path/to/out
```

The old fork rule, from a worktree of `264cf57`, into the same output directory:

```bash
git worktree add /path/to/matam_264cf57 264cf57
(cd /path/to/matam_264cf57 && git submodule update --init --recursive && python build.py)
SKIP_UPSTREAM=1 FORK_REPO=/path/to/matam_264cf57 FORK_LABEL=oldrule \
    eval_and_release/paired_end_assembly/run_compare.sh /path/to/out
```

Summary tables, with evaluation against the true references:

```bash
eval_and_release/paired_end_assembly/summarise.py /path/to/out up oldrule fork
```

The pair-outcome breakdown (section 3.2) and the merge analysis (section 4.1)
were computed from a run's `workdir/*.components.csv` and filtered SAM, with
species taken from the read-name prefix (e.g. `SFL_X96963-…`).
