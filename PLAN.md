# MATAM modernisation plan
# Repo: xvazquezc/matam  Branch: master

---

## Current repository state

- Upstream base: `v1.6.2` (`e7cd9b2`)
- Fork remote: `git@github.com:xvazquezc/matam.git`
- All the work below is on `origin/master` (2026-11-11). Later changes
  (paired-end rework, index records, `--stop_after`, `--version`, input
  validation) are tracked in `TESTING_PLAN.md`.

## Completed work (8 commits on top of v1.6.2)

| Commit | Message |
|--------|---------|
| e7b1799 | fix(#124): move module-level binary assertions to point of use |
| 821999e | feat: add gzip support for fastq input files |
| 9f5732f | feat: add paired-end reads support via --forward/--reverse |
| d5c230f | feat: migrate to SortMeRNA v4+ API |
| d24feae | perf: replace custom Python fastx utilities with seqkit and cd-hit-est |
| f619e65 | feat: migrate to SortMeRNA v5+ API and Python >=3.14 (environment.yml only) |
| c5f5a06 | improve db prep |
| d371b1d | perf(ovgraphbuild): dynamically schedule pair-processing chunks |

All of these commits are on the fork.

---

## SortMeRNA v5+ migration — completed

### Phase A — environment update (`f619e65`)

| Was | Now |
|-----|-----|
| `python >=3` | `python >=3.14` |
| `sortmerna>=4.0` | `sortmerna>=5.0` |
| `quast==5.0.2` | *(removed — no Python 3.14 build; optional test dep only)* |

### Phases B and C — CLI and work-directory update (`c5f5a06`)

- Replaced deprecated `-a` with `--threads` in both assembly mappings and
  abundance mapping.
- Replaced `--index 1` with explicit index-only `--task 5` in
  `index_ref_db.py`, `matam_db_preprocessing.py`, and `compute_abundance.py`.
- Added `matam_assembly.py --workdir`; it defaults to the output directory.
- Kept immutable reference indexes under `<database-prefix>.smr/idx` while
  putting SortMeRNA run state under the selected working directory.
- Confirmed both complete and clustered FASTA references are indexed during
  database preparation.

### Phase D — index compatibility documentation

The README warns about the SortMeRNA 6 BBHash index-format break. Users upgrading from SortMeRNA 4/5
must rebuild both `.smr` indexes.

---

## `ovgraphbuild` parallelisation and sharding — completed locally (`d371b1d`)

- Added `--threads`, `--pair-shard-index`, and `--pair-shard-count`.
- MATAM passes `--cpu` to `ovgraphbuild --threads`.
- Split the triangular read-pair space into pair-count-balanced contiguous
  ranges.
- Dynamically schedules eight deterministic microchunks per worker to reduce
  imbalance caused by concentrated shared-reference comparisons.
- Caps graph workers at 24 while accepting and reporting larger requested CPU
  counts.
- Uses thread-local statistics and output fragments, then merges fragments in
  range order to preserve byte-identical node and edge output.
- Emits shard range, pair-count, requested/effective thread, thread-limit, and
  chunk-count metadata.
- Cleans temporary fragments and reports output/worker failures as normal
  non-zero exits.

### Validation

- Test suite: `25 passed, 6 skipped`.
- CSV output is byte-identical to the previous serial implementation for every
  tested thread count.
- ASQG output is byte-identical at 1 and 24 threads.
- Eight external shards cover all pairs once and reconstruct the serial edge
  order exactly.
- Larger test: 11,650 reads, 88,105 SAM records, and 67,855,425 candidate
  read pairs.

| Threads | Wall time | Speedup |
|--------:|----------:|--------:|
| 1 | 15.00 s | 1.00x |
| 2 | 7.89 s | 1.90x |
| 4 | 4.53 s | 3.31x |
| 8 | 2.95 s | 5.08x |
| 16 | 1.80 s | 8.33x |
| 24 | 1.65 s | 9.09x |

Peak memory remained approximately 130 MB. Scaling flattens between 16 and 24
workers, so the 24-worker cap is retained.

---

## Production-size benchmark (synthetic, 2026-11-11)

The slow production SAM wasn't available, so two production-size filtered SAMs
were made with real MATAM steps (`--stop_after alignments_filtering`):
simulated ART 2x150 reads, mapped against a 51,483-reference SILVA 144 subset
(10 families, 95% clustering).

- *mixed*: a 34-taxon community (controls, strain pairs, novel pairs);
- *skewed*: the same plus 5 Enterobacteriaceae (E. coli/Shigella,
  Salmonella, Klebsiella, Enterobacter, Citrobacter) at 2500x. That is a
  dense, conserved part of SILVA, so most of those reads share references.

AMD Ryzen AI 9 HX 370 (12 cores, 24 threads). Serial = 1 thread on shards
0 and 40 of 64, averaged and scaled by 64.

| | mixed | skewed |
|---|---:|---:|
| Mapped reads (mates) | 169,624 | 180,320 |
| SAM records | 929,301 | 1,014,978 |
| Candidate read pairs | 14.39 G | 16.26 G |
| Compatible (edges) | 77.6 M (0.54%) | 149.1 M (0.92%) |
| Neither compatible nor incompatible | 13.61 G (94.6%) | 14.72 G (90.5%) |
| Wall, 24 threads, all pairs | 613 s | 1,034 s |
| Serial estimate | ~5,450 s (1.5 h) | ~9,700 s (2.7 h) |
| Speedup at 24 threads | ~8.9x | ~9.4x |
| Serial throughput | ~2.6 M pairs/s | ~1.7 M pairs/s |
| Peak memory | 0.99 GB | 1.08 GB |

Shard costs vary with the taxa they hold, because ART read names are grouped
by taxon after the QNAME sort: shards 0 and 40 of *mixed* took 118 s and 53 s.

**What this settles, and what it doesn't**

- At the production read count (170k), a 24-thread job builds the graph in
  10-17 min on these data. Time grows with the square of the read count.
- The production run (~100 h serial for 170k reads, i.e. ~40 k pairs/s) is
  40-65x slower per pair. These data don't reproduce that. Likely causes, not
  verified: up to 10 alignments per read against full SILVA (5.5 here), and
  far more pairs sharing references and overlapping. Only the production SAM
  can confirm.
- With the ~9x thread speed-up, a 100 h serial graph is ~11 h in one 24-thread
  job, within a 48 h walltime. Splitting the graph into separate jobs is
  therefore not needed for that case. It would be needed only for graphs that
  stay above the walltime at 24 threads (e.g. 2x the reads: ~44 h).
- 90-95% of the candidate pairs end as "neither compatible nor incompatible":
  the two reads have no reference in common (pairs that share a reference
  but don't overlap on it are counted as incompatible). They never give an
  edge. The all-pairs loop still tests every pair (up to
  10x10 alignment comparisons each).

## Remaining work

1. **Done (2026-11-12, not committed yet): candidate-pair pruning in
   `ovgraphbuild`.**
   - **What it does:** by default each read is compared only with later reads
     that share a reference, found through a reference → reads index built
     once. They are compared in ascending order.
   - **Statistics:** the skipped pairs are added in bulk to the "neither" and
     "incompatible" counts. Edges, nodes and statistics are byte-identical to
     the all-pairs loop, which `--all-pairs` keeps.
   - **Threads and shards:** both work as before. The new
     `compared_pair_count` is written to stdout and to the shard metadata.
   - **Tests:** `tests/test_ovgraphbuild_pair_pruning.py` (12). The ASan and
     UBSan build is clean.
   - **Benchmark** (SILVA 144 subset sets above; 24 threads, plus 1 thread on
     shards 0 and 40 of 64; outputs identical to the baseline):

     | Set | Pairs compared | 24 threads, wall | CPU | 1 thread, shards 0 / 40 |
     |---|---:|---:|---:|---:|
     | mixed, 161,680 mates | 455 M of 13.07 G (3.5%) | 544 → 388 s (1.40×) | 9,811 → 6,188 s (1.59×) | 46.6 → 26.3 s / 47.5 → 22.6 s |
     | skewed, 161,846 mates | 729 M of 13.10 G (5.6%) | 697 → 543 s (1.28×) | 14,273 → 10,590 s (1.35×) | 40.4 → 18.8 s / 85.7 → 61.6 s |

     The index adds about 35 MB to peak memory.
   - **Why the gain is modest:** about 95% of the pairs are skipped, but they
     were only 25-40% of the CPU time. The rest is spent on the pairs that
     share a reference:
     - the 52 M (mixed) or 83 M (skewed) compatible pairs, each re-aligned
       de novo with dynamic programming (`align2ReadSequences`);
     - the overlap walk.

     gprof (`-O2 -pg -fno-inline`, skewed shard 40) puts about a third of the
     samples in the seqan DP code and about a third in seqan string iterators.
     Its attribution of anonymous-namespace symbols is unreliable, and
     `-fno-inline` inflates the iterator cost.
   - **Done (2026-11-20): shards and thread chunks balanced by expected
     work, not pair count.**
     - Work of read i: 1, plus one per (later read, shared reference),
       computed from the reference index with binary searches.
     - `--all-pairs` still balances by pair count. Shard metadata gains
       `assigned_work` and `total_work`. Outputs are unchanged.
     - Benchmark (same sets; 16 shards run at once at 1 thread each):

       | Set | 24 threads, wall | Busy threads | Slowest shard (CPU) | Shard max/mean |
       |---|---:|---:|---:|---:|
       | mixed | 402 → 304 s (1.32×) | 15.7 → 21.7 | 867 → 457 s | 3.04 → 1.45 |
       | skewed | 531 → 481 s (1.10×) | 19.9 → 22.1 | 996 → 711 s | 2.24 → 1.44 |

       Compared pairs per shard went from 13.6-75.0 M to 20.9-34.0 M
       (mixed) and from 11.9-78.0 M to 42.8-50.7 M (skewed).
     - **What is left of the imbalance:** shard CPU now follows the number
       of compatible pairs (edges), not the comparisons. A least-squares
       fit over the 16 shards (no intercept, R² 0.96-0.97) gives 2.6-2.9 s
       per M comparisons and about 72 s per M compatible pairs (the de novo
       DP): the compatible pairs are 74-76% of the CPU.
     - Tests: `test_ovgraphbuild_pair_pruning.py` (lopsided input: by pair
       count, 4 shards compare `[0, 0, 0, 19900]` pairs; now each about a
       quarter), and `--all-pairs` keeps pair-count shards.
   - **Possible follow-ups:**
     - weight the balancing by an estimate of compatible pairs: count the
       later reads that overlap by >= `minOverlap` on a shared reference
       (interval counting on sorted starts). The same count could prune the
       pairs that share a reference but can't overlap;
     - ~~make the de novo check cheaper~~: done exactly, see below.
   - **Done (2026-11-20, not committed yet): perfect overlaps skip the full
     de novo alignment.**
     - **Profile first** (cycle counters on 4 single-thread shards): the
       SeqAn DP of `align2ReadSequences` was 84-88% of the pair-comparison
       time; `computeAlignmentStats` and the output were about 1% each (the
       gprof split above was wrong). 94-97% of re-aligned pairs became edges.
     - **Why it is exact:** the de novo alignment scores +1/-1/-1 with free
       end gaps. Every alignment ending at a last-row/column cell E scores at
       most L, the length of the diagonal from E back to the first
       row/column, and only the ungapped all-match diagonal reaches L. So when
       the best score is reached at one end cell only and equals L there,
       SeqAn has a single optimal alignment to return, with no tie-breaking.
     - **How:** `perfectOverlap.cpp` computes the score only, 16 cells of an
       anti-diagonal at a time (GCC vector extensions; plain SSE2 on the
       conda build). When the condition holds it fills the statistics and
       overlap exactly as `computeAlignmentStats` would; otherwise (ties,
       mismatches, indels) SeqAn runs as before. `--full-realignment` keeps
       SeqAn for every pair; `--check-perfect-overlaps` (testing) runs both
       and stops on a difference. Re-aligned and shortcut counts are in the
       `-v` output and the shard metadata.
     - **Checks:** 2.3 M random pairs (repeats, homopolymers, IUPAC,
       substitutions, indels; 0 differences, also under ASan/UBSan); every
       shortcut pair of both benchmark sets compared with SeqAn (0
       differences); a deliberately broken tie test is caught by the unit
       harness and by `test_tied_overlaps_take_the_seqan_route`.
     - **Benchmark** (same sets; edges, nodes, merged shards and statistics
       identical):

       | Set | Shortcut pairs | 24 threads, wall | CPU | 16 shards, max CPU |
       |---|---:|---:|---:|---:|
       | mixed | 52.20 M of 52.94 M (98.6%) | 304 → 36 s (8.4×) | 6,637 → 744 s (8.9×) | 484 → 57 s |
       | skewed | 82.94 M of 85.41 M (97.1%) | 480 → 69 s (6.9×) | 10,657 → 1,504 s (7.1×) | 736 → 113 s |

       The SeqAn alignments left (0.7 M and 2.5 M pairs) are about 4% and 12%
       of the CPU; the rest is the alignment-pair comparison loop.
     - Tests: `tests/test_ovgraphbuild_perfect_overlap.py` (32).
2. **Done (2026-11-23): production benchmark.**
   - **Data:** the MarkerMAG-nf run on Gadi (`ferrari/dnasip`, MATAM 1.6.0,
     8 CPUs, 24 h walltime) whose 50-100% subsamples timed out. Filtered SAMs
     and the clustered SILVA 138.2 NR99 reference (422,709 sequences) copied
     here. Reads are MarkerMAG's `.1`/`.2`-suffixed reads, so no mates are
     detected (`Mate` column all 0).
   - **The earlier premise was wrong.** The slow graphs were not 170k reads at
     40-65× the synthetic cost per pair. They had 0.68-1.97 M reads, and the
     serial rate was about 2.9 M pairs/s on Gadi, like the synthetic sets.
   - **Results** (this machine, 16 threads, `nice`; Gadi = 1.6.0 serial):

     | Set | Reads (records) | Pairs | Compared (shared ref) | Edges | Gadi 1.6.0 | Here, 16 threads | Peak RSS |
     |---|---:|---:|---:|---:|---:|---:|---:|
     | Heavy 100% | 679,088 (2.88 M) | 230.6 G | 11.1 G (4.8%) | 294.4 M | 79,956 s (22.2 h) | 951 s | 4.1 GB |
     | Bulk 50% | 984,506 (4.88 M) | 484.6 G | 14.0 G (2.9%) | 353.8 M | 134,380 s (37.3 h) | 1,390 s | 6.2 GB |
     | Bulk 100% | 1,969,143 (9.75 M) | 1.94 T | 55.8 G (2.9%) | 1.417 G | timed out | 5,878 s (1 h 38 min) | 11.4 GB |

     Heavy and Bulk 50% edges are byte-identical to the 1.6.0 outputs (md5),
     and so are the nodes once the new `Mate` column is removed.
   - **Serial rate** (Heavy, 1 thread, 4 shards of 64, mean × 64):

     | Mode | Serial (here) |
     |---|---:|
     | `--all-pairs --full-realignment` (≈ 1.6.0) | 47,870 s |
     | `--full-realignment` (pruning only) | 27,800 s |
     | default (pruning + shortcut) | 10,370 s (4.6× less than ≈ 1.6.0) |

     The ≈ 1.6.0 mode against the Gadi time gives Gadi `normalbw` ≈ 1.67×
     slower per core than this machine. The shortcut covers 86.7% of the
     re-aligned pairs (97-99% on the synthetic sets) and gives 2.7× on real
     data. What is left: about 25% SeqAn on the other 13%, 75% the comparison
     loop over the 4.8% of pairs that share a reference.
   - **Gadi estimates for the graph step** (×1.67, from the serial rates):
     Heavy 100% 4.8 h serial, Bulk 50% 7.1 h, Bulk 100% about 30 h; with 8
     CPUs roughly 0.7 h, 1 h and 4-5 h.
   - **New bottlenecks.**
     - `componentsearch` (graph compaction, single-threaded): Heavy 100%
       took 703 s here (20 min on Gadi) and **14.6 GB peak, about 50 bytes
       per edge**. Bulk 100% (1.417 G edges) would need about 70 GB and over
       1 h on Gadi. Fixed by the rewrite in item 5.
     - SortMeRNA mapping: 2.3 h (Heavy) and 7.5 h (Bulk 100%) per subsample
       on Gadi at 8 CPUs. "Map once" (item 4) does this once per sample.
3. ~~Add MarkerMAG-nf scatter/merge orchestration only if a graph is still
   over the walltime at 24 threads after 1 and 2.~~ Not needed: the largest
   production graph (Bulk 100%) is an estimated 4-5 h on 8 Gadi CPUs. The
   limit was `componentsearch` memory (see items 2 and 5).
4. **Done (2026-11-20, not committed yet): map once, subsample afterwards.**
   - `--read_alignments SAM` takes the read mapping from SortMeRNA alignments
     made earlier. It writes the same files the mapping would
     (`scripts/read_alignments.py`), plus a completion record so that resuming
     works.
     - The SAM's `@PG` line must show the same database, `--best` and `-e`.
     - The clustered index is not needed.
   - Every run publishes `marker_alignments.sam.gz`.
   - MarkerMAG-nf uses this for all its subsamples (see its
     `docs/matam-performance-plan.md`).
   - **Not exact.** SortMeRNA's minimum score grows with the total read length
     mapped together, and `--num_alignments` keeps the first alignments that
     reach it. Minimum scores were 65, 68 and 72 for 1%, 10% and 100% of the
     mixed set, and 84 at an emulated raw-read size.
     - After filtering, 1–3% of a subsample's reads had different alignments.
     - On the same reads, the filtered SAM and the scaffolds are
       byte-identical.
     - On 1–50% subsamples, the assemblies differed less from the
       per-subsample ones than two MATAM seeds differ from each other.
       Recovery of the true sequences was unchanged within the seed range.
   - Tests: `tests/test_read_alignments.py` (14) and a SortMeRNA contract test
     for the threshold effect. The docstring of the older subsetting test was
     corrected: its property holds only away from the threshold.

---

5. **Done (2026-11-24, not committed yet): `componentsearch` rewrite, same
   output.**
   - `componentsearch` is the upstream submodule
     (`bonsai-team/ComponentSearch2`, f705924). It held the read graph as
     `std::map` adjacency sets and copied it into the meta graph, about 50
     bytes per edge.
   - The rewrite keeps the read graph once, as a compressed adjacency array
     (8 bytes per edge, filled in two passes over the edges file). Meta nodes
     hold read indices only, and every step changes the meta graph in place.
     A Fenwick tree picks the k-th remaining read for each BFS start, so the
     `rand()` sequence and every tie are as upstream.
   - **Same output for the same seed, quirks included.** Upstream passes `-N`,
     not `-E`, to the edge filter, and that filter counts links by read index,
     not by meta node. Both are kept; whether to fix them is a separate
     decision, to be judged on assemblies.
   - One behaviour change: an edge to a read missing from the nodes file is an
     error (exit 1) instead of undefined behaviour.
   - **Validation:**
     - 600 random graphs (sparse, trees, stars, interval, dense, duplicates
       and self-loops; random seeds, `-N`, `-E`, `-o`) against the upstream
       binary: stdout and the three outputs byte-identical, also with
       ASan/UBSan.
     - Production graphs, `-r 42 -N 1 -E 1`: components, meta nodes, meta
       edges, stdout and stderr byte-identical.

       | Graph | Edges | Upstream | Rewrite |
       |---|---:|---:|---:|
       | Heavy 100% | 294.4 M | 9 min 11 s, 14.8 GB | 2 min 4 s, 2.3 GB |
       | Bulk 50% | 353.8 M | 18 min 8 s, 17.8 GB | 3 min 26 s, 2.8 GB |
       | Bulk 100% | 1.417 G | not run (~70 GB) | 5 min 26 s, 11.1 GB |

       Bulk 100% gives 1,969,143 reads, 21,041 meta nodes and 13,041 meta
       edges after filtering.
   - Tests: `tests/test_componentsearch.py` (15; expected outputs from the
     upstream binary). Suite 365 passed.
   - The submodule now points to the fork `xvazquezc/ComponentSearch2`
     (`.gitmodules`), where the rewrite goes on top of upstream f705924 and
     the fork's CI workflow (`make`).

## Note on conda-forge vs bioconda AVX-512

Verified with `objdump -d | grep -c "zmm\|evex"` across environments:

| Build | zmm/evex count | Runs? |
|-------|---------------|-------|
| bioconda `4.3.7=hdbdd923_1` | 2771 | **SIGILL crash** |
| conda-forge `7.0.0=py314hd80a2e3_0` | 0 | OK |
| conda-forge `7.0.0=py312h8ced662_0` | 0 | OK |

conda-forge sortmerna >=5.0 contains no AVX-512 instructions.
bioconda `4.3.7=hdbdd923_1` does and crashes on AVX2-only hardware.
Always use the conda-forge channel for sortmerna >=5.0.

---

## Test setup

Conda env `matam_test` uses Python 3.14 + sortmerna 7.0.0 (conda-forge).
MATAM C++ binaries compile with `python build.py`; the current `ovgraphbuild`
validation binary was also rebuilt directly with CMake.

Reference DB prepared end-to-end with the current MATAM code:
```
/db/markermag/SILVA_138.2_SSURef_NR99_tax_silva.10_matam/
```

`-d` prefix to pass to `matam_assembly.py`:
```
/db/markermag/SILVA_138.2_SSURef_NR99_tax_silva.10_matam/SILVA_138.2_SSURef_NR99_tax_silva.10_NR95
```

Both `.complete.smr/idx` and `.clustered.smr/idx` exist and were built
successfully with SortMeRNA 7.0.0. SortMeRNA run state should remain under
`/tmp` because `/home` has limited free space.

Example test command:
```bash
MATAMDIR=/home/xabi/Desktop/github/matam
conda activate matam_test && export PATH="$MATAMDIR/bin:$PATH"
python $MATAMDIR/scripts/matam_assembly.py \
  -d /db/markermag/SILVA_138.2_SSURef_NR99_tax_silva.10_matam/SILVA_138.2_SSURef_NR99_tax_silva.10_NR95 \
  -i $MATAMDIR/examples/16sp_simulated_dataset/16sp.art_HS25_pe_100bp_50x.fq \
  --cpu 24 --max_memory 10000 \
  -o /tmp/matam_test_out/output \
  --workdir /tmp/matam_test_out/run \
  --perform_taxonomic_assignment -v
```
