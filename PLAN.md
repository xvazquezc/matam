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
2. Benchmark the slow production filtered SAM with its clustered reference
   (`ovgraphbuild` with `--threads 24`, plus 1 thread on a few
   `--pair-shard-*` slices for the serial rate). Record the read count,
   alignments per read, compatible and "neither" fractions, and pairs/s.
3. Add MarkerMAG-nf scatter/merge orchestration only if a graph is still over
   the walltime at 24 threads after 1 and 2.
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
