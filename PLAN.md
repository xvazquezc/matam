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

Commits through `c5f5a06` are on the fork. `d371b1d` is committed locally but
has not yet been pushed.

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

The README warning for the SortMeRNA 6 BBHash index-format break is present in
the working tree but is not committed yet. Users upgrading from SortMeRNA 4/5
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
  no alignment of one read overlaps one of the other on a common reference.
  They never give an edge. The all-pairs loop still tests every pair (up to
  10x10 alignment comparisons each).

## Remaining work

1. Candidate-pair pruning in `ovgraphbuild` (likely the biggest gain). Per
   read, compare only reads that share a reference (better: overlap on one),
   in ascending order. Edges then stay byte-identical, and the skipped pairs
   are added to the "incompatible" and "neither" counts so the statistics
   don't change either. This works with the existing threads and pair shards.
   Measure it on the synthetic sets above, then on the production SAM.
2. Benchmark the slow production filtered SAM with its clustered reference
   (`ovgraphbuild` with `--threads 24`, plus 1 thread on a few
   `--pair-shard-*` slices for the serial rate). Record the read count,
   alignments per read, compatible and "neither" fractions, and pairs/s.
3. Add MarkerMAG-nf scatter/merge orchestration only if a graph is still over
   the walltime at 24 threads after 1 and 2.

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
