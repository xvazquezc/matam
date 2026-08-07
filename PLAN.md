# MATAM modernisation plan
# Repo: xvazquezc/matam  Branch: master

---

## Current repository state

- Upstream base: `v1.6.2` (`e7cd9b2`)
- Fork remote: `git@github.com:xvazquezc/matam.git`
- `origin/master`: `c5f5a06`
- Local `master`: `d371b1d` (one commit ahead of `origin/master`)
- Uncommitted documentation: `README.md` and this `PLAN.md`

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

## Remaining work

1. Commit the README index-compatibility and work-directory documentation.
2. Push `d371b1d` and the documentation commit to `origin/master`.
3. Benchmark the slow production filtered SAM with its matching clustered
   reference; the current 67.9-million-pair dataset validates correctness and
   scaling but is not the pathological workload.
4. Add MarkerMAG-nf scatter/merge orchestration only if external graph shards
   are still needed after production benchmarking.

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
