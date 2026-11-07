# MATAM bug fixes and unit-test backlog

Tracks the bugs and unit tests identified during the test-coverage review
(2026-11-06). Update the status of each item as it is worked on.

**Status:** `[ ]` open · `[~]` in progress · `[x]` done (fixed + tested) · `[-]` won't fix

**Testing rule:** run all tests in the `testmatam` conda environment, e.g.
`conda run -n testmatam python -m pytest tests/ -q`. Installing test
dependencies into `testmatam` is allowed.

## Existing coverage (before this backlog)

`binary_utils`, `compute_abundance.abundance_calculation`, `krona`, `rdp`,
`components_assembly` pair extraction, SortMeRNA paired contract and
subsampling equivalence, `ovgraphbuild` threading/sharding, and one end-to-end
functional run (`tests/test_functional.py`).

---

## Tier 0: upstream design issue

### U1. [x] Read identity for paired reads with identical mate IDs
Everything from the SortMeRNA SAM onward identifies a read by its SAM QNAME
(the FASTQ ID up to the first whitespace). Modern Illumina (Casava 1.8+)
mates share that ID (`@ID 1:N:0` / `@ID 2:N:0`), so R1 and R2 get the same
QNAME.

**Verified 2026-11-06 with SortMeRNA 7.1.0** in `testmatam`, using MATAM's
step-1 options on simulated pairs from `ABA_X81660`:
- Casava IDs → **identical QNAME for both mates**, FLAG only `0`/`16`, no
  pairing bits (0x1/0x40/0x80), RNEXT `*`, PNEXT `0`. Nothing in the SAM
  distinguishes the mates except record order, and MATAM's
  `sort -k1,1V -k12,12Vr` destroys that order. FASTX and BLAST outputs
  carry the same shared ID.
- R2 with 12 mismatches (AS 240 vs R1 AS 300) → `filter_score_multialign
  -t 0.9`, in both geometric and straight mode, **drops R2 entirely**
  (240 < 270).
- The same reads renamed to `ID/1` / `ID/2` → SortMeRNA keeps the suffixes
  in the QNAME, and both mates survive the filter.

Steps affected (from reading the code, beyond the filter):
- `filter_score_multialign.py`: groups consecutive alignments by QNAME, so
  the two mates compete, and the weaker mate's alignments are dropped when
  they score below `threshold × best` (or the geometric chain).
- `ovgraphbuild`: a node is a run of consecutive SAM records with the same
  QNAME, so both mates collapse into one node and their alignments are
  treated as alternative placements of one read. See also #8.
- componentsearch → `read_metanode_component` and `compute_lca_from_tab`:
  keyed by QNAME (one row per merged pair).
- `components_assembly.extract_*_paired_reads_by_component`: an ID →
  component dict; both mates look up the same entry, so the
  "same-component pair" check never rejects anything.
- Abundance: the SortMeRNA BLAST output of reads vs scaffolds goes through
  `get_best_matches_from_blast` (grouped by read ID) and then
  `abundance_calculation` (set of read IDs). With shared IDs the mates
  compete for best hits and a pair counts once; with `/1`, `/2` it counts
  twice. The abundance scale therefore depends on the read naming.

**Decision (user, 2026-11-06):** pair-aware processing throughout, keyed on
(QNAME, mate). **No renaming of input reads.** Mate identity must be
recovered inside the pipeline. A proper end-to-end test with
realistic-named PE input is required. (The earlier "Option A: rename IDs at
input" is rejected.)

**Mate recovery: feasibility verified 2026-11-06** (SortMeRNA 7.1.0, the
5,825 ART example pairs converted to Casava naming, 16 references,
`--threads 4 --num_alignments 10`):
- SortMeRNA 7.1 has **no option** that sets SAM pairing flags. `--paired`,
  `--paired_in` and `--out2` only change FASTX routing.
- `--out2` writes the aligned R1/R2 FASTX in step: 5,825 / 5,825 records,
  0 ID mismatches line by line.
- Matching SAM SEQ (reverse-complemented when FLAG has 0x10) against the
  R1/R2 FASTX assigns **all 41,986 records**: 21,144 R1, 20,842 R2,
  0 ambiguous, 0 unmatched.
- Each QNAME's records are contiguous, with R1 records before R2.
- Approach: after mapping, MATAM tags each record with
  0x1 + 0x40/0x80 (plus mate fields), then every later step keys on
  (QNAME, mate).

Scripts used: `make_pairs.py`, `to_casava.py` and `check_mates.py` in the
session scratchpad.

With mate-keyed nodes, #2 stays latent (each mate is in one component).

**Implemented 2026-11-06** (plan: `~/.claude/plans/witty-hatching-corbato.md`).
Paired reads stay in separate R1/R2 files everywhere (user requirement):
- step-1 SortMeRNA uses `--paired_in --out2` (`<aligned>_fwd/_rev.fq.gz`);
  tagging, component extraction and abundance mapping all use the R1/R2 pair.
  No interleaved file is produced or read any more (the interleaved
  extractor and the tagger's interleaved mode were removed).
- MATAM publishes the selected pairs as `<out_dir>/marker_reads_R1.fq.gz` /
  `marker_reads_R2.fq.gz`, a stable `--filter_only` output contract.
- **MarkerMAG-nf** (`nextflow/bin/markermag_matam_filter.py`,
  `MarkerMAG/matam_16s.py`, `MarkerMAG/read_pairing.copy_fastq_pair`)
  consumes that pair directly and checks it stays in step, instead of
  globbing and splitting an interleaved file. Its downstream pairing was
  already flag-based (`read_pairing.sam_read_id_and_mate`).
  Tests: `nextflow/test/test_python314_and_paired_reads.py::FilterOutputTests`.
  The real wrapper was run against real MATAM on CASAVA pairs: OK.

**C++ verification (ovgraphbuild + `seqanObjectsFormating.h`), 2026-11-06.**
Builds of the old (`HEAD`) and new sources were made in the scratchpad:
- Clean build: 0 compiler warnings in `ovgraphbuild/src` (old and new); all
  113 warnings come from the seqan headers, identically in both.
- ASan + UBSan (`-fsanitize=address,undefined`): all ovgraphbuild tests and
  an 8-thread run on the real mate-tagged e2e SAM give no memory errors and
  no UB. LeakSanitizer reports 33 bytes per reference name: seqan's skiplist
  `Map` (`refNameToIndexMap`) never frees its keys. Identical in the old
  build, so pre-existing, bounded, and not addressed.
- Old vs new on the real SAM: unflagged legacy-named input gives
  byte-identical `edges.csv` and stdout statistics, with nodes identical
  apart from the new `Mate` column. New on CASAVA + mate flags vs old on
  legacy names gives byte-identical edges and 0/2,914 node mismatches. The
  old binary on CASAVA + flags gives 1,457 nodes (the bug).
- `--asqg` (VT/ED counts = nodes/edges) and `--pair-shard-*` (the union of
  3 shards = the full edge set) work on mate-tagged input.
- Permanent test:
  `test_ovgraphbuild_mates.py::test_flagged_shared_names_build_the_same_graph_as_unique_legacy_names`.
  The new ovgraphbuild tests fail against the old binary, as intended.
- **Memory/CPU evaluation (2026-11-06).** Synthetic MATAM-like SAMs
  (33-char Illumina names, 150 bp reads with qualities, 3 alignments per
  read; 20k / 100k / 400k reads), old vs new on identical input, with
  `/usr/bin/time -v`:
  - The first version checked contiguity with an `unordered_set<string>`
    filled during loading. It cost **+2.0% peak RSS** at load (≈ 69 B/read:
    +26 MB at 400k reads) and **+8.5% CPU** in the unchanged pairwise loop
    (single-thread medians 32.07 s vs 29.55 s, 3 reps; ~+11% in
    16-thread runs at 100k reads). The no-check variant matched old, so the
    per-read allocations interleaved with the records slowed the loop
    (worse locality).
  - **Replaced** by `exitIfReadsAreSplit`: after loading, sort one index
    vector by (QNAME, mate) and check neighbours. Result: CPU 29.32 s
    (= old, within noise) and load peak 1,310.4 MB vs old 1,310.6 MB at 400k
    reads. Edges are byte-identical in every variant; ASan/UBSan clean; all
    tests pass.
- `scripts/sam_mates.py` (new): sets 0x1 + 0x40/0x80 (+0x8 when the mate has
  no alignment) by SEQ match, streaming the SAM and the read pairs in
  lockstep; inconsistencies are hard errors. Used right after mapping and in
  the abundance step.
- `filter_score_multialign.py`: groups by (QNAME, mate) and ranks by `AS:i`
  internally (no external score sort).
- `ovgraphbuild`: node = (QNAME, mate); nodes CSV `Id;Label;Mate;Specie`
  (ComponentSearch copies it through unchanged); non-contiguous reads are
  rejected. Also fixes #8 and the unmapped-then-mapped grouping bug (#12).
- `scripts/pe_tables.py` (new): read-component and complete-taxonomy tables
  keyed by (read, mate), replacing the shell `cut`/`join` pipelines.
  `compute_lca_from_tab.py -g 1,2` gives one vote per mate.
- `components_assembly.py`: pair extraction looks up (ID, 1)/(ID, 2) and
  rejects out-of-sync mates. SGA appends `/1`, `/2` itself for identical
  names (verified).
- Abundance: mapping writes SAM (not BLAST), tagged by `sam_mates`, best
  matches per (read, mate) at ≥ 0.99 × best; each mate counts as one read
  (user decision).
- Input read count reported in mates (#11).

**Tests:** `test_sam_mates.py`, `test_filter_score_multialign.py`,
`test_ovgraphbuild_mates.py`, `test_pe_tables_and_lca.py`,
`test_sga_paired_contract.py` (mate-keyed extraction),
`test_compute_abundance.py` (per-mate), `test_sortmerna_pe_contract.py`
(SortMeRNA assumptions, `--out2`), and **`test_functional_pe.py`**, the end-to-end test.
It builds a 16-reference DB, runs MATAM on CASAVA (gz), SRA-style and legacy
`/1`, `/2` copies of the same pairs, and requires identical nodes, edges,
components, SGA pairs, scaffolds and abundances (also with
`--coverage_threshold`). It also checks the `--filter_only` R1/R2 output.
Before the change it failed: 1,457 vs 2,914 nodes.

## Tier 1: confirmed bugs (each reproduced with a probe)

### 1. [x] Component LCA carries over to components without an LCA entry
- **Where:** `scripts/components_assembly.py`, `concat_components_fasta_with_lca`
- **Bug:** `component_lca = 'NULL'` was set once before the loop, so a
  component missing from the LCA table inherited the previous component's
  lineage.
- **Expected:** a missing component is written as `lca=NULL`; contigs are
  renumbered 1..N across components; empty sequences are skipped.
- **Impact:** wrong label only in the intermediate contigs FASTA in
  `<out_dir>/workdir/`; nothing parses `lca=` downstream, so final results
  are unaffected.
- **Fix:** look up `component_lca_dict.get(component_id, 'NULL')` inside the loop.
- **Test:** `tests/test_sga_paired_contract.py::test_component_without_lca_is_labelled_null`
- **Status:** fixed and tested 2026-11-06; committed in `bc37861`.

### 2. [x] LCA grouping leaks rows across factors
- **Context:** the component LCA only feeds the `lca=` labels in the
  intermediate contigs FASTA (the graph-stats consumer is commented out).
  Final taxonomy is RDP on the scaffolds, independent of this step.
- **Where:** `scripts/compute_lca_from_tab.py`, `read_tab_file_handle_sorted`
- **Bug:** when the same group ID (read) ends one factor (component) and
  starts the next, its rows are dropped from the first factor and added to
  the second one's group.
- **Repro:** rows `r1 A`, `r2 A`, `r2 B`, `r3 B` (factor = column 3) give
  A = `[[r1]]`, B = `[[r2(A), r2(B)], [r3]]`. Expected: A = `[[r1], [r2]]`,
  B = `[[r2], [r3]]`.
- **CLI effect** (`-t 4 -f 3 -g 1`; inputs and outputs verified 2026-11-06):
  - the second factor is **written under the first factor's ID**, because
    the ID is read from the factor's first row (the leaked A row). Output is
    `A Bac;Firm` / `A Bac` (`-m 1.0`) or `A Bac;Firm` / `A Bac;Prot`
    (`-m 0.51`), instead of `A Bac;Firm` / `B Bac;Prot`.
    `extract_lca_by_component` keeps the last value, so A gets the wrong LCA
    and B gets none (it becomes `NULL` after fix #1);
  - if a factor consists only of the shared read, an empty list is yielded
    and the CLI crashes with `IndexError`. `matam_assembly.py` runs this step
    with `logged_call` (not `check_call`), so the run carries on with a
    truncated LCA table and only exits non-zero at the end.
- **Impact:** latent. It needs the same read ID listed in two consecutive
  components, as the last read of one and the first of the next. With unique
  mate IDs (the example data uses `/1`, `/2`) each read is in one component.
  Open question: mates with identical IDs (Casava 1.8 style
  `@ID 1:N:0` / `@ID 2:N:0`) may make it reachable. Not verified how
  SortMeRNA/ovgraphbuild name nodes in that case.
- **Fix idea:** check for a factor change first. On a change, flush the
  current group into the factor list, yield it, and reset both lists;
  otherwise flush on a group change only. Never yield an empty factor list
  (empty input should yield nothing).
- **Fixed 2026-11-06** together with multi-column `-g`.
  **Tests:** `test_pe_tables_and_lca.py::test_group_continuing_into_the_next_factor_stays_with_its_factor`,
  `::test_cli_labels_each_component_with_its_own_reads`,
  `::test_cli_handles_a_component_made_of_one_shared_read`, `::test_empty_input_yields_no_factor`.

### 3. [x] `is_phred33` rejects valid input and crashes on empty files
- **Where:** `scripts/fastq_utils.py`, `is_phred33` (and `read_fastq_file_handle`)
- **Bugs:**
  - empty FASTQ raises `ValueError: min() iterable argument is empty`,
    because `read_fastq_file_handle` yields `('', '', '')` on empty input;
  - uniform high quality (`IIII`) returns False (requires min < 59);
  - any quality character ≥ `K` (Q42+) returns False.
- **Impact:** valid Phred+33 input (simulated reads, HiFi, some newer
  Illumina) is rejected. Called from `matam_assembly.py` input validation.
- **Fix idea:** decide based on the absence of characters below `;`/`@`
  (Phred+64 marker) rather than requiring a low minimum; handle empty input.
- **Fixed 2026-11-06:** a character below `;` (59) proves Phred+33; otherwise
  the file is accepted while every character is within `!`..`J`, and rejected
  only with the Phred+64 profile (nothing below `;`, something above `J`) or a
  character outside `!`..`~`. `read_fastq_file_handle` yields nothing for an
  empty file; an empty file is not rejected. Still only checked with
  `--verbose` (Tier 3 item). **Tests:** `test_input_parsers.py` (FASTQ section).

### 4. [x] `parse_cigar` crashes on repeated adjacent operations
- **Where:** `scripts/compute_contigs_compatibility.py`, `SamAlignment.parse_cigar`
  (live path, used with `--contigs_binning`). Same copy in
  `sample_sam_by_coverage.py`, `compute_assembly_stats.py`,
  `filter_sam_by_pid.py`, `compute_ref_coverage_histogram.py`.
- **Bug:** `cigar_tab[-1][1] += count` assigns into a tuple, so `5M5M` raises
  `TypeError`.
- **Impact:** low likelihood (aligners rarely emit adjacent repeated ops).
- **Fix idea:** replace the last entry with a new tuple, or share one
  `parse_cigar` helper.
- **Fixed 2026-11-06** in the 4 merging copies (the merged entry is rebuilt as
  a new tuple); `compute_assembly_stats.py` doesn't merge and was fine. Sharing
  one helper stays a Tier 3 item. **Tests:** `test_input_parsers.py`
  (every copy, parametrised).

### 5. [x] `scaffold_contigs` raises an undefined exception
- **Where:** `scripts/scaffold_contigs.py`, `iter_read_bases`
- **Bug:** `raise ParsingError(...)` but `ParsingError` is never defined, so
  any unexpected mpileup character gives a `NameError`.
- **Note:** `#` (reverse-strand deletion) only appears with samtools
  `--reverse-del`, which MATAM doesn't pass, so default output is fine.
- **Fix idea:** raise `ValueError` with the offending character and position.
- **Fixed 2026-11-06:** `ValueError` naming the character, its index and the
  read-bases field; the stray `print` to stdout (the scaffolds output when
  `-o -`) was removed. **Tests:** `test_input_parsers.py` (mpileup section).

### 6. [x] `runner.logged_call` masks failures in shell pipelines
- **Where:** `scripts/runner.py`, `logged_call` / `logged_check_call`
- **Bug:** commands run via `/bin/sh` without `pipefail`, so only the last
  command's exit code counts. `false | true` and
  `gzip -cd /nonexistent.gz | wc -l` both return 0.
- **Impact:** many MATAM steps are pipelines (e.g. LCA labelling in
  `matam_assembly.py`, `get_best_matches` in `compute_abundance.py`), so a
  failed step can pass silently and leave empty outputs. Highest real-run impact.
- **Related:** in verbose mode the read loop stops once the process exits,
  so the tail of the output can be lost.
- **Fix idea:** run with `executable='/bin/bash'` and prefix `set -o pipefail;`.
- **Partial, 2026-11-06:** opt-in `pipefail=True` added to `logged_call` /
  `logged_check_call` and used for the pipelines touched by U1 (alignment
  filtering, coverage sampling, LCA). Not yet the default: existing pipelines
  use `grep -v`, which exits 1 when nothing matches, so each one needs
  checking before switching. **Tests:** `test_runner_and_coverage.py`.
- **Done 2026-11-06:** all 30 `logged_call`/`logged_check_call` call sites
  audited: none has an inner command that may legitimately fail (the
  `grep -v "^@"` went with U1), so **pipefail is now the default**
  (`pipefail=False` to opt out). Verbose mode reads the output until EOF with
  an incremental UTF-8 decoder, so the tail is no longer lost; the unclosed
  `/dev/null` handle was replaced by `subprocess.DEVNULL`.
  **Tests:** `test_runner_and_coverage.py` (default pipefail, `check_call`
  exit, verbose tail with a deliberately slow log handler).
- **Not covered:** `matam_db_preprocessing.py` runs its own pipelines with
  `subprocess.call(shell=True)` and sums the return codes, so failures inside
  its pipelines are still masked there.

### 7. [x] `sample_sam_by_coverage.py` uses removed `np.int`
- **Where:** `scripts/sample_sam_by_coverage.py`, `compute_ref_coverage`
- **Bug:** `np.int` was removed in NumPy 1.24, so the script crashes. The
  env pins Python 3.14, which only has NumPy 2.x builds.
- **Impact:** live when `--coverage_threshold N` (N > 0) is used. It runs as
  `sample_sam_cov_bin` in the alignment-filtering step of `matam_assembly.py`.
  The default is `0` (off). *Corrected 2026-11-06: an earlier note said the
  script was never run; that grep matched only the filename, not the variable.*
- **Fix idea:** `dtype=int`. Add a unit test for `compute_ref_coverage`.
  NumPy is not yet installed in `testmatam`.
- **Note:** `tests/test_functional.py` passes `--coverage_threshold 2000`, so
  it would hit this crash. It's skipped locally only because the SILVA DB is
  absent.
- **Fixed 2026-11-06** (`dtype=int`). **Tests:**
  `test_runner_and_coverage.py::test_reference_coverage_counts_each_aligned_reference_position`,
  and `test_functional_pe.py::test_coverage_subsampling_keeps_reads_whole_and_naming_independent` (end to end).

### 9. [x] Reference headers containing tabs corrupt the SortMeRNA SAM (minor)
- SortMeRNA uses the reference header up to the first space as RNAME, so
  tabs end up inside it and shift every SAM column. Seen with
  `examples/.../16sp.fasta` (`>ABA_X81660<TAB>1459bp<TAB>...`).
- `matam_db_preprocessing.py` doesn't reduce headers to their first token.
  SILVA uses spaces, so it's unaffected; custom DBs could break.
- **Fix idea:** sanitise headers (whitespace → space) during DB preprocessing.
- **Real impact (2026-11-06):** a DB built from the raw `16sp.fasta` makes
  MATAM crash at scaffolding (`invalid literal for int(): 'Klebsiella'`, the
  shifted POS column). Not cosmetic.
- **Fixed 2026-11-06:** the cleaning step turns whitespace runs in headers into
  single spaces and trims trailing ones, so RNAME = reference ID and the
  description is kept. Full MATAM run on the raw-header DB: 8 scaffolds, no
  error. **Tests:** `test_db_preprocessing.py` (headers, IDs agree across the
  complete/taxo/clustered files, SortMeRNA SAM columns).

### 8. [x] Use-after-move in ovgraphbuild read grouping (never triggered)
- **Where:** `ovgraphbuild/src/ovgraphbuild.cpp`, `readSamFile` (~line 211-214,
  unchanged since the 2016 initial commit)
- **Suspicion:** `readBamRecordBuffer.push_back(std::move(bamRecord))` is
  followed by `oldQName = bamRecord.qName`. If SeqAn's string move empties
  `qName`, every alignment starts a new node, so nodes are alignments rather
  than reads and the mapped-read stats are inflated.
- **Verify:** not possible locally (binary not built, `ovgraphbuild/lib/seqan`
  submodule not checked out, ovgraphbuild tests skip). Test: a SAM with one
  read aligned twice should give 1 node in `.nodes.csv`.
- **Fix idea:** assign `oldQName` before the move.
- **Verified 2026-11-06** after building: 41,986 SAM records gave 11,650
  nodes (one per mate), so SeqAn's move leaves `qName` intact and the bug
  never triggered. Fixed anyway (the key is taken before the move).
  **Test:** `test_ovgraphbuild_mates.py::test_all_alignments_of_one_read_form_one_node`.

### 10. [x] Coverage subsampling ignored `--seed` (non-reproducible runs)
- `sample_sam_by_coverage.py` called `random.seed()` (system randomness), and
  ties at one position depended on `sort`'s whole-line fallback, so read
  names. Found when the coverage-sampled e2e runs differed by one node.
- **Fixed 2026-11-06:** the sampler takes `--seed`, `matam_assembly.py` passes
  `--seed` to it, and the position sort is stable (`sort -s`).
- Note: ComponentSearch also seeds from the clock (`time(0)`) unless
  `--seed` is given. The e2e test passes `--seed 42`.

### 11. [x] "Identified as marker" mixed pairs and mates (200%)
- The input count was pairs (R1 lines / 4) and the selection was mates.
- **Fixed 2026-11-06:** the input is counted as R1 + R2 mates.
  **Test:** `test_functional_pe.py::test_read_counts_are_in_mates`.

### 12. [x] ovgraphbuild merged a mapped record into the previous read after an unmapped one
- In `readSamFile`, a mapped record that followed an unmapped record of the
  same read didn't start a new buffer, so it joined the previous read's node.
- **Fixed 2026-11-06.** **Test:**
  `test_ovgraphbuild_mates.py::test_mapped_record_after_an_unmapped_one_starts_its_own_node`.

### F1. [x] Follow-up: over half of the pairs are discarded before SGA
- With correct mate identity, mates often fall into different graph
  components. The "same-component pair" rule then discards the pair:
  3,250 of 5,825 pairs (56%) on the full example, and 804 of 1,457 on the e2e
  subsample. This was already the case with `/1`, `/2` naming before U1. SGA
  therefore assembles each component from under half of its reads.
- **Origin (checked 2026-11-06):** the same-component rule was added in
  `82c990d` ("move to only support separated pe-reads"). Upstream MATAM
  (fork point `e7cd9b2`, bonsai-team, 2025-03-28) had no pair rule:
  `extract_reads_by_component` put every read of a component into that
  component's file and SGA always assembled it as single-end reads (the
  wrapper's `--paired_end` branch was a stub that was never passed). The fork
  kept that until `82c990d`. `d0a4a6f` has the rule with
  the same outcome as now for `/1`, `/2` names (distinct node keys, identical
  graph). With CASAVA names the mates collapsed, so the check never fired.
  Breakdown on the CASAVA run (every 4th pair, `--seed 42`): of 1,457 pairs,
  653 same component, 706 different components, 94 with one mate outside every
  component, 4 with neither. Now 1,306 mates are assembled; the original rule
  would have assembled all 2,812 mates that are in a component.
- **Upstream vs fork assembly (2026-11-06)**, full example (5,825 pairs),
  same references, defaults, seeds 1–3. Upstream = bioconda `matam` 1.6.1
  (scripts identical to `e7cd9b2`, SortMeRNA 2.1b) in a scratchpad env;
  fork = working tree. Evaluated with `evaluate_assembly.py` vs the 16 refs.
  Means: reads given to SGA 11,377 vs 6,193; contigs 468 vs 271 (63.2 vs
  39.7 kb); ref coverage by contigs 97.9% vs 94.6%; final scaffolds 16 vs
  16.7, 23.6 vs 19.5 kb; **ref coverage by final assembly 97.9% vs 80.9%**;
  refs ≥ 90% covered 16/16 vs 9.7/16. The fork misses *S. flexneri* in every
  seed (contigs cover it at 94–95%, but no scaffold survives). Error rate
  ~0% for both. Confounded by tool versions (SortMeRNA 2.1b vs 7.1, etc.);
  an ablation with the upstream assembly rule in the fork would isolate F1.
  Scripts: `scratchpad/compare/run_compare.sh`, `summarise.py`.
- **Why (2026-11-07):** components are short (median span 204 bp on the
  reference, IQR 180–280) and barely longer than a fragment (median 149 bp
  here), so 41% of pairs straddle a boundary even though their mates overlap.
  Real libraries (300–550 bp inserts) would lose more.
- **Rejected: merging mate-linked components.** 38/98 components are
  multi-species (conserved regions); merging on ≥1 linking pair gives 3 groups
  (largest 6,407 reads, 9 species), and even ≥20 pairs leaves a 4,895-read,
  8-species group. Chimera risk.
- **Fixed 2026-11-07:** `_split_pairs_by_component` sends each pair to every
  component holding one of its mates (spanning pairs go to both; a pair with
  one mate in a component goes there; components are never merged). SGA only
  uses pairing in `preprocess --pe-mode 1`, so the partner mate is extra
  same-molecule sequence across the boundary. Log line:
  `Read pairs by component: W within one component, S spanning two (given to
  both), O with one mate in a component, U in no component (not assembled)`.
  Means of 3 seeds, full example (upstream / old fork rule / new routing):
  reads to SGA 11,377 / 6,193 / 16,589; components with contigs 104 / 79 /
  102; ref coverage by contigs 97.9% / 94.6% / 98.2%; **final ref coverage
  97.9% / 80.9% / 98.3%**; refs ≥ 90% covered 16 / 9.7 / 16 (S. flexneri
  recovered at 97–99%); final N50 1,492 / 1,402 / 1,499; error ~0% for all;
  wall time 20.3 / 8.5 / 11.6 s. **Tests:** `test_sga_paired_contract.py`
  (spanning pair goes to both, one-mate pairs, routing log line, CASAVA);
  `test_functional_pe.py` compares the routing counts across namings.

---

## Tier 2: characterization tests (no bugs found, high impact on results)

- [ ] **`compute_contigs_compatibility`** (`--contigs_binning`):
  `return_compatibility_status` overlap boundaries, the early exit in
  `compute_contigs_compatibility_matrix`, the secondary-incompatibility rule
  in `compute_bin_list`, and an empty input. Edge case to decide on: a
  1-base overlap is reported as "no overlap" (inclusive end coordinate).
- [ ] **`compute_lca`**: strict LCA (`min_proportion=1.0`) vs `--quorum`
  mode, tie-breaking, and the `Root` fallback.
- [ ] **`scaffold_contigs.find_called_base`**: majority base, `*` deletions,
  insertion threshold (≥ coverage/2), and ties (settled alphabetically, so
  `T` beats `A`).
- [ ] **Abundance round-trip**: `complete_fasta_with_abundance` →
  `get_abundance_from_fasta`; a missing ID becomes 0.
- [ ] **`components_assembly.assemble_component`** with a fake assembler
  (monkeypatch `AssemblerFactory.get`): re-runs with read correction only
  when coverage exceeds the threshold, returns None on empty output, and
  turns `SystemExit` into `Exception`. Also `estimate_coverage` /
  `nucleotidic_number` on `.gz` input and zero-length contigs.
- [ ] **`matam_assembly` helpers**: `FastaStats` on empty input (returns -1),
  `fastq_line_count` on `.gz`, the ambiguous-match case in
  `find_sortmerna_output`, and `format_seq`.
- [ ] **`matam_assembly.parse_arguments`**: range checks, `--debug` turning
  on verbose and keep-tmp, the coverage threshold being cleared when read
  correction isn't `auto`, and `--workdir` defaulting to `--out_dir`. Needs
  `parse_arguments(argv=None)` or monkeypatching `sys.argv`. Also fix the
  quorum error message (says `[0.51,1]`, checks `[0,1]`).
- [ ] **`AssemblerFactory.get`** with an unknown name (`KeyError`), and SGA
  `build_command_line` flags for each `read_correction` mode.

---

## Tier 3: structural changes that make testing cheaper

- [ ] **Merge the 19 copies of `read_fasta_file_handle`** (plus repeated
  `format_seq`, `parse_cigar`, `read_tab_file_handle_sorted`) into shared
  helpers (`fasta_utils.py`) and test them once. The copies behave
  differently on edge cases (some yield `('', '')` on empty input;
  `matam_assembly.read_tab_file_handle_sorted` yields `[[]]`).
- [ ] **Extract FASTQ input validation into a function**: the line count
  being a multiple of 4, matching read counts in the pair, and the Phred+33
  check run only with `--verbose` in `matam_assembly.main`. Moving them into
  `validate_paired_fastq(fwd, rev)` would make them testable, and they should
  probably run every time.

---

## Progress log

| Date | Item | Change |
|------|------|--------|
| 2026-11-06 | — | Backlog created from the test-coverage review |
| 2026-11-06 | #1 | Fixed `concat_components_fasta_with_lca`; added regression test (uncommitted) |
| 2026-11-06 | #2 | Ran the CLI: factor ID mislabelled and `IndexError` on single-read factors; item updated |
| 2026-11-06 | U1, #7, #8 | Traced read identity upstream (U1); #7 corrected to live with `--coverage_threshold`; added suspected ovgraphbuild bug #8 |
| 2026-11-06 | U1 | Installed SortMeRNA 7.1.0 in `testmatam`; verified shared QNAME/no mate flags, weaker mate dropped by the filter, and `/1`, `/2` suffixes preserved |
| 2026-11-06 | U1, #7, #9 | Decision: pair-aware throughout, no input renaming. Verified mate recovery via `--out2` + SEQ matching (100% of 41,986 records). Functional test hits #7; added #9 |
| 2026-11-06 | U1, #2, #6, #7, #8, #10–#12 | Pair-aware pipeline implemented (no input renaming, mate = counting unit); PE end-to-end naming-invariance test added; suite 84 passed, 6 skipped (SILVA functional test) |
| 2026-11-06 | U1 | Paired format end to end: MATAM `--out2`, R1/R2 through tagging/assembly/abundance, `marker_reads_R1/R2` published; MarkerMAG-nf hand-off consumes the pair. MATAM suite 87 passed, 6 skipped; MarkerMAG-nf 13 passed |
| 2026-11-06 | U1, #8 | C++ verification: clean build without new warnings, ASan/UBSan clean (seqan Map leak pre-existing), old vs new differential byte-identical, ASQG/shards OK; naming-equivalence graph test added |
| 2026-11-06 | U1 | Memory/CPU eval of ovgraphbuild: hash-set contiguity check cost +2% RSS, +8.5% CPU; replaced with a post-load sort check (no measurable cost) |
| 2026-11-06 | #3, #4, #5, #6, #9 | Fixed; `test_input_parsers.py`, `test_db_preprocessing.py`, runner tests (17 new checks fail on the previous code). Pipefail default after auditing all call sites. Suite 115 passed, 6 skipped |
| 2026-11-07 | F1 | Cause: components ≈ fragment length. Merging mate-linked components rejected (chains species). Pairs now go to every component holding a mate: final ref coverage 98.3% vs upstream 97.9% (old rule 80.9%). Suite 116 passed, 6 skipped |
| 2026-11-07 | U1, F1 | Documented in `eval_and_release/paired_end_assembly/README.md`, with portable `run_compare.sh` / `summarise.py`; README links to it |
