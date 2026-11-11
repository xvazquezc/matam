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

### 13. [x] Rerun or `--resume_from` into an existing output directory fails (SortMeRNA 7 run state)
- SortMeRNA 7 keeps run state in `<workdir>/kvdb` and aborts with
  `kvdb resume aborted: inputs mismatch` when the next run's inputs differ.
  MATAM reuses `sortmerna_reads_mapping/` and `sortmerna_scaffolding/`, so
  resuming from any step up to `scaffolding` fails at scaffolding (the contigs
  are rebuilt). A rerun with different reads into the same directory fails at
  step 1. Resuming from `abundance_calculation` or `taxonomic_assignment` works
  (abundance uses a fresh temp dir).
- Found 2026-11-10 by a smoke test. Dates from the SortMeRNA v5+ migration
  (`c5f5a06`), not the paired-end rework.
- **Fix idea:** remove `<workdir>/kvdb` before each SortMeRNA mapping (the
  state isn't reusable across runs); add a resume test.
- **Fixed 2026-11-10** (user's spec: check each step for completeness; remove
  incomplete output before continuing). SortMeRNA 7 has no force/overwrite
  option (only `--kvdb` and checkpointing via `--flush-delay`).
  New `scripts/sortmerna_step.py`, used for the read mapping and the
  scaffolding mapping:
  - after a successful run, `<aligned>.matam_complete.json` records the
    command, the inputs (size, mtime) and the required outputs (size); run
    state (`kvdb`, `readb`) is then removed;
  - a matching record means the step is complete and its output is reused;
    otherwise (no record, missing/truncated output, changed inputs or command)
    the step's SortMeRNA outputs and run state are removed and it reruns;
  - `--threads` and `-v` are not part of the record; paths inside the run are
    stored relative, so a moved or copied run directory stays resumable;
  - when `--resume_from` skips the read mapping, its output must be complete,
    or MATAM stops and asks to resume from `reads_mapping`; output from an
    older MATAM without a record is used with a warning.
  Verified: resume from every step (on runs from the old and new code) gives a
  byte-identical final assembly; in place, resuming from `reads_mapping`
  reuses the mapping, from `scaffolding` reuses the contig mapping; a
  SortMeRNA run killed mid-mapping leaves `kvdb` behind and the rerun clears
  it and completes; a rerun into the same directory with other reads (same
  file names) redoes both mappings.
  **Tests:** `test_sortmerna_step.py` (13, fake SortMeRNA that aborts on stale
  `kvdb` like the real one), `test_functional_pe.py` resume tests (5; all fail
  on the previous code, 2 with `kvdb resume aborted`).

### 14. [x] `--contigs_binning` fails with current samtools (fixed, untested; no test planned)
- `compute_contigs_compatibility.py` writes placeholder `@SQ ... LN:9999999999`
  headers (upstream code); current samtools rejects lengths above the BAM limit
  (`output format is not compatible with very large references`), so
  `samtools sort` fails.
- Found 2026-11-10 by a smoke test. Pre-existing upstream code exposed by
  modern samtools.
- **Fixed 2026-11-11:** the placeholder headers are written as `LN:2147483647`
  (the BAM maximum). **Untested:** no unit or end-to-end test was added, and
  none is planned (MarkerMAG-nf doesn't use the option).
- **Evaluated for MarkerMAG-nf, 2026-11-11** (before the fix). With
  `LN:2147483647` in a scratch copy, all 21 binning runs complete. Setup:
  SILVA 144 subset DB (10 families, 51,483 refs, 95% clustering); simulated
  community (ART HS25 2x150): 6 controls in the DB, 8 strain pairs (DB
  sequence + copy with 0.3/0.6/1.0/1.5% random substitutions, 1:1 and 3:1),
  6 novel pairs (95-98.4% apart, removed with their >=98.5% relatives). Run
  through MarkerMAG-nf's own `markermag_matam_assemble_pct.py` at 1-100%,
  then its post-processing (vsearch 0.99, longest per cluster, >=1200 bp).
  3 replicates (read and MATAM seeds), default / binning:
  - MarkerMAG-nf doesn't pass the option. If it did, its wrapper reads
    `workdir/scaffolds.NR.min_500bp.fa`, which binning renames
    (`scaffolds.contigs_binning.NR...`), so every depth comes back empty.
  - Final 16S 16-18 / 18-25; exact (>=99.9% to a truth) 10-13 / 5-8.
  - Strain pairs >=1% apart with both members recovered (>=99%): 0-1/4 /
    3-4/4. Strain truths exact: 3-6/16 / 1/16. Mosaic sequences (alleles of
    both siblings) 8-10 / 12-17.
  - Novel truths recovered 3-4/12 / 1-3/12: binning splits the mosaics into
    clean pieces of ~700 bp, below the 1200 bp linking minimum.
  - High-coverage control (900x): binning turns error contigs into their own
    scaffold, and longest-per-cluster can pick it (99.4% instead of 100%).
  - Default mode also makes mosaics whenever related organisms share a
    reference, including a 95%-apart novel pair (31 vs 5 alleles).
  Verdict: a trade-off, not an improvement; keep it off in MarkerMAG-nf.
  Scripts: `scratchpad/binning/` (session 5b25223e).

### 15. [x] SortMeRNA reference index: reuse unchecked, rebuild unsafe, indexing not optional
- **Does mapping reuse the prebuilt index?** Yes. With `--idx-dir <db>.<kind>.smr/idx`
  SortMeRNA 7.1 logs `Found 4 non-empty index files. Skipping indexing.` But it
  checks nothing else:
  - the index files are named after a hash of the fasta's *file name* only (not
    its path, content, `-L`, `-m` or `--max_pos`);
  - a fasta rewritten under the same name is used with the old index: same size
    gives silently wrong seeds; a shorter fasta segfaulted (rc 139);
  - an index built with another `-L`/`-m` is used as is. `-m` changes results:
    a 600-ref index split into 9 parts (`-m 1`) gave 9,883 vs 5,801 blast lines
    with `--num_alignments 10` (alignments are selected per part);
  - index files are created empty at the start of a build and filled at the end,
    so an index left by a killed build, or one still being written by another
    run, passes the "non-empty" test.
- **Without a prebuilt index** (before this fix): if `<db>.<kind>.smr` is
  missing, mapping fails (`filesystem error: cannot create directory`: SortMeRNA
  creates `idx/` but not its parent). If `.smr` exists, SortMeRNA builds in the
  DB directory with no lock: 2 builders plus 1 mapping run started together
  each built the full index (3x the work, 15 min each on 30k SILVA seqs).
- `matam_db_preprocessing.py --max_memory` was ignored (no `-m` since the port).
- `--max_pos`: SortMeRNA 7.1's help says the default is 1000, but the effective
  default is 10000, as in 2.1b: an index built without `--max_pos` is
  byte-identical to `--max_pos 10000` and differs from `--max_pos 1000`. The
  port didn't change the index. A real 1000 cap (30k SILVA seqs, run side by
  side): 201 vs 267 s, 0.96 vs 1.25 GB RSS, 515 vs 579 MB. Best hit unchanged
  for 94.9% of 10k reads, better with 10000 for 3.4%, better with 1000 for
  1.6%. Kept at 10000 (same as 2.1b and the existing indexes).
- **Cost on a full SILVA SSU NR99 fasta** (`SILVA_144_SSURef_NR99_tax_silva_trunc.fasta`,
  U→T, 905,628 seqs, 1.37 Gbp; clustered with vsearch at 95%: 288,470
  centroids, 446 Mbp, 1 h 41 min on 8 threads, 5.7 GB). `sortmerna --task 5
  -m 10000` is single-threaded. Random subsets, with `/usr/bin/time -v`, on a
  quiet machine: 30k complete 233-267 s, 1.25 GB, 579 MB; 30k centroids 238 s,
  1.6 GB, 679 MB; 100k complete 846 s, 2.94 GB, 1.64 GB on disk. Time grows
  slightly faster than linearly (30k to 100k: ×3.6 for ×3.3 seqs). The same
  runs under load were 2-4.5× slower.
  Extrapolated (exponent 1.07): complete index ≈ 2-2.5 h and ≈ 15 GB on disk;
  clustered ≈ 40-45 min and ≈ 6.5 GB. RSS grew linearly with the sequence count
  (≈ 24 MB per 1,000 seqs); not extrapolated to the full DB, since `-m` splits
  large indexes into parts (unknown at -m 10000). A full complete-DB run
  reached 5.2 GB after 63 min under load, then I stopped it (it would have
  taken 4-5 h at that load).
- **Fixed 2026-11-10.** New `scripts/sortmerna_index.py`:
  - one definition of the indexing parameters: `-L 18 --max_pos 10000 -m
    <--max_memory>`, passed explicitly. The DB stores them in
    `<db>.<kind>.smr/matam_index_params.json`, and every index of that DB,
    prebuilt or built at mapping time, uses them;
  - indexes are built in a temp dir and renamed into place with a record
    (`idx/matam_index.json`: fasta name, size, mtime, SHA-256, parameters,
    SortMeRNA version, index file sizes), under a lock file;
  - before each mapping MATAM checks the record (SHA-256 only when the size
    matches but the mtime doesn't). Stale, incomplete or other-parameter
    indexes are rebuilt;
  - an index without a record (existing DBs) is adopted if SortMeRNA's
    `<hash>.idx_stats.yaml` matches: same fasta name, size, number of sequences
    and total length, `L`, `max_positions` (≤ 10000 and not exactly 1000), a
    partitioning `-m` gives, intact files. The record is then added (0.26 s on
    30k seqs). A DB without a params file is treated as indexed by older MATAM
    with SortMeRNA's defaults (`-m 3072`, any number of parts), and keeps them
    when rebuilt. Not detectable: a hand-picked cap below 10000 other than 1000,
    and a same-size edit that keeps sequence count and total length;
  - mapping passes `--idx-dir ... --index 0` (SortMeRNA never builds itself),
    as a run option outside the completion record (#13). The index is resolved
    only when the step actually runs;
  - no usable index: built once in the DB dir if writable (concurrent runs wait
    on the lock), else privately in the step's SortMeRNA workdir
    (`sortmerna_*/idx`, reused by reruns). Each case is logged;
  - `matam_db_preprocessing.py --skip_index` (default: build), and `-m` restored
    there and in `index_ref_db.py`, which uses the same helper.
- **Default stays "build"**: every run needs both indexes, and the build is
  single-threaded and long on SILVA (see the progress log). Skipping it moves
  that cost into the first assembly run, or into every run on a read-only DB.
  `--skip_index` suits `--filter_only` users (only the clustered index is ever
  built), and DBs indexed later or elsewhere.
- **Not fixed (separate):** a read-only DB still fails at scaffolding:
  `samtools view -T <db>.complete.fasta` writes `<db>.complete.fasta.fai` next
  to the fasta (concurrent runs on a writable DB race on it too).
- **Tests:** `test_sortmerna_index.py` (23, real SortMeRNA: record, stale by
  size/checksum, copied fasta, other parameters, truncated files, adoption and
  its rejections, failed build, writable/read-only DB, 3 concurrent runs build
  once); `test_sortmerna_step.py` (`prepare` hook); `test_db_preprocessing.py`
  (record with `-m` = `--max_memory`, `--skip_index`); `test_functional_pe.py`
  (prebuilt index used; unindexed writable DB indexed once; unindexed read-only
  DB gives the same assembly as the prebuilt one, with private indexes).

### 16. [x] DB preprocessing always re-clusters (no way to reuse a clustered set)
- **Problem:** the 95% vsearch clustering of SILVA NR99 takes ~1 h 40 min on 8
  threads (288,470 centroids from 905,628 sequences), and was redone on every
  preprocessing, even with a clustered set at hand. (SILVA's own NR99 99%
  clustering is not a substitute: the clustered db is MATAM's deliberately less
  redundant recruitment set.)
- **Fix:** `matam_db_preprocessing.py --clustered_ref FASTA`. Only its IDs are
  used: the sequences are selected from the cleaned complete db (`seqkit grep`),
  so they get the same cleaning and their IDs are always in the complete db and
  taxonomy. IDs missing from the cleaned db are reported and left out; none found
  is an error. Mutually exclusive with `--clustering_id_threshold` and
  `--by_domain`. The db prefix is named after the clustered file.
- **Tests:** `test_db_preprocessing.py`: selection keeps the cleaned sequences
  in db order with ID-only headers, no vsearch run, missing IDs reported, tmp
  files removed; no matching ID is an error; conflicting options rejected.
- [x] **Are the dropped centroids orphan taxa?**
  - **2026-11-10, full SILVA, taxon level:** 1,811 of 288,470 centroids were
    dropped, all for >5 consecutive non-ACGT bases. That day's taxon labels are
    shifted by one: SILVA 144 has an extra rank, domain;kingdom;phylum;…, so
    "genus" was really the family level.
  - **2026-11-11, member level, on a random 10k-sequence SILVA subset** (full
    run stopped at the user's request). Clustered at 95% with `--uc`: 6,978
    centroids.
    - Cleaning removed 29 sequences, 23 of them centroids. N-rich sequences
      seldom reach 95% to anything, so they become centroids, mostly
      singletons. The 23 dropped clusters had only 7 members.
    - 2 members (both eukaryotes) have no kept centroid at >=95%. Their
      lineages are still represented, and their best kept centroids are at
      91.9% and 94.8%.
    - 4 prokaryote genera and 1 family lost their only centroid, but all of
      them also have no sequence left in the complete db. Cleaning removes
      them whatever the clustering, so `--clustered_ref` added no taxon loss
      here.
  - **Verdict:** orphaning from external centroids is marginal. The report
    (#17) flags it per run.
- **Also:** `--by_kingdom`/`--kingdoms` renamed `--by_domain`/`--domains`. The
  old names are still accepted, with a deprecation warning. The default db
  suffix is now `_bd`, was `_bk`. The domain is matched on the first taxonomy
  rank instead of anywhere in the header. On SILVA 144 this made no
  difference, but a name such as "Some bacteria endosymbiont" would have put a
  eukaryote into both sets.

### 17. [x] Safeguard: sequences lost before the clustered db are reported, with their nearest remaining relative
- **Problem:** sequences can be lost on the way to the clustered db with only a
  count in the log, or nothing at all:
  - removed by cleaning: a run of more than `--max_consecutive_n` (default 5)
    non-ACGT characters, N or any IUPAC code, or the optional length filters;
  - `--clustered_ref` centroids not in the cleaned db;
  - with `--by_domain`, sequences outside `--domains`, lost silently.

  Inside preprocessing, vsearch clusters the cleaned db, so no centroid is lost
  after clustering; there, the loss is the sequences removed before it.
- **Criterion: a close relative must remain, not a taxon name.** The first
  version checked taxon names. On the 10k subset names and similarity
  disagreed: Rathayibacter kept a namesake while its nearest remaining
  reference was 79.5%; Actinobacillus--other lost its name while Basfia
  remained at 94.6%.
  - Each discarded sequence (other domains excepted) is aligned to the
    remaining references, i.e. the cleaned complete db, with `vsearch
    --usearch_global` and the clustering's semi-global setting. At least 90%
    of the query must be covered: without that, free end gaps gave 1-column
    "100%" hits.
  - A close relative is left if the best identity is >= `--close_relative_id`
    (default: the clustering identity; 0 skips the search).
  - vsearch counts N and IUPAC codes as matching compatible bases.
  - The clustered db is no reference for this: centroids are, by
    construction, below the clustering identity to each other.
- **Files:**
  - `<db>.discarded.tsv`: id, reason, taxonomy, nearest remaining reference,
    identity, its lineage, and whether a close relative is left.
  - `<db>.discarded_taxa.tsv`: per reason and taxon, the level (1 = domain),
    lineage flag, discarded count, count without a close relative, and counts
    left by name in the complete and clustered db.
  - Levels are used because rank counts differ: SILVA 144 adds a kingdom rank;
    eukaryotes have 2-20 fields and are reported at domain level only. The
    organism name (last field) is not a taxon.
  - The log gives one line per reason, a WARNING when sequences leave no close
    relative, with the worst examples and their nearest identity.
- **10k SILVA subset** (1.1% of SILVA, so relatives are much sparser than in
  the full db):
  - All 6 non-centroid cleaning removals have a close relative left.
  - 20 of the 23 dropped centroids have none (nearest 79.5-94.6%).
  - Centroids are the longest sequences of their clusters; shorter members
    often cover <90% of them, which counts as no close relative over the
    whole sequence.
- **Tests:**
  - `test_discarded_refs_report.py` (9 unit tests).
  - `test_db_preprocessing.py`:
    - unrelated vs near-duplicate removed sequence (nearest ABA_X81660,
      >=99%);
    - `--close_relative_id 1.0` and `0`;
    - discarded `--clustered_ref` centroids, including one not in the input;
    - other domains.

### 18. [x] Pipelines silently index an unindexed database in every job (MarkerMAG-nf)
- Without a usable prebuilt index, `sortmerna_index.ensure` builds one: in the
  db directory if writable (shared, under a `flock`), otherwise privately in
  the run. A MarkerMAG-nf run on an unindexed SILVA db (read-only, or where
  `flock` fails, e.g. some Lustre mounts) therefore indexes it in every
  MATAM_FILTER and MATAM_ASSEMBLE task (samples x 7 depths, ~2.5-3 h each),
  with only a log warning.
- **Fixed 2026-11-11** (first as an opt-in `--require_index`, then made the
  default at the user's request):
  - `matam_assembly.py` uses only the db's prebuilt indexes by default
    (`sortmerna_index.prebuilt`: record check, or adoption of a record-less
    index). It never builds one and takes no lock. Every index the run needs
    (clustered; complete unless `--filter_only`) is checked before any step.
    Otherwise it exits 1, giving the reason and the `index_ref_db.py -i <db>`
    command to run. `--no_require_index` restores building on first use (in
    the db if writable, else in the run). `sortmerna_index.ensure` defaults
    to `require=True`.
  - `build.py` links `index_ref_db.py` into `bin/`. `--skip_index` help, log
    and README now say an index must be added before assembly.
  - MarkerMAG-nf: the wrappers never pass `--no_require_index`.
    `resolve_matam_db_prefix` stops the workflow before any task when either
    index has neither a MATAM record nor SortMeRNA index statistics (legacy,
    adopted by MATAM).
- **Tests:**
  - `test_sortmerna_index.py`: by default the prebuilt index is used without
    a lock; a missing, stale or other-parameter index is never built; a
    record-less index is adopted; a missing fasta is reported. Build-path
    tests pass `require=False`.
  - `test_functional_pe.py`: by default runs on prebuilt indexes, stops
    before any step on an unindexed db (full and `--filter_only`), and checks
    the scaffolding index before the read mapping; `index_ref_db.py` makes a
    db usable; the build-on-first-use tests pass `--no_require_index`.
  - MarkerMAG-nf: `test_python314_and_paired_reads.py` (neither wrapper
    passes `--no_require_index`).
  - Checked by hand: Nextflow 26.10.5 `-preview` against 5 db states (indexed
    as dir and prefix, legacy, unindexed, complete index missing); the real
    wrappers on the indexed SILVA-subset db (OK) and an unindexed db (filter
    exits 1 at once, nothing indexed).
- **Not addressed:** with `--no_require_index`, a failing `flock` on the
  run's own work directory (the private-index fallback) is not caught.

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
| 2026-11-10 | #13, #14 | Smoke test of options outside the e2e test: full run with taxonomy/Krona, read correction yes/auto, `--optimize_components`, `--straight_mode`, `--workdir`, gz input without `-v` all OK; resume (#13) and `--contigs_binning` (#14) fail |
| 2026-11-10 | #13 | Completeness records for SortMeRNA steps (`sortmerna_step.py`): complete output reused, incomplete/outdated output and run state removed before rerunning; resume and rerun verified end to end. README documents resuming. Suite 134 passed, 6 skipped |
| 2026-11-10 | #15 | SortMeRNA index: mapping reuses the prebuilt index but never checks it (stale index used silently or segfaults; `-m` changes results; concurrent builds race). New `sortmerna_index.py`: one set of indexing params (`-L 18 --max_pos 10000 -m <max_memory>`) per DB, record-checked indexes, adoption of record-less ones, locked atomic builds in the DB dir or the run. `--skip_index` (default: build). SILVA NR99 indexing extrapolated ≈ 2-2.5 h (complete) + 45 min (clustered). Suite 164 passed, 6 skipped |
| 2026-11-10 | #16 | `--clustered_ref` reuses a clustered set instead of re-clustering. Full SILVA NR99 with the 95% vsearch centroids: preprocessing (no index) in 32 s instead of ~1 h 40 min for vsearch; 286,659 of 288,470 centroids selected, 1,811 reported as removed by the N filter (centroids made from the uncleaned fasta). Suite 168 passed, 6 skipped |
| 2026-11-10 | #16 | `--clustered_ref` help/README no longer tied to 95%. `--by_kingdom` renamed to `--by_domain` (old names deprecated aliases; domain anchored to the first rank). Taxon-level orphan check of the 1,811 dropped centroids done; member-level check stopped and scheduled for next session. Suite 170 passed, 6 skipped |
| 2026-11-11 | #17 | Report of the reference sequences that don't reach the clustered db (cleaning, discarded `--clustered_ref` centroids, `--by_domain` other domains) with taxonomy and lost taxa per rank: `<db>.discarded.tsv`, `<db>.discarded_taxa.tsv`, log summary. Full SILVA re-clustering with `--uc` started for the real-data check and the scheduled orphan test |
| 2026-11-11 | #16, #17 | Report uses taxonomy levels, not rank names: SILVA 144 adds a kingdom rank, which shifted the 11-10 labels. Eukaryotes are reported at domain level only. Lineages lost only from the clustered db are told apart from those gone from the complete db. Orphan test on a 10k SILVA subset: 23 dropped centroids, 7 members, 2 orphans (euk, best kept centroid 91.9%/94.8%); the 4 genera + 1 family that lost their centroid are all gone from the complete db too. Suite 181 passed, 6 skipped |
| 2026-11-11 | #17 | Lost taxa judged by sequence similarity, not names: each discarded sequence's nearest remaining reference (vsearch, >=90% query coverage), close relative if >= `--close_relative_id` (default clustering identity). 10k subset: 0/6 cleaning removals and 20/23 dropped centroids without a close relative (sparse subset). Suite 183 passed, 6 skipped |
| 2026-11-11 | #14 | `--contigs_binning` evaluated for MarkerMAG-nf on a simulated SILVA-subset community (3 replicates): resolves >=1%-apart strains (3-4/4 vs 0-1/4) but halves exact sequences, adds mosaics and error duplicates, and loses novel 16S under 1200 bp. MarkerMAG-nf's wrapper would get no scaffolds with it (file renamed). Kept off; #14 fix (`LN:2147483647`) verified in a scratch copy only |
| 2026-11-11 | #14 | Fixed: `compute_contigs_compatibility.py` writes `LN:2147483647` instead of `LN:9999999999`. Untested by design: no test added or planned |
| 2026-11-11 | #18 | `--require_index` (prebuilt indexes only, checked before any step, no lock); `index_ref_db.py` in `bin/`; MarkerMAG-nf passes it and checks the db indexes before any task. MATAM suite 194 passed, 6 skipped; MarkerMAG-nf 14 passed; Nextflow `-preview` on 5 db states |
| 2026-11-11 | #18 | Inverted: prebuilt indexes are required by default (`sortmerna_index.ensure(require=True)`); `--no_require_index` builds them on first use; MarkerMAG-nf wrappers rely on the default. MATAM suite 194 passed, 6 skipped; MarkerMAG-nf 14 passed; Nextflow `-preview` and real wrappers re-checked |
