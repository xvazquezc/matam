MATAM
=========
*Mapping-Assisted Targeted-Assembly for Metagenomics*

MATAM is a software dedicated to the fast and accurate targeted assembly of short reads sequenced from a genomic marker of interest.
It has been applied to the assembly of 16S rRNA markers and is validated on simulated, synthetic and genuine metagenomes.

The related article of this method is available [here](https://doi.org/10.1093/bioinformatics/btx644).


# Table of contents
- [MATAM](#matam)
- [Table of contents](#table-of-contents)
- [1. Hardware requirements](#1-hardware-requirements)
- [2. Installation](#2-installation)
- [3. Run MATAM](#3-run-matam)
  - [3.1 Input data](#31-input-data)
  - [3.2 Database preparation (clustering \& indexation)](#32-database-preparation-clustering--indexation)
    - [3.2.1 SILVA database](#321-silva-database)
    - [3.2.2 Custom database](#322-custom-database)
  - [3.3 De-novo assembly](#33-de-novo-assembly)
  - [3.4 Example with the provided dataset](#34-example-with-the-provided-dataset)
- [4. Samples comparaison](#4-samples-comparaison)
- [5. Release versioning](#5-release-versioning)


# <a id="hardware-requirements"></a>1. Hardware requirements

We recommend running MATAM with at least 10Gb of free RAM. You can try
running MATAM with less RAM by setting  --max_memory  to a lower value
(eg. --max_memory 4000 for 4Gb). However, this could degrade the performance of the
software.

Some steps of MATAM are highly parallelized. So you can get a
significant speed up  by setting the --cpu option according to your
hardware configuration.


# <a id="installation"></a>2. Installation

MATAM is installed from source. Install Conda or Miniconda first, then clone
this repository, create its Conda environment, and compile the two native
executables:

```bash
git clone --recursive https://github.com/xvazquezc/matam.git
cd matam
conda env create --file environment.yml
conda activate matam
./build.py
export PATH="$(pwd)/bin:$PATH"
```

To install the dependencies into an existing environment instead, replace the
environment creation and activation commands above with:

```bash
conda env update --name <environment-name> --file environment.yml
conda activate <environment-name>
```

The `export` command affects the current shell only. Add the absolute path to
this repository's `bin` directory to your shell or workflow configuration if
MATAM must remain available after the shell exits. Verify the installation
with:

```bash
matam_assembly.py --help
```

`build.py` also initializes missing Git submodules, so an existing non-recursive
clone can be built with the same command.

`matam_assembly.py --version` (also `matam_db_preprocessing.py` and
`index_ref_db.py`) prints the version, which is also written at the start of
every run log. On the commit tagged `v<VERSION>` it is the number in the
`VERSION` file. On any other commit the short commit hash is added
(`<VERSION>+<hash>`). While `VERSION` is empty, the version is just the hash.

# <a id="run-matam"></a>3. Run MATAM

MATAM accepts paired reads only. Supply forward and reverse FastQ files with
`--forward` and `--reverse`; single-end and pre-interleaved input are not
supported.

Paired reads stay in separate R1/R2 files throughout. The read pairs selected by
the mapping step (pairs with at least one mate matching the database, both mates
kept) are written to the output directory as `marker_reads_R1.fq.gz` and
`marker_reads_R2.fq.gz`, with the original read names, in step. With
`--filter_only`, MATAM stops after this step, so these two files are its output.

## <a id="input-data"></a>3.1 Input data

The reads used to reconstruct the markers must be provided under FASTQ format (plain or gzipped),
with **Phred+33 encoded** quality scores.

The two files must list the mates of each pair in the same order. **Read identifiers must be
unique per pair**, and both naming conventions are supported as is, so there is no need to rename your reads:

- identical identifiers for both mates, as written by Illumina (`@ID 1:N:0:...` / `@ID 2:N:0:...`) or by SRA dumps (`@ID` / `@ID`);
- legacy mate suffixes (`@ID/1` / `@ID/2`).

MATAM works on mates: after the read mapping it records which mate each alignment comes from
(standard SAM flags 0x40/0x80), so the two mates of a pair are filtered, placed in the overlap
graph (one node per mate), and counted for abundance as two distinct reads, whatever their
naming. The results do not depend on the naming convention.

Each read pair is assembled in every graph component that holds one of its mates.
How pairs are handled from mapping to assembly, and how the results compare with
upstream MATAM, is documented in
[evaluation/paired_end_assembly](evaluation/paired_end_assembly/README.md).

## <a id="database-preparation"></a>3.2 Database preparation (clustering & indexation)

SortMeRNA 6 introduced an index format that is incompatible with indexes created by SortMeRNA 4 or 5. When upgrading, delete the `<prefix>.complete.smr` and `<prefix>.clustered.smr` directories, then rerun `index_ref_db.py` or `matam_db_preprocessing.py` before assembly.

### <a id="silva-database"></a>3.2.1 SILVA database

MATAM has no bundled database (the SILVA 128 database of upstream MATAM is no
longer downloadable). Build one from the SSU Ref NR99 sequences of the current
[SILVA release](https://www.arb-silva.de/):

```bash
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
gunzip SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
matam_db_preprocessing.py -i SILVA_138.2_SSURef_NR99_tax_silva.fasta -d $DBDIR \
    --out_db_name SILVA_138.2_NR95 --cpu 8 --max_memory 10000 -v
```

where `$DBDIR` is the directory used to store the database. The database prefix,
given to `matam_assembly.py -d`, is then `$DBDIR/SILVA_138.2_NR95`. Clustering and
indexing take a few hours for SILVA NR99 (details below).

### <a id="custom-database"></a>3.2.2 Custom database

Any reference fasta with SILVA-formatted headers (`>ID taxonomy`) can be prepared
the same way:

`matam_db_preprocessing.py -i ref_db.fasta -d $DBDIR --cpu 4 --max_memory 10000 -v`

where `$DBDIR` is the directory used to store the database.

This also builds the SortMeRNA indexes of the database (`<prefix>.complete.smr/idx`,
`<prefix>.clustered.smr/idx`), with `-m` set by `--max_memory`. Indexing is single-threaded:
for a full SILVA SSU NR99 database, expect about 2-2.5 h and 15 GB for the complete index
and 45 min and 6.5 GB for the clustered one. Each index carries a record of the fasta and
indexing parameters it was built from (`matam_index.json`).

`matam_assembly.py` uses only these prebuilt indexes. If an index the run needs is missing,
incomplete, or doesn't match its fasta or the database's parameters, it stops with an error
before any step. This stops, for example, every job of a pipeline from spending hours on
indexing. It takes no file lock, so it also works on filesystems without lock support. To
index an existing database (e.g. one made with `--skip_index`), run
`index_ref_db.py -i $DBDIR/<prefix> --max_memory 10000`.

With `--no_require_index`, `matam_assembly.py` instead builds each missing or out-of-date
index the first time it needs it. It builds it in the database directory if that is writable
(once, shared by later runs; concurrent runs wait for it), otherwise inside the run, in
every run.

The clustered database (used for read mapping and the overlap graph) is made by clustering
the references with vsearch at `--clustering_id_threshold` identity (default 0.95), which
takes about 1 h 40 min on 8 threads for SILVA NR99 at 95%. `--by_domain` clusters each domain
(first rank of the taxonomy) separately. If you already have a clustered set, at whatever
identity you chose, e.g. the clustered database of an earlier preprocessing, pass it with
`--clustered_ref clustered.fasta` to skip the clustering: its sequence IDs select the
clustered sequences from the cleaned input database, and the database prefix is named after
that file. Make it from the cleaned database: centroids removed by the cleaning (N or length
filters) leave their clusters without a representative, and are reported. SILVA's NR99 is
only made non-redundant at 99%; mapping against it unclustered is slower and spreads each
read's best alignments over near-identical references.

Every reference sequence that does not reach the clustered database is reported in
`<prefix>.discarded.tsv`. That covers sequences removed by cleaning (a run of more than
`--max_consecutive_n` ambiguous bases, or the optional length filters), `--clustered_ref`
centroids that are not in the cleaned database, and with `--by_domain` sequences of other
domains. Except for other domains, each one is aligned to the remaining references (the
cleaned complete database) and its nearest reference is reported, with the identity and
taxonomy. It has a close relative left if that identity is at least `--close_relative_id`
(default: the clustering identity) over at least 90% of the sequence; `--close_relative_id 0`
skips the search. `<prefix>.discarded_taxa.tsv` summarises this by taxon and taxonomy level
(1 = domain; eukaryotes at the domain level only). The log warns about discarded sequences
that leave no close relative.

## <a id="de-novo-assembly"></a>3.3 De-novo assembly

When your database is ready, then you will be able to reconstruct your markers:
* Assembly only  
  In this mode, MATAM will reconstruct the full length sequences present in the sample.  
  `matam_assembly.py -d $DBDIR/prefix --forward reads_R1.fastq --reverse reads_R2.fastq --cpu 4 --max_memory 10000 -v`

* Assembly and taxonomic assignment  
  In this mode, MATAM additionnaly provides a taxonomic classification of the sequences found, together with their abundance. Note that the
  classification is done with RDP with the default training model "16srrna". So this mode may be not suitable for other phylogenetic markers.  
  `matam_assembly.py -d $DBDIR/prefix --forward reads_R1.fastq --reverse reads_R2.fastq --cpu 4 --max_memory 10000 -v --perform_taxonomic_assignment`
    The taxonomic assignment is done with [RDP classifier](https://rdp.cme.msu.edu/) and the training model used by default is "16srrna"

SortMeRNA run state is stored in the MATAM output directory by default. Use `--workdir WORKDIR` to store it elsewhere; reference indexes stay under the database `.smr/idx` directories (or, with `--no_require_index` and a read-only database without prebuilt indexes, in the SortMeRNA work directories of the run).

A run can be resumed with `--resume_from STEP`, or rerun into the same output directory.
After each SortMeRNA mapping, MATAM records the command, inputs and outputs in a
`*.matam_complete.json` file next to the outputs. A complete mapping of the same inputs is
reused; an incomplete one (e.g. an interrupted run) or one made from other inputs is
removed, together with SortMeRNA's run state, and redone. Resuming after the read mapping
requires its output to be complete; otherwise MATAM asks to resume from `reads_mapping`.

`--stop_after STEP` ends the run after that step (same step names as `--resume_from`). The
files the next steps read are kept, so the run can be continued later with
`--resume_from` the next step. `--stop_after scaffolding` gives the final assembly
(`final_assembly.fa`, the scaffolds without abundances) and skips the abundance calculation,
which maps every selected read to the scaffolds. With it, a run that assembles no scaffolds
ends normally with an empty assembly.

MATAM requests gzip-compressed SortMeRNA read-mapping reports to reduce
temporary storage. It automatically accepts SortMeRNA's `.fq` or `.fastq`
FASTQ naming, with or without a `.gz` suffix.

where `$DBDIR` is the database directory and `prefix` is the common prefix used to name the database files.
For example, for the database built in [3.2.1](#silva-database), the prefix is SILVA_138.2_NR95.


## <a id="example-with-the-provided-dataset"></a>3.4 Example with the provided dataset

1. The example dataset, [16 bacterial species simulated dataset](examples/16sp_simulated_dataset/16sp.art_HS25_pe_100bp_50x.fq),
   is a single interleaved file. Split it into R1/R2 files:

   ```bash
   paste - - - - - - - - < examples/16sp_simulated_dataset/16sp.art_HS25_pe_100bp_50x.fq \
       | awk -F'\t' '{print $1"\n"$2"\n"$3"\n"$4 > "16sp_R1.fastq"; print $5"\n"$6"\n"$7"\n"$8 > "16sp_R2.fastq"}'
   ```

2. Build the SILVA database ([3.2.1](#silva-database)).

3. De-novo assembly and taxonomic assignment:

   `matam_assembly.py -d $DBDIR/SILVA_138.2_NR95 --forward 16sp_R1.fastq --reverse 16sp_R2.fastq --cpu 4 --max_memory 10000 -v --perform_taxonomic_assignment`


# <a id="samples-comparaison"></a>4. Samples comparaison

We provide a script to compare the abundances of different samples. Available only if the `--perform_taxonomic_assignment` was used when running MATAM.

`matam_compare_samples.py -s samples_to_compare.tsv -t contingency_table.tsv -c comparaison_table.tsv`

The `samples_to_compare.tsv` file is a tabulated file listing the FASTA & RDP files of each sample to compare (see example below).  
The `contingency_table.tsv` file will report the abundance for each sequence.  
The `comparaison_table.tsv` file will report a comparaison by "taxonomic path" of the abundance for the samples.


Example
```
sample1 <tab> $WORKDIR/matam_sample1/final_assembly.fa <tab> $WORKDIR/matam_sample1/rdp.tab
sample2 <tab> $WORKDIR/matam_sample2/final_assembly.fa <tab> $WORKDIR/matam_sample2/rdp.tab
```

The first column is the ID of the sample and it must be unique among the file.  

# <a id="release-versioning"></a>5. Release versioning

MATAM releases will be following the Semantic Versioning 2.0.0 rules described here: http://semver.org/spec/v2.0.0.html
