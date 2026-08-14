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
  - [2.1 Current development version](#21-current-development-version)
  - [2.2 Published Conda package](#22-published-conda-package)
  - [2.3 Published Docker image](#23-published-docker-image)
- [3. Run MATAM](#3-run-matam)
  - [3.1 Input data](#31-input-data)
  - [3.2 Database preparation (clustering \& indexation)](#32-database-preparation-clustering--indexation)
    - [3.2.1 Provided database](#321-provided-database)
    - [3.2.2 Custom database](#322-custom-database)
  - [3.3 De-novo assembly](#33-de-novo-assembly)
  - [3.4 Example with default database and provided dataset](#34-example-with-default-database-and-provided-dataset)
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

The current development version must be installed from this repository. The
published Bioconda package and Docker image contain the upstream MATAM release,
not the changes maintained in this fork.

## <a id="matam-from-source-code"></a>2.1 Current development version

Install Conda or Miniconda first. Then clone this repository, create its Conda
environment, and compile the two native executables:

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

## <a id="matam-with-conda"></a>2.2 Published Conda package

[![Anaconda-Server Badge](https://anaconda.org/bioconda/matam/badges/version.svg)](https://anaconda.org/bioconda/matam)

The upstream release, without this fork's development changes, is available as
`bioconda::matam`.

```bash
conda install -c conda-forge -c bioconda bioconda::matam
```

## <a id="matam-in-docker"></a>2.3 Published Docker image

The published BioContainer also contains the upstream release. To retrieve it,
run:

`docker pull quay.io/biocontainers/matam:{release_version}`

Then all the commands used in this README will be available as:

`docker run -v host_directory:/workdir quay.io/biocontainers/matam:{release_version} CMD`

Noticed that you have to specify a docker volume to share data between the host and the container and use this workdir for your analysis. Otherwise your data will be lost when exiting the container.

Finally, if you prefer an interactive session with the container, run:

`docker run -it quay.io/biocontainers/matam:{release_version}`


# <a id="run-matam"></a>3. Run MATAM

## <a id="input-data"></a>3.1 Input data

The reads used to reconstruct the markers must be provided under FASTQ format.  
Be aware than the **the sequence identifiers have to be unique** and the **the score must be Phred+33 encoded**.

## <a id="database-preparation"></a>3.2 Database preparation (clustering & indexation)

SortMeRNA 6 introduced an index format that is incompatible with indexes created by SortMeRNA 4 or 5. When upgrading, delete the `<prefix>.complete.smr` and `<prefix>.clustered.smr` directories, then rerun `index_ref_db.py` or `matam_db_preprocessing.py` before assembly.

### <a id="provided-database"></a>3.2.1 Provided database

By default, MATAM provides a SSU rRNA reference database where the clusterisation step has already been done (i.e. the sequences sharing 95% of identity have been clustered with [Sumaclust](https://git.metabarcoding.org/obitools/sumaclust/wikis/home)).  
The  [FASTA](https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz) file used for this database comes from [SILVA 128 release](https://www.arb-silva.de/documentation/release-128/).

To use the default SSU rRNA reference database, run the following command:

`index_default_ssu_rrna_db.py -d $DBDIR --max_memory 10000`

where `$DBDIR` is the directory used to store the database.

### <a id="custom-database"></a>3.2.2 Custom database

If the provided database does not fulfill your needs, you can prepare a custom database of your own by running the following command:

`matam_db_preprocessing.py -i ref_db.fasta -d $DBDIR --cpu 4 --max_memory 10000 -v`

where `$DBDIR` is the directory used to store the database.

## <a id="de-novo-assembly"></a>3.3 De-novo assembly

When your database is ready, then you will be able to reconstruct your markers:
* Assembly only  
  In this mode, MATAM will reconstruct the full length sequences present in the sample.  
  `matam_assembly.py -d $DBDIR/prefix -i reads.fastq --cpu 4 --max_memory 10000 -v`

* Assembly and taxonomic assignment  
  In this mode, MATAM additionnaly provides a taxonomic classification of the sequences found, together with their abundance. Note that the
  classification is done with RDP with the default training model "16srrna". So this mode may be not suitable for other phylogenetic markers.  
  `matam_assembly.py -d $DBDIR/prefix -i reads.fastq --cpu 4 --max_memory 10000 -v --perform_taxonomic_assignment`  
    The taxonomic assignment is done with [RDP classifier](https://rdp.cme.msu.edu/) and the training model used by default is "16srrna"

SortMeRNA run state is stored in the MATAM output directory by default. Use `--workdir WORKDIR` to store it elsewhere; reference indexes remain under the database `.smr/idx` directories.

where `$DBDIR` is the database directory and `prefix` is the common prefix used to name the database files.
For example, with the default database, the prefix is SILVA_128_SSURef_NR95.  


## <a id="example-with-default-database-and-provided-dataset"></a>3.4 Example with default database and provided dataset

1. Retrieve the example dataset: [16 bacterial species simulated dataset](examples/16sp_simulated_dataset/16sp.art_HS25_pe_100bp_50x.fq)

  `wget https://raw.githubusercontent.com/bonsai-team/matam/master/examples/16sp_simulated_dataset/16sp.art_HS25_pe_100bp_50x.fq`

2. Getting and indexing default SSU rRNA reference database

  `index_default_ssu_rrna_db.py -d $DBDIR --max_memory 10000`

3. De-novo assembly

  `matam_assembly.py -d $DBDIR/SILVA_128_SSURef_NR95 -i 16sp.art_HS25_pe_100bp_50x.fq --cpu 4 --max_memory 10000 -v --perform_taxonomic_assignment`


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
