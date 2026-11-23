<h1> <p align="center"> Welcome to MOPP2 - Multiomics Processing Pipeline 2! </p> </h1>

<h3> <p align="center"> How can MOPP2 help me? </p> </h3>

MOPP2 processes paired metagenomic (**metaG**), metatranscriptomic (**metaT**) and
metaRibo-Seq (**metaRS**, ribosome-protected fragments) data from the same samples.
Every sample gets its **own reference database**, built from the genomes actually
detected in that sample, and every library is mapped to its sample's database. The
result is a set of gene- and genome-level count, TPM and RPKM tables for each
omic.

There are two typical use cases:

(1) To run every step from one config file, run `mopp2 workflow --help`

(2) To run a single step, run `mopp2 trim --help`, `mopp2 map --help`, etc.

```
mopp2 metadata      create / validate the sample metadata
mopp2 trim          Trim Galore, per-omic parameters
mopp2 host-deplete  remove host reads                 (host-associated samples only)
mopp2 subdb         per-sample reference databases    (--mode curated | reuse | globdb)
mopp2 map           bowtie2 -> per-sample BAMs
mopp2 annotate      prodigal + eggNOG-mapper
mopp2 blacklist     ORFs lying inside rRNA / tRNA genes
mopp2 saf           per-sample SAF without those ORFs
mopp2 features      featureCounts -> count / TPM / RPKM tables
```

<h3> <p align="center"> What do I need to run MOPP2? </p> </h3>

1. An input folder with your raw sequencing data (`*.fastq.gz` or `*.fq.gz`).
2. A tab-delimited metadata file describing the sequencing data, following this
   [template](tests/test_data/metadata_template.tsv). `mopp2 metadata` can generate
   it from the file names. Identifiers may only contain letters, digits, `_` and `-`.
   Single-end libraries (always the case for metaRS) have only an `r1` row.
3. A config file (`mopp2 init-config`), which points every step after `trim` at
   your reads and databases.
4. Reference databases (see [Databases](#databases)): a GlobDB sylph sketch plus
   the GlobDB genomes (or your own genome pool), the eggNOG-mapper database, and a
   bowtie2 index of a tRNA database. A host genome bowtie2 index is only needed for
   host-associated samples.

MOPP2 runs on Linux. The step scripts use GNU coreutils and GNU grep (`grep -P`).

***
<h2> <p align="center"> Installation </p> </h2>

To install the most up-to-date version of mopp2, run:
```bash
git clone https://github.com/<your-account>/mopp2.git
cd mopp2
conda env create -f environment.yml      # or: mamba env create -f environment.yml
conda activate mopp2
pip install .
```

Run the following commands to make sure the installation is complete:
```bash
mopp2 --help
pytest tests          # runs a real Trim Galore job on synthetic reads
```

If you did `pip install -e .`, which is a local (editable) installation, you need
to keep the source code for the package to run.

### Building the conda package

The repository also contains a conda recipe. To build mopp2 as a conda package and
install it, with all its dependencies, into a new environment:
```bash
conda install -n base conda-build          # once
conda build conda-recipe -c conda-forge -c bioconda
conda create -n mopp2 -c local -c conda-forge -c bioconda mopp2
conda activate mopp2
```
`conda build` runs the tests in `tests/` against the built package.

### Optional: CheckM2

`mopp2 subdb --mode globdb` uses `checkm2` to score the project's own MAGs, if it
is on `PATH`. Without it, every genome gets `quality_score=70`. CheckM2 is not in
`environment.yml` because of its size. Install it in the same env if you have MAGs.

***
<h2> <p align="center"> Databases </p> </h2>

None of the databases ships with mopp2. Point the config at them:

| config variable | used by | what |
|---|---|---|
| `SYLPH_DB_R1` | `subdb --mode globdb` | sylph sketch of GlobDB (we use `globdb_r226_sylph_c200.syldb`, built with `sylph sketch -c 200`) |
| `GLOBDB_ROOT` | `subdb --mode globdb` | the directory the sketch was built in, holding the GlobDB genome fasta files. sylph's `Genome_file` paths are relative to it. |
| `EGGNOG_DATA_DIR` | `annotate` | eggNOG-mapper data directory with `eggnog.db` and `eggnog_proteins.dmnd` (we use emapper DB 5.0.2; `download_eggnog_data.py`) |
| `TRNA_IDX` | `blacklist` | bowtie2 index of tRNA sequences, e.g. all bacterial tRNAs from GtRNAdb, built with `bowtie2-build` |
| `HUMAN_IDX` | `host-deplete` | bowtie2 index of the host genome, e.g. GRCh38.p14 |

***
<h2> <p align="center"> Documentation </p> </h2>

<h4> <p align="center"> test runs </p> </h4>

Use case 1: `mopp2 workflow`

```bash
mopp2 metadata -i rawdata/ -o project/metadata.tsv -p     # 1. metadata
mopp2 init-config -o project/project.config.sh --project project/   # 2. config
$EDITOR project/project.config.sh                          #    fill in databases + libraries

mopp2 workflow project/project.config.sh \
   --subdb-mode globdb \             # how the per-sample reference databases are built
   -i rawdata/ \                     # raw fastq.gz (enables trim)
   -m project/metadata.tsv \         # metadata (enables trim)
   -t 4                              # libraries trimmed at the same time
```

Use case 2: `mopp2 <module>`

If you prefer to run each step separately:
```bash
# step 1a: auto-generate metadata
mopp2 metadata -i rawdata/ -o project/metadata.tsv -p

# step 1b: validate metadata
mopp2 metadata -i project/metadata.tsv -p

# step 2: trim
mopp2 trim -i rawdata/ -m project/metadata.tsv -o project/2.trim -t 4

# step 3: write a config and fill it in
mopp2 init-config -o project/project.config.sh --project project/

# step 3 (optional): host depletion, for stool / skin / plant tissue
mopp2 host-deplete project/project.config.sh

# step 4: per-sample reference databases
mopp2 subdb project/project.config.sh --mode globdb

# step 5: map every library to its sample's database
mopp2 map project/project.config.sh

# step 6: gene prediction + functional annotation
mopp2 annotate project/project.config.sh

# step 7: structural-RNA ORF blacklist, then the SAFs without those ORFs
mopp2 blacklist project/project.config.sh
mopp2 saf project/project.config.sh

# step 8: feature tables
mopp2 features project/project.config.sh
```

`map` and `annotate` are independent, so they can run at the same time.

<h2> <p align="center"> </p> </h2>

<h4> <p align="center"> mopp2 workflow </p> </h4>

usage: `mopp2 workflow <CONFIG> --subdb-mode <curated|reuse|globdb> [-i <Input Directory> -m <Metadata (tsv)>] [--host-depletion] [--from STEP] [--to STEP]`

Runs the whole chain from one config:

```
trim -> [host-deplete] -> subdb -> map -> annotate -> blacklist -> saf -> features
```

- `trim` runs only when `-i` and `-m` are given. It writes to `--trim-dir`, which
  defaults to `TRIM_DIR` from the config, else `<PROJECT>/2.trim`.
- `host-deplete` runs only with `--host-depletion`.
- Every step skips work that is already done. If a step fails, fix it and resume
  with `--from <step>`. `--to <step>` stops early.

<h2> <p align="center"> </p> </h2>

<h4> <p align="center"> mopp2 metadata </p> </h4>

usage: `mopp2 metadata -i <Input Directory/Input Path> [-o <Output Path>] [-p] [-m]`

This module creates and/or validates the metadata.

| sample_name | identifier | omic | strand |
|---|---|---|---|
| S1_metaG_S121_L004_R1_001.fastq.gz | S1 | metaG | r1 |
| S1_metaG_S121_L004_R2_001.fastq.gz | S1 | metaG | r2 |
| S1_metaT_S37_L004_R1_001.fastq.gz | S1 | metaT | r1 |
| S1_metaT_S37_L004_R2_001.fastq.gz | S1 | metaT | r2 |
| S1_metaRS_S13_L004_R1_001.fastq.gz | S1 | metaRS | r1 |

If the input is a **directory**, mopp2 lists every file ending in `.fastq.gz` or
`.fq.gz` that has a `metaG` / `metaT` / `metaRS` token in its name. The identifier
is the part before that token, and the strand is the last `_R1` / `_R2` after it
(`r1` if there is none). The result is written to `-o` (default
`<dir>/metadata_auto.tsv`) and validated. Files whose names cannot be parsed are
reported and left out, so check the result and edit it if needed.

If the input is a **metadata file**, it is validated. `-o` is ignored.

`--paired-end` requires an r1 and an r2 for every metaG / metaT library.
`--multiomics` requires metaG, metaT and metaRS for every sample. metaRS is always
single-end. A log file is written next to the metadata.

<h2> <p align="center"> </p> </h2>

<h4> <p align="center"> mopp2 trim </p> </h4>

usage: `mopp2 trim -i <Input Directory> -o <Output Directory> -m <Metadata (tsv)> [-t <Num Jobs>]`

`mopp2 trim` trims every library in the metadata with Trim Galore. The omic
decides the parameters:

| omic | mode | Trim Galore options |
|---|---|---|
| metaG, metaT | paired-end | `--paired --length 20` |
| metaRS | single-end | `--length 18 --max_length 75` |

The metaRS limits keep ribosome-protected fragments (~18-75 nt after adapter
removal) and drop reads that were never cut, which are not footprints. FastQC runs
on the trimmed reads unless `--no-fastqc` is given.

**R1 and R2 are kept as separate files**, because `subdb` and `map` use them
as pairs. The output names are what the config expects by default:

```
2.trim/
├── <id>_metaG_R1.fq.gz   <id>_metaG_R2.fq.gz
├── <id>_metaT_R1.fq.gz   <id>_metaT_R2.fq.gz
├── <id>_metaRS.fq.gz
├── trim_summary.tsv                 one row per library: trimmed / skipped / failed
├── mopp2_trim_<time>.log
└── trimming_reports/                Trim Galore reports, FastQC, one log per library
```

A library whose outputs already exist is skipped, so an interrupted run can simply
be started again (`--force` re-trims). `-t` sets how many libraries are trimmed at
once and `--cores-per-job` passes `--cores` to each Trim Galore job.
`--min-length`, `--metars-min-length` and `--metars-max-length` change the length
limits. `--dry-run` prints the Trim Galore commands without running them.

example:
```bash
mopp2 trim \
   -i rawdata/ \
   -o project/2.trim \
   -m project/metadata.tsv \
   -t 4
```

rRNA removal from metaT (e.g. SortMeRNA) is not part of mopp2. Point `METAT_DIR`
at its clean reads, or at `2.trim` with `METAT_R1_SUFFIX=_metaT_R1.fq.gz` to map
the trimmed reads as they are.

<h2> <p align="center"> </p> </h2>

<h4> <p align="center"> config file </p> </h4>

usage: `mopp2 init-config -o <Config Path> [--project <Project Directory>]`

Every step after `trim` takes one argument, a bash config file. `init-config`
writes the commented template ([`mopp2/scripts/config.template.sh`](mopp2/scripts/config.template.sh)).
Fill in `PROJECT`, the database paths and one section per library type, and delete
the sections your project does not have. Reads are found as
`${DIR}/<sample><SUFFIX>`, for example `${METAG_DIR}/<sample>_metaG_R1.fq.gz`.
`STRAND` is the featureCounts `-s` value (metaT is usually reverse-stranded, `2`).

By default the project directory fills up like this:

```
<PROJECT>/
├── 2.trim/              mopp2 trim
├── 3.host_depletion/    mopp2 host-deplete
├── 5.subdb/             mopp2 subdb, annotate, blacklist, saf
├── 6.mapping/           mopp2 map
└── 7.features/          mopp2 features
```

<h2> <p align="center"> </p> </h2>

<h4> <p align="center"> mopp2 host-deplete </p> </h4>

usage: `mopp2 host-deplete <CONFIG>`

Maps `TRIM_DIR` against `HUMAN_IDX` with bowtie2 and keeps the reads that do not
align (`HF_OUT/filtered`). Only use it for host-associated samples, and then point
the library directories at `READS_DIR`.

metaRS is left alone while `SKIP_METARS_FILTER=1` (the default). Short (~28 bp)
reads align to a 3 Gb host genome by chance. In one stool sample, 11.3% of metaRS
reads looked human (78% of them multi-mapping) against 0.67% of metaG reads.

`mopp2 drop-trimmed <CONFIG> [--watch]` deletes a sample's trimmed metaG once
its host-depleted copy has been verified. It is meant for when there is no disk
space to keep both copies.

<h2> <p align="center"> </p> </h2>

<h4> <p align="center"> mopp2 subdb </p> </h4>

usage: `mopp2 subdb <CONFIG> --mode <curated|reuse|globdb>`

Builds one reference database per sample, so reads are never mapped against
genomes that sample does not contain. Choose the mode by what you already have:

| `--mode` | you have | what it does |
|---|---|---|
| `curated` | a SynCom: two hand-curated genome pools (`MAGS_POOL`, `GLOBDB_POOL`) | sylph per sample against the pool |
| `reuse` | a dropout set: another project's `unified_db.syldb` (`SUBDB_SRC`) | sylph per sample against that database |
| `globdb` | a real community: only reads | sylph against GlobDB finds the pool; optionally adds the project's MAGs (`REFINEMENT_ROOT`), isolates (`ISOLATE_LIST`) and a fungus (`FUNGAL_FASTA`) |

In `globdb` mode, GlobDB genomes detected in any sample are pooled with the MAGs,
and the pool is dereplicated with skani at 99% ANI. Genomes with 95-99% ANI are kept
as separate strains. Each sample then keeps the genomes detected in it. They are
selected either with Eff_cov cutoffs (`MIN_EFF_COV`, `MIN_EFF_COV_R2`) or by
sylph's own confidence tier (`SYLPH_TIER_MAX=T1..T4`). The project's own genomes
are included in a sample only where sylph calls them present
(`KNOWN_MIN_ANI_LO`). `EXCLUDE_REFS` keeps chosen genomes out entirely, for
example genomes hit only by a spiked-in plasmid.

Outputs: `SUBDB_ROOT/unified_db.tsv` and `SUBDB_ROOT/samples/<sample>/{references.tsv,ref.fna,bowtie2 index}`.

<h2> <p align="center"> </p> </h2>

<h4> <p align="center"> mopp2 map </p> </h4>

usage: `mopp2 map <CONFIG>`

Maps every library to its own sample's database with bowtie2, then sorts and
indexes the BAMs. metaG and metaT are mapped paired-end on default settings.
metaRS is mapped single-end, twice:

| variant | bowtie2 | |
|---|---|---|
| `metaRS_strict` | `--sensitive -L 18 -N 1 --score-min L,0,-0.01` | about zero mismatches |
| `metaRS_mopp2` | `--sensitive -L 18 -N 1 --score-min L,0,-0.3` | about one mismatch per footprint (used downstream) |

`METARS_BAM_SUBDIR` picks which one `features` counts.

Outputs: `MAPPING_OUT/<sample>/{metaG,metaT,metaRS_strict,metaRS_mopp2}/<sample>.bam` and `MAPPING_OUT/manifest.tsv`.

<h2> <p align="center"> </p> </h2>

<h4> <p align="center"> mopp2 annotate </p> </h4>

usage: `mopp2 annotate <CONFIG>`

Predicts genes on the unified genome pool with prodigal (`-p meta`) and annotates
them with eggNOG-mapper. A fungal genome can be annotated from JGI gene models
instead (`FUNGAL_REF_ID`, `FUNGAL_GFF`, `FUNGAL_PROT`). `SKIP_EGGNOG=1` stops after
prodigal. Everything `blacklist`, `saf` and `features` need is still written, and
a later run with `SKIP_EGGNOG=0` adds the functional layer (`gene_metadata.tsv`)
without redoing anything else.

<h2> <p align="center"> </p> </h2>

<h4> <p align="center"> mopp2 blacklist & mopp2 saf </p> </h4>

usage: `mopp2 blacklist <CONFIG>` then `mopp2 saf <CONFIG> [--cds <FILE>] [--mode auto|filter|rebuild]`

Prodigal predicts spurious ORFs inside rRNA and tRNA genes, and those ORFs draw a
large share of metaRS reads. `blacklist` removes an ORF when more than half of it
lies inside an rRNA locus found by barrnap (`MIN_OVERLAP_FRAC`), or when it has an
aligned block of at least 50 bp against the tRNA database (`MIN_TRNA_BLOCK`). Do not
lower `MIN_TRNA_BLOCK`: at ~20 bp the false-positive rate is ~98%. `saf` then
writes each sample's `ref.saf` without those ORFs, so featureCounts never sees
them. `features` refuses to count metaRS without the blacklist.

<h2> <p align="center"> </p> </h2>

<h4> <p align="center"> mopp2 features </p> </h4>

usage: `mopp2 features <CONFIG>`

Counts reads per gene with featureCounts (`-O --fraction -M --fraction`, strandedness
from `*_STRAND`) and builds the tables:

```
7.features/
├── fc_per_sample/<sample>.<lib>.counts(.summary)
└── tables/
    ├── gene_features.<lib>.{count,tpm}.tsv
    ├── genome_features.<lib>.{count,tpm,rpkm}.tsv
    ├── gene_annotation.tsv
    ├── genome_annotation.tsv
    └── counting_summary.tsv
```

Gene TPM is normalized by gene length. Genome-level counts are gene counts summed
per genome, and genome TPM and RPKM are normalized by genome length (from
`unified_db.tsv`). `genome_annotation.tsv` also lists each genome's gene count and
total CDS length, if you need a different normalization.

***
<h2> <p align="center"> Running the scripts without the mopp2 command </h2>

Each step is a bash script in `mopp2/scripts/` (`mopp2 scripts-dir` prints the
installed location) and can be run as `bash <script>.sh <CONFIG>`. Run that way,
the scripts `conda activate` the environments named in the config
(`MAP_CONDA_ENV`, `FEAT_CONDA_ENV`, `CONDA_ENV`). Through the `mopp2` command they
use the active environment instead. Pass `mopp2 --conda-switch <step> ...` to keep
the old behaviour.

The design notes behind each step, with the measurements that motivated them, are
in [`docs/`](docs/).

***
<h2> <p align="center"> FAQs </p> </h2>

1. `mopp2 subdb --mode globdb` says every genome is missing.

   `GLOBDB_ROOT` must be the directory the sylph sketch was built in (e.g. `SylDB/`),
   not the `globdb_r226_genome_fasta/` directory inside it.

2. `mopp2 subdb --mode globdb` runs out of memory.

   Profiling against the GlobDB r226 sketch (34 GB) needs ~40 GB of RAM. All samples
   are profiled in one sylph call, so the sketch is loaded only once.

3. `mopp2 features` refuses to run on metaRS.

   Run `mopp2 blacklist` and `mopp2 saf` first. Without the blacklist, 38-68% of
   metaRS CDS-assigned reads land on rRNA ORFs, unevenly across samples. Set
   `ALLOW_UNFILTERED_METARS=1` only if you know why you want that.

4. `trim_galore` / `bowtie2` / ... not found.

   Activate the environment first (`conda activate mopp2`). If you run the scripts
   with `bash` directly, check `MAP_CONDA_ENV` / `FEAT_CONDA_ENV` in the config.
