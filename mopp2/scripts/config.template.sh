# config.template.sh - one config drives every mopp2 step.
#
#   mopp2 init-config -o myproject.config.sh --project /path/to/project
#   $EDITOR myproject.config.sh
#
#   mopp2 trim          -i rawdata/ -m metadata.tsv -o /path/to/project/2.trim
#   mopp2 host-deplete  myproject.config.sh        # host-associated only
#   mopp2 subdb         myproject.config.sh --mode curated|reuse|globdb
#   mopp2 map           myproject.config.sh
#   mopp2 annotate      myproject.config.sh
#   mopp2 blacklist     myproject.config.sh
#   mopp2 saf           myproject.config.sh
#   mopp2 features      myproject.config.sh
#
# or all of it at once:  mopp2 workflow myproject.config.sh --subdb-mode globdb
#
# Every variable below except PROJECT has a working default or is optional.
# Sections marked [3c] and [host] are only read by those steps.

# ============================================================================
# REQUIRED - paths
# ============================================================================
PROJECT=/path/to/project

SUBDB_ROOT=${PROJECT}/5.subdb        # output of build_subdb + build_annotation
OUT_ROOT=${SUBDB_ROOT}               # same directory; the name build_subdb/5a use
MAPPING_OUT=${PROJECT}/6.mapping     # output of 4_run_mapping.sh
FEATURES_OUT=${PROJECT}/7.features   # output of 6_build_features.sh

# ============================================================================
# [host] host depletion - 3_host_depletion.sh
# ============================================================================
# Delete this section for a syncom or an environmental cohort.
# Reads go 2.trim -> 3.host_depletion/filtered, and every later step reads the
# filtered copy, so point the library DIRs below at $READS_DIR.
TRIM_DIR=${PROJECT}/2.trim
HF_OUT=${PROJECT}/3.host_depletion
READS_DIR=${HF_OUT}/filtered
HUMAN_IDX=/path/to/db/human_genome/GCA_000001405.29_GRCh38.p14_genomic   # bowtie2 index prefix
# Leave short-read RPF libraries out of it: at ~28 bp the host "hits" are mostly
# chance matches to repetitive regions (measured 11.28% vs metaG's 0.67% on the
# same sample, 78% of them multi-mapping). With this on, point METARS_DIR at the
# unfiltered reads. See README.
SKIP_METARS_FILTER=1
# Cohort size. Only an orchestrator needs it, to tell a finished host-depletion
# run from a half-finished one; it cannot count $TRIM_DIR, because
# drop_trimmed_metaG.sh removes those files as it verifies them.
N_SAMPLES=

# ============================================================================
# conda environments
# ============================================================================
# Only read when a step script is run directly with `bash <script>.sh`.
# Through the `mopp2` command every tool comes from the active (mopp2) env and
# these are ignored; pass `mopp2 --conda-switch ...` to honour them.
MAP_CONDA_ENV=base                   # bowtie2 (+ samtools on PATH)
FEAT_CONDA_ENV=features              # featureCounts + samtools + python3
# 5b additionally needs barrnap + bowtie2 + python3 on PATH.

# ============================================================================
# resources
# ============================================================================
THREADS=24
SORT_THREADS=4
SORT_MEM=2G
TMP_ROOT=                            # empty => a sort_tmp beside each output

# ============================================================================
# libraries - delete any section the project does not have.
# Reads are discovered as ${DIR}/<sample><SUFFIX>; the sample list comes from
# the sub-DB, so no sample file is needed.
#
# STRAND is the featureCounts -s value: 0 unstranded, 1 forward, 2 reverse.
# ============================================================================

# ---- metaG (paired) ----
METAG_DIR=${PROJECT}/2.trim
METAG_R1_SUFFIX=_metaG_R1.fq.gz
METAG_R2_SUFFIX=_metaG_R2.fq.gz
METAG_STRAND=0

# ---- metaT ----
# `mopp2 trim` writes <sample>_metaT_R1/_R2.fq.gz into 2.trim. rRNA removal
# (e.g. SortMeRNA) is not part of mopp2; point METAT_DIR at its clean reads, or
# at 2.trim with METAT_R1_SUFFIX=_metaT_R1.fq.gz to map the trimmed reads as-is.
METAT_DIR=${PROJECT}/3.sortmerna/clean_reads
METAT_PAIRED=1                       # 0 for single-end, then set METAT_SUFFIX
METAT_R1_SUFFIX=_metaT_clean_R1.fq.gz
METAT_R2_SUFFIX=_metaT_clean_R2.fq.gz
#METAT_SUFFIX=_metaT_clean.fq.gz     # used when METAT_PAIRED=0
METAT_STRAND=2

# ---- metaRS (single-end ribosome-protected fragments) ----
METARS_DIR=${PROJECT}/2.trim
METARS_SUFFIX=_metaRS.fq.gz
METARS_STRAND=0
METARS_BAM_SUBDIR=metaRS_mopp2       # which mapping variant featureCounts reads
#METARS_DB_OVERRIDE=                 # map metaRS against a different sub-DB sample

# ============================================================================
# [3c] discovering the genome pool from GlobDB - 3c_build_subdb_from_globdb.sh
# ============================================================================
# Only for 3c. 3a takes MAGS_POOL + GLOBDB_POOL instead; 3b takes SUBDB_SRC.
#
# Leave REFINEMENT_ROOT empty when the cohort has no MAGs of its own; set it to
# the metawrap refinement root to add them to the pool, which changes nothing
# else in the script.
REFINEMENT_ROOT=
BIN_SUBDIR=metawrap_50_10_bins

# The sketch to profile against, and where the round-1 table lands. STAGE 0
# builds the table if it is not already there; hand it an existing one from
# another run to skip straight to pool selection.
SYLPH_DB_R1=/path/to/db/SylDB/globdb_r226_sylph_c200.syldb
SYLPH_TSV_R1=${OUT_ROOT}/sylph_round1/results.tsv
# sylph reports Genome_file relative to the directory the sketch was built in,
# e.g. "globdb_r226_genome_fasta/R226CHUNK013/GCF_000742135.fa.gz". GLOBDB_ROOT
# must be the directory those paths are relative to - SylDB itself, NOT the
# genome_fasta directory inside it. Getting it wrong makes every genome miss.
GLOBDB_ROOT=/path/to/db/SylDB/

MIN_EFF_COV=0.1        # round 1: which genomes enter the pool at all
MIN_EFF_COV_R2=0.05    # round 2: which genomes enter one sample's sub-DB
ANI_DEDUP=99.0         # collapse near-duplicates above this
ANI_STRAIN_BAND=95.0   # keep 95-99% as separate strains
ALIGN_FRAC_MIN=30      # skani minimum alignment fraction (%)
SYLPH_MIN_ANI=95       # below 95 sylph profile is documented to be inaccurate
SYLPH_C=200

# [v4] Select by sylph's own confidence instead of the two Eff_cov cutoffs
# above. Empty keeps the cutoffs; T1..T4 replaces both, for the pool and for
# each sample:
#   T1  Eff_lambda HIGH (>= ~3x; ANI measured, not corrected)
#   T2  Eff_cov >= 1   and ANI 5% bound >= 95
#   T3  Eff_cov >= 0.1 and ANI 5% bound >= 95
#   T4  everything sylph reported
SYLPH_TIER_MAX=
# [v4] ref_ids kept out of every sub-DB, whitespace-separated - e.g. GlobDB
# genomes hit only by reads from a plasmid added to the samples.
EXCLUDE_REFS=
# [v4] The project's own isolate genomes, one per line: path[<TAB>ref_id]
# (ref_id defaults to the basename without extension). Like FUNGAL_FASTA below,
# each goes into a sample's sub-DB only where STAGE 0b calls it present
# (KNOWN_MIN_ANI_LO), and there a GlobDB genome of the same species (skani ANI
# >= ISOLATE_MIN_ANI) gives way to it.
ISOLATE_LIST=
ISOLATE_MIN_ANI=95
# [v4] Read sketches; STAGE 0 writes them and STAGE 0b profiles them again.
#READ_SKETCH_DIR=$(dirname ${SYLPH_TSV_R1})/sketches

# [v4] One fungal genome GlobDB does not contain. FUNGAL_REF_ID in the 5a
# section must be this file's basename without extension - 3c checks.
FUNGAL_FASTA=
# [v4] The project's own genomes (isolates, fungus) are profiled against
# themselves and called present per sample where the ANI 5% bound reaches this.
# Tiers do not apply to them (see README).
KNOWN_MIN_ANI_LO=98

# ============================================================================
# 5a_build_annotation.sh
# ============================================================================
EGGNOG_DATA_DIR=/path/to/db/emapperdb-5.0.2   # must hold eggnog.db + eggnog_proteins.dmnd
PARALLEL_PRODIGAL=24
EGGNOG_BLOCK_SIZE=2.0
EGGNOG_SENS_MODE=iterate
EGGNOG_TAX_SCOPE=auto

# SKIP_EGGNOG=1 stops 5a after prodigal. gene_inventory.tsv, cds_features.tsv
# and the SAFs are still written, which is everything 5b, 5c and 6 need - only
# the functional layer (gene_metadata.tsv) is left out, and a later rerun with
# SKIP_EGGNOG=0 adds it without recomputing anything else. Worth it whenever the
# pool is large: a GlobDB-derived gut pool can reach millions of proteins, where
# emapper runs for days and the count tables should not have to wait.
SKIP_EGGNOG=0

# A fungal genome annotated from a JGI GeneCatalog rather than by prodigal.
# Leave FUNGAL_REF_ID empty for an all-prokaryote pool: 5a then skips STAGE 2b
# and STAGE 4b, and every other stage already filters on source != "Fungus".
# With 3c, set FUNGAL_FASTA in the [3c] section as well; 3c puts the genome in
# the sub-DBs, 5a annotates it.
FUNGAL_REF_ID=
FUNGAL_GFF=
FUNGAL_PROT=

# ============================================================================
# 5b_make_structural_rna_blacklist.sh
# ============================================================================
# ORF nucleotide sequences. The default is fine unless this project REUSES
# ANOTHER PROJECT'S ANNOTATION, in which case point it at the source project -
# the blacklist belongs to the annotation, not to the project.
#ORF_FFN=${SUBDB_ROOT}/unified_orfs.ffn
#ORF_FFN=/path/to/source_project/5.subdb/unified_orfs.ffn

# Genome inventory and per-gene coordinates. Defaults are correct for a sub-DB
# built by this pipeline; override alongside ORF_FFN when reusing an annotation.
#UNIFIED_TSV=${SUBDB_ROOT}/unified_db.tsv
#GENE_INV=${SUBDB_ROOT}/gene_inventory.tsv

# bowtie2 index of the tRNA database.
TRNA_IDX=/path/to/db/tRNA_bt2/tRNAdb   # bowtie2 index of a tRNA set, e.g. GtRNAdb bacteria

# An ORF must lie at least this far inside an rRNA locus to be called spurious.
MIN_OVERLAP_FRAC=0.5

# Minimum aligned block against the tRNA database. DO NOT RELAX THIS. At ~20 bp
# the false-positive rate is ~98% and it removes >10k real genes.
MIN_TRNA_BLOCK=50

# ============================================================================
# 6_build_features.sh
# ============================================================================
# Where the blacklist lives. The default is where 5b writes it.
#ORF_BLACKLIST=${SUBDB_ROOT}/structural_rna_orf_blacklist.txt

# 6 refuses to run when the manifest has metaRS libraries and no blacklist is
# found, because 38-68% of metaRS CDS-assigned reads would be rRNA on spurious
# ORFs, unevenly across samples. Only set this if you know why you want that.
#ALLOW_UNFILTERED_METARS=1
