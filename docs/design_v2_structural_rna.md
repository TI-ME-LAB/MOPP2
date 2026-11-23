# mopp2_v2

All-English rework of the mopp2 chain. The structural-RNA ORF blacklist now comes
from the annotation instead of from the finished feature tables, which removes a
circular dependency that could block a run at step 6.

The scripts in the mopp2 root are untouched. Nothing here has been run yet.

## Pipeline

```
3a_build_subdb.sh           sylph profiling -> per-sample genome set + bowtie2 indexes
3b_build_subdb_round2.sh    second-round sub-DB build
        |
4_run_mapping.sh            bowtie2 -> per-sample BAMs        [parallel with 5*]
        |
5a_build_annotation.sh      prodigal / JGI / eggNOG
        |                   -> gene_inventory.tsv, unified_orfs.ffn
        |                   -> samples/*/ref.saf   (PROVISIONAL)
5b_make_structural_rna_blacklist.sh
        |                   -> structural_rna_orf_blacklist.{tsv,txt}
5c_gen_saf.sh               -> samples/*/ref.saf   (spurious ORFs removed)
        |
6_build_features.sh         featureCounts + count/TPM/RPKM matrices
                            -> $FEATURES_OUT/tables/   (final, no cleanup follows)
```

5a/5b/5c are three sub-steps of one stage: all three write into `5.subdb/`.

Every script stands alone — there is no shared library. `gen_manifest`, which
builds the manifest of (sample, library, BAM) rows that drives counting, is
duplicated verbatim in `4_run_mapping.sh` and `6_build_features.sh`. Step 6
cannot rely on step 4 having written `manifest.tsv`: the older `4.run_mapping.sh`
never wrote one, so a run resumed from those BAMs arrives at step 6 with no
manifest. **If you change one copy, change the other.**

## If you are stuck at step 6

You ran build_subdb, build_annotation and mapping, and now `build_features` wants
a structural-RNA blacklist you do not have. The old `7.clean_structural_rna.sh`
could not give you one: it read `${FEATURES_OUT}/tables/gene_annotation.tsv`,
which does not exist until the feature tables are built. The blacklist needed the
tables and the tables needed the blacklist.

That circle is gone. `5b` derives the blacklist from the annotation alone
(`unified_db.tsv` + `gene_inventory.tsv` + `unified_orfs.ffn`), all of which you
already have. Three commands, same interface:

```bash
cp config.template.sh myproject.config.sh    # fill in PROJECT and the read paths
bash 5b_make_structural_rna_blacklist.sh myproject.config.sh
bash 5c_gen_saf.sh                       myproject.config.sh
bash 6_build_features.sh                 myproject.config.sh
```

Nothing from steps 3–5 is recomputed and your BAMs are not touched.

Two traps these scripts handle for you:

- **`cds_features.tsv` may not exist.** Only `12model.build_annotation.sh` writes
  it; the older `5.build_annotation.sh` builds `ref.saf` straight from
  `gene_inventory.tsv`. `5c` detects this and switches to filtering your existing
  SAFs in place, which also preserves multi-exon fungal features exactly as
  annotated. Originals are kept as `ref.saf.unfiltered`, and re-running always
  filters from that copy, so it never compounds.
- **featureCounts output is cached.** If `fc_per_sample/` holds results from a run
  against the unfiltered SAF, step 6 would skip those samples and the spurious
  ORFs would survive in your tables. Step 6 compares mtimes and recounts anything
  older than its SAF.

## Configuration

Copy `config.template.sh` per project. Required: `PROJECT`, `SUBDB_ROOT`,
`MAPPING_OUT`, `FEATURES_OUT`, and one section per library. Everything else has a
working default.

| Variable | Used by | Default |
|---|---|---|
| `ORF_FFN` | 5b | `$SUBDB_ROOT/unified_orfs.ffn` |
| `UNIFIED_TSV`, `GENE_INV` | 5b | `$SUBDB_ROOT/{unified_db,gene_inventory}.tsv` |
| `ORF_BLACKLIST` | 5c, 6 | `$SUBDB_ROOT/structural_rna_orf_blacklist.txt` |
| `TRNA_IDX` | 5b | `/media/ExtraDisk/db/tRNA_bt2/tRNAdb` |
| `MIN_OVERLAP_FRAC` | 5b | `0.5` |
| `MIN_TRNA_BLOCK` | 5b | `50` — do not relax, see below |
| `METARS_BAM_SUBDIR` | 6 | `metaRS_mopp2` |
| `ALLOW_UNFILTERED_METARS` | 6 | unset — escape hatch |

A project **reusing another project's annotation** sets `ORF_FFN` (and usually
`UNIFIED_TSV` / `GENE_INV`) to the source project: the blacklist belongs to the
annotation, not to the project.

---

# Why any of this exists

Prodigal occasionally predicts a CDS **inside an rRNA or tRNA gene**. rRNA is
present at very high copy number, so in metaRS (ribosome profiling, very short
single-end reads) the rRNA reads pile onto these spurious ORFs and featureCounts
records them as translation signal.

Measured share of CDS-assigned metaRS reads landing on spurious ORFs:

| Project | share | per-sample range |
|---|---|---|
| 12.model `7.features_noSpike` | 38.8% | 5% – 69% |
| 13.dropout `7.features` | 62.0% | 33% – 89% |
| 13.dropout `7.features_mags16` | 66.1% | 19% – 81% |
| 13.dropout `7.features_mags16b` | 67.8% | 46% – 84% |

**The damage is not uniform inflation — it is that each sample is inflated by a
different amount** (19% – 89%). Compare samples directly and this purely
technical bias reads as biology. In 12.model the metaRS TPM of *Klebsiella
michiganensis* moves by −30% to **+297%** once corrected; S11 nearly quadruples.

metaT (0.1%) and metaG (0.1%) are barely affected in magnitude, but all three
libraries must use the **same** feature set: a spurious ORF is not a gene, and
leaving it in metaG while removing it from metaRS would put the numerator and
denominator of TE on different gene sets.

## Why the blacklist comes before the SAF

metaRS measures ribosome-protected fragments and metaT is rRNA-depleted, so both
quantify **mRNA**. The CDS feature set therefore *means* "mRNA regions". An ORF
predicted inside a 16S gene is not a contaminated count to subtract afterwards —
it is a feature that should never have been in the SAF.

1. **Every table is correct as produced.** The old cleanup pass had to recompute
   TPM/RPKM (the original denominators included the spurious ORFs — each column
   summing to exactly 1e6 proves it) and re-aggregate the genome level from the
   gene level. That whole error-prone step is gone.
2. **The numbers differ slightly.** featureCounts runs with `-O --fraction
   -M --fraction`, so a read straddling a spurious ORF and a real gene
   contributes 0.5 to each. Deleting the spurious row afterwards discards that
   0.5; filtering the SAF returns it to the real gene. In 13.dropout, 8 of 102
   blacklisted ORFs overlap a real gene (3–181 bp) — small but real. The `-M`
   path is not quantified.

Expect metaRS `pct_assigned` in `counting_summary.tsv` to **drop**. Reads
formerly Assigned to spurious ORFs now land in `Unassigned_NoFeatures` — the
correct classification: they are rRNA, not translation.

## Why the blacklist does not depend on sylph

Whether an ORF sits inside a 16S gene is a property of the genome and its
annotation, not of which genomes sylph selected.

The old script built its barrnap input from `samples/*/ref.fna`, the
sylph-selected subset. In 12.model that was 22 of 24 genomes, leaving **129
ORF-bearing contigs never scanned** — ORFs on `ArtspOAP107_contigs` and
`MarspOAE513_scaffolds.trim` could never be flagged even though they were in the
SAF. The tRNA arm meanwhile scanned `unified_orfs.ffn`, which *is* pool-wide, so
the two halves of the evidence had different scope. Two projects sharing one
annotation therefore disagreed: 12.model 133 entries, 13.dropout 102, with 10
unique to 13.dropout.

`5b` now scans every genome in `unified_db.tsv`, matching the scope of
`gene_inventory.tsv` and `unified_orfs.ffn` exactly.

## Sanity check on the blacklist

A complete genome such as Klebsiella should have ~25 rRNA loci (7–8 operons ×
16S/23S/5S); a MAG usually 1–3, because multi-copy regions collapse during
assembly. **The blacklist should be hundreds of entries, not tens of thousands.**
5b warns if it exceeds 2% of all ORFs.

## Things not to change, and why

**`MIN_TRNA_BLOCK=50`.** The first attempt filtered ORFs against SILVA with
bowtie2 `--local --score-min G,20,3`, which for a 255 bp sequence only demands
about 19 bp of perfect match. Of 10,644 hits, **83% had an aligned block of only
20–29 bp**, the false-positive rate was ~98%, and it removed 10,715 real genes.

> When filtering by homology on alignment score, check aligned block length and
> coverage. "Did it hit" is not enough.

**Do not run tRNAscan-SE on the union.** Very slow, and its stderr once wrote
**166 GB** of warnings in an hour and filled a 19 TB disk. Every long-running
tool in 5b has stderr discarded for that reason.

**Do not add `set -u` to 5b.** conda's deactivate scripts dereference unset
variables and would abort the run.

**Do not bowtie2 against SortMeRNA's SILVA FASTA.**
`/media/ExtraDisk/db/sortmerna_db/smr_v4.3_*.fasta` uses an RNA alphabet (U, not
T); U is 19.9% of it, and bowtie2 treats U as an ambiguous non-ACGT base, so ~20%
of every reference is dead. Same 200 reads: blastn (megablast) 96% hit rate,
bowtie2 `--local` 4%, bowtie2 `--very-sensitive-local` 4%.
Detect it with `bowtie2-inspect -s <idx> | awk '/^Sequence-/{n+=$NF} END{print n}'`
— 103,610,311 against a FASTA holding 133,533,894 bases; the 22% gap is the U.
Convert with `awk '/^>/{print;next}{gsub(/U/,"T");print}'` before indexing.

> When handed someone else's reference database, check the alphabet first. An
> index base count noticeably below the FASTA base count is the signal.

## Why there is no sequence-level rRNA depletion step

Annotation-based removal is sufficient and more reliable. Every read from S01/S06
mapping to the bacterial sub-DB was classified as rRNA or not, and cross-tabulated
by where it landed:

| Landing site | share in S01 | fate of the rRNA reads |
|---|---|---|
| rRNA locus, no ORF | 72.1% | NoFeatures; never enters counts |
| intergenic | 14.8% | same, not counted |
| real CDS | 8.7% | only 54 rRNA reads leak in (0.2% of all rRNA) |
| spurious ORF | 4.4% | the only contamination route — closed by the blacklist |

About 30% of reads mapping to bacteria are rRNA, but annotation catches 88–89% of
them and only 0.1–0.2% leak into real CDS. **The single contamination route into
the count table is the spurious ORFs.**

> Coordinates beat sequences. Of the 82,162 reads landing on rRNA loci, SILVA
> recognised only 28,412 (35%) — the rest are too short to identify. But the fact
> that they lie inside a 16S gene is certain. For very short reads, "where it
> landed" is not limited by alignment sensitivity.

## What the blacklist does not cover

It only handles rRNA/tRNA loci. Amplification artefacts and structural RNA piling
onto **ordinary CDS** — such as ampE in K. michiganensis: 36k reads, 96% of them
on one position with a single sequence — are invisible to barrnap and survive
this pipeline intact. Screening for those is a separate gene-level QC step on the
finished tables (`8.footprint_qc.py` in the mopp2 root), not part of this chain.


## Origin of each file

| File | Origin |
|---|---|
| `3a_build_subdb.sh` | `12model.build_subdb.sh`, unchanged |
| `3b_build_subdb_round2.sh` | `build_subdb_round2.sh`, unchanged |
| `4_run_mapping.sh` | `run_mapping2.sh`, unchanged |
| `5a_build_annotation.sh` | `12model.build_annotation.sh` + note that its SAF is provisional |
| `5b_make_structural_rna_blacklist.sh` | new; steps 1–4 of `7.clean_structural_rna.sh`, rescoped pool-wide |
| `5c_gen_saf.sh` | `gen_saf_from_cds.sh` + blacklist, filter mode, provenance |
| `6_build_features.sh` | `build_features2.sh` + blacklist handling and guards |
| `config.template.sh` | new; general replacement for the per-project configs |
| — | `manifest_lib.sh` inlined into steps 4 and 6 |

`7.clean_structural_rna.sh` has no counterpart: steps 1–4 became `5b`, and step 5
(rebuilding clean tables) is unnecessary once the SAF is filtered.

## Status

Nothing here has been run, and no `structural_rna_orf_blacklist.txt` exists under
any `5.subdb*/` yet. The two blacklists on disk
(`13.dropout/7.features*_clean/` and `12.model/8.TE/*_v2`) came from the older
sylph-scoped code, disagree with each other, and should not be reused.
