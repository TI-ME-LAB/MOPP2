# mopp2_v3

v2 plus what the RA cohort needed. v2's README is still the reference for *why*
the chain is shaped the way it is — the structural-RNA blacklist, the ordering of
5b/5c/6, the traps in 5b — none of which changed. This file covers only what v3
adds.

Every addition is off by default: with a v2 config, v3 runs exactly as v2 did.

## Pipeline

```
3_host_depletion.sh          bowtie2 vs host, keep the unaligned    [new, optional]
drop_trimmed_metaG.sh        delete trimmed metaG once verified     [new, utility]
        |
3a_build_subdb.sh            pool already curated by hand      (a syncom)
3b_build_subdb_round2.sh     reuse another project's sketch    (a dropout set)
3c_build_subdb_from_globdb.sh  discover the pool from GlobDB   [new]
        |
4_run_mapping.sh             bowtie2 -> per-sample BAMs      [parallel with 5a]
        |
5a_build_annotation.sh       prodigal / JGI / eggNOG
5b_make_structural_rna_blacklist.sh
5c_gen_saf.sh
        |
6_build_features.sh          featureCounts + count/TPM/RPKM matrices
```

Stage 3 now has three entry points. Pick by what you already have:

| | you have | script |
|---|---|---|
| a syncom | two hand-curated pool directories | `3a` |
| a dropout set | another project's `unified_db.syldb` | `3b` |
| a real community | only reads | **`3c`** |

## What v3 adds

### 3c: discovering the genome pool from the reads

`3a` needs `MAGS_POOL` + `GLOBDB_POOL`, which only a syncom has. mopp2's original
`3.build_subdb.sh` instead wanted cohort MAGs plus a round-1 sylph table computed
somewhere else. A cohort that has been sequenced and nothing more has neither.

`3c` is that original script with two changes:

- **STAGE 0 (new)** sketches every metaG sample and profiles all the sketches
  against a pre-built GlobDB/GTDB sketch in **one** `sylph profile` call, writing
  the round-1 table the rest of the script expects. One call, not one per sample:
  the GlobDB r226 sketch is 34 GB, and profiling sample by sample would load it
  once per sample. STAGE 2 reads the table as a union over samples, so they do
  not need to stay separate. Loading it costs ~40 GB of RAM — check the machine
  before starting. Skipped entirely if `$SYLPH_TSV_R1` already exists.
- **STAGE 1 is a no-op** when `REFINEMENT_ROOT` is empty. `mag_inventory.tsv` is
  written header-only and every later stage already handled a MAG-free pool.
  Setting `REFINEMENT_ROOT` turns the MAG arm back on with no other change, so a
  cohort assembled later does not need a different script.

- **The GlobDB lookup no longer scans per genome.** sylph reports `Genome_file`
  as a path relative to the directory its sketch was built in, e.g.
  `globdb_r226_genome_fasta/R226CHUNK013/GCF_000742135.fa.gz`, so `GLOBDB_ROOT`
  must be the directory those paths are relative to — `SylDB` itself, not the
  `globdb_r226_genome_fasta` inside it. Point it one level too deep and every
  genome misses, and the original fell back to one `find` over the whole tree
  **per genome**: GlobDB r226 holds 306,260 files, so a few thousand genomes
  would have taken days. It now resolves directly, and only if something really
  is missing does it index the tree once and join on basename. It also **fails**
  when fewer genomes were linked than selected, instead of carrying on with a
  quietly smaller pool.

- **No second sylph round when the pool has no MAGs.** Round 2 exists because a
  cohort's own MAGs are not in GlobDB: round 1 cannot say which sample holds
  which MAG, so the merged MAG+GlobDB pool has to be profiled again. With no
  MAGs there is no new database to profile against — the pool *is* the round-1
  hits — so STAGE 7 now splits the round-1 table per sample instead, and only
  re-profiles when `REFINEMENT_ROOT` is set.

  Measured on a 48-sample GlobDB profile (10_SI) to check what that gives up:

  | | |
  |---|---|
  | genomes detected in round 1 | 1163 |
  | pairs at >=99% ANI | **2** |
  | pool after dedup | 1161 |

  So round 2 would have profiled against a database differing from round 1 by two
  genomes. `sylph profile` is already non-redundant within a sample — no genome
  appears twice — which is the property a second round would otherwise supply.

  Dedup itself still earns its place, for a reason that is about one sample
  rather than across samples. In that dataset one sample held **all three**
  members of a >=99% cluster at once (Adjusted_ANI 97.6–98.5, Eff_cov
  0.16–0.24 — the real organism is not in GlobDB, so sylph spread weak evidence
  over several distant relatives). Three near-identical genomes in one sub-DB
  means bowtie2 multi-maps, and featureCounts runs `-M --fraction`, so one
  organism's reads get split three ways across three sets of gene ids. STAGE 7
  therefore rewrites collapsed ids to the representative STAGE 5 kept and merges
  the duplicate rows, keeping the better-covered one and summing abundance.

  `MIN_EFF_COV_R2` still applies; set it to 0 to take every row sylph reported
  (its own floor is around Eff_cov 0.02).

  Output moved from `sylph_round2/` to `sample_references/`, since in the no-MAG
  path nothing is profiled twice.

Nothing else in it is touched.

### 5a: the fungal genome is optional

12.model carried one (*Trichoderma harzianum*, gene models from JGI rather than
prodigal) and 5a required `FUNGAL_REF_ID`/`FUNGAL_GFF`/`FUNGAL_PROT`. They are
now optional: an empty `FUNGAL_REF_ID` skips STAGE 2b (inject JGI gene models)
and STAGE 4b (JGI KEGG/GO/KOG/IPR supplement). Every other stage already filtered
on `source != "Fungus"` and was left alone.

### 5a: SKIP_EGGNOG

Step 6 reads `gene_inventory.tsv`, `unified_db.tsv`, the SAFs and the blacklist.
It never reads eggNOG output — that only produces `gene_metadata.tsv`, the
functional layer. So the count tables do not have to wait for emapper, and on a
pool discovered from GlobDB they should not: a human gut pool can reach millions
of proteins, where emapper runs for days.

`SKIP_EGGNOG=1` stops 5a after prodigal; 5b, 5c and 6 then run to completion.
Rerunning 5a later with `SKIP_EGGNOG=0` adds the functional layer and recomputes
nothing else — STAGE 3 skips itself when the annotations already exist.

### 3_host_depletion.sh

For stool, skin or plant tissue, every downstream number — sylph abundance,
bowtie2 mapping rate — is diluted by host reads by a per-sample amount. bowtie2
against `$HUMAN_IDX`; reads that fail to align are what continue.

Leave the step out for a syncom or an environmental cohort.

**`SKIP_METARS_FILTER=1` leaves metaRS alone, and for short RPF reads that is the
right default.** Measured on one RA stool sample against GRCh38.p14:

| library | read length | "human" | aligned 1x | aligned >1x |
|---|---|---|---|---|
| metaG | 85 bp | 0.67% | — | — |
| metaRS | 27.8 bp | **11.28%** | 2.51% | **8.77%** |

One sample cannot be 0.67% host in its DNA and 11.28% host in its RPFs, and a
real host read would align uniquely rather than 78%-multi-mapping. The 28 bp
reads are landing in repetitive regions of a 3.1 Gb genome by chance; filtering
would have thrown away 8.7M real bacterial reads from that sample alone. Point
`METARS_DIR` at the unfiltered reads when this is on, and do not delete them.

This is the third time the same mistake has been available in this codebase:

| | permissive match | false positives |
|---|---|---|
| ORFs vs SILVA, `--score-min G,20,3` | ~20 bp block | ~98%, removed 10,715 real genes |
| bowtie2 vs SortMeRNA SILVA FASTA | RNA alphabet | 4% hit rate vs blastn's 96% |
| metaRS vs GRCh38, 28 bp | whole read | 17x metaG's rate, 78% multi-mapping |

> Short reads plus a permissive alignment produce hits that mean nothing. Check
> aligned block length, uniqueness, and whether the number is plausible against
> an independent measurement of the same thing.

Host depletion for *deposition* is a different problem: the requirement applies
to every library, but it cannot be met this way for 28 bp data. That belongs on
the path to submission, not in this chain.

### drop_trimmed_metaG.sh

A disk-pressure utility, not a pipeline stage. It deletes a sample's trimmed
metaG once the host-depleted copy is verified: the bowtie2 log reached
`overall alignment rate`, both outputs pass `gzip -t`, and R1/R2 hold equally
many reads. `--watch` polls, so space is freed while host depletion is still
running. Deletions are logged to `dropped_trimmed_metaG.tsv`, and what it removes
is reproducible from the raw reads by re-running trim-galore.

Use it only when there is no room to keep both copies.

## Origin of each file

| File | Origin |
|---|---|
| `3_host_depletion.sh` | new (RA) |
| `3a_build_subdb.sh` | v2, unchanged |
| `3b_build_subdb_round2.sh` | v2, unchanged |
| `3c_build_subdb_from_globdb.sh` | mopp2 `3.build_subdb.sh` + STAGE 0, MAG arm optional |
| `4_run_mapping.sh` | v2, unchanged |
| `5a_build_annotation.sh` | v2 + optional fungal stages + `SKIP_EGGNOG` |
| `5b_make_structural_rna_blacklist.sh` | v2, unchanged |
| `5c_gen_saf.sh` | v2, unchanged |
| `6_build_features.sh` | v2, unchanged |
| `drop_trimmed_metaG.sh` | new (RA) |
| `config.template.sh` | v2 + the `[host]`, `[3c]` and 5a sections |

## Status

v3 was assembled from the RA run (`/media/ExtraDisk/ongoing/RA/scripts/`), which
is where these scripts are actually being exercised — it is in progress, so
nothing here has been through a complete end-to-end run yet. Two filenames differ
there for historical reasons: RA calls `3c_build_subdb_from_globdb.sh`
`3_build_subdb.sh`, `3_host_depletion.sh` `3_human_filter.sh`, and
`drop_trimmed_metaG.sh` `3b_drop_trimmed_metaG.sh`. The contents are the same
apart from the header comments. RA also carries two orchestrators that are
specific to it and are not part of v3: `run_all_RA.sh` and `chain_subdb.sh`.
