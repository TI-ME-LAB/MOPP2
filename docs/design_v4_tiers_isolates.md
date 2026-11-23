# mopp2_v4

v3 plus what Exp3_as (FunSyn, soil metaG + metaRS, a Trichoderma inoculant and
a spiked-in plasmid) needed. Only `3c_build_subdb_from_globdb.sh` and
`config.template.sh` changed; every other script is v3's. v3's README
(`../mopp2_v3/README.md`) still covers 3c's design, host depletion and SKIP_EGGNOG,
and v2's covers the structural-RNA chain.

Every addition is off by default: with a v3 config, v4's 3c selects the same
genomes v3 did. The fixes in the last section change only what a rerun does.

## 3c: selecting by sylph's confidence

v3 chose genomes with two Eff_cov cutoffs, 0.1 for the pool and 0.05 per sample.
On Exp3_as that kept 265 of the 270 round-1 rows: in practice everything sylph
reported. sylph itself says how sure it is, and `SYLPH_TIER_MAX` uses that
instead:

| tier | rule | Exp3_as rows |
|---|---|---|
| T1 | `Eff_lambda == HIGH` (>= ~3x, ANI measured, not coverage-corrected) and ANI >= 95 | 102 |
| T2 | `Eff_cov >= 1` and ANI 5% bound >= 95 | 32 |
| T3 | `0.1 <= Eff_cov < 1` and ANI 5% bound >= 95 | 109 |
| T4 | everything else sylph reported (Eff_cov < 0.1, or the 5% bound < 95) | 27 |

A tier replaces both cutoffs and applies to the pool (STAGE 3) and to each
sample (STAGE 7) alike. What the tiers buy on Exp3_as:

| | species per sample | sylph sequence abundance kept |
|---|---|---|
| T1 | 10-15 | 89-98% |
| <= T2 | 13-20 | 96-99% |
| <= T3 | 23-39 | 99.8-100% |

T3 nearly doubles the species count for 1-3% more abundance, and that is where
low-coverage traces from other samples land (27 rows at <= 2% of the same
species' coverage elsewhere, almost all T3/T4). Exp3_as uses T2.

In practice <= T2 is an Eff_cov >= 1x cutoff: on Exp3_as the two select the
same 134 rows. The ANI interval only removed 11 rows, all below 0.62x. It
matters below 1x, which is where the project's own genomes are judged (below).

The tier rule is written twice, in STAGE 2&3 and STAGE 7, so each heredoc stands
alone. Change both together. T1's ANI check never fires on a `sylph profile`
table, whose `--minimum-ani 95` already removed such rows. A `sylph query` table
keeps them down to 90: in Exp3_as, a related *Arthrobacter* at 7.6x showed up
there at HIGH with 91.2% ANI.

## 3c: EXCLUDE_REFS

ref_ids that never enter any sub-DB. Exp3_as needed it for two GlobDB
*Escherichia* genomes that every sample hit at 0.7-2.4x and 95.6/97.8% ANI. They
were fragments of a plasmid added to the samples, not an organism. Left in, the
plasmid's reads would be counted as *E. coli* genes. The two references are
95.3% ANI to each other, so the 99% dedup would not have merged them either.

## 3c: the project's own genomes (STAGE 0b, 3b)

A SynCom's members are sequenced genomes, and GlobDB holds at best their
relatives (and no fungi at all). `ISOLATE_LIST` names the isolate genomes and
`FUNGAL_FASTA` a fungus. STAGE 0b profiles all of them together against the read
sketches and calls each one present or absent per sample. Present genomes go
into that sample's sub-DB. Isolates as source `Isolate`; the fungus as `Fungus`,
which 5a's JGI gene-model stage (`FUNGAL_REF_ID`/`FUNGAL_GFF`/`FUNGAL_PROT`,
added in v3 for 12.model) annotates with no further change.

**Why not put every member in every sub-DB.** Exp3_as gave the live SynCom to
Group_4 and an autoclaved one to Group_8, so it is tempting to give all eight
samples the full SynCom. The sub-DBs are also used for metaRS, and ~28 bp reads
find a place to land on a genome that is not there. A genome only enters a
sample that carries it.

**Why not by tier.** The tiers gate on Eff_cov, which measures depth, not
certainty. A 40 Mb genome gives far more k-mers than a 4 Mb one at the same
depth, and a known genome is checked against itself, not a relative. Presence
is decided on the ANI 5% bound (`KNOWN_MIN_ANI_LO`, default 98; v4's first
draft called it `FUNGAL_MIN_ANI_LO`, which is still read). For the fungus:

| Exp3_as sample | Eff_cov | k-mers matched | ANI (5-95%) | present |
|---|---|---|---|---|
| Group_8_2 | 11.1 | 168,737 | 99.92 (HIGH) | yes |
| Group_4_1 | 1.31 | 120,629 | 99.79-99.83 | yes |
| Group_4_2 | 0.15 | 19,701 | 99.11-99.37 | yes |
| Group_8_3 | 0.27 | 30,850 | 98.93-99.11 | yes |
| Group_8_4 | 0.07 | 3,660 | 95.66-96.63 | no |

At 0.15x, Group_4_2 is a T3 by the bacterial rule, yet it matched ~20k k-mers
with a 0.3-point interval at 99.2: no doubt about presence. What separates
present from absent is whether as many k-mers turn up as that exact genome would
give at that depth. sylph's coverage-adjusted ANI reports exactly that.
Group_8_4 found about a third of the expected k-mers: something related was
there, but not this strain. The isolates behave the same way. Under T2,
*Brevibacillus* OAP136, *Burkholderia* OAS925, *Marmoricola* OAE513 and
*Niastella* OAS944 were left out of Group_4, where they sat at 98.1-99.7% (5%
bound) but below 1x.

**Same-species GlobDB genomes (STAGE 3b).** Where a sample carries an isolate,
a GlobDB genome of the same species (skani ANI >= `ISOLATE_MIN_ANI`, default
95) that the selection also picked there is a stand-in for it and is dropped
from that sample:

| Exp3_as, Group_4 | sample vs GlobDB stand-in | sample vs isolate |
|---|---|---|
| *Chitinophaga* OAE865 (for *C. polysaccharea_A*) | 97.6 | 100.00 |
| *Lysobacter* OAE881 (for *Lysobacter_J soli*) | 97.5 | 100.00 |
| *Bradyrhizobium* OAE829 (for sp018130165, 95.35% to it) | 96.2-96.4 | 99.65-99.96 |

Where the isolate is not called present, the GlobDB genome is a background strain
like any other and stays on the tier rule. GlobDB genomes left with no sample
leave the pool, and so do isolates no sample carries.

Details that matter if this is extended:

- **Isolates are never collapsed in dedup.** Each was called present on its own,
  and a >= 99% GlobDB relative still in the pool is there for samples without
  the isolate.
- **Their rows come from STAGE 0b, never from the table being split.** With MAGs
  in the pool, round 2 profiles the isolates too (they are in
  `unified_db.syldb`, which keeps their k-mers from being reassigned to
  relatives), but presence is still STAGE 0b's call. The fungus stays out of
  `unified_db.syldb` and joins `unified_db.tsv` after dedup.
- **Their `Taxonomic_abundance` is `NA`.** They were profiled without the
  background, so sylph's figures mean nothing next to the GlobDB rows'.
- **`FUNGAL_REF_ID` must be the FASTA's basename without extension.** STAGE 8
  names contigs `<ref_id>__<scaffold>`, and 5a writes gene coordinates against
  the same name. 3c refuses to start if the two disagree. If no sample carries
  the fungus, 3c leaves it out and warns. Clear `FUNGAL_REF_ID` before 5a then,
  or 5a injects gene models for a genome no sub-DB contains.
- `ISOLATE_LIST` lines are `path[<TAB>ref_id]`. The ref_id defaults to the
  basename without extension. Isolates enter `unified_db.tsv` as source
  `Isolate`. 5a treats every source except `Fungus` as a prokaryote for
  prodigal, so nothing downstream changes.
- Read sketches live in `READ_SKETCH_DIR` (default: `sketches/` beside
  `SYLPH_TSV_R1`), so STAGE 0b reuses STAGE 0's instead of re-sketching the reads.

## 3c: fixes that change reruns

Tiers make it natural to rerun with a different selection, and three v3
behaviours broke on exactly that:

- **`globdb_pool/` kept every link from earlier runs.** STAGE 5 takes the whole
  directory, so genomes from a looser selection stayed in, and the
  linked-vs-selected count check then failed. It is now relinked from scratch.
- **`skani_all.tsv` and `unified_db.syldb` were reused whenever they existed.**
  After the pool grew, the new genomes had no skani pairs and escaped dedup. Each
  is now reused only when its input list is byte-identical to the last run's.
- **STAGE 8 re-applied `MIN_EFF_COV_R2`.** That was redundant in v3 (STAGE 7 had
  applied the same cutoff) but would undo a tier. It is gone, and genomes that
  are selected but missing from `unified_db.tsv` are now counted, not dropped
  silently.

Two more that affect a first run as well:

- **STAGE 9's summary was empty.** `build_one_sample` changed into the sample
  directory in the main shell, so STAGE 9's relative paths missed and
  `summary.tsv` was written into the last sample. Each build now runs in a
  subshell.
- **3c defaulted to a conda env** that has no sylph, skani or samtools,
  so its tool check could not pass. It now defaults to `base`, like 3a and 5a.

## Status

Exercised on Exp3_as (`.../FunSyn.com/Exp3_as/as/exp3_as.config.sh`), no-MAG path.
The MAG path (round 2) has the v4 changes but has not been run with them.
