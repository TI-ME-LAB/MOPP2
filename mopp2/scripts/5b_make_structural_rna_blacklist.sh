#!/bin/bash
# 5b_make_structural_rna_blacklist.sh - identify prodigal ORFs that are actually
# rRNA / tRNA, so they can be kept out of the CDS feature set.
#
#   Usage: bash 5b_make_structural_rna_blacklist.sh <mapping.config.sh>
#
# PIPELINE POSITION
#   3.build_subdb.sh      sylph selects genomes, builds samples/<s>/ref.fna
#   5.build_annotation.sh prodigal -> gene_inventory.tsv, cds_features.tsv, unified_orfs.ffn
#   5b (this script)      -> structural_rna_orf_blacklist.txt
#   5c_gen_saf.sh  builds ref.saf, dropping the blacklisted gene_ids
#   6_build_features.sh     featureCounts + matrices
#
# It depends on the annotation, not on sylph. Whether an ORF sits inside a 16S
# gene is a property of the genome and its annotation; it does not change when
# sylph picks a different genome subset for a different project. So the scan
# covers every genome in unified_db.tsv - exactly the scope of gene_inventory.tsv
# and unified_orfs.ffn - and the result is reusable by any project sharing that
# annotation.
#
#   (The earlier version built its union from samples/*/ref.fna, i.e. the
#   sylph-selected subset. In 12.model that was 22 of 24 genomes, so ORFs on
#   ArtspOAP107_contigs and MarspOAE513 could never be flagged as rRNA even
#   though they were present in the SAF. It also made two projects sharing one
#   annotation disagree about which ORFs are real.)
#
# WHY THIS RUNS BEFORE THE SAF
#   metaRS (ribosome-protected fragments) and metaT (rRNA-depleted) both measure
#   mRNA. The CDS feature set therefore *means* "mRNA regions". A prodigal ORF
#   predicted inside a 16S/23S gene is not a contaminated count to subtract
#   later - it is a feature that should never have entered the SAF at all.
#   Filtering here means featureCounts never sees it and every downstream table
#   is correct as produced: no TPM/RPKM recomputation, no genome-level re-rollup.
#
#   It also matters numerically. featureCounts runs with `-O --fraction
#   -M --fraction`, so a read straddling a fake ORF and a real gene is split
#   0.5/0.5. Deleting the fake ORF's row afterwards discards that 0.5 instead of
#   returning it to the real gene; filtering the SAF up front does not.
#
# INPUTS
#   ${SUBDB_ROOT}/unified_db.tsv        ref_id -> genome FASTA path (plain or .gz)
#   ${SUBDB_ROOT}/gene_inventory.tsv    ORF coordinates (gene_id/ref_id/contig/start/end)
#   ${ORF_FFN}                          ORF nucleotide sequences, for tRNA alignment
#
# OUTPUTS
#   ${SUBDB_ROOT}/structural_rna_orf_blacklist.tsv   full table with evidence
#   ${SUBDB_ROOT}/structural_rna_orf_blacklist.txt   gene_id list (consumed by 5c)
#   ${SUBDB_ROOT}/_structural_rna_work/              intermediates, resumable
#
# NOTE: stderr of every long-running tool is discarded on purpose. tRNAscan-SE
# once wrote 166 GB of warnings and filled a 19 TB disk; barrnap and bowtie2 are
# throttled here for the same reason.
# NOTE: do NOT add `set -u` - conda's deactivate scripts dereference unset
# variables and would abort the run.
set -eo pipefail

CONFIG=${1:?usage: bash 5b_make_structural_rna_blacklist.sh <mapping.config.sh>}
source "$CONFIG"

: "${SUBDB_ROOT:?config is missing SUBDB_ROOT}"

# ORF nucleotide sequences. Defaults to the sub-DB being annotated. Projects that
# reuse another project's annotation (e.g. 13.dropout reusing 12.model) must set
# ORF_FFN explicitly in their config to point at the source project - and should
# point SUBDB_ROOT-derived inputs there too, since the blacklist belongs to the
# annotation, not to the project.
ORF_FFN=${ORF_FFN:-${SUBDB_ROOT}/unified_orfs.ffn}
GENE_INV=${GENE_INV:-${SUBDB_ROOT}/gene_inventory.tsv}
UNIFIED_TSV=${UNIFIED_TSV:-${SUBDB_ROOT}/unified_db.tsv}
TRNA_IDX=${TRNA_IDX:?config is missing TRNA_IDX (bowtie2 index of a tRNA database)}
THREADS=${THREADS:-24}

# An ORF must lie at least this far inside an rRNA locus to be called spurious.
MIN_OVERLAP_FRAC=${MIN_OVERLAP_FRAC:-0.5}
# Minimum aligned block against the tRNA database. Do not relax this: at ~20 bp
# the false-positive rate is ~98% and it kills >10k real genes.
MIN_TRNA_BLOCK=${MIN_TRNA_BLOCK:-50}

WORK=${SUBDB_ROOT}/_structural_rna_work
BLACKLIST_TSV=${SUBDB_ROOT}/structural_rna_orf_blacklist.tsv
BLACKLIST_TXT=${SUBDB_ROOT}/structural_rna_orf_blacklist.txt

mkdir -p "$WORK"
echo "[$(date '+%F %T')] config=$CONFIG"
echo "  SUBDB_ROOT=$SUBDB_ROOT"
echo "  UNIFIED_TSV=$UNIFIED_TSV"
echo "  GENE_INV=$GENE_INV"
echo "  ORF_FFN=$ORF_FFN"

[[ -s "$ORF_FFN"     ]] || { echo "ERROR: ORF sequences not found: $ORF_FFN"; echo "  (projects reusing another annotation must set ORF_FFN in the config)"; exit 1; }
[[ -s "$GENE_INV"    ]] || { echo "ERROR: gene inventory not found: $GENE_INV (run 5.build_annotation.sh first)"; exit 1; }
[[ -s "$UNIFIED_TSV" ]] || { echo "ERROR: unified_db.tsv not found: $UNIFIED_TSV"; exit 1; }

for t in barrnap bowtie2 python3; do
    command -v "$t" >/dev/null 2>&1 || { echo "ERROR: $t not found in PATH (env ${CONDA_DEFAULT_ENV:-none})"; exit 1; }
done

DF=$(df -BG --output=avail "$SUBDB_ROOT" | tail -1 | tr -dc 0-9)
echo "  free disk: ${DF}G"
if [[ "$DF" -lt 20 ]]; then echo "ERROR: need at least 20G free"; exit 1; fi

# ---- 1. union reference over the whole genome pool --------------------------
# Contigs are renamed <ref_id>__<original_contig>, the same convention
# build_subdb.sh uses, so they join directly against the contig column of
# gene_inventory.tsv.
UNION=${WORK}/union.fna
if [[ ! -s "$UNION" ]]; then
    echo "[1/4] building union over all genomes in unified_db.tsv ..."
    python3 - "$UNIFIED_TSV" "$UNION" <<'PY'
import csv, gzip, sys, os
unified, out = sys.argv[1], sys.argv[2]
n_ref = n_contig = 0; missing = []
with open(out, "w") as o:
    for r in csv.DictReader(open(unified), delimiter='\t'):
        rid, path = r['ref_id'], r['path']
        if not os.path.isfile(path):
            missing.append(rid); continue
        opener = gzip.open if path.endswith('.gz') else open
        with opener(path, 'rt') as fh:
            for line in fh:
                if line[0] == '>':
                    o.write(f">{rid}__{line[1:].split()[0]}\n"); n_contig += 1
                else:
                    o.write(line)
        n_ref += 1
print(f"      genomes: {n_ref}, contigs: {n_contig}")
if missing:
    print(f"      WARNING: {len(missing)} genome file(s) missing: {', '.join(missing[:5])}")
    sys.exit(1)
PY
else
    echo "[1/4] union.fna exists, skip"
fi

# ---- 2. rRNA loci: run both the bacterial and eukaryotic models, take the union
for kd in bac euk; do
    GFF=${WORK}/barrnap_${kd}.gff
    if [[ ! -s "$GFF" ]]; then
        echo "[2/4] barrnap --kingdom $kd ..."
        barrnap --kingdom "$kd" --threads "$THREADS" --quiet "$UNION" 2>/dev/null > "$GFF"
    fi
    echo "      $kd rRNA loci: $(grep -vc '^#' "$GFF" || true)"
done

# ---- 3. tRNA: align ORFs to the tRNA database, require a long aligned block --
# tRNAscan-SE is not used here: it is very slow on a large union and its stderr
# is a disk hazard. A strict bowtie2 alignment gives the same answer in minutes.
TRNA_HITS=${WORK}/trna_orf_strict.tsv
if [[ ! -s "$TRNA_HITS" ]]; then
    echo "[3/4] tRNA-derived ORFs (aligned block >= ${MIN_TRNA_BLOCK} bp) ..."
    bowtie2 -x "$TRNA_IDX" -f -U "$ORF_FFN" --local -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 \
        --score-min G,20,3 --no-unal --threads "$THREADS" 2>/dev/null \
    | awk -v m="$MIN_TRNA_BLOCK" '$1!~/^@/{
          n=0; c=$6; while(match(c,/[0-9]+M/)){ n+=substr(c,RSTART,RLENGTH-1); c=substr(c,RSTART+RLENGTH) }
          if(n>=m) print $1"\t"n }' | sort -u -k1,1 > "$TRNA_HITS"
fi
echo "      tRNA-derived ORFs: $(wc -l < "$TRNA_HITS")"

# ---- 4. intersect with ORF coordinates and emit the blacklist ---------------
echo "[4/4] writing blacklist ..."
python3 - "$WORK" "$GENE_INV" "$BLACKLIST_TSV" "$BLACKLIST_TXT" "$MIN_OVERLAP_FRAC" <<'PY'
import sys, csv
from collections import defaultdict
work, gene_inv, out_tsv, out_txt, minfrac = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5])

# rRNA loci per contig, sorted so the scan below can break early
loci=defaultdict(list)
for kd in ("bac","euk"):
    for line in open(f"{work}/barrnap_{kd}.gff"):
        if line.startswith('#'): continue
        p=line.rstrip('\n').split('\t')
        if len(p)>=5 and int(p[3])<int(p[4]):
            loci[p[0]].append((int(p[3])-1,int(p[4]),p[8].split(';')[0].replace('Name=','')))
for c in loci: loci[c].sort()

trna=set(l.split('\t')[0] for l in open(f"{work}/trna_orf_strict.tsv"))

# Every contig carrying an ORF should have been scanned by barrnap. If not, the
# union and the annotation are out of sync and the rRNA arm is silently blind.
seen_contigs=set()

hits=[]; n=0
for row in csv.DictReader(open(gene_inv), delimiter='\t'):
    n+=1; g=row['gene_id']; seen_contigs.add(row['contig'])
    s,e=int(row['start'])-1,int(row['end']); L=e-s
    best=0; bt=None
    for ls,le,t in loci.get(row['contig'],()):
        if le<=s: continue
        if ls>=e: break
        ov=min(e,le)-max(s,ls)
        if ov>best: best,bt=ov,t
    # rRNA takes precedence: positional evidence is stronger than a sequence
    # alignment for features this short.
    if L>0 and best/L>=minfrac: hits.append((g,row['ref_id'],'rRNA:'+str(bt),best,L,best/L))
    elif g in trna:             hits.append((g,row['ref_id'],'tRNA',0,L,0.0))

with open(out_tsv,"w") as f:
    f.write("gene_id\tref_id\trna_type\toverlap_bp\torf_len\tfrac\n")
    for h in sorted(hits): f.write("%s\t%s\t%s\t%d\t%d\t%.2f\n"%h)
with open(out_txt,"w") as f:
    for h in sorted(hits): f.write(h[0]+"\n")

print(f"      blacklisted: {len(hits)} / {n} ORFs ({len(hits)/n*100:.2f}%)")
by=defaultdict(int)
for h in hits: by[h[2].split(':')[0]]+=1
print("      by type: " + ", ".join(f"{k}={v}" for k,v in sorted(by.items())))

# Scope check: the two arms of the evidence must cover the same genomes.
scanned=set()
for kd in ("bac","euk"):
    for line in open(f"{work}/barrnap_{kd}.gff"):
        if not line.startswith('#'): scanned.add(line.split('\t')[0])
refs_with_orfs={c.split('__',1)[0] for c in seen_contigs}
refs_scanned={c.split('__',1)[0] for c in scanned}
blind=refs_with_orfs-refs_scanned
if blind:
    print(f"      NOTE: {len(blind)} genome(s) have ORFs but no rRNA hit at all: {', '.join(sorted(blind)[:5])}")
    print(f"            (fine if they genuinely carry no rRNA; suspicious if many)")

# Sanity bound. A complete genome carries ~25 rRNA loci (7-8 operons x 16S/23S/5S);
# a MAG usually 1-3 because multi-copy regions collapse during assembly. The
# expected blacklist size is hundreds, not tens of thousands.
if len(hits) > n*0.02:
    print(f"      WARNING: >2% of ORFs blacklisted - inspect the barrnap output before trusting this.")
PY

echo
echo "[$(date '+%F %T')] done."
echo "  $BLACKLIST_TSV"
echo "  $BLACKLIST_TXT"
echo "  next: bash 5c_gen_saf.sh --subdb $SUBDB_ROOT --blacklist $BLACKLIST_TXT"
echo "  ($WORK/union.fna can be deleted; it rebuilds in minutes)"
