#!/bin/bash
# 5a_build_annotation.sh
#
# Functional annotation of the 12.model unified sub-DB.
#
#   * 22 prokaryotic genomes: prodigal -p meta  (ORF prediction)
#   * 1 fungal genome (Trichoderma harzianum M10): NOT prodigal.
#       Gene models + proteins taken from the JGI MycoCosm GeneCatalog.
#   * All proteins (bacteria + fungus) -> eggNOG-mapper (uniform schema).
#   * JGI native KEGG/GO/KOG/InterPro kept as a fungal supplement.
#
# Outputs (under $OUT_ROOT):
#   unified_orfs/<ref_id>/...                per-bacterial-ref prodigal
#   unified_orfs.faa / .gff / .ffn           merged (bacteria + fungus proteins/gff)
#   gene_inventory.tsv                       per-gene metadata (all genes)
#   cds_features.tsv                         per-CDS-interval (multi-exon for fungus) -> SAF source
#   unified_annotations/unified.emapper.*    eggNOG-mapper output
#   gene_metadata.tsv                        gene_inventory + eggNOG columns (all genes)
#   fungal_jgi_annotation.tsv                JGI native KEGG/GO/KOG/IPR for fungal genes
#   samples/<s>/ref.saf                      per-sample SAF for featureCounts
#
# NOTE: the ref.saf written by STAGE 5 is PROVISIONAL. It still contains prodigal
# ORFs predicted inside rRNA/tRNA genes, which must not be counted as CDS for
# metaRS or metaT. 5b_make_structural_rna_blacklist.sh identifies them and
# 5c_gen_saf.sh rewrites every ref.saf without them.
# 5c is the authority on SAF content - always run it after this script.
# See README.md.
#
# Usage:
#   bash 5a_build_annotation.sh configs/12model.config.sh

set -uo pipefail

CONFIG=${1:?config file required}
source "$CONFIG"

: "${OUT_ROOT:?}"
: "${EGGNOG_DATA_DIR:?}"
# RA is an all-prokaryote cohort: there is no fungal genome in the pool, so
# FUNGAL_REF_ID is empty and STAGE 2b / STAGE 4b are skipped. Every other
# stage already filters on source != "Fungus", so nothing else changes.
: "${FUNGAL_REF_ID:=}"
: "${FUNGAL_GFF:=}"
: "${FUNGAL_PROT:=}"
: "${THREADS:=16}"
: "${PARALLEL_PRODIGAL:=${THREADS}}"
: "${EGGNOG_BLOCK_SIZE:=2.0}"
: "${EGGNOG_TAX_SCOPE:=auto}"
# SKIP_EGGNOG=1 stops after prodigal and still writes gene_inventory.tsv,
# cds_features.tsv and the SAFs - which is everything 5b, 5c and 6 need. Only
# gene_metadata.tsv (the functional layer) is left out, and a later rerun with
# SKIP_EGGNOG=0 adds it without recomputing anything else. RA's GlobDB-derived
# pool can reach millions of proteins, where emapper takes days; the count
# tables should not have to wait for it.
: "${SKIP_EGGNOG:=0}"
: "${EGGNOG_SENS_MODE:=iterate}"

# ---- conda activation -------------------------------------------------------
CONDA_ENV=${CONDA_ENV:-base}
if [[ -z "${MOPP2_NO_CONDA_SWITCH:-}" ]] && { [[ -z "${CONDA_DEFAULT_ENV:-}" ]] || [[ "${CONDA_DEFAULT_ENV}" != "$CONDA_ENV" ]]; }; then
    set +u
    for CONDA_BASE in "$HOME/miniconda3" "$HOME/anaconda3" /opt/miniconda3 /opt/anaconda3; do
        if [[ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
            source "$CONDA_BASE/etc/profile.d/conda.sh"; conda activate "$CONDA_ENV"; break
        fi
    done
    set -u
fi

for tool in prodigal emapper.py python3; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: $tool not in PATH (env ${CONDA_DEFAULT_ENV:-none})"; exit 1; }
done
for f in eggnog.db eggnog_proteins.dmnd; do
    [[ -f "$EGGNOG_DATA_DIR/$f" ]] || { echo "ERROR: $EGGNOG_DATA_DIR/$f missing"; exit 1; }
done

UNIFIED_TSV="${OUT_ROOT}/unified_db.tsv"
[[ -f "$UNIFIED_TSV" ]] || { echo "ERROR: $UNIFIED_TSV not found (run 12model.build_subdb.sh first)"; exit 1; }

ORF_ROOT="${OUT_ROOT}/unified_orfs"
EMAP_ROOT="${OUT_ROOT}/unified_annotations"
MERGED_FAA="${OUT_ROOT}/unified_orfs.faa"
MERGED_GFF="${OUT_ROOT}/unified_orfs.gff"
MERGED_FFN="${OUT_ROOT}/unified_orfs.ffn"
GENE_INV="${OUT_ROOT}/gene_inventory.tsv"
CDS_FEAT="${OUT_ROOT}/cds_features.tsv"
GENE_META="${OUT_ROOT}/gene_metadata.tsv"
FUN_SUPP="${OUT_ROOT}/fungal_jgi_annotation.tsv"
# Defined unconditionally: STAGE 2c tests it even when STAGE 2b is skipped.
FUN_CDS_FEAT="${OUT_ROOT}/.fungal_cds_features.tsv"

mkdir -p "$ORF_ROOT" "$EMAP_ROOT"
cd "$OUT_ROOT"
exec > >(tee annotation.log) 2>&1

echo "[$(date '+%F %T')] === 12.model annotation ==="
echo "  Unified refs:  $(tail -n +2 "$UNIFIED_TSV" | wc -l)"
echo "  Fungal ref:    ${FUNGAL_REF_ID:-<none: all-prokaryote pool>}"
echo "  eggNOG data:   $EGGNOG_DATA_DIR"
echo "  prodigal:      $(prodigal -v 2>&1 | head -1)"
echo "  emapper:       $(emapper.py --version 2>&1 | head -1)"

# ============================================================================
# STAGE 1: prodigal -p meta on prokaryotic refs (source != Fungus)
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 1: prodigal ORF prediction (prokaryotes)"

run_prodigal_one() {
    local rid="$1" rpath="$2"
    local outdir="${ORF_ROOT}/${rid}"
    local faa="${outdir}/${rid}.faa" gff="${outdir}/${rid}.gff" ffn="${outdir}/${rid}.ffn"
    local log="${outdir}/${rid}.prodigal.log" done_flag="${outdir}/.done"
    mkdir -p "$outdir"
    [[ -f "$done_flag" && -s "$faa" && -s "$gff" ]] && return 0
    local CAT="cat"; [[ "$rpath" == *.gz ]] && CAT="zcat"
    if ! $CAT "$rpath" | prodigal -p meta -f gff -q -o "$gff" -a "$faa" -d "$ffn" >/dev/null 2>"$log"; then
        echo "    [$rid] prodigal FAILED (see $log)"; rm -f "$faa" "$gff" "$ffn" "$done_flag"; return 1
    fi
    touch "$done_flag"
}
export -f run_prodigal_one
export ORF_ROOT

JOBLIST=$(mktemp)
awk -F'\t' -v fr="$FUNGAL_REF_ID" 'NR>1 && $2!="Fungus" && $1!=fr {print $1"\t"$4}' "$UNIFIED_TSV" > "$JOBLIST"
N_REFS=$(wc -l < "$JOBLIST")
echo "  prokaryotic refs: $N_REFS ($PARALLEL_PRODIGAL in parallel)"
< "$JOBLIST" xargs -L1 -P "$PARALLEL_PRODIGAL" bash -c 'run_prodigal_one "$1" "$2"' _ || true
rm -f "$JOBLIST"
N_DONE=$(find "$ORF_ROOT" -maxdepth 2 -name .done | wc -l)
echo "  prodigal done: $N_DONE / $N_REFS"

# ============================================================================
# STAGE 2: merge prokaryotic ORFs -> gene_inventory + merged faa/gff/ffn
#          (GFF is the source of truth for contig names + coords)
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 2: merge prokaryotic ORFs"

: > "$MERGED_FAA"; : > "$MERGED_GFF"; : > "$MERGED_FFN"
printf 'gene_id\tref_id\tsource\tcontig\tstart\tend\tstrand\tlength_nt\tlength_aa\tpartial\tstart_type\n' > "$GENE_INV"

python3 - "$UNIFIED_TSV" "$ORF_ROOT" "$MERGED_FAA" "$MERGED_GFF" "$MERGED_FFN" "$GENE_INV" "$FUNGAL_REF_ID" <<'PY'
import sys, os, re
unified_tsv, orf_root, out_faa, out_gff, out_ffn, gene_inv, fungal_ref = sys.argv[1:8]

refs = {}
with open(unified_tsv) as f:
    hdr = f.readline().rstrip('\n').split('\t'); ci = {n: hdr.index(n) for n in ['ref_id','source']}
    for line in f:
        p = line.rstrip('\n').split('\t')
        if not p or not p[0]: continue
        if p[ci['source']] == 'Fungus' or p[ci['ref_id']] == fungal_ref: continue
        refs[p[ci['ref_id']]] = p[ci['source']]

attr_re = re.compile(r'(\w+)=([^;]+)')
id_re = re.compile(r'ID=([^;\s]+)')
n_used = gene_total = 0

with open(out_faa,'w') as faa_out, open(out_gff,'w') as gff_out, \
     open(out_ffn,'w') as ffn_out, open(gene_inv,'a') as inv:
    gff_out.write("##gff-version 3\n")
    for rid, source in refs.items():
        rdir = os.path.join(orf_root, rid)
        if not os.path.isfile(os.path.join(rdir, '.done')): continue
        faa_in, gff_in, ffn_in = (os.path.join(rdir, f"{rid}.{e}") for e in ('faa','gff','ffn'))
        if not (os.path.isfile(faa_in) and os.path.isfile(gff_in)): continue
        n_used += 1
        id_meta = {}
        with open(gff_in) as gin:
            for line in gin:
                if line.startswith('#') or not line.strip(): continue
                c = line.rstrip('\n').split('\t')
                if len(c) < 9 or c[2] != 'CDS': continue
                seqid, start, end, strand, attrs = c[0], c[3], c[4], c[6], c[8]
                d = dict(attr_re.findall(attrs)); pid = d.get('ID')
                if not pid: continue
                gidx = pid.rsplit('_',1)[-1]; nc = f"{rid}__{seqid}"
                id_meta[pid] = {'gene_id': f"{nc}_{gidx}", 'new_contig': nc,
                                'start': start, 'end': end, 'strand': strand,
                                'partial': d.get('partial','00'), 'start_type': d.get('start_type','NA')}
        with open(faa_in) as fin:
            for line in fin:
                if not line.startswith('>'): faa_out.write(line); continue
                m = id_re.search(line); meta = id_meta.get(m.group(1)) if m else None
                if not meta: faa_out.write(line); continue
                tail = line.rstrip('\n').split(' # ',1); tail = tail[1] if len(tail)>1 else ''
                faa_out.write(f">{meta['gene_id']} # {tail}\n"); gene_total += 1
                inv.write(f"{meta['gene_id']}\t{rid}\t{source}\t{meta['new_contig']}\t"
                          f"{meta['start']}\t{meta['end']}\t{meta['strand']}\t"
                          f"{int(meta['end'])-int(meta['start'])+1}\t0\t{meta['partial']}\t{meta['start_type']}\n")
        if os.path.isfile(ffn_in):
            with open(ffn_in) as fin:
                for line in fin:
                    if not line.startswith('>'): ffn_out.write(line); continue
                    m = id_re.search(line); meta = id_meta.get(m.group(1)) if m else None
                    if meta:
                        tail = line.rstrip('\n').split(' # ',1); tail = tail[1] if len(tail)>1 else ''
                        ffn_out.write(f">{meta['gene_id']} # {tail}\n")
                    else: ffn_out.write(line)
        with open(gff_in) as gin:
            for line in gin:
                if line.startswith('#') or not line.strip(): continue
                c = line.rstrip('\n').split('\t')
                if len(c) < 9 or c[2] != 'CDS': continue
                d = dict(attr_re.findall(c[8])); meta = id_meta.get(d.get('ID'))
                if not meta: continue
                c[0] = meta['new_contig']
                c[8] = f"ID={meta['gene_id']};" + ";".join(f"{k}={v}" for k,v in d.items() if k!='ID')
                gff_out.write("\t".join(c) + "\n")
sys.stderr.write(f"  prok refs merged: {n_used}, genes: {gene_total}\n")
PY

# fill length_aa for prokaryotic genes from merged FAA
python3 - "$MERGED_FAA" "$GENE_INV" <<'PY'
import sys, os
faa, inv = sys.argv[1], sys.argv[2]
aa = {}; gid=None; buf=[]
def flush():
    if gid is None: return
    s=''.join(buf); s=s[:-1] if s.endswith('*') else s; aa[gid]=len(s)
for line in open(faa):
    if line.startswith('>'): flush(); gid=line[1:].split()[0]; buf=[]
    else: buf.append(line.rstrip('\n'))
flush()
tmp=inv+'.tmp'
with open(inv) as fin, open(tmp,'w') as fout:
    h=fin.readline(); fout.write(h); cols=h.rstrip('\n').split('\t')
    ig=cols.index('gene_id'); ia=cols.index('length_aa')
    for line in fin:
        c=line.rstrip('\n').split('\t')
        if len(c)>ia: c[ia]=str(aa.get(c[ig],0))
        fout.write("\t".join(c)+"\n")
os.replace(tmp,inv)
PY
echo "  prokaryotic genes: $(tail -n +2 "$GENE_INV" | wc -l)"

# ============================================================================
# STAGE 2b: inject fungal gene models from JGI GeneCatalog
#   - append proteins to merged FAA (header -> gene_id)
#   - append CDS rows to merged GFF (seqid -> ref__scaffold, ID -> gene_id)
#   - append one gene row per protein to gene_inventory (span + summed CDS len)
#   - write fungal_cds_features.tsv (one row per CDS exon) for SAF
# ============================================================================
echo
if [[ -z "$FUNGAL_REF_ID" ]]; then
echo "[$(date '+%F %T')] STAGE 2b: no fungal genome in this pool, skipping"
else
echo "[$(date '+%F %T')] STAGE 2b: injecting fungal gene models (JGI GeneCatalog)"

python3 - "$FUNGAL_REF_ID" "$FUNGAL_GFF" "$FUNGAL_PROT" "$MERGED_FAA" "$MERGED_GFF" "$GENE_INV" "$FUN_CDS_FEAT" <<'PY'
import sys, os, re, gzip
ref, gff, prot, out_faa, out_gff, gene_inv, cds_feat = sys.argv[1:8]

def op(p): return gzip.open(p,'rt') if p.endswith('.gz') else open(p)

# ---- parse GeneCatalog GFF: proteinId -> list of CDS exons ----
# columns: seqid JGI CDS start end . strand phase  name "..."; proteinId N; ...
pid_re = re.compile(r'proteinId\s+(\d+)')
genes = {}   # pid -> {'scaffold','strand','exons':[(start,end)]}
with open(gff) as f:
    for line in f:
        if not line.strip() or line.startswith('#'): continue
        c = line.rstrip('\n').split('\t')
        if len(c) < 9 or c[2] != 'CDS': continue
        m = pid_re.search(c[8])
        if not m: continue
        pid = m.group(1)
        g = genes.setdefault(pid, {'scaffold': c[0], 'strand': c[6], 'exons': []})
        g['exons'].append((int(c[3]), int(c[4])))

# ---- proteins: pid -> aa length; also stream to append with new headers ----
def prot_id(h):
    # >jgi|TriharM10_1|100473|CE100472_2802
    parts = h[1:].split('|')
    return parts[2] if len(parts) >= 3 else h[1:].split()[0]

aa_len = {}
with op(prot) as f, open(out_faa,'a') as faa_out:
    pid=None; buf=[]
    def flush():
        if pid is None: return
        s=''.join(buf); s=s[:-1] if s.endswith('*') else s
        aa_len[pid]=len(s)
    for line in f:
        if line.startswith('>'):
            # write previous
            if pid is not None and pid in genes:
                faa_out.write(f">{ref}__gene{pid}\n")
                seq=''.join(buf)
                for i in range(0,len(seq),70): faa_out.write(seq[i:i+70]+"\n")
            flush(); pid = prot_id(line); buf=[]
        else:
            buf.append(line.rstrip('\n'))
    # last
    if pid is not None:
        s=''.join(buf); s2=s[:-1] if s.endswith('*') else s; aa_len[pid]=len(s2)
        if pid in genes:
            faa_out.write(f">{ref}__gene{pid}\n")
            for i in range(0,len(s),70): faa_out.write(s[i:i+70]+"\n")

# ---- append gene rows (gene_inventory), CDS rows (gff), exon rows (cds_feat)
n_gene=n_exon=0
with open(out_gff,'a') as gff_out, open(gene_inv,'a') as inv, open(cds_feat,'w') as cf:
    for pid, g in genes.items():
        gene_id = f"{ref}__gene{pid}"
        contig  = f"{ref}__{g['scaffold']}"
        strand  = g['strand']
        exons   = sorted(g['exons'])
        gstart  = exons[0][0]; gend = exons[-1][1]
        cds_len = sum(e-s+1 for s,e in exons)
        inv.write(f"{gene_id}\t{ref}\tFungus\t{contig}\t{gstart}\t{gend}\t{strand}\t"
                  f"{cds_len}\t{aa_len.get(pid,0)}\tNA\tNA\n")
        n_gene += 1
        for s,e in exons:
            gff_out.write(f"{contig}\tJGI\tCDS\t{s}\t{e}\t.\t{strand}\t.\tID={gene_id};proteinId={pid}\n")
            cf.write(f"{gene_id}\t{ref}\t{contig}\t{s}\t{e}\t{strand}\n")
            n_exon += 1
sys.stderr.write(f"  fungal genes: {n_gene}, CDS exons: {n_exon}\n")
PY
echo "  fungal genes injected (total gene_inventory now: $(tail -n +2 "$GENE_INV" | wc -l))"
fi

# ============================================================================
# STAGE 2c: build cds_features.tsv (SAF source): prok = 1 row/gene, fungus = exons
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 2c: building cds_features.tsv"
{
    echo -e "gene_id\tref_id\tcontig\tstart\tend\tstrand"
    # prokaryotic: 1 interval per gene from gene_inventory
    awk -F'\t' 'NR>1 && $3!="Fungus" {print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7}' "$GENE_INV"
    # fungal: exon-level
    [[ -s "$FUN_CDS_FEAT" ]] && cat "$FUN_CDS_FEAT"
} > "$CDS_FEAT"
rm -f "$FUN_CDS_FEAT"
echo "  cds_features rows: $(($(wc -l < "$CDS_FEAT")-1))"

# ============================================================================
# STAGE 2d: sanity check contig names vs a sample ref.fna
# ============================================================================
python3 - "$OUT_ROOT" "$GENE_INV" <<'PY'
import os, sys, csv
root, gi = sys.argv[1], sys.argv[2]
sd = os.path.join(root,'samples')
if not os.path.isdir(sd): sys.exit(0)
test=None
for e in sorted(os.listdir(sd)):
    if os.path.isfile(os.path.join(sd,e,'ref.fna')) and os.path.isfile(os.path.join(sd,e,'references.tsv')):
        test=e; break
if not test: print("  (no sample ref.fna to check)"); sys.exit(0)
fna=set()
for line in open(os.path.join(sd,test,'ref.fna')):
    if line.startswith('>'): fna.add(line[1:].split()[0])
used=set()
with open(os.path.join(sd,test,'references.tsv')) as f:
    h=f.readline().rstrip('\n').split('\t'); ir=h.index('ref_id')
    for line in f:
        c=line.rstrip('\n').split('\t')
        if len(c)>ir: used.add(c[ir])
n=m=0; miss=None
with open(gi) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r['ref_id'] not in used: continue
        n+=1
        if r['contig'] in fna: m+=1
        elif miss is None: miss=(r['ref_id'],r['contig'])
        if n>=1000: break
if n and m==n: print(f"  OK: sample {test}, {n} gene contigs all match ref.fna")
elif n: print(f"  *** WARN: sample {test} {m}/{n} match; first mismatch {miss}"); sys.exit(2)
else: print(f"  (sample {test} uses no annotated refs)")
PY
[[ $? -ne 0 ]] && { echo "  Sanity check FAILED — aborting before eggNOG."; exit 1; }

# ============================================================================
# STAGE 3: eggNOG-mapper on merged proteome (bacteria + fungus)
# ============================================================================
echo
if [[ "$SKIP_EGGNOG" == "1" ]]; then
echo "[$(date '+%F %T')] STAGE 3+4: SKIP_EGGNOG=1, skipping eggNOG-mapper and gene_metadata.tsv"
echo "  rerun this script with SKIP_EGGNOG=0 to add them later"
else
echo "[$(date '+%F %T')] STAGE 3: eggNOG-mapper"
EMAP_ANN="${EMAP_ROOT}/unified.emapper.annotations"
if [[ -s "$EMAP_ANN" ]]; then
    echo "  Exists, skipping. (rm -rf $EMAP_ROOT to rerun)"
else
    [[ -s "$MERGED_FAA" ]] || { echo "ERROR: empty $MERGED_FAA"; exit 1; }
    case "$EGGNOG_SENS_MODE" in
        iterate) SENS="--dmnd_iterate yes" ;;
        *)       SENS="--sensmode $EGGNOG_SENS_MODE" ;;
    esac
    mkdir -p "$EMAP_ROOT/tmp"
    ( cd "$EMAP_ROOT"
      emapper.py -i "$MERGED_FAA" -o unified --output_dir "$EMAP_ROOT" \
        --data_dir "$EGGNOG_DATA_DIR" --itype proteins -m diamond $SENS \
        --dmnd_db "$EGGNOG_DATA_DIR/eggnog_proteins.dmnd" \
        --tax_scope "$EGGNOG_TAX_SCOPE" --block_size "$EGGNOG_BLOCK_SIZE" \
        --cpu "$THREADS" --temp_dir "$EMAP_ROOT/tmp" --override 2>&1 | tee emapper.log )
    [[ -s "$EMAP_ANN" ]] || { echo "ERROR: no annotations produced"; exit 1; }
fi
echo "  annotated genes: $(grep -vc '^#' "$EMAP_ANN" 2>/dev/null || echo 0)"

# ============================================================================
# STAGE 4: gene_metadata.tsv = gene_inventory + eggNOG columns (all genes)
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 4: gene_metadata.tsv"
python3 - "$GENE_INV" "$EMAP_ANN" "$GENE_META" <<'PY'
import sys
gi, ann, out = sys.argv[1:4]
emap={}; cols=[]
for line in open(ann):
    if line.startswith('##'): continue
    if line.startswith('#'): cols=line.lstrip('#').rstrip('\n').split('\t'); continue
    if not cols: continue
    r=line.rstrip('\n').split('\t')
    if len(r)<len(cols): r+=['']*(len(cols)-len(r))
    emap[r[0]]=r
WANT=['seed_ortholog','evalue','score','eggNOG_OGs','max_annot_lvl','COG_category',
      'Description','Preferred_name','GOs','EC','KEGG_ko','KEGG_Pathway','KEGG_Module',
      'KEGG_Reaction','KEGG_rclass','BRITE','KEGG_TC','CAZy','BiGG_Reaction','PFAMs']
idx={c:cols.index(c) for c in WANT if c in cols}
with open(gi) as fin, open(out,'w') as fout:
    h=fin.readline().rstrip('\n').split('\t')
    fout.write("\t".join(h+[c for c in WANT if c in idx])+"\n")
    ig=h.index('gene_id'); na=nt=0
    for line in fin:
        c=line.rstrip('\n').split('\t'); nt+=1
        row=emap.get(c[ig])
        extra=[row[idx[k]] for k in WANT if k in idx] if row else ['' for k in WANT if k in idx]
        if row: na+=1
        fout.write("\t".join(c+extra)+"\n")
sys.stderr.write(f"  genes: {nt}, eggNOG-annotated: {na} ({100*na/max(nt,1):.1f}%)\n")
PY
fi

# ============================================================================
# STAGE 4b: fungal JGI native annotation supplement (KEGG/GO/KOG/IPR)
# ============================================================================
echo
if [[ -z "$FUNGAL_REF_ID" ]]; then
echo "[$(date '+%F %T')] STAGE 4b: no fungal genome in this pool, skipping"
else
echo "[$(date '+%F %T')] STAGE 4b: fungal JGI native annotation supplement"
python3 - "$FUNGAL_REF_ID" "$FUN_SUPP" "${FUNGAL_KEGG:-}" "${FUNGAL_GO:-}" "${FUNGAL_KOG:-}" "${FUNGAL_IPR:-}" <<'PY'
import sys, os, gzip
ref, out, kegg, go, kog, ipr = sys.argv[1:7]
def op(p): return gzip.open(p,'rt') if p.endswith('.gz') else open(p)
def rows(p, pid_col=0):
    if not p or not os.path.isfile(p): return
    with op(p) as f:
        hdr=f.readline().rstrip('\n').split('\t')
        for line in f:
            c=line.rstrip('\n').split('\t')
            if len(c)<=pid_col: continue
            yield hdr, c
data={}   # pid -> dict of sets
def add(pid, key, val):
    if not val or val in ('\\N',''): return
    data.setdefault(pid,{}).setdefault(key,set()).add(val)
# KEGG: proteinId ecNum ... pathway ...
for hdr,c in rows(kegg):
    pid=c[0]; d=dict(zip(hdr,c))
    add(pid,'KEGG_EC', d.get('ecNum','')); add(pid,'KEGG_pathway', d.get('pathway',''))
# GO: proteinId ... goAcc
for hdr,c in rows(go):
    pid=c[0]; d=dict(zip(hdr,c)); add(pid,'GO', d.get('goAcc',''))
# KOG: transcriptId proteinId kogid kogdefline kogClass ...
for hdr,c in rows(kog, pid_col=1):
    d=dict(zip(hdr,c)); pid=d.get('proteinId','')
    add(pid,'KOG', d.get('kogid','')); add(pid,'KOG_class', (d.get('kogClass','') or '').strip())
# IPR: proteinId iprId iprDesc ...
for hdr,c in rows(ipr):
    pid=c[0]; d=dict(zip(hdr,c)); add(pid,'InterPro', d.get('iprId',''))
KEYS=['KEGG_EC','KEGG_pathway','GO','KOG','KOG_class','InterPro']
with open(out,'w') as o:
    o.write("gene_id\tproteinId\t"+"\t".join(KEYS)+"\n")
    for pid in sorted(data, key=lambda x:int(x) if x.isdigit() else x):
        vals=["; ".join(sorted(data[pid].get(k,[]))) for k in KEYS]
        o.write(f"{ref}__gene{pid}\t{pid}\t"+"\t".join(vals)+"\n")
sys.stderr.write(f"  fungal genes with JGI annotation: {len(data)} -> {os.path.basename(out)}\n")
PY
fi

# ============================================================================
# STAGE 5: per-sample SAF for featureCounts (from cds_features.tsv)
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 5: per-sample SAF"
mapfile -t SAMPLES < <(for d in "$OUT_ROOT"/samples/*/; do s=$(basename "$d"); [[ -f "${d}references.tsv" ]] && echo "$s"; done)
echo "  samples with references.tsv: ${#SAMPLES[@]}"

python3 - "$OUT_ROOT" "$CDS_FEAT" "${SAMPLES[@]}" <<'PY'
import sys, os
root, cds = sys.argv[1], sys.argv[2]; samples=sys.argv[3:]
ref_feats={}
with open(cds) as f:
    h=f.readline().rstrip('\n').split('\t'); ci={n:h.index(n) for n in ['gene_id','ref_id','contig','start','end','strand']}
    for line in f:
        c=line.rstrip('\n').split('\t')
        if len(c)<max(ci.values())+1: continue
        ref_feats.setdefault(c[ci['ref_id']],[]).append(
            (c[ci['gene_id']],c[ci['contig']],c[ci['start']],c[ci['end']],c[ci['strand']]))
for s in samples:
    sdir=os.path.join(root,'samples',s); used=[]
    with open(os.path.join(sdir,'references.tsv')) as f:
        h=f.readline().rstrip('\n').split('\t'); ir=h.index('ref_id')
        for line in f:
            c=line.rstrip('\n').split('\t')
            if len(c)>ir: used.append(c[ir])
    with open(os.path.join(sdir,'ref.saf'),'w') as out:
        out.write("GeneID\tChr\tStart\tEnd\tStrand\n"); ng=nr=0
        for rid in used:
            fe=ref_feats.get(rid)
            if not fe: continue
            nr+=1
            for gid,contig,st,en,strand in fe:
                out.write(f"{gid}\t{contig}\t{st}\t{en}\t{strand}\n"); ng+=1
    sys.stderr.write(f"  [{s}] refs={len(used)} with_feats={nr} rows={ng}\n")
PY

# ============================================================================
# STAGE 6: Summary
# ============================================================================
echo
echo "[$(date '+%F %T')] === Summary ==="
{
    printf 'metric\tvalue\n'
    printf 'n_unified_refs\t%s\n'   "$(tail -n +2 "$UNIFIED_TSV" | wc -l)"
    printf 'n_genes_total\t%s\n'    "$(tail -n +2 "$GENE_INV" | wc -l)"
    printf 'n_genes_fungal\t%s\n'   "$(awk -F'\t' 'NR>1 && $3=="Fungus"' "$GENE_INV" | wc -l)"
    printf 'n_genes_annotated\t%s\n' "$(grep -vc '^#' "$EMAP_ANN" 2>/dev/null || echo 0)"
    printf 'n_samples_with_saf\t%s\n' "$(ls "$OUT_ROOT"/samples/*/ref.saf 2>/dev/null | wc -l)"
} | column -t -s $'\t' | tee annotation_summary.tsv

echo
echo "[$(date '+%F %T')] ANNOTATION COMPLETE."
echo "  Merged proteome:   $MERGED_FAA"
echo "  Gene inventory:    $GENE_INV"
echo "  CDS features:      $CDS_FEAT"
if [[ "$SKIP_EGGNOG" == "1" ]]; then
  echo "  eggNOG annotation: <skipped>"
else
  echo "  eggNOG annotation: $EMAP_ANN"
  echo "  Gene metadata:     $GENE_META"
fi
[[ -n "$FUNGAL_REF_ID" ]] && echo "  Fungal JGI suppl.: $FUN_SUPP"
echo "  Per-sample SAFs:   $OUT_ROOT/samples/<s>/ref.saf"
