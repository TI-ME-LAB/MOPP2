#!/bin/bash
# 6_build_features.sh - manifest-driven feature tables (featureCounts + matrices)
#
#   Usage: bash 6_build_features.sh <mapping_config.sh>
#
# Runs after 5b_make_structural_rna_blacklist.sh and 5c_gen_saf.sh.
# By that point the SAF already contains only real protein-coding features, so
# every table this script writes is correct as produced. There is no cleanup
# pass afterwards - the old ${FEATURES_OUT}_clean/ stage no longer exists.
#
# CONSUMES
#   $MAPPING_OUT/manifest.tsv                  (from run_mapping2.sh; regenerated if absent)
#   $MAPPING_OUT/<read_id>/<variant>/<read_id>.bam
#   $SUBDB_ROOT/samples/<db_sample>/ref.saf    (from 5c_gen_saf.sh)
#   $SUBDB_ROOT/gene_inventory.tsv, unified_db.tsv
#   $SUBDB_ROOT/structural_rna_orf_blacklist.txt   (from step 1; see below)
#
# PRODUCES under $FEATURES_OUT
#   fc_per_sample/<read_id>.<lib>.counts(.summary)
#   tables/gene_features.<lib>.{count,tpm}.tsv
#   tables/genome_features.<lib>.{count,tpm,rpkm}.tsv
#   tables/gene_annotation.tsv, genome_annotation.tsv, counting_summary.tsv
#
# THE BLACKLIST IS APPLIED TWICE, ON PURPOSE
#   - In the SAF (step 2): featureCounts never sees the fake ORFs, so its
#     `-O --fraction` overlap arithmetic is right and reads shared with a real
#     neighbouring gene are credited in full.
#   - Here, to the gene universe: gene_inventory.tsv is the raw prodigal record
#     and still lists them, so gene_annotation.tsv and the count matrices would
#     otherwise disagree with the SAF. Filtering here also makes the script
#     produce correct tables from an older, unfiltered SAF - degraded, because
#     the fractional overlaps are already lost, but not wrong in the denominators.
#
# metaRS counts use the $METARS_BAM_SUBDIR variant (default metaRS_mopp2).

set -uo pipefail
CONFIG=${1:?config file required}
source "$CONFIG"
# ---------------------------------------------------------------------------
# gen_manifest() - emit a TSV to stdout describing every (reads, library, DB) job:
#   read_id  libtype  db_sample  paired  strand  r1  r2
#
# Normally read_id == db_sample (a sample's reads map to its own sub-DB). The
# METARS_DB_OVERRIDE array lets a metaRS sample borrow another sample's DB,
# e.g. syncom_1_day0 metaRS -> REF_1 sub-DB.
#
# Config vars consumed (set in the project config):
#   SUBDB_ROOT
#   METAG_DIR METAG_R1_SUFFIX METAG_R2_SUFFIX METAG_STRAND
#   METAT_DIR METAT_PAIRED METAT_R1_SUFFIX METAT_R2_SUFFIX METAT_SUFFIX METAT_STRAND
#   METARS_DIR METARS_SUFFIX METARS_STRAND
#   METARS_DB_OVERRIDE=( "read_id:db_sample" ... )   # optional
#
# A library is disabled simply by leaving its *_DIR empty/unset.
#
# NOTE: this function is duplicated in 4_run_mapping.sh and 6_build_features.sh so
# that each script stands alone. 6 cannot rely on 4 having written the manifest -
# the older 4.run_mapping.sh never wrote one, so a run resumed from those BAMs
# arrives here with no manifest.tsv. Keep the two copies in sync.
# ---------------------------------------------------------------------------
gen_manifest() {
    local samples_dir="${SUBDB_ROOT}/samples"
    local dbs=() d s
    for d in "$samples_dir"/*/; do
        [[ -d "$d" ]] || continue
        s=$(basename "$d")
        [[ -f "${d}ref.1.bt2" || -f "${d}ref.1.bt2l" ]] && dbs+=("$s")
    done

    printf 'read_id\tlibtype\tdb_sample\tpaired\tstrand\tr1\tr2\n'

    local r1 r2 t1 t2 u
    for s in "${dbs[@]}"; do
        # ---- metaG (paired) ----
        if [[ -n "${METAG_DIR:-}" ]]; then
            r1="${METAG_DIR}/${s}${METAG_R1_SUFFIX}"; r2="${METAG_DIR}/${s}${METAG_R2_SUFFIX}"
            [[ -f "$r1" && -f "$r2" ]] && \
                printf '%s\tmetaG\t%s\t1\t%s\t%s\t%s\n' "$s" "$s" "${METAG_STRAND:-0}" "$r1" "$r2"
        fi
        # ---- metaT (paired or single) ----
        if [[ -n "${METAT_DIR:-}" ]]; then
            if [[ "${METAT_PAIRED:-1}" == "1" ]]; then
                t1="${METAT_DIR}/${s}${METAT_R1_SUFFIX}"; t2="${METAT_DIR}/${s}${METAT_R2_SUFFIX}"
                [[ -f "$t1" && -f "$t2" ]] && \
                    printf '%s\tmetaT\t%s\t1\t%s\t%s\t%s\n' "$s" "$s" "${METAT_STRAND:-2}" "$t1" "$t2"
            else
                t1="${METAT_DIR}/${s}${METAT_SUFFIX}"
                [[ -f "$t1" ]] && \
                    printf '%s\tmetaT\t%s\t0\t%s\t%s\t\n' "$s" "$s" "${METAT_STRAND:-2}" "$t1"
            fi
        fi
        # ---- metaRS (single-end RPF, own DB) ----
        if [[ -n "${METARS_DIR:-}" ]]; then
            u="${METARS_DIR}/${s}${METARS_SUFFIX}"
            [[ -f "$u" ]] && \
                printf '%s\tmetaRS\t%s\t0\t%s\t%s\t\n' "$s" "$s" "${METARS_STRAND:-0}" "$u"
        fi
    done

    # ---- metaRS overrides: read_id borrows another sample's DB ----
    if [[ -n "${METARS_DB_OVERRIDE+x}" ]]; then
        local entry read_id db
        for entry in "${METARS_DB_OVERRIDE[@]}"; do
            read_id="${entry%%:*}"; db="${entry##*:}"
            u="${METARS_DIR}/${read_id}${METARS_SUFFIX}"
            [[ -f "$u" ]] || { echo "  WARN(manifest): override reads missing: $u" >&2; continue; }
            if [[ -f "${SUBDB_ROOT}/samples/${db}/ref.1.bt2" || -f "${SUBDB_ROOT}/samples/${db}/ref.1.bt2l" ]]; then
                printf '%s\tmetaRS\t%s\t0\t%s\t%s\t\n' "$read_id" "$db" "${METARS_STRAND:-0}" "$u"
            else
                echo "  WARN(manifest): override DB missing index: samples/${db}" >&2
            fi
        done
    fi
}

: "${SUBDB_ROOT:?}"
: "${MAPPING_OUT:?}"
: "${FEATURES_OUT:?}"
: "${THREADS:=16}"
: "${METARS_BAM_SUBDIR:=metaRS_mopp2}"

CONDA_ENV=${FEAT_CONDA_ENV:-${CONDA_ENV:-features}}
if [[ -z "${MOPP2_NO_CONDA_SWITCH:-}" ]] && { [[ -z "${CONDA_DEFAULT_ENV:-}" ]] || [[ "${CONDA_DEFAULT_ENV}" != "$CONDA_ENV" ]]; }; then
    set +u
    for CB in "$HOME/miniconda3" "$HOME/anaconda3" /opt/miniconda3 /opt/anaconda3; do
        [[ -f "$CB/etc/profile.d/conda.sh" ]] && { source "$CB/etc/profile.d/conda.sh"; conda activate "$CONDA_ENV"; break; }
    done
    set -u
fi
for t in featureCounts python3 samtools; do
    command -v "$t" >/dev/null 2>&1 || { echo "ERROR: $t not found (env ${CONDA_DEFAULT_ENV:-none})"; exit 1; }
done

GENE_INV="${SUBDB_ROOT}/gene_inventory.tsv"
UNIFIED_TSV="${SUBDB_ROOT}/unified_db.tsv"
[[ -s "$GENE_INV" ]] || { echo "ERROR: $GENE_INV missing (run annotation first)"; exit 1; }

# Blacklist location: overridable in the config, otherwise the sub-DB default.
ORF_BLACKLIST=${ORF_BLACKLIST:-${SUBDB_ROOT}/structural_rna_orf_blacklist.txt}
[[ -s "$ORF_BLACKLIST" ]] || ORF_BLACKLIST=""

MANIFEST="${MAPPING_OUT}/manifest.tsv"
[[ -s "$MANIFEST" ]] || { echo "  (manifest missing - regenerating)"; gen_manifest > "$MANIFEST"; }

# ---- refuse to silently produce contaminated metaRS tables ------------------
# metaRS measures translation on mRNA. Without the blacklist, 38-68% of its
# CDS-assigned reads are rRNA landing on spurious ORFs, and the fraction varies
# 19-89% between samples - which reads as biology in any cross-sample comparison.
if [[ -z "$ORF_BLACKLIST" ]] && grep -q $'\tmetaRS\t' "$MANIFEST"; then
    if [[ "${ALLOW_UNFILTERED_METARS:-0}" != "1" ]]; then
        echo "ERROR: the manifest contains metaRS libraries but no structural-RNA blacklist was found."
        echo "  expected: ${SUBDB_ROOT}/structural_rna_orf_blacklist.txt"
        echo "  fix:      bash 5b_make_structural_rna_blacklist.sh $CONFIG"
        echo "            bash 5c_gen_saf.sh $SUBDB_ROOT <cds_features.tsv> <blacklist.txt>"
        echo "  override: ALLOW_UNFILTERED_METARS=1 (only if you know why)"
        exit 1
    fi
    echo "WARNING: proceeding without a blacklist because ALLOW_UNFILTERED_METARS=1"
fi

mkdir -p "$FEATURES_OUT/fc_per_sample" "$FEATURES_OUT/tables"
cd "$FEATURES_OUT"
exec > >(tee "${FEATURES_OUT}/features.log") 2>&1

echo "[$(date '+%F %T')] === Feature tables (manifest-driven) ==="
echo "  Sub-DB:     $SUBDB_ROOT"
echo "  Mapping:    $MAPPING_OUT"
echo "  Output:     $FEATURES_OUT"
echo "  metaRS BAM: $METARS_BAM_SUBDIR"
echo "  Blacklist:  ${ORF_BLACKLIST:-none} ${ORF_BLACKLIST:+($(wc -l < "$ORF_BLACKLIST") gene_ids)}"
echo "  featureCounts: $(featureCounts -v 2>&1 | grep -oE 'v[0-9.]+' | head -1)"

# ---- stale-SAF check -------------------------------------------------------
# If blacklisted gene_ids are still present in a SAF, that SAF predates step 2.
# The tables below will still come out with correct denominators, but the
# fractional overlap counts were already lost inside featureCounts.
if [[ -n "$ORF_BLACKLIST" ]]; then
    STALE=0
    for SAFF in "$SUBDB_ROOT"/samples/*/ref.saf; do
        [[ -s "$SAFF" ]] || continue
        if cut -f1 "$SAFF" | grep -qxFf "$ORF_BLACKLIST"; then STALE=$((STALE+1)); fi
    done
    if [[ "$STALE" -gt 0 ]]; then
        echo "  WARNING: $STALE SAF file(s) still contain blacklisted gene_ids."
        echo "           Re-run 5c_gen_saf.sh with the blacklist, then delete"
        echo "           $FEATURES_OUT/fc_per_sample/ so featureCounts re-runs."
    fi
fi

# ============================================================================
# STAGE 1: featureCounts per (read_id, library)
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 1: featureCounts..."
while IFS=$'\t' read -r READ LIB DB PAIRED STRAND R1 R2; do
    [[ "$READ" == "read_id" ]] && continue
    SAF="${SUBDB_ROOT}/samples/${DB}/ref.saf"
    case "$LIB" in
        metaG)  VARIANT=metaG ;;
        metaT)  VARIANT=metaT ;;
        metaRS) VARIANT="$METARS_BAM_SUBDIR" ;;
        *) continue ;;
    esac
    BAM="${MAPPING_OUT}/${READ}/${VARIANT}/${READ}.bam"
    OUT="${FEATURES_OUT}/fc_per_sample/${READ}.${LIB}.counts"
    # featureCounts output is cached, but it is only valid for the SAF it was
    # produced from. After 5c rewrites the SAFs, anything older is stale and must
    # be recounted - otherwise the spurious ORFs silently survive in the tables.
    if [[ -s "$OUT" && -s "${OUT}.summary" ]]; then
        if [[ "$SAF" -nt "$OUT" ]]; then
            echo "  [$READ/$LIB] SAF is newer than the counts - recounting"
        else
            echo "  [$READ/$LIB] exists, skip"; continue
        fi
    fi
    if [[ ! -s "$SAF" ]]; then echo "  [$READ/$LIB] SAF missing: $SAF (skip)"; continue; fi
    if [[ ! -s "$BAM" ]]; then echo "  [$READ/$LIB] BAM missing: $BAM (skip)"; continue; fi
    FC_ARGS=(-a "$SAF" -F SAF -T "$THREADS" -Q 0 -O --fraction -M --fraction -s "$STRAND" -o "$OUT")
    [[ "$PAIRED" == "1" ]] && FC_ARGS+=(-p --countReadPairs)
    echo "  [$READ/$LIB] db=$DB -s $STRAND paired=$PAIRED"
    featureCounts "${FC_ARGS[@]}" "$BAM" > "${OUT}.log" 2>&1 || echo "  [$READ/$LIB] FAILED (see ${OUT}.log)"
done < "$MANIFEST"

# ============================================================================
# STAGE 2: annotation tables (gene_annotation, genome_annotation)
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 2: annotation tables..."
python3 - "$FEATURES_OUT" "$SUBDB_ROOT" "$ORF_BLACKLIST" <<'PY'
import os, sys, csv
from collections import defaultdict
out_root, subdb, blacklist_path = sys.argv[1], sys.argv[2], sys.argv[3]

excluded=set()
if blacklist_path:
    with open(blacklist_path) as f:
        excluded={ln.strip() for ln in f if ln.strip()}

unified={}
with open(os.path.join(subdb,'unified_db.tsv')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        unified[r['ref_id']]={'source':r.get('source','NA'),'folder':r.get('folder','NA'),'length':r.get('length','NA')}

# gene_annotation.tsv is gene_inventory.tsv plus a folder column, minus the
# structural-RNA ORFs. It must describe exactly the gene set of the count
# matrices, so the same exclusion is applied in STAGE 3.
gi=os.path.join(subdb,'gene_inventory.tsv')
ga=os.path.join(out_root,'tables','gene_annotation.tsv')
seen=defaultdict(lambda:{'n':0,'bp':0}); n_drop=0
with open(gi) as fin, open(ga,'w') as fout:
    rdr=csv.DictReader(fin,delimiter='\t'); cols=list(rdr.fieldnames)
    if 'folder' not in cols: cols.insert(cols.index('source')+1 if 'source' in cols else len(cols),'folder')
    fout.write('\t'.join(cols)+'\n')
    for r in rdr:
        if r['gene_id'] in excluded: n_drop+=1; continue
        rid=r.get('ref_id','')
        r['folder']=unified.get(rid,{}).get('folder','NA')
        fout.write('\t'.join(r.get(c,'') for c in cols)+'\n')
        seen[rid]['n']+=1
        try: seen[rid]['bp']+=int(r.get('length_nt') or 0)
        except ValueError: pass

gna=os.path.join(out_root,'tables','genome_annotation.tsv')
with open(gna,'w') as out:
    out.write('ref_id\tsource\tfolder\tgenome_length_bp\tn_genes\ttotal_cds_bp\n')
    for r in sorted(set(unified)|set(seen)):
        u=unified.get(r,{}); s=seen.get(r,{'n':0,'bp':0})
        out.write(f"{r}\t{u.get('source','NA')}\t{u.get('folder','NA')}\t{u.get('length','NA')}\t{s['n']}\t{s['bp']}\n")
print(f"  gene_annotation + genome_annotation written (excluded {n_drop} structural-RNA ORFs)", file=sys.stderr)
PY

# ============================================================================
# STAGE 3: abundance matrices (per library: gene count/tpm, genome count/tpm/rpkm)
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 3: abundance matrices..."
python3 - "$FEATURES_OUT" "$SUBDB_ROOT" "$MANIFEST" "$ORF_BLACKLIST" <<'PY'
import os, sys, csv
from collections import defaultdict, OrderedDict
out_root, subdb, manifest, blacklist_path = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

excluded=set()
if blacklist_path:
    with open(blacklist_path) as f:
        excluded={ln.strip() for ln in f if ln.strip()}

# ---- gene_inventory: gene_id -> ref_id, length_nt (structural RNA removed) ----
gene2ref={}; gene_len={}; ref2genes=defaultdict(list)
with open(os.path.join(subdb,'gene_inventory.tsv')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        g=r['gene_id']
        if g in excluded: continue
        gene2ref[g]=r['ref_id']
        try: gene_len[g]=int(r.get('length_nt') or 0)
        except ValueError: gene_len[g]=0
        ref2genes[r['ref_id']].append(g)

# ---- unified_db: ref_id -> genome length ----
ref_len={}
with open(os.path.join(subdb,'unified_db.tsv')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        try: ref_len[r['ref_id']]=int(r.get('length') or 0)
        except ValueError: ref_len[r['ref_id']]=0

# ---- manifest: per lib -> ordered read_ids, and read_id -> db_sample ----
lib_reads=OrderedDict()          # lib -> [read_id,...]
readdb={}                        # (read_id,lib) -> db_sample
with open(manifest) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        lib=r['libtype']; rid=r['read_id']
        lib_reads.setdefault(lib,[])
        if rid not in lib_reads[lib]: lib_reads[lib].append(rid)
        readdb[(rid,lib)]=r['db_sample']

def load_fc(rid, lib):
    """Per-gene counts from one featureCounts output.

    Blacklisted genes are dropped here as well as from the gene universe. With a
    filtered SAF they are absent anyway; with an older unfiltered SAF this keeps
    them out of the TPM denominator, the genome-level rollup and the RPKM
    library size, all of which would otherwise be inflated by rRNA reads.
    """
    fp=os.path.join(out_root,'fc_per_sample',f"{rid}.{lib}.counts")
    d={}
    if not os.path.isfile(fp): return d
    with open(fp) as f:
        for line in f:
            if line.startswith('#') or line.startswith('Geneid'): continue
            p=line.rstrip('\n').split('\t')
            if p[0] in excluded: continue
            try: d[p[0]]=float(p[-1])
            except (ValueError,IndexError): pass
    return d

tdir=os.path.join(out_root,'tables')
for lib, reads in lib_reads.items():
    fc={rid: load_fc(rid,lib) for rid in reads}
    # ---- gene-level count ----
    all_genes=sorted(gene2ref)
    with open(os.path.join(tdir,f'gene_features.{lib}.count.tsv'),'w') as out:
        out.write('gene_id\t'+'\t'.join(reads)+'\n')
        for g in all_genes:
            out.write(g+'\t'+'\t'.join(f"{fc[rid].get(g,0.0):.3f}" for rid in reads)+'\n')
    # ---- gene-level tpm ----
    scale={}
    for rid in reads:
        denom=sum(fc[rid][g]/(gene_len.get(g,0)/1000.0) for g in fc[rid] if gene_len.get(g,0)>0)
        scale[rid]=1e6/denom if denom>0 else 0.0
    with open(os.path.join(tdir,f'gene_features.{lib}.tpm.tsv'),'w') as out:
        out.write('gene_id\t'+'\t'.join(reads)+'\n')
        for g in all_genes:
            ln=gene_len.get(g,0)
            row=[((fc[rid].get(g,0.0)/(ln/1000.0))*scale[rid] if ln>0 else 0.0) for rid in reads]
            out.write(g+'\t'+'\t'.join(f"{v:.4f}" for v in row)+'\n')
    # ---- genome-level: aggregate gene counts by ref ----
    gcount=defaultdict(lambda: defaultdict(float))
    for rid in reads:
        for g,c in fc[rid].items():
            r=gene2ref.get(g)
            if r is not None: gcount[r][rid]+=c
    all_refs=sorted(ref_len)
    with open(os.path.join(tdir,f'genome_features.{lib}.count.tsv'),'w') as out:
        out.write('ref_id\t'+'\t'.join(reads)+'\n')
        for r in all_refs:
            out.write(r+'\t'+'\t'.join(f"{gcount[r].get(rid,0.0):.3f}" for rid in reads)+'\n')
    # genome tpm (by genome length)
    with open(os.path.join(tdir,f'genome_features.{lib}.tpm.tsv'),'w') as out:
        out.write('ref_id\t'+'\t'.join(reads)+'\n')
        gsc={}
        for rid in reads:
            denom=sum(gcount[r].get(rid,0.0)/(ref_len[r]/1000.0) for r in all_refs if ref_len[r]>0)
            gsc[rid]=1e6/denom if denom>0 else 0.0
        for r in all_refs:
            bp=ref_len[r]
            row=[((gcount[r].get(rid,0.0)/(bp/1000.0))*gsc[rid] if bp>0 else 0.0) for rid in reads]
            out.write(r+'\t'+'\t'.join(f"{v:.4f}" for v in row)+'\n')
    # genome rpkm (by genome length + library size)
    libsize={rid: sum(fc[rid].values()) for rid in reads}
    with open(os.path.join(tdir,f'genome_features.{lib}.rpkm.tsv'),'w') as out:
        out.write('ref_id\t'+'\t'.join(reads)+'\n')
        for r in all_refs:
            bp=ref_len[r]
            row=[(gcount[r].get(rid,0.0)*1e9/(bp*libsize[rid]) if bp>0 and libsize[rid]>0 else 0.0) for rid in reads]
            out.write(r+'\t'+'\t'.join(f"{v:.4f}" for v in row)+'\n')
    print(f"  {lib}: {len(reads)} samples, {len(all_genes)} genes, gene+genome tables written", file=sys.stderr)
PY

# ============================================================================
# STAGE 4: counting summary
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 4: counting summary..."
python3 - "$FEATURES_OUT" "$MANIFEST" <<'PY'
import os, sys, csv
out_root, manifest = sys.argv[1], sys.argv[2]
rows=[]
with open(manifest) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        rows.append((r['read_id'], r['libtype']))
def psum(p):
    if not os.path.isfile(p): return None
    d={}
    with open(p) as f:
        next(f,None)
        for line in f:
            x=line.rstrip('\n').split('\t')
            if len(x)>=2:
                try: d[x[0]]=float(x[1])
                except ValueError: pass
    return d
op=os.path.join(out_root,'tables','counting_summary.tsv')
with open(op,'w') as out:
    out.write('read_id\tlibrary\tassigned\tno_features\tambiguity\tmultimapping\tunmapped\ttotal\tpct_assigned\n')
    for rid,lib in rows:
        s=psum(os.path.join(out_root,'fc_per_sample',f'{rid}.{lib}.counts.summary'))
        if s is None: out.write(f"{rid}\t{lib}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"); continue
        a=s.get('Assigned',0); tot=sum(s.values()); pct=100*a/tot if tot>0 else 0
        out.write(f"{rid}\t{lib}\t{a:.0f}\t{s.get('Unassigned_NoFeatures',0):.0f}\t{s.get('Unassigned_Ambiguity',0):.0f}\t"
                  f"{s.get('Unassigned_MultiMapping',0):.0f}\t{s.get('Unassigned_Unmapped',0):.0f}\t{tot:.0f}\t{pct:.2f}\n")
print(f"  -> {op}", file=sys.stderr)
PY

echo
echo "=== counting_summary.tsv (head) ==="
# With the blacklist applied at the SAF, reads that used to be Assigned to fake
# ORFs now fall into no_features. A lower pct_assigned for metaRS is the expected
# and correct outcome, not a regression.
column -t -s $'\t' < "${FEATURES_OUT}/tables/counting_summary.tsv" 2>/dev/null | head -20

echo
echo "[$(date '+%F %T')] FEATURE TABLES COMPLETE."
echo "  Tables: ${FEATURES_OUT}/tables/"
echo "  These tables are final - there is no separate cleanup step."
