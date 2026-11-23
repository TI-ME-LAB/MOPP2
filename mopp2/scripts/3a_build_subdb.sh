#!/bin/bash
# 12model.build_subdb.sh
#
# Per-sample sub-database construction for the 12.model project.
#
# The genome pools (mags_pool + globdb_pool) are ALREADY curated by hand, so
# this script skips the MAG-collection / GlobDB-selection stages of the original
# 3.build_subdb.sh and starts from the two existing pool directories.
#
# Workflow:
#   Stage 1: Inventory both pools -> unified_db.tsv (ref_id, source, folder, path, length)
#            source = MAG | Fungus | GlobDB. Handles mixed extensions.
#   Stage 2: skani all-vs-all (report high-ANI pairs; drop only if DEDUP_POOL=1)
#   Stage 3: Sketch unified database for sylph
#   Stage 4: Round-2 sylph: each metaG sample (paired) vs unified database
#   Stage 5: Per-sample sub-DB (ref.fna with renamed contigs + bowtie2 index)
#   Stage 6: Summary
#
# Usage:
#   bash 12model.build_subdb.sh 12model.config.sh

set -uo pipefail

CONFIG=${1:?config file required}
source "$CONFIG"

: "${OUT_ROOT:?}"
: "${MAGS_POOL:?}"
: "${GLOBDB_POOL:?}"
: "${METAG_DIR:?}"
: "${METAG_R1_SUFFIX:=_metaG_R1.fq.gz}"
: "${METAG_R2_SUFFIX:=_metaG_R2.fq.gz}"
: "${FUNGAL_REF_ID:=}"
: "${MIN_EFF_COV_R2:=0.05}"
: "${SYLPH_C:=200}"
: "${ANI_DEDUP:=99.0}"
: "${ANI_STRAIN_BAND:=95.0}"
: "${ALIGN_FRAC_MIN:=30}"
: "${DEDUP_POOL:=0}"
: "${THREADS:=16}"

# ---- conda activation -------------------------------------------------------
CONDA_ENV=${CONDA_ENV:-base}
if [[ -z "${MOPP2_NO_CONDA_SWITCH:-}" ]] && { [[ -z "${CONDA_DEFAULT_ENV:-}" ]] || [[ "${CONDA_DEFAULT_ENV}" != "$CONDA_ENV" ]]; }; then
    set +u
    for CONDA_BASE in "$HOME/miniconda3" "$HOME/anaconda3" /opt/miniconda3 /opt/anaconda3; do
        if [[ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
            source "$CONDA_BASE/etc/profile.d/conda.sh"
            conda activate "$CONDA_ENV"
            break
        fi
    done
    set -u
fi

for tool in sylph skani bowtie2 bowtie2-build samtools; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: $tool not found (env: ${CONDA_DEFAULT_ENV:-none})"; exit 1
    fi
done

mkdir -p "$OUT_ROOT"; cd "$OUT_ROOT"
exec > >(tee build_subdb.log) 2>&1

echo "[$(date '+%F %T')] === 12.model sub-DB build ==="
echo "  mags_pool:   $MAGS_POOL"
echo "  globdb_pool: $GLOBDB_POOL"
echo "  metaG dir:   $METAG_DIR ($METAG_R1_SUFFIX / $METAG_R2_SUFFIX)"
echo "  fungal ref:  ${FUNGAL_REF_ID:-<none>}"
echo "  dedup pool:  $DEDUP_POOL (0 = keep all, report only)"

# ============================================================================
# STAGE 1: Inventory both pools -> unified_db.tsv
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 1: Inventorying pools..."

python3 - "$MAGS_POOL" "$GLOBDB_POOL" "$FUNGAL_REF_ID" > unified_db.tsv <<'PY'
import sys, os, gzip

mags_pool, globdb_pool, fungal_ref = sys.argv[1], sys.argv[2], sys.argv[3]

EXTS = ('.fa.gz', '.fna.gz', '.fasta.gz', '.fa', '.fna', '.fasta', '.gz')
def ref_id(fn):
    for e in EXTS:
        if fn.endswith(e):
            return fn[:-len(e)]
    return fn

def glen(path):
    op = gzip.open if path.endswith('.gz') else open
    n = 0
    try:
        with op(path, 'rt') as f:
            for line in f:
                if not line.startswith('>'):
                    n += len(line.rstrip())
    except Exception:
        return 0
    return n

rows = []
seen = set()
def scan(pool, default_source):
    if not os.path.isdir(pool):
        return
    for fn in sorted(os.listdir(pool)):
        path = os.path.join(pool, fn)
        if not os.path.isfile(path):
            continue
        if not fn.endswith(EXTS):
            continue
        rid = ref_id(fn)
        if rid in seen:
            sys.stderr.write(f"  WARN: duplicate ref_id {rid} ({path}) skipped\n")
            continue
        seen.add(rid)
        src = default_source
        if fungal_ref and rid == fungal_ref:
            src = 'Fungus'
        folder = os.path.basename(pool.rstrip('/'))
        rows.append((rid, src, folder, path, glen(path)))

scan(mags_pool, 'MAG')
scan(globdb_pool, 'GlobDB')

print("ref_id\tsource\tfolder\tpath\tlength")
for rid, src, folder, path, length in sorted(rows, key=lambda r: (r[1], r[0])):
    print(f"{rid}\t{src}\t{folder}\t{path}\t{length}")
PY

N_TOTAL=$(tail -n +2 unified_db.tsv | wc -l)
N_MAG=$(awk -F'\t' 'NR>1 && $2=="MAG"' unified_db.tsv | wc -l)
N_FUN=$(awk -F'\t' 'NR>1 && $2=="Fungus"' unified_db.tsv | wc -l)
N_GDB=$(awk -F'\t' 'NR>1 && $2=="GlobDB"' unified_db.tsv | wc -l)
echo "  Unified pool: $N_TOTAL genomes (MAG:$N_MAG Fungus:$N_FUN GlobDB:$N_GDB)"

# ============================================================================
# STAGE 2: skani all-vs-all (report high-ANI pairs; optional dedup)
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 2: skani ANI (report; dedup=$DEDUP_POOL)..."

awk -F'\t' 'NR>1 {print $4}' unified_db.tsv > all_genomes.list
if [[ ! -s skani_all.tsv ]]; then
    skani triangle -t "$THREADS" --sparse -l all_genomes.list \
        --min-af "$ALIGN_FRAC_MIN" -s "$ANI_STRAIN_BAND" \
        > skani_all.tsv 2> skani.log || echo "  (skani returned non-zero; check skani.log)"
fi
echo "  skani pairs at ANI >= ${ANI_STRAIN_BAND}%: $(($(wc -l < skani_all.tsv)-1))"

python3 - "$OUT_ROOT" "$ANI_DEDUP" "$DEDUP_POOL" <<'PY'
import sys, os, csv
root, ani_dedup, do_dedup = sys.argv[1], float(sys.argv[2]), sys.argv[3] == '1'

cand = {}
with open(os.path.join(root, 'unified_db.tsv')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        cand[r['ref_id']] = {'length': int(r['length']), 'source': r['source']}

def to_id(p):
    bn = os.path.basename(p)
    for e in ('.fa.gz','.fna.gz','.fasta.gz','.fa','.fna','.fasta','.gz'):
        if bn.endswith(e): return bn[:-len(e)]
    return bn

high = []
with open(os.path.join(root, 'skani_all.tsv')) as f:
    rdr = csv.reader(f, delimiter='\t')
    hdr = next(rdr, None)
    if hdr:
        ir, iq, ia = hdr.index('Ref_file'), hdr.index('Query_file'), hdr.index('ANI')
        for row in rdr:
            if len(row) <= max(ir, iq, ia): continue
            a = float(row[ia]); ra, rb = to_id(row[ir]), to_id(row[iq])
            if ra == rb or ra not in cand or rb not in cand: continue
            if a >= ani_dedup:
                high.append((ra, rb, a))

with open(os.path.join(root, 'high_ani_pairs.tsv'), 'w') as out:
    out.write("ref_a\tref_b\tani\n")
    for ra, rb, a in sorted(high, key=lambda x: -x[2]):
        out.write(f"{ra}\t{rb}\t{a:.2f}\n")

drop = set()
if do_dedup and high:
    parent = {r: r for r in cand}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for ra, rb, _ in high:
        parent[find(ra)] = find(rb)
    clusters = {}
    for r in cand:
        clusters.setdefault(find(r), []).append(r)
    for members in clusters.values():
        if len(members) == 1: continue
        keep = sorted(members, key=lambda m: (-cand[m]['length'], m))[0]
        for m in members:
            if m != keep: drop.add(m)

# Rewrite unified_db.tsv keeping only survivors (identity if not dedup)
inp = os.path.join(root, 'unified_db.tsv')
rows = list(csv.reader(open(inp), delimiter='\t'))
hdr, body = rows[0], rows[1:]
with open(inp, 'w') as out:
    out.write('\t'.join(hdr) + '\n')
    for row in body:
        if row and row[0] not in drop:
            out.write('\t'.join(row) + '\n')

sys.stderr.write(f"  high-ANI (>= {ani_dedup}%) pairs: {len(high)} -> high_ani_pairs.tsv\n")
if do_dedup:
    sys.stderr.write(f"  dedup dropped: {len(drop)} genomes\n")
elif high:
    sys.stderr.write("  NOTE: high-ANI pairs present but kept (DEDUP_POOL=0). Review high_ani_pairs.tsv.\n")
PY

# ============================================================================
# STAGE 3: Sketch unified database
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 3: Sketching unified database..."
awk -F'\t' 'NR>1 {print $4}' unified_db.tsv > unified_paths.txt
if [[ ! -f unified_db.syldb ]]; then
    sylph sketch -t "$THREADS" -c "$SYLPH_C" -o unified_db -l unified_paths.txt \
        2>&1 | tee sylph_sketch.log | tail -3
else
    echo "  Reusing unified_db.syldb"
fi

# ============================================================================
# STAGE 4: Round-2 sylph per sample (paired metaG vs unified database)
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 4: Round-2 sylph per sample..."
mapfile -t SAMPLES < <(
    ls "$METAG_DIR"/*"$METAG_R1_SUFFIX" 2>/dev/null \
        | xargs -n1 basename | sed "s/${METAG_R1_SUFFIX}\$//"
)
echo "  metaG samples: ${#SAMPLES[@]}"

mkdir -p sylph_round2
for s in "${SAMPLES[@]}"; do
    out="sylph_round2/${s}.tsv"
    r1="${METAG_DIR}/${s}${METAG_R1_SUFFIX}"
    r2="${METAG_DIR}/${s}${METAG_R2_SUFFIX}"
    [[ -s "$out" ]] && { echo "  [$s] exists, skip"; continue; }
    [[ -f "$r1" && -f "$r2" ]] || { echo "  [$s] fastq missing, skip"; continue; }
    echo "  [$s] sylph profile..."
    sylph profile unified_db.syldb -1 "$r1" -2 "$r2" -t "$THREADS" \
        > "$out" 2>> sylph_round2.log
done

# ============================================================================
# STAGE 5: Per-sample sub-DB (ref.fna + bowtie2 index)
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 5: Per-sample sub-DB construction..."

build_one() {
    local SAMPLE=$1
    local SDIR="${OUT_ROOT}/samples/${SAMPLE}"
    local SYLPH="${OUT_ROOT}/sylph_round2/${SAMPLE}.tsv"
    [[ -s "$SYLPH" ]] || { echo "  [$SAMPLE] no round2 sylph, skip"; return; }
    mkdir -p "$SDIR"; cd "$SDIR"

    python3 - "$SAMPLE" "$OUT_ROOT" "$MIN_EFF_COV_R2" > references.tsv <<'PY'
import sys, os, csv
sample, root, mincov = sys.argv[1], sys.argv[2], float(sys.argv[3])
unified = {}
with open(os.path.join(root, 'unified_db.tsv')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        unified[r['ref_id']] = r
hits = []
with open(os.path.join(root, 'sylph_round2', f'{sample}.tsv')) as f:
    rdr = csv.reader(f, delimiter='\t')
    hdr = next(rdr)
    ci = {n: hdr.index(n) for n in ['Genome_file','Taxonomic_abundance','Adjusted_ANI','Eff_cov']}
    for row in rdr:
        try:
            if float(row[ci['Eff_cov']]) < mincov: continue
        except (ValueError, IndexError):
            continue
        gid = os.path.basename(row[ci['Genome_file']])
        for e in ('.fa.gz','.fna.gz','.fasta.gz','.fa','.fna','.fasta','.gz'):
            if gid.endswith(e): gid = gid[:-len(e)]; break
        if gid in unified:
            hits.append((gid, float(row[ci['Taxonomic_abundance']]),
                         float(row[ci['Adjusted_ANI']]), float(row[ci['Eff_cov']])))
print("ref_id\tsource\tfolder\tpath\ttax_abund\tani\teff_cov")
for gid, ta, ani, ec in sorted(hits, key=lambda x: -x[1]):
    u = unified[gid]
    print(f"{gid}\t{u['source']}\t{u['folder']}\t{u['path']}\t{ta:.4f}\t{ani:.2f}\t{ec:.3f}")
PY

    local N=$(tail -n +2 references.tsv | wc -l)
    echo "  [$SAMPLE] refs: $N"

    : > ref.fna
    echo -e "new_contig\tref_id\tsource\toriginal_contig" > ref.contig_map.tsv
    while IFS=$'\t' read -r rid src folder rpath _; do
        [[ -f "$rpath" ]] || { echo "    MISSING: $rpath"; continue; }
        local CAT="cat"; [[ "$rpath" == *.gz ]] && CAT="zcat"
        $CAT "$rpath" | awk -v rid="$rid" -v src="$src" -v mf="ref.contig_map.tsv" '
            /^>/ { orig=substr($1,2); new=rid"__"orig;
                   print ">"new; print new"\t"rid"\t"src"\t"orig >> mf; next }
            { print }' >> ref.fna
    done < <(tail -n +2 references.tsv)

    local nc=$(grep -c "^>" ref.fna || true)
    local bp=$(awk '!/^>/ {n+=length($0)} END{print n+0}' ref.fna)
    echo "  [$SAMPLE] ref.fna: $nc contigs, $(awk -v b=$bp 'BEGIN{printf "%.1f", b/1e6}') Mb"

    if [[ "$nc" -gt 0 ]]; then
        bowtie2-build --threads "$THREADS" ref.fna ref > bt2_build.log 2>&1
        samtools faidx ref.fna 2>/dev/null
    else
        echo "  [$SAMPLE] WARNING: 0 references passed threshold; no index built."
    fi
    cd "$OUT_ROOT"
}

for s in "${SAMPLES[@]}"; do build_one "$s"; done

# ============================================================================
# STAGE 6: Summary
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 6: Summary"
{
    echo -e "sample\tn_refs\tn_MAG\tn_Fungus\tn_GlobDB\tref_Mb"
    for s in "${SAMPLES[@]}"; do
        d="samples/${s}"; [[ -f "$d/references.tsv" ]] || continue
        n=$(tail -n +2 "$d/references.tsv" | wc -l)
        nm=$(awk -F'\t' 'NR>1 && $2=="MAG"' "$d/references.tsv" | wc -l)
        nf=$(awk -F'\t' 'NR>1 && $2=="Fungus"' "$d/references.tsv" | wc -l)
        ng=$(awk -F'\t' 'NR>1 && $2=="GlobDB"' "$d/references.tsv" | wc -l)
        bp=0; [[ -f "$d/ref.fna" ]] && bp=$(awk '!/^>/ {n+=length($0)} END{printf "%.1f", n/1e6}' "$d/ref.fna")
        echo -e "${s}\t${n}\t${nm}\t${nf}\t${ng}\t${bp}"
    done
} | column -t -s $'\t' | tee subdb_summary.tsv

echo
echo "[$(date '+%F %T')] SUB-DB BUILD COMPLETE."
