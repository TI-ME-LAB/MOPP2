#!/bin/bash
# build_subdb_round2.sh
#
# Reuse an EXISTING unified sub-DB sketch (built by 12model.build_subdb.sh in
# some source project) to construct per-sample sub-DBs for a DIFFERENT set of
# metaG samples. Only round-2 sylph + per-sample assembly run here — the genome
# pool, skani and sketch are taken as-is from $SUBDB_SRC.
#
# Reads from $SUBDB_SRC:  unified_db.syldb, unified_db.tsv
# Writes under $OUT_DIR:   sylph_round2/<s>.tsv, samples/<s>/{references.tsv,ref.fna,ref*.bt2,...}
#
# Usage:
#   bash build_subdb_round2.sh <config.sh>

set -uo pipefail
CONFIG=${1:?config file required}
source "$CONFIG"

: "${SUBDB_SRC:?SUBDB_SRC (dir with unified_db.syldb + unified_db.tsv) required}"
: "${OUT_DIR:?OUT_DIR (target for round2 + samples) required}"
: "${METAG_DIR:?}"
: "${METAG_R1_SUFFIX:=_metaG_R1.fq.gz}"
: "${METAG_R2_SUFFIX:=_metaG_R2.fq.gz}"
: "${MIN_EFF_COV_R2:=0.05}"
: "${THREADS:=16}"

CONDA_ENV=${CONDA_ENV:-base}
if [[ -z "${MOPP2_NO_CONDA_SWITCH:-}" ]] && { [[ -z "${CONDA_DEFAULT_ENV:-}" ]] || [[ "${CONDA_DEFAULT_ENV}" != "$CONDA_ENV" ]]; }; then
    set +u
    for CB in "$HOME/miniconda3" "$HOME/anaconda3" /opt/miniconda3 /opt/anaconda3; do
        [[ -f "$CB/etc/profile.d/conda.sh" ]] && { source "$CB/etc/profile.d/conda.sh"; conda activate "$CONDA_ENV"; break; }
    done
    set -u
fi
for tool in sylph bowtie2 bowtie2-build samtools; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: $tool not found (env ${CONDA_DEFAULT_ENV:-none})"; exit 1; }
done

SYLDB="${SUBDB_SRC}/unified_db.syldb"
UNIFIED="${SUBDB_SRC}/unified_db.tsv"
[[ -f "$SYLDB"   ]] || { echo "ERROR: $SYLDB not found (build the source sub-DB first)"; exit 1; }
[[ -f "$UNIFIED" ]] || { echo "ERROR: $UNIFIED not found"; exit 1; }

mkdir -p "$OUT_DIR"; cd "$OUT_DIR"
exec > >(tee build_subdb_round2.log) 2>&1

echo "[$(date '+%F %T')] === Round-2 sub-DB build (reusing sketch) ==="
echo "  Source sketch: $SYLDB"
echo "  Unified DB:    $UNIFIED ($(tail -n +2 "$UNIFIED" | wc -l) genomes)"
echo "  metaG dir:     $METAG_DIR ($METAG_R1_SUFFIX / $METAG_R2_SUFFIX)"
echo "  Output:        $OUT_DIR"

mapfile -t SAMPLES < <(
    ls "$METAG_DIR"/*"$METAG_R1_SUFFIX" 2>/dev/null \
        | xargs -n1 basename | sed "s/${METAG_R1_SUFFIX}\$//"
)
echo "  metaG samples: ${#SAMPLES[@]}"

# ---- Stage A: round-2 sylph -------------------------------------------------
echo
echo "[$(date '+%F %T')] Stage A: round-2 sylph per sample..."
mkdir -p sylph_round2
for s in "${SAMPLES[@]}"; do
    out="sylph_round2/${s}.tsv"
    r1="${METAG_DIR}/${s}${METAG_R1_SUFFIX}"
    r2="${METAG_DIR}/${s}${METAG_R2_SUFFIX}"
    [[ -s "$out" ]] && { echo "  [$s] exists, skip"; continue; }
    [[ -f "$r1" && -f "$r2" ]] || { echo "  [$s] fastq missing, skip"; continue; }
    echo "  [$s] sylph profile..."
    sylph profile "$SYLDB" -1 "$r1" -2 "$r2" -t "$THREADS" > "$out" 2>> sylph_round2.log
done

# ---- Stage B: per-sample sub-DB --------------------------------------------
echo
echo "[$(date '+%F %T')] Stage B: per-sample sub-DB construction..."
build_one() {
    local SAMPLE=$1
    local SDIR="${OUT_DIR}/samples/${SAMPLE}"
    local SYLPH="${OUT_DIR}/sylph_round2/${SAMPLE}.tsv"
    [[ -s "$SYLPH" ]] || { echo "  [$SAMPLE] no round2 sylph, skip"; return; }
    mkdir -p "$SDIR"; cd "$SDIR"

    python3 - "$SAMPLE" "$UNIFIED" "$SYLPH" "$MIN_EFF_COV_R2" > references.tsv <<'PY'
import sys, os, csv
sample, unified_tsv, sylph_path, mincov = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
unified = {}
with open(unified_tsv) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        unified[r['ref_id']] = r
hits = []
with open(sylph_path) as f:
    rdr = csv.reader(f, delimiter='\t'); hdr = next(rdr)
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
        echo "  [$SAMPLE] WARNING: 0 refs passed threshold; no index."
    fi
    cd "$OUT_DIR"
}
for s in "${SAMPLES[@]}"; do build_one "$s"; done

# ---- Summary ----------------------------------------------------------------
echo
echo "[$(date '+%F %T')] Summary"
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
echo "[$(date '+%F %T')] ROUND-2 SUB-DB BUILD COMPLETE."
echo "  Note: functional annotation (ORFs/eggNOG/cds_features) is genome-based and"
echo "  identical to the source project; per-sample SAFs for these samples can be"
echo "  generated from \$SUBDB_SRC/cds_features.tsv + each samples/<s>/references.tsv."
