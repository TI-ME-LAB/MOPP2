#!/bin/bash
# run_mapping2.sh — manifest-driven mapping (supports metaRS DB overrides)
#
# Maps every job in the generated manifest to its sub-DB bowtie2 index.
#   metaG / metaT : paired or single, bowtie2 default sensitive.
#   metaRS        : single-end RPF, mapped TWICE (strict + mopp2 score-min).
#
# Output BAMs:
#   $MAPPING_OUT/<read_id>/metaG/<read_id>.bam
#   $MAPPING_OUT/<read_id>/metaT/<read_id>.bam
#   $MAPPING_OUT/<read_id>/metaRS_strict/<read_id>.bam
#   $MAPPING_OUT/<read_id>/metaRS_mopp2/<read_id>.bam
#
# Usage:
#   bash run_mapping2.sh <mapping_config.sh>

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
: "${THREADS:=16}"
: "${SORT_THREADS:=4}"
: "${SORT_MEM:=2G}"
: "${TMP_ROOT:=}"
: "${METARS_STRICT_SCORE_MIN:=L,0,-0.01}"
: "${METARS_MOPP2_SCORE_MIN:=L,0,-0.3}"

# ---- conda env (bowtie2); samtools expected on PATH ----
CONDA_ENV=${MAP_CONDA_ENV:-${CONDA_ENV:-base}}
if [[ -z "${MOPP2_NO_CONDA_SWITCH:-}" ]] && { [[ -z "${CONDA_DEFAULT_ENV:-}" ]] || [[ "${CONDA_DEFAULT_ENV}" != "$CONDA_ENV" ]]; }; then
    set +u
    for CB in "$HOME/miniconda3" "$HOME/anaconda3" /opt/miniconda3 /opt/anaconda3; do
        [[ -f "$CB/etc/profile.d/conda.sh" ]] && { source "$CB/etc/profile.d/conda.sh"; conda activate "$CONDA_ENV"; break; }
    done
    set -u
fi
for t in bowtie2 samtools; do
    command -v "$t" >/dev/null 2>&1 || { echo "ERROR: $t not found (env ${CONDA_DEFAULT_ENV:-none})"; exit 1; }
done

mkdir -p "$MAPPING_OUT"
MANIFEST="${MAPPING_OUT}/manifest.tsv"
gen_manifest > "$MANIFEST"
exec > >(tee "${MAPPING_OUT}/mapping.log") 2>&1

N_JOBS=$(($(wc -l < "$MANIFEST")-1))
echo "[$(date '+%F %T')] === Mapping (manifest-driven) ==="
echo "  Sub-DB root: $SUBDB_ROOT"
echo "  Output:      $MAPPING_OUT"
echo "  Manifest:    $MANIFEST ($N_JOBS jobs)"
echo "  bowtie2:     $(bowtie2 --version 2>&1 | head -1)"
echo "  samtools:    $(samtools --version 2>&1 | head -1)"
echo "  Libraries:   $(tail -n +2 "$MANIFEST" | cut -f2 | sort | uniq -c | tr '\n' ' ')"

# ----------------------------------------------------------------------------
map_variant() {
    # read_id libtype db_sample paired strand r1 r2 variant score_min
    local READ=$1 LIB=$2 DB=$3 PAIRED=$4 R1=$6 R2=$7 VARIANT=$8 SCORE_MIN=$9
    local IDX="${SUBDB_ROOT}/samples/${DB}/ref"
    local OUTDIR="${MAPPING_OUT}/${READ}/${VARIANT}"
    local BAM="${OUTDIR}/${READ}.bam" LOG="${OUTDIR}/${READ}.bowtie2.log"
    local SORTLOG="${OUTDIR}/${READ}.sort.log" FLAG="${OUTDIR}/${READ}.flagstat"
    mkdir -p "$OUTDIR"
    if [[ -s "$BAM" && -s "${BAM}.bai" ]]; then echo "    [$READ/$VARIANT] exists, skip"; return 0; fi
    if [[ ! -f "${IDX}.1.bt2" && ! -f "${IDX}.1.bt2l" ]]; then echo "    [$READ/$VARIANT] DB index missing: $IDX"; return 1; fi

    local SORT_TMP
    [[ -n "$TMP_ROOT" ]] && SORT_TMP="${TMP_ROOT}/sort_${READ}_${VARIANT}_$$" || SORT_TMP="${OUTDIR}/sort_tmp"
    mkdir -p "$SORT_TMP"

    local ARGS=(-x "$IDX" --no-unal --threads "$THREADS")
    if [[ "$PAIRED" == "1" ]]; then
        ARGS+=(-1 "$R1" -2 "$R2")
    else
        ARGS+=(--sensitive -L 18 -N 1 --score-min "$SCORE_MIN" -U "$R1")
    fi

    echo "    [$READ/$VARIANT] db=$DB start $(date '+%T')"
    bowtie2 "${ARGS[@]}" 2> "$LOG" \
        | samtools sort -@ "$SORT_THREADS" -m "$SORT_MEM" -T "${SORT_TMP}/${READ}" -o "$BAM" - 2> "$SORTLOG"
    local RC=("${PIPESTATUS[@]}")
    rm -rf "$SORT_TMP"
    if [[ ${RC[0]} -ne 0 ]]; then echo "    [$READ/$VARIANT] bowtie2 FAILED"; tail -3 "$LOG"|sed 's/^/      /'; rm -f "$BAM"; return 1; fi
    if [[ ${RC[1]} -ne 0 ]]; then echo "    [$READ/$VARIANT] sort FAILED"; tail -3 "$SORTLOG"|sed 's/^/      /'; rm -f "$BAM"; return 1; fi
    samtools index -@ "$SORT_THREADS" "$BAM" 2>>"$SORTLOG" || { echo "    [$READ/$VARIANT] index FAILED"; return 1; }
    samtools flagstat "$BAM" > "$FLAG" 2>/dev/null || true
    local RATE=$(grep "overall alignment rate" "$LOG" 2>/dev/null | awk '{print $1}')
    echo "    [$READ/$VARIANT] done ${RATE:-?} $(date '+%T')"
}

# ----------------------------------------------------------------------------
while IFS=$'\t' read -r READ LIB DB PAIRED STRAND R1 R2; do
    [[ "$READ" == "read_id" ]] && continue
    case "$LIB" in
        metaG|metaT)
            map_variant "$READ" "$LIB" "$DB" "$PAIRED" "$STRAND" "$R1" "$R2" "$LIB" "default" || true
            ;;
        metaRS)
            map_variant "$READ" "$LIB" "$DB" "$PAIRED" "$STRAND" "$R1" "$R2" "metaRS_strict" "$METARS_STRICT_SCORE_MIN" || true
            map_variant "$READ" "$LIB" "$DB" "$PAIRED" "$STRAND" "$R1" "$R2" "metaRS_mopp2"  "$METARS_MOPP2_SCORE_MIN"  || true
            ;;
        *) echo "  WARN: unknown libtype $LIB for $READ" ;;
    esac
done < "$MANIFEST"

# ---- summary ----
echo
echo "[$(date '+%F %T')] === Mapping summary ==="
{
    echo -e "read_id\tvariant\talign_rate\tn_mapped"
    for sdir in "$MAPPING_OUT"/*/; do
        rid=$(basename "$sdir")
        [[ "$rid" == "manifest.tsv" ]] && continue
        for vdir in "$sdir"*/; do
            [[ -d "$vdir" ]] || continue
            v=$(basename "$vdir"); log="${vdir}${rid}.bowtie2.log"; flag="${vdir}${rid}.flagstat"
            [[ -f "$log" ]] || continue
            rate=$(grep "overall alignment rate" "$log" 2>/dev/null | awk '{print $1}')
            mapped="NA"; [[ -f "$flag" ]] && mapped=$(awk '/[[:space:]]mapped \(/ && !/primary/ {print $1; exit}' "$flag")
            echo -e "${rid}\t${v}\t${rate:-NA}\t${mapped:-NA}"
        done
    done
} | column -t -s $'\t' | tee "${MAPPING_OUT}/summary.tsv"

echo
echo "[$(date '+%F %T')] MAPPING COMPLETE."
