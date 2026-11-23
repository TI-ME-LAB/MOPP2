#!/bin/bash
# 3_host_depletion.sh - remove host reads before sub-DB construction and mapping.
#
# For a host-associated cohort (stool, skin, plant tissue) every downstream
# number - sylph abundance, bowtie2 mapping rate - is diluted by host reads by a
# per-sample amount unless they are taken out here. bowtie2 against the host
# index in $HUMAN_IDX; the reads that FAIL to align are kept.
#
# Not needed for a syncom or an environmental cohort: leave this step out.
#
#   metaG  paired  -> --un-conc-gz : a pair is kept only if it does not align
#                                    concordantly to the host
#   metaRS single  -> --un-gz       : skipped when SKIP_METARS_FILTER=1
#
# End-to-end (not --local) on purpose: metaRS reads are ~28 bp, and a full-length
# 28 bp match to a 3.1 Gb genome is specific (~4e-8 by chance), while --local
# would call a 20 bp partial match human.
#
# The SAM is discarded - only the unaligned reads are wanted.
#
#   Usage: bash 3_host_depletion.sh <config.sh>
set -uo pipefail
CONFIG=${1:?config file required}
source "$CONFIG"

: "${TRIM_DIR:?}"
: "${HF_OUT:?}"
: "${HUMAN_IDX:?}"
: "${THREADS:=24}"
# SKIP_METARS_FILTER=1 leaves metaRS alone. Worth it when the cohort turns out
# to be nearly host-free: RA's first sample was 0.67% human in metaG, and metaRS
# host reads would have to spuriously align to a bacterial CDS to reach a count
# table at all. Set METARS_DIR to the unfiltered reads when this is on.
: "${SKIP_METARS_FILTER:=0}"

CONDA_ENV=${MAP_CONDA_ENV:-base}
if [[ -z "${MOPP2_NO_CONDA_SWITCH:-}" ]] && [[ "${CONDA_DEFAULT_ENV:-}" != "$CONDA_ENV" ]]; then
    set +u
    for CB in "$HOME/miniconda3" "$HOME/anaconda3" /opt/miniconda3 /opt/anaconda3; do
        [[ -f "$CB/etc/profile.d/conda.sh" ]] && { source "$CB/etc/profile.d/conda.sh"; conda activate "$CONDA_ENV"; break; }
    done
    set -u
fi
for t in bowtie2 samtools; do command -v $t >/dev/null || { echo "ERROR: $t not found"; exit 1; }; done
[[ -f "${HUMAN_IDX}.1.bt2l" || -f "${HUMAN_IDX}.1.bt2" ]] || { echo "ERROR: no bowtie2 index at $HUMAN_IDX"; exit 1; }

OUT="${HF_OUT}/filtered"; LOGS="${HF_OUT}/logs"
mkdir -p "$OUT" "$LOGS"

# sample list = whatever metaG R1 files exist
mapfile -t SAMPLES < <(ls "$TRIM_DIR"/*"${METAG_R1_SUFFIX}" 2>/dev/null | xargs -n1 basename | sed "s/${METAG_R1_SUFFIX}\$//")
echo "[$(date '+%F %T')] samples: ${#SAMPLES[@]} -> ${SAMPLES[*]}"

SUMMARY="${HF_OUT}/human_filter_summary.tsv"
[[ -s "$SUMMARY" ]] || echo -e "sample\tlibrary\ttotal\tkept\tpct_human" > "$SUMMARY"

# A bowtie2 log only ends with "overall alignment rate" if that run finished.
# Testing the output file instead would accept a truncated one from a killed run.
finished() { grep -q 'overall alignment rate' "$1" 2>/dev/null; }

record() {  # sample lib bowtie2_log
    local s=$1 lib=$2 log=$3
    local total rate kept
    total=$(grep -oP '^\d+(?= reads; of these)' "$log" | head -1)
    rate=$(grep -oP '^[\d.]+(?=% overall alignment rate)' "$log" | head -1)
    [[ -z "$total" ]] && return
    kept=$(awk -v t="$total" -v r="${rate:-0}" 'BEGIN{printf "%d", t*(100-r)/100}')
    echo -e "${s}\t${lib}\t${total}\t${kept}\t${rate}" >> "$SUMMARY"
}

for s in "${SAMPLES[@]}"; do
    # ---- metaG (paired) ----
    r1="${TRIM_DIR}/${s}${METAG_R1_SUFFIX}"
    r2="${TRIM_DIR}/${s}${METAG_R2_SUFFIX}"
    o1="${OUT}/${s}${METAG_R1_SUFFIX}"
    o2="${OUT}/${s}${METAG_R2_SUFFIX}"
    if [[ -s "$o1" && -s "$o2" ]] && finished "${LOGS}/${s}.metaG.log"; then
        echo "[$(date '+%F %T')] [$s] metaG already filtered, skip"
    elif [[ -f "$r1" && -f "$r2" ]]; then
        if [[ -e "$o1" || -e "$o2" ]]; then
            echo "[$(date '+%F %T')] [$s] metaG output exists but its run did not finish, redoing"
            rm -f "$o1" "$o2" "${OUT}/${s}_metaG_R1.fq.gz" "${OUT}/${s}_metaG_R2.fq.gz"
        fi
        echo "[$(date '+%F %T')] [$s] metaG human filter..."
        bowtie2 -x "$HUMAN_IDX" -1 "$r1" -2 "$r2" -p "$THREADS" \
            --un-conc-gz "${OUT}/${s}_metaG_R%.fq.gz" \
            -S /dev/null 2> "${LOGS}/${s}.metaG.log"
        # bowtie2 writes _R1/_R2; rename only if the suffix differs
        [[ -f "${OUT}/${s}_metaG_R1.fq.gz" && "${OUT}/${s}_metaG_R1.fq.gz" != "$o1" ]] && mv "${OUT}/${s}_metaG_R1.fq.gz" "$o1"
        [[ -f "${OUT}/${s}_metaG_R2.fq.gz" && "${OUT}/${s}_metaG_R2.fq.gz" != "$o2" ]] && mv "${OUT}/${s}_metaG_R2.fq.gz" "$o2"
        record "$s" metaG "${LOGS}/${s}.metaG.log"
    else
        echo "[$(date '+%F %T')] [$s] no metaG, skip"
    fi

    # ---- metaRS (single) ----
    rs="${TRIM_DIR}/${s}${METARS_SUFFIX}"
    os="${OUT}/${s}${METARS_SUFFIX}"
    if [[ "$SKIP_METARS_FILTER" == "1" ]]; then
        [[ -f "$rs" ]] && echo "[$(date '+%F %T')] [$s] metaRS: SKIP_METARS_FILTER=1, left unfiltered"
    elif [[ -s "$os" ]] && finished "${LOGS}/${s}.metaRS.log"; then
        echo "[$(date '+%F %T')] [$s] metaRS already filtered, skip"
    elif [[ -f "$rs" ]]; then
        if [[ -e "$os" ]]; then
            echo "[$(date '+%F %T')] [$s] metaRS output exists but its run did not finish, redoing"
            rm -f "$os"
        fi
        echo "[$(date '+%F %T')] [$s] metaRS human filter..."
        bowtie2 -x "$HUMAN_IDX" -U "$rs" -p "$THREADS" \
            --un-gz "$os" -S /dev/null 2> "${LOGS}/${s}.metaRS.log"
        record "$s" metaRS "${LOGS}/${s}.metaRS.log"
    else
        echo "[$(date '+%F %T')] [$s] no metaRS (expected for 130_*), skip"
    fi
done

echo
echo "[$(date '+%F %T')] human filter complete. Summary:"
column -t -s $'\t' "$SUMMARY"
