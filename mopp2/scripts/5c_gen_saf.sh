#!/bin/bash
# 5c_gen_saf.sh - write per-sample ref.saf with structural-RNA ORFs removed.
#
#   Usage:
#     5c_gen_saf.sh <mapping.config.sh> [--cds <FILE>] [--mode auto|filter|rebuild] [--no-backup]
#     5c_gen_saf.sh --subdb <DIR> --blacklist <FILE> [same options]
#
#   The config form takes SUBDB_ROOT from the config and the blacklist from
#   ORF_BLACKLIST (default $SUBDB_ROOT/structural_rna_orf_blacklist.txt), so the
#   three steps share one interface:
#     bash 5b_make_structural_rna_blacklist.sh myproject.config.sh
#     bash 5c_gen_saf.sh                       myproject.config.sh
#     bash 6_build_features.sh                 myproject.config.sh
#
# Runs after 5b_make_structural_rna_blacklist.sh and before 6_build_features.sh.
#
# WHY THIS EXISTS
#   metaRS measures ribosome-protected fragments and metaT is rRNA-depleted, so
#   both quantify mRNA. The SAF handed to featureCounts therefore means "mRNA
#   regions". A prodigal ORF predicted inside an rRNA/tRNA gene is not a
#   contaminated count to subtract later - it is a feature that should never have
#   been in the SAF. Removing it here means featureCounts never sees it, so every
#   downstream table is correct as produced.
#
#   It also matters numerically: featureCounts runs with `-O --fraction
#   -M --fraction`, so a read straddling a spurious ORF and a real gene is split
#   0.5/0.5. Deleting the spurious row from the count table afterwards discards
#   that 0.5; removing it from the SAF gives it back to the real gene.
#
# TWO MODES, because annotation runs differ in what they leave behind
#
#   filter   Rewrite the existing samples/*/ref.saf, dropping blacklisted GeneIDs.
#            Needs nothing but the SAFs themselves, so it works no matter which
#            annotation script produced them, and it preserves multi-exon
#            (fungal) features exactly as annotated.
#            THIS IS THE MODE TO USE IF YOU HAVE ALREADY RUN build_subdb AND
#            build_annotation AND JUST NEED TO GET PAST 6_build_features.sh.
#
#   rebuild  Regenerate each ref.saf from a feature table (cds_features.tsv, or
#            gene_inventory.tsv as a fallback) joined against
#            samples/<s>/references.tsv. Use this when a project reuses another
#            project's annotation, e.g.
#              --subdb /path/13.dropout/5.subdb \
#              --cds   /path/12.model/5.subdb/cds_features.tsv
#
#   auto     rebuild when a cds_features.tsv is available, otherwise filter the
#            existing SAFs, otherwise rebuild from gene_inventory.tsv.
#
# Note that only 12model.build_annotation.sh writes cds_features.tsv; the older
# 5.build_annotation.sh does not, and builds ref.saf straight from
# gene_inventory.tsv. That is exactly why filter mode exists.
#
# IDEMPOTENT: filter mode saves the untouched SAF as ref.saf.unfiltered on first
# run and always filters from that copy, so re-running never compounds.
#
# featureCounts groups rows sharing a GeneID, so several rows per gene is
# intentional and is preserved.

set -uo pipefail

SUBDB_ROOT=""; BLACKLIST=""; CDS_FEAT=""; MODE="auto"; BACKUP=1

# A leading argument that is not an option is a config file: source it and take
# the paths from there.
if [[ $# -gt 0 && "$1" != -* ]]; then
    CONFIG=$1; shift
    [[ -f "$CONFIG" ]] || { echo "ERROR: no such config: $CONFIG"; exit 1; }
    source "$CONFIG"
    SUBDB_ROOT=${SUBDB_ROOT:?config is missing SUBDB_ROOT}
    BLACKLIST=${ORF_BLACKLIST:-${SUBDB_ROOT}/structural_rna_orf_blacklist.txt}
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        --subdb)     SUBDB_ROOT=$2; shift 2 ;;
        --blacklist) BLACKLIST=$2;  shift 2 ;;
        --cds)       CDS_FEAT=$2;   shift 2 ;;
        --mode)      MODE=$2;       shift 2 ;;
        --no-backup) BACKUP=0;      shift ;;
        -h|--help)   sed -n '2,52p' "$0"; exit 0 ;;
        *) echo "ERROR: unknown argument: $1"; exit 1 ;;
    esac
done

[[ -n "$SUBDB_ROOT" ]] || { echo "ERROR: --subdb is required"; exit 1; }
[[ -d "$SUBDB_ROOT" ]] || { echo "ERROR: no such directory: $SUBDB_ROOT"; exit 1; }
[[ -n "$BLACKLIST"  ]] || { echo "ERROR: --blacklist is required (or pass a config)"; exit 1; }
[[ -s "$BLACKLIST"  ]] || { echo "ERROR: blacklist not found: $BLACKLIST"; echo "  run 5b_make_structural_rna_blacklist.sh first"; exit 1; }
case "$MODE" in auto|filter|rebuild) ;; *) echo "ERROR: --mode must be auto, filter or rebuild"; exit 1 ;; esac

N_BL=$(wc -l < "$BLACKLIST")
echo "SUBDB_ROOT: $SUBDB_ROOT"
echo "Blacklist:  $BLACKLIST ($N_BL gene_ids)"

# ---- resolve the mode --------------------------------------------------------
N_SAF=$(ls "$SUBDB_ROOT"/samples/*/ref.saf 2>/dev/null | wc -l)
if [[ -z "$CDS_FEAT" && -s "${SUBDB_ROOT}/cds_features.tsv" ]]; then
    CDS_FEAT="${SUBDB_ROOT}/cds_features.tsv"
fi
if [[ "$MODE" == "auto" ]]; then
    if   [[ -n "$CDS_FEAT" ]];  then MODE=rebuild
    elif [[ "$N_SAF" -gt 0 ]];  then MODE=filter
    else                             MODE=rebuild
    fi
    echo "Mode:       $MODE (auto)"
else
    echo "Mode:       $MODE"
fi

if [[ "$MODE" == "rebuild" ]]; then
    if [[ -z "$CDS_FEAT" ]]; then
        # The older 5.build_annotation.sh never writes cds_features.tsv. Its
        # gene_inventory.tsv carries the same columns, one row per gene, which is
        # equivalent for prokaryotes; multi-exon fungal genes collapse to a single
        # interval, so prefer filter mode when fungal genomes are involved.
        CDS_FEAT="${SUBDB_ROOT}/gene_inventory.tsv"
        echo "            no cds_features.tsv - falling back to gene_inventory.tsv"
        echo "            (multi-exon fungal features collapse to one interval;"
        echo "             use --mode filter to keep them intact)"
    fi
    [[ -s "$CDS_FEAT" ]] || { echo "ERROR: feature table not found: $CDS_FEAT"; exit 1; }
    echo "Features:   $CDS_FEAT"
fi

# ==============================================================================
# filter mode
# ==============================================================================
if [[ "$MODE" == "filter" ]]; then
    [[ "$N_SAF" -gt 0 ]] || { echo "ERROR: no samples/*/ref.saf to filter under $SUBDB_ROOT"; exit 1; }
    echo "SAF files:  $N_SAF"
    python3 - "$SUBDB_ROOT" "$BLACKLIST" "$BACKUP" <<'PY'
import glob, os, shutil, sys
root, blacklist_path, backup = sys.argv[1], sys.argv[2], sys.argv[3] == "1"
excluded = {ln.strip() for ln in open(blacklist_path) if ln.strip()}

tot_in = tot_out = 0; matched = set()
for saf in sorted(glob.glob(os.path.join(root, 'samples', '*', 'ref.saf'))):
    pristine = saf + '.unfiltered'
    # Always filter from the untouched original, so re-running is a no-op rather
    # than filtering an already-filtered file.
    if backup and not os.path.exists(pristine):
        shutil.copy2(saf, pristine)
    src = pristine if os.path.exists(pristine) else saf

    kept, dropped = [], 0
    with open(src) as f:
        header = f.readline()
        for line in f:
            gid = line.split('\t', 1)[0]
            if gid in excluded:
                dropped += 1; matched.add(gid)
            else:
                kept.append(line)
    with open(saf, 'w') as out:
        out.write(header); out.writelines(kept)
    tot_in += len(kept) + dropped; tot_out += len(kept)
    sys.stderr.write(f"  [{os.path.basename(os.path.dirname(saf))}] "
                     f"rows {len(kept)+dropped} -> {len(kept)}  (removed {dropped})\n")

sys.stderr.write(f"  total rows {tot_in} -> {tot_out} (removed {tot_in-tot_out})\n")
unmatched = len(excluded) - len(matched)
if unmatched:
    sys.stderr.write(f"  NOTE: {unmatched} of {len(excluded)} blacklisted gene_ids were not "
                     f"present in any SAF (genomes not selected for these samples)\n")
if tot_in == tot_out:
    sys.stderr.write("  WARNING: nothing was removed. Either the SAFs were already "
                     "filtered, or the blacklist was built against a different annotation.\n")
PY
fi

# ==============================================================================
# rebuild mode
# ==============================================================================
if [[ "$MODE" == "rebuild" ]]; then
    mapfile -t SAMPLES < <(for d in "$SUBDB_ROOT"/samples/*/; do s=$(basename "$d"); [[ -f "${d}references.tsv" ]] && echo "$s"; done)
    echo "Samples with references.tsv: ${#SAMPLES[@]}"
    [[ ${#SAMPLES[@]} -gt 0 ]] || { echo "ERROR: no samples/*/references.tsv under $SUBDB_ROOT"; exit 1; }

    python3 - "$SUBDB_ROOT" "$CDS_FEAT" "$BLACKLIST" "$BACKUP" "${SAMPLES[@]}" <<'PY'
import sys, os, shutil
root, cds, blacklist_path, backup = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4] == "1"
samples = sys.argv[5:]
excluded = {ln.strip() for ln in open(blacklist_path) if ln.strip()}

ref_feats = {}
n_total = n_dropped = 0
with open(cds) as f:
    h = f.readline().rstrip('\n').split('\t')
    ci = {n: h.index(n) for n in ['gene_id','ref_id','contig','start','end','strand']}
    for line in f:
        c = line.rstrip('\n').split('\t')
        if len(c) < max(ci.values())+1: continue
        n_total += 1
        gid = c[ci['gene_id']]
        if gid in excluded:
            n_dropped += 1; continue
        ref_feats.setdefault(c[ci['ref_id']], []).append(
            (gid, c[ci['contig']], c[ci['start']], c[ci['end']], c[ci['strand']]))
sys.stderr.write(f"  feature rows: {n_total}, excluded: {n_dropped}, kept: {n_total-n_dropped}\n")
if n_dropped == 0:
    sys.stderr.write("  WARNING: the blacklist matched nothing in this feature table - "
                     "it was probably built against a different annotation.\n")

for s in samples:
    sdir = os.path.join(root, 'samples', s); used = []
    with open(os.path.join(sdir, 'references.tsv')) as f:
        hh = f.readline().rstrip('\n').split('\t'); ir = hh.index('ref_id')
        for line in f:
            c = line.rstrip('\n').split('\t')
            if len(c) > ir: used.append(c[ir])
    saf = os.path.join(sdir, 'ref.saf')
    if backup and os.path.exists(saf) and not os.path.exists(saf + '.unfiltered'):
        shutil.copy2(saf, saf + '.unfiltered')
    with open(saf, 'w') as out:
        out.write("GeneID\tChr\tStart\tEnd\tStrand\n"); ng = nr = 0
        for rid in used:
            fe = ref_feats.get(rid)
            if not fe: continue
            nr += 1
            for gid, contig, st, en, strand in fe:
                out.write(f"{gid}\t{contig}\t{st}\t{en}\t{strand}\n"); ng += 1
    sys.stderr.write(f"  [{s}] refs={len(used)} with_feats={nr} rows={ng}\n")
PY
fi

# Provenance, so a stale SAF is detectable later (6_build_features.sh reports it).
{
    echo -e "generated_at\t$(date '+%F %T')"
    echo -e "mode\t$MODE"
    echo -e "features\t${CDS_FEAT:-<existing ref.saf>}"
    echo -e "blacklist\t$BLACKLIST"
    echo -e "blacklist_n\t$N_BL"
} > "${SUBDB_ROOT}/saf_provenance.tsv"

echo
echo "Done. Provenance: ${SUBDB_ROOT}/saf_provenance.tsv"
[[ "$BACKUP" -eq 1 ]] && echo "Originals kept as samples/*/ref.saf.unfiltered"
echo "Next: bash 6_build_features.sh <mapping.config.sh>"
echo "      (delete \$FEATURES_OUT/fc_per_sample/ first if featureCounts already ran)"
