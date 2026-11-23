#!/bin/bash
# 3c_build_subdb_from_globdb.sh - stage-3 entry point for a project that has
# neither a curated genome pool nor its own MAGs: the pool is discovered from
# the reads by profiling them against a pre-built GlobDB/GTDB sketch.
#
#   Which stage-3 script to use:
#     3a  the genome pool is already curated by hand          (a syncom)
#     3b  another project's unified sketch is being reused     (a dropout set)
#     3c  neither - discover the pool from GlobDB              (this script)
#
# Two-round sylph pipeline for per-sample sub-database construction.
#
# WHAT DIFFERS FROM mopp2's build_subdb_v6.sh, WHICH THIS IS DERIVED FROM
#   That script took the genome pool from two places: MAGs assembled and binned
#   from the cohort itself ($REFINEMENT_ROOT), plus GlobDB genomes selected from
#   a round-1 sylph profile that had already been run elsewhere ($SYLPH_TSV_R1).
#   A cohort with no assembly, no binning and no earlier sylph run has neither,
#   so:
#
#     * STAGE 0 (new) runs round-1 sylph itself, profiling every metaG sample
#       against the pre-built GlobDB r226 sketch, and writes the combined
#       results.tsv the original expected to be handed. It is skipped if
#       $SYLPH_TSV_R1 already exists.
#     * STAGE 1 is a no-op when $REFINEMENT_ROOT is empty: the pool is GlobDB
#       only. Everything downstream already handles a MAG-free pool, since it
#       only ever reads mag_inventory.tsv, which is then header-only. Setting
#       $REFINEMENT_ROOT turns the MAG arm back on with no other change, so a
#       cohort that is assembled later does not need a different script.
#
#   Nothing else is changed.
#
# WHAT v4 ADDS (Exp3_as). Each is off by default; a v3 config runs as v3 did.
#   * SYLPH_TIER_MAX selects genomes by sylph's own confidence instead of an
#     Eff_cov cutoff, for the pool (STAGE 3) and per sample (STAGE 7) alike:
#       T1  Eff_lambda == HIGH          (>= ~3x; ANI measured, not corrected)
#       T2  Eff_cov >= 1  and ANI 5% bound >= 95
#       T3  Eff_cov >= 0.1 and ANI 5% bound >= 95
#       T4  anything else sylph reported
#     Empty keeps MIN_EFF_COV / MIN_EFF_COV_R2.
#   * EXCLUDE_REFS keeps named ref_ids out of every sub-DB - e.g. GlobDB
#     genomes hit only by reads from a spiked-in plasmid.
#   * ISOLATE_LIST: the project's own isolate genomes (a SynCom's members).
#     Each is called present per sample from its own profile (STAGE 0b, same
#     rule as the fungus below), and where a sample carries it, a GlobDB genome
#     of the same species (skani ANI >= ISOLATE_MIN_ANI) gives way to it in that
#     sample (STAGE 3b). Isolates no sample carries stay out.
#   * FUNGAL_FASTA adds one fungal genome that GlobDB does not contain. Its
#     presence, like each isolate's, is decided per sample by sylph against the
#     genome itself (STAGE 0b), on ANI rather than tier: a 40 Mb genome gives
#     tens of thousands of matching k-mers below 1x, so Eff_cov says how deep
#     it was sequenced, not whether it is there. Present samples get it as source "Fungus", which 5a's JGI
#     gene-model stage (FUNGAL_REF_ID/FUNGAL_GFF/FUNGAL_PROT) already expects.
#   * Reruns no longer inherit the previous pool: globdb_pool/ is relinked from
#     scratch and skani/sylph outputs are reused only if their input list is
#     unchanged. Switching tiers otherwise left stale genomes in the pool.
#
# Workflow:
#   Stage 0: Round-1 sylph: every metaG sample vs the GlobDB r226 sketch  [new]
#   Stage 0b: Own genomes (isolates, fungus): presence per sample  (ISOLATE_LIST/FUNGAL_FASTA)
#   Stage 1: Collect all MAGs + CheckM quality   (skipped when no REFINEMENT_ROOT)
#   Stage 2: Read first-round sylph results
#   Stage 3: Select GlobDB genomes referenced by any sample (eff_cov >= 0.1, or tier)
#   Stage 3b: Same-species GlobDB genomes give way to a carried isolate (ISOLATE_LIST)
#   Stage 4: Build unified pool and skani all-vs-all
#   Stage 5: Dedup unified pool at ANI >= ANI_DEDUP, preserving strain band 95-99%
#            -> unified database
#   Stage 6: Sketch unified database for second-round sylph
#   Stage 7: Second-round sylph: each sample's metaG reads vs unified database
#            -> per-sample reference list (eff_cov >= 0.05)
#   Stage 8: Per-sample sub-DB construction (fasta + bowtie2 index)
#
# Usage:
#   bash 3c_build_subdb_from_globdb.sh <config.sh>

set -uo pipefail

CONFIG=${1:?config file required}
source "$CONFIG"

: "${REFINEMENT_ROOT:=}"               # empty => no MAG pool, GlobDB only
: "${SYLPH_TSV_R1:?}"                  # first-round sylph results; built by STAGE 0 if absent
: "${SYLPH_DB_R1:=}"                   # sylph sketch to profile against in STAGE 0
: "${GLOBDB_ROOT:?}"
: "${METAG_DIR:?}"                     # for second-round sylph
: "${METAG_R1_SUFFIX:=_metaG_R1.fq.gz}"
: "${METAG_R2_SUFFIX:=_metaG_R2.fq.gz}"
: "${OUT_ROOT:?}"
: "${BIN_SUBDIR:=metawrap_50_10_bins}"
: "${MIN_EFF_COV:=0.1}"
: "${MIN_EFF_COV_R2:=0.05}"
: "${ANI_DEDUP:=99.0}"
: "${ANI_STRAIN_BAND:=95.0}"
: "${ALIGN_FRAC_MIN:=30}"
: "${SYLPH_MIN_ANI:=95}"
: "${SYLPH_C:=200}"
: "${THREADS:=16}"
: "${SYLPH_TIER_MAX:=}"                # empty => Eff_cov cutoffs; T1..T4 => tier selection
: "${EXCLUDE_REFS:=}"                  # whitespace-separated ref_ids kept out of every sub-DB
: "${FUNGAL_FASTA:=}"                  # fungal genome outside GlobDB; empty => none
: "${FUNGAL_REF_ID:=}"                 # must be FUNGAL_FASTA's basename minus extension
: "${KNOWN_MIN_ANI_LO:=${FUNGAL_MIN_ANI_LO:-98}}"  # own genome present if ANI 5% bound (ANI when HIGH) >= this
: "${READ_SKETCH_DIR:=$(dirname "$SYLPH_TSV_R1")/sketches}"
: "${ISOLATE_LIST:=}"                  # the project's own isolate genomes; empty => none
: "${ISOLATE_MIN_ANI:=95}"             # same species as a selected GlobDB genome => replace it

case "$SYLPH_TIER_MAX" in
    ""|T1|T2|T3|T4) ;;
    *) echo "ERROR: SYLPH_TIER_MAX must be empty or one of T1 T2 T3 T4 (got '$SYLPH_TIER_MAX')"; exit 1 ;;
esac
if [[ -n "$SYLPH_TIER_MAX" ]]; then
    SEL_R1="tier <= $SYLPH_TIER_MAX"; SEL_R2=$SEL_R1
else
    SEL_R1="eff_cov >= $MIN_EFF_COV"; SEL_R2="eff_cov >= $MIN_EFF_COV_R2"
fi
if [[ -n "$FUNGAL_FASTA" ]]; then
    [[ -f "$FUNGAL_FASTA" ]] || { echo "ERROR: FUNGAL_FASTA not found: $FUNGAL_FASTA"; exit 1; }
    FUNGAL_FASTA=$(readlink -f "$FUNGAL_FASTA")
    fstem=$(basename "$FUNGAL_FASTA"); fstem=${fstem%.gz}
    for e in .fasta .fna .fa; do fstem=${fstem%"$e"}; done
    # 5a reads FUNGAL_REF_ID from the same config to find this genome's rows,
    # and STAGE 8 names its contigs <ref_id>__<scaffold>; the two must agree.
    if [[ "$FUNGAL_REF_ID" != "$fstem" ]]; then
        echo "ERROR: FUNGAL_REF_ID ('$FUNGAL_REF_ID') must be '$fstem', the basename of FUNGAL_FASTA"
        exit 1
    fi
fi

# Optional: conda env activation. Defaults to base, as in 3a and 5a; the env
# must provide sylph, skani and samtools for the tool check below.
CONDA_ENV=${CONDA_ENV:-base}
if [[ -z "${MOPP2_NO_CONDA_SWITCH:-}" ]] && { [[ -z "${CONDA_DEFAULT_ENV:-}" ]] || [[ "${CONDA_DEFAULT_ENV}" != "$CONDA_ENV" ]]; }; then
    set +u  
    for CONDA_BASE in /opt/miniconda3 /opt/anaconda3 "$HOME/miniconda3" "$HOME/anaconda3"; do
        if [[ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
            source "$CONDA_BASE/etc/profile.d/conda.sh"
            conda activate "$CONDA_ENV"
            break
        fi
    done
    set -u 
fi

for tool in sylph skani bowtie2 samtools bowtie2-build; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: $tool not found (env: ${CONDA_DEFAULT_ENV:-none})"
        exit 1
    fi
done

# checkm2 optional; warning if missing
if ! command -v checkm2 >/dev/null 2>&1; then
    echo "WARNING: checkm2 not found. GlobDB will use default quality_score=70."
    echo "  Install with: conda install -c bioconda checkm2"
fi

mkdir -p "$OUT_ROOT"; cd "$OUT_ROOT"
exec > >(tee pipeline.log) 2>&1

echo "[$(date '+%F %T')] === Sub-DB v6 (two-round sylph) ==="
echo "  Refinement root: ${REFINEMENT_ROOT:-<none: GlobDB-only pool>}"
echo "  First sylph:     $SYLPH_TSV_R1"
echo "  metaG dir:       $METAG_DIR"
echo "  Dedup ANI:       $ANI_DEDUP%"
echo "  Strain band:     $ANI_STRAIN_BAND%-$ANI_DEDUP% (preserved)"
echo "  Pool selection:  $SEL_R1 (for selecting GlobDB to include)"
echo "  Per sample:      $SEL_R2 (for per-sample reference selection)"
echo "  Excluded refs:   ${EXCLUDE_REFS:-<none>}"
echo "  Isolates:        ${ISOLATE_LIST:-<none>}${ISOLATE_LIST:+ (same species as GlobDB at ANI >= $ISOLATE_MIN_ANI)}"
[[ -n "$FUNGAL_FASTA$ISOLATE_LIST" ]] && \
echo "  Own genomes:     present per sample if ANI 5% bound >= $KNOWN_MIN_ANI_LO"
echo "  Fungal genome:   ${FUNGAL_FASTA:-<none>}"

# ============================================================================
# STAGE 0: Round-1 sylph (RA addition)
# ============================================================================
# Profile every metaG sample against the pre-built GlobDB r226 sketch to find
# out which genomes are in this cohort at all. The original script was handed
# this table; RA has to make it.
#
# The reads are sketched first and all 12 sketches are then profiled in ONE
# `sylph profile` call. Profiling sample by sample would load the 34 GB GlobDB
# sketch 12 times over; this loads it once. STAGE 2 reads the combined table as
# a union over samples, so nothing downstream needs them kept apart.
SKDIR="$READ_SKETCH_DIR"

# Sketch every metaG sample that has no sketch in $SKDIR yet. STAGE 0 and
# STAGE 0b both profile these; whichever runs first makes them.
sketch_reads() {
    mkdir -p "$SKDIR"
    local s n
    local -a samples=() r1=() r2=() names=()
    mapfile -t samples < <(
        ls "$METAG_DIR"/*"$METAG_R1_SUFFIX" 2>/dev/null \
            | xargs -n1 basename | sed "s/${METAG_R1_SUFFIX}\$//"
    )
    echo "  samples: ${#samples[@]}"
    [[ ${#samples[@]} -gt 0 ]] || { echo "ERROR: no metaG reads under $METAG_DIR"; exit 1; }

    for s in "${samples[@]}"; do
        if [[ -s "${SKDIR}/${s}.paired.sylsp" ]]; then
            echo "  [$s] sketch exists, skip"
        else
            r1+=("${METAG_DIR}/${s}${METAG_R1_SUFFIX}")
            r2+=("${METAG_DIR}/${s}${METAG_R2_SUFFIX}")
            names+=("$s")
        fi
    done
    if [[ ${#names[@]} -gt 0 ]]; then
        echo "  sketching ${#names[@]} samples (c=$SYLPH_C)..."
        sylph sketch -1 "${r1[@]}" -2 "${r2[@]}" -S "${names[@]}" \
            -c "$SYLPH_C" -t "$THREADS" -d "$SKDIR" \
            > "${SKDIR}/sylph_sketch.log" 2>&1 \
            || { echo "ERROR: sylph sketch failed, see ${SKDIR}/sylph_sketch.log"; exit 1; }
    fi
    n=$(ls "$SKDIR"/*.sylsp 2>/dev/null | wc -l)
    echo "  read sketches: $n"
    [[ "$n" -eq "${#samples[@]}" ]] || { echo "ERROR: expected ${#samples[@]} sketches in $SKDIR, have $n"; exit 1; }
}

if [[ -s "$SYLPH_TSV_R1" ]]; then
    echo
    echo "[$(date '+%F %T')] STAGE 0: reusing $SYLPH_TSV_R1"
else
    echo
    echo "[$(date '+%F %T')] STAGE 0: round-1 sylph vs GlobDB..."
    [[ -f "$SYLPH_DB_R1" ]] || { echo "ERROR: SYLPH_DB_R1 not found: $SYLPH_DB_R1"; exit 1; }
    R1_DIR=$(dirname "$SYLPH_TSV_R1"); mkdir -p "$R1_DIR"
    sketch_reads
    R1_SAMPLES=("$SKDIR"/*.sylsp)

    # ---- one profile call over all sketches ----
    echo "  profiling all sketches against $(basename "$SYLPH_DB_R1")..."
    sylph profile "$SYLPH_DB_R1" "$SKDIR"/*.sylsp \
        -t "$THREADS" -s 4 --minimum-ani "$SYLPH_MIN_ANI" \
        -o "${SYLPH_TSV_R1}.part" 2> "${R1_DIR}/sylph_round1.log" \
        || { echo "ERROR: sylph profile failed, see ${R1_DIR}/sylph_round1.log"; rm -f "${SYLPH_TSV_R1}.part"; exit 1; }
    mv "${SYLPH_TSV_R1}.part" "$SYLPH_TSV_R1"
    echo "  round-1 table: $(( $(wc -l < "$SYLPH_TSV_R1") - 1 )) rows over ${#R1_SAMPLES[@]} samples"
fi

# ============================================================================
# STAGE 0b: Which of the project's own genomes each sample carries (v4)
# ============================================================================
# The project's own genomes - its isolates (ISOLATE_LIST) and a fungus GlobDB
# lacks (FUNGAL_FASTA) - are profiled together against the read sketches, and
# each is called present per sample on a rule of its own, not on a tier. The
# tiers gate on Eff_cov, which measures depth: in Exp3_as a 40 Mb genome at
# 0.15x still matched ~20k k-mers with an ANI interval of 99.1-99.4. Presence is
# whether as many k-mers turn up as that exact genome would give at that depth,
# which is what sylph's coverage-adjusted ANI reports, so the rule is on its 5%
# bound (KNOWN_MIN_ANI_LO).
#
# A genome not called present stays out of that sample's sub-DB even when the
# experiment put it in every pot: one that is not there only gives metaRS's
# ~28 bp reads a wrong place to land.
KNOWN_DIR="${OUT_ROOT}/known"
KNOWN_PROFILE="${KNOWN_DIR}/profile.tsv"
KNOWN_PRESENCE="${KNOWN_DIR}/presence.tsv"
mkdir -p isolates_pool "$KNOWN_DIR"
find isolates_pool -maxdepth 1 -type l -delete
if [[ -n "$ISOLATE_LIST" ]]; then
    [[ -s "$ISOLATE_LIST" ]] || { echo "ERROR: ISOLATE_LIST missing or empty: $ISOLATE_LIST"; exit 1; }
    while IFS=$'\t' read -r ipath iid _; do
        [[ -z "$ipath" || "$ipath" == \#* ]] && continue
        [[ -f "$ipath" ]] || { echo "ERROR: isolate genome not found: $ipath"; exit 1; }
        if [[ -z "${iid:-}" ]]; then
            iid=$(basename "$ipath"); iid=${iid%.gz}
            for e in .fasta .fna .fa; do iid=${iid%"$e"}; done
        fi
        ext=.fa; [[ "$ipath" == *.gz ]] && ext=.fa.gz
        ln -sf "$(readlink -f "$ipath")" "isolates_pool/${iid}${ext}"
    done < "$ISOLATE_LIST"
fi
printf 'sample\tref_id\tpresent\tAdjusted_ANI\tANI_5-95_percentile\tEff_cov\tEff_lambda\tContainment_ind\n' > "$KNOWN_PRESENCE"
N_FUN_PRESENT=0
if [[ -n "$FUNGAL_FASTA" || -n "$ISOLATE_LIST" ]]; then
    echo
    echo "[$(date '+%F %T')] STAGE 0b: presence of the project's own genomes..."
    {
        ls isolates_pool/* 2>/dev/null
        [[ -n "$FUNGAL_FASTA" ]] && echo "$FUNGAL_FASTA"
    } > "${KNOWN_DIR}/genomes.list.new"
    echo "  genomes: $(wc -l < "${KNOWN_DIR}/genomes.list.new")"
    if [[ -s "$KNOWN_PROFILE" ]] && cmp -s "${KNOWN_DIR}/genomes.list.new" "${KNOWN_DIR}/genomes.list"; then
        rm -f "${KNOWN_DIR}/genomes.list.new"
        echo "  reusing $KNOWN_PROFILE (same genomes)"
    else
        mv "${KNOWN_DIR}/genomes.list.new" "${KNOWN_DIR}/genomes.list"
        sketch_reads
        sylph sketch -l "${KNOWN_DIR}/genomes.list" -c "$SYLPH_C" -t "$THREADS" \
            -o "${KNOWN_DIR}/known" > "${KNOWN_DIR}/sketch.log" 2>&1 \
            || { echo "ERROR: sylph sketch of the known genomes failed, see ${KNOWN_DIR}/sketch.log"; exit 1; }
        sylph profile "${KNOWN_DIR}/known.syldb" "$SKDIR"/*.sylsp \
            -t "$THREADS" --minimum-ani "$SYLPH_MIN_ANI" \
            -o "${KNOWN_PROFILE}.part" 2> "${KNOWN_DIR}/profile.log" \
            || { echo "ERROR: sylph profile of the known genomes failed, see ${KNOWN_DIR}/profile.log"; rm -f "${KNOWN_PROFILE}.part"; exit 1; }
        mv "${KNOWN_PROFILE}.part" "$KNOWN_PROFILE"
    fi
    # Decided afresh every run, so changing KNOWN_MIN_ANI_LO needs no re-profile.
    python3 - "$KNOWN_PROFILE" "${KNOWN_DIR}/genomes.list" "$KNOWN_MIN_ANI_LO" "$SKDIR" >> "$KNOWN_PRESENCE" <<'PY'
import sys, os, csv, glob
prof, glist, min_lo, skdir = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]

def sample_of(p):
    b = os.path.basename(p)
    for e in ('.paired.sylsp', '.sylsp'):
        if b.endswith(e):
            return b[:-len(e)]
    return b

def gid(p):
    b = os.path.basename(p)
    for e in ('.fa.gz', '.fna.gz', '.fasta.gz', '.fa', '.fna', '.fasta'):
        if b.endswith(e):
            return b[:-len(e)]
    return b

samples = sorted(sample_of(p) for p in glob.glob(os.path.join(skdir, '*.sylsp')))
genomes = [gid(l.strip()) for l in open(glist) if l.strip()]
rows = {}
with open(prof) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        rows[(sample_of(r['Sample_file']), gid(r['Genome_file']))] = r

for g in genomes:
    hit = []
    for s in samples:
        r = rows.get((s, g))
        if r is None:                      # nothing passed sylph's --minimum-ani
            print(f"{s}\t{g}\tno\tNA\tNA\tNA\tNA\tNA")
            continue
        if r['Eff_lambda'] == 'HIGH':      # measured directly, no interval reported
            lo = float(r['Adjusted_ANI'])
        else:
            try:
                lo = float(r['ANI_5-95_percentile'].split('-')[0])
            except ValueError:
                lo = float('nan')
        present = 'yes' if lo >= min_lo else 'no'
        if present == 'yes':
            hit.append(s)
        print(f"{s}\t{g}\t{present}\t{r['Adjusted_ANI']}\t{r['ANI_5-95_percentile']}\t"
              f"{r['Eff_cov']}\t{r['Eff_lambda']}\t{r['Containment_ind']}")
    sys.stderr.write(f"  {g:<32} {len(hit)}/{len(samples)}  {' '.join(hit)}\n")
PY
    [[ -n "$FUNGAL_FASTA" ]] && \
        N_FUN_PRESENT=$(awk -F'\t' -v f="$FUNGAL_REF_ID" 'NR>1 && $2==f && $3=="yes"' "$KNOWN_PRESENCE" | wc -l)
fi

# ============================================================================
# STAGE 1: Collect MAGs + CheckM quality
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 1: Collecting MAGs..."
mkdir -p mags_pool

echo -e "ref_id\tfolder\toriginal_path\tcompleteness\tcontamination\tquality_score" > mag_inventory.tsv

# RA sets REFINEMENT_ROOT empty: no bins to collect, mag_inventory.tsv stays
# header-only and the pool is GlobDB alone.
for sample_dir in ${REFINEMENT_ROOT:+"$REFINEMENT_ROOT"/*/}; do
    folder=$(basename "$sample_dir")
    bins_dir="${sample_dir}${BIN_SUBDIR}"
    [[ -d "$bins_dir" ]] || continue
    stats_file="${sample_dir}${BIN_SUBDIR}.stats"

    for fa in "$bins_dir"/*.fa; do
        [[ -f "$fa" ]] || continue
        bin_basename=$(basename "$fa")
        bin_stem="${bin_basename%.fa}"
        ref_id="${folder}__${bin_stem}"
        ln -sf "$(readlink -f "$fa")" "mags_pool/${ref_id}.fa"

        comp="NA"; cont="NA"; qs="NA"
        if [[ -f "$stats_file" ]]; then
            line=$(awk -v b="$bin_stem" -F'\t' '$1==b' "$stats_file" || true)
            if [[ -n "$line" ]]; then
                comp=$(echo "$line" | cut -f2)
                cont=$(echo "$line" | cut -f3)
                qs=$(awk -v c="$comp" -v x="$cont" 'BEGIN{printf "%.2f", c - 5*x}')
            fi
        fi
        echo -e "${ref_id}\t${folder}\t${fa}\t${comp}\t${cont}\t${qs}" >> mag_inventory.tsv
    done
done

N_MAG=$(tail -n +2 mag_inventory.tsv | wc -l)
echo "  MAGs collected: $N_MAG"

# ============================================================================
# STAGE 2 & 3: Parse first-round sylph -> determine which GlobDB to include
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 2&3: Reading first-round sylph + selecting GlobDB..."

python3 - "$SYLPH_TSV_R1" "$MIN_EFF_COV" "$SYLPH_TIER_MAX" "$EXCLUDE_REFS" selected_pairs.tsv \
    > globdb_to_include.txt <<'PY'
import sys, os, csv
tsv, mincov, tier_max, exclude = sys.argv[1], float(sys.argv[2]), sys.argv[3], set(sys.argv[4].split())
pairs_out = open(sys.argv[5], 'w')       # which sample selected which genome, for STAGE 3b
pairs_out.write('Sample_file\tref_id\n')

def gid(p):
    b = os.path.basename(p)
    for e in ('.fa.gz', '.fna.gz', '.fa', '.fna', '.fasta', '.gz'):
        if b.endswith(e):
            return b[:-len(e)]
    return b

# Same rule as STAGE 7; see SYLPH_TIER_MAX in the config.
def tier(row, ci):
    if row[ci['Eff_lambda']] == 'HIGH':
        # profile's --minimum-ani already keeps these >= 95; a query table does not
        return 1 if float(row[ci['Adjusted_ANI']]) >= 95 else 4
    cov = float(row[ci['Eff_cov']])
    try:
        lo = float(row[ci['ANI_5-95_percentile']].split('-')[0])
    except ValueError:
        lo = float('nan')
    if cov < 0.1 or not lo >= 95:
        return 4
    return 2 if cov >= 1 else 3

genomes, excluded = set(), set()
with open(tsv) as f:
    rdr = csv.reader(f, delimiter='\t')
    hdr = next(rdr)
    ci = {n: hdr.index(n) for n in ['Genome_file', 'Eff_cov', 'Eff_lambda', 'ANI_5-95_percentile', 'Adjusted_ANI']}
    for row in rdr:
        try:
            keep = (tier(row, ci) <= int(tier_max[1:])) if tier_max else \
                   (float(row[ci['Eff_cov']]) >= mincov)
        except (ValueError, IndexError):
            continue
        if not keep:
            continue
        if gid(row[ci['Genome_file']]) in exclude:
            excluded.add(gid(row[ci['Genome_file']]))
            continue
        genomes.add(row[ci['Genome_file']])
        pairs_out.write(f"{row[hdr.index('Sample_file')]}\t{gid(row[ci['Genome_file']])}\n")

pairs_out.close()
for g in sorted(genomes):
    print(g)
if excluded:
    sys.stderr.write(f"  excluded by EXCLUDE_REFS: {' '.join(sorted(excluded))}\n")
missing = exclude - excluded
if missing:
    sys.stderr.write(f"  EXCLUDE_REFS not selected anyway: {' '.join(sorted(missing))}\n")
PY

N_GDB_NEED=$(wc -l < globdb_to_include.txt)
echo "  Unique GlobDB to include ($SEL_R1): $N_GDB_NEED"

# Relink from scratch. Links left by an earlier run with a looser selection
# would otherwise stay in the pool, since STAGE 5 takes every file in here.
mkdir -p globdb_pool
find globdb_pool -maxdepth 1 -type l -delete

# sylph reports Genome_file as a path relative to the directory the sketch was
# built in, e.g. "globdb_r226_genome_fasta/R226CHUNK013/GCF_000742135.fa.gz", so
# GLOBDB_ROOT has to be the directory that path is relative to. Get it wrong and
# every genome misses, which used to mean one `find` over the whole tree PER
# GENOME - 306k files in GlobDB r226, so thousands of genomes would take days.
# Resolve directly first; only if something actually misses, index the tree once.
link_one() {   # link_one <resolved_path> <rel>
    local src=$1 rel=$2 bn gid
    bn=$(basename "$rel")
    gid="${bn%.gz}"; gid="${gid%.fa}"; gid="${gid%.fna}"; gid="${gid%.fasta}"
    ln -sf "$(readlink -f "$src")" "globdb_pool/${gid}.fa.gz"
}

: > .globdb_missing
while IFS= read -r rel; do
    [[ -z "$rel" ]] && continue
    if   [[ -f "$GLOBDB_ROOT/$rel" ]]; then link_one "$GLOBDB_ROOT/$rel" "$rel"
    elif [[ -f "$rel"             ]]; then link_one "$rel"               "$rel"
    else echo "$rel" >> .globdb_missing
    fi
done < globdb_to_include.txt

N_MISS=$(wc -l < .globdb_missing)
if [[ "$N_MISS" -gt 0 ]]; then
    echo "  $N_MISS genomes did not resolve under GLOBDB_ROOT=$GLOBDB_ROOT"
    echo "  indexing the GlobDB tree once to find them..."
    find "$GLOBDB_ROOT" \( -name '*.fa.gz' -o -name '*.fna.gz' -o -name '*.fa' -o -name '*.fna' \) \
        -printf '%f\t%p\n' 2>/dev/null > .globdb_index
    echo "  indexed $(wc -l < .globdb_index) files"
    # One pass over the index, not one scan per missing genome: a per-genome
    # lookup would be O(missing x 306k) and turn the safety net into its own
    # bottleneck.
    awk -F'\t' 'NR==FNR {p[$1]=$2; next}
                { n=split($0,a,"/"); b=a[n]; if (b in p) print $0"\t"p[b] }' \
        .globdb_index .globdb_missing > .globdb_found
    while IFS=$'\t' read -r rel src; do
        [[ -n "$src" && -f "$src" ]] && link_one "$src" "$rel"
    done < .globdb_found
    echo "  recovered $(wc -l < .globdb_found) of $N_MISS by basename"
fi
rm -f .globdb_missing .globdb_index .globdb_found

N_GDB=$(ls globdb_pool/*.fa.gz 2>/dev/null | wc -l)
echo "  GlobDB linked: $N_GDB / $N_GDB_NEED"
# Losing genomes here is silent damage - they simply never enter the pool, and
# the sub-DBs are quietly smaller than the profile said they should be.
if [[ "$N_GDB" -ne "$N_GDB_NEED" ]]; then
    echo "ERROR: $((N_GDB_NEED - N_GDB)) of $N_GDB_NEED GlobDB genomes could not be located."
    echo "       Check GLOBDB_ROOT: sylph's Genome_file paths are relative to it."
    exit 1
fi

# ============================================================================
# STAGE 3b: Same-species GlobDB genomes give way to the isolate (v4)
# ============================================================================
# Where a sample carries one of the project's isolates (STAGE 0b), a GlobDB
# genome of the same species (skani ANI >= ISOLATE_MIN_ANI) that the selection
# also picked there is only a stand-in for it: in Exp3_as the samples were
# 99.6-100% ANI to the isolates against 96.2-97.6% to stand-ins such as
# Chitinophaga polysaccharea_A. It is dropped from that sample. Where the
# isolate is not present, the GlobDB genome is a background strain like any
# other and stays on the selection rule. GlobDB genomes left with no sample, and
# isolates that no sample carries, leave the pool.
SAME_SP="${OUT_ROOT}/isolate_same_species.tsv"
printf 'globdb\tisolate\tANI\tAlign_fraction_ref\tAlign_fraction_query\n' > "$SAME_SP"
if [[ -n "$ISOLATE_LIST" ]]; then
    echo
    echo "[$(date '+%F %T')] STAGE 3b: isolates vs same-species GlobDB genomes..."
    ls isolates_pool/* > isolates.list
    ls globdb_pool/*.fa.gz > .globdb_pool.list
    echo "  isolates: $(wc -l < isolates.list), selected GlobDB: $(wc -l < .globdb_pool.list)"
    skani dist --ql isolates.list --rl .globdb_pool.list -t "$THREADS" \
        --min-af "$ALIGN_FRAC_MIN" -s "$ISOLATE_MIN_ANI" \
        > isolate_vs_globdb.tsv 2> skani_isolates.log \
        || { echo "ERROR: skani isolates vs GlobDB failed, see skani_isolates.log"; exit 1; }
    rm -f .globdb_pool.list isolates.list

    python3 - isolate_vs_globdb.tsv "$ISOLATE_MIN_ANI" selected_pairs.tsv "$KNOWN_PRESENCE" \
        "$SAME_SP" > .prune.tsv <<'PY'
import sys, os, csv
tsv, min_ani, pairs, presence, same_out = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5]

def to_id(p):
    b = os.path.basename(p)
    for e in ('.fa.gz', '.fa'):
        if b.endswith(e):
            return b[:-len(e)]
    return b

# skani's -s screen is approximate, so the threshold is applied again here.
# A GlobDB genome close to two isolates goes to the closer one.
best = {}
with open(tsv) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        a = float(r['ANI'])
        if a < min_ani:
            continue
        g = to_id(r['Ref_file'])
        if g not in best or a > best[g][1]:
            best[g] = (to_id(r['Query_file']), a, r['Align_fraction_ref'], r['Align_fraction_query'])
with open(same_out, 'a') as o:
    for g, (i, a, afr, afq) in sorted(best.items()):
        o.write(f"{g}\t{i}\t{a:.2f}\t{afr}\t{afq}\n")

present, samples, isolates = set(), set(), set()
with open(presence) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        samples.add(r['sample'])
        isolates.add(r['ref_id'])
        if r['present'] == 'yes':
            present.add((r['sample'], r['ref_id']))
by_len = sorted(samples, key=len, reverse=True)
def sample_of(path):
    b = os.path.basename(path)
    for e in ('.paired.sylsp', '.sylsp'):
        if b.endswith(e):
            b = b[:-len(e)]
            break
    return b if b in samples else next((s for s in by_len if b.startswith(s)), b)

needed, shadowed = set(), 0
with open(pairs) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        s, g = sample_of(r['Sample_file']), r['ref_id']
        i = best.get(g, (None,))[0]
        if i and (s, i) in present:
            shadowed += 1
        else:
            needed.add(g)
for g, (i, a, _, _) in sorted(best.items()):
    fate = 'kept for samples without the isolate' if g in needed else 'dropped'
    sys.stderr.write(f"  {g:<18} ~ {i:<26} ANI {a:6.2f}  {fate}\n")
sys.stderr.write(f"  {shadowed} sample-genome selections handed to the isolate\n")
for g in sorted(set(best) - needed):
    print(f"globdb\t{g}")
carried = {i for (_, i) in present}
in_pool = {to_id(f) for f in os.listdir('isolates_pool')}   # the fungus is in presence too
for i in sorted((isolates & in_pool) - carried):
    print(f"isolate\t{i}")
PY
    while IFS=$'\t' read -r kind id; do
        if [[ "$kind" == globdb ]]; then rm -f "globdb_pool/${id}.fa.gz"
        else rm -f "isolates_pool/${id}.fa" "isolates_pool/${id}.fa.gz"; fi
    done < .prune.tsv
    N_GONE=$(awk -F'\t' '$1=="isolate"' .prune.tsv | wc -l)
    rm -f .prune.tsv
    N_GDB=$(ls globdb_pool/*.fa.gz 2>/dev/null | wc -l)
    echo "  isolates in no sample, left out: $N_GONE; GlobDB left: $N_GDB"
fi
N_ISO=$(ls isolates_pool/ 2>/dev/null | wc -l)

# ============================================================================
# STAGE 4: Global skani on (all MAG + isolates + selected GlobDB)
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 4: Global skani ANI..."
echo "  Total genomes: $((N_MAG + N_ISO + N_GDB))"

{
    ls mags_pool/*.fa 2>/dev/null
    ls isolates_pool/* 2>/dev/null
    ls globdb_pool/*.fa.gz 2>/dev/null
} > all_genomes.list.new
# Reuse only for the same pool: after a change of selection, old pairs would
# leave the new genomes undeduplicated.
if [[ -s skani_all.tsv ]] && cmp -s all_genomes.list.new all_genomes.list; then
    rm -f all_genomes.list.new
    echo "  Reusing skani_all.tsv (same genome list)"
else
    mv all_genomes.list.new all_genomes.list

    skani triangle -t "$THREADS" --sparse -l all_genomes.list \
        --min-af "$ALIGN_FRAC_MIN" \
        -s "$ANI_STRAIN_BAND" \
        > skani_all.tsv 2> skani.log
    echo "  Skani done: $(wc -l < skani_all.tsv) pairs at ANI >= ${ANI_STRAIN_BAND}%"
fi

# ============================================================================
# STAGE 5: Dedup unified pool at ANI_DEDUP, using genome length as quality proxy
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 5: Unified pool dedup (length-based ranking)..."

python3 - "$OUT_ROOT" "$ANI_DEDUP" > unified_db.tsv <<'PY'
import sys, os, csv, gzip

root, ani_dedup = sys.argv[1], float(sys.argv[2])

def genome_length(path):
    """Sum of all sequence lengths in a fasta (.fa or .fa.gz)"""
    opener = gzip.open if path.endswith('.gz') else open
    n = 0
    try:
        with opener(path, 'rt') as f:
            for line in f:
                if not line.startswith('>'):
                    n += len(line.rstrip())
    except Exception:
        return 0
    return n

# Load MAGs (own MetaWRAP MAGs)
mag_inv = {}
with open(os.path.join(root, 'mag_inventory.tsv')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        mag_inv[r['ref_id']] = r

# Build candidates with length
candidates = {}
print(f"  Computing genome lengths...", file=sys.stderr)
for mid, info in mag_inv.items():
    path = os.path.join(root, 'mags_pool', f'{mid}.fa')
    candidates[mid] = {
        'ref_id': mid,
        'source': 'MAG',
        'folder': info['folder'],
        'path': path,
        'length': genome_length(path),
    }

for f in os.listdir(os.path.join(root, 'globdb_pool')):
    if not f.endswith('.fa.gz'): continue
    gid = f[:-6]
    path = os.path.join(root, 'globdb_pool', f)
    candidates[gid] = {
        'ref_id': gid,
        'source': 'GlobDB',
        'folder': 'GlobDB',
        'path': path,
        'length': genome_length(path),
    }

# The project's own isolates (STAGE 3b). Only those that replaced a GlobDB
# genome are in here.
iso_dir = os.path.join(root, 'isolates_pool')
for f in (os.listdir(iso_dir) if os.path.isdir(iso_dir) else []):
    iid = f[:-6] if f.endswith('.fa.gz') else f[:-3] if f.endswith('.fa') else None
    if iid is None: continue
    path = os.path.join(iso_dir, f)
    candidates[iid] = {
        'ref_id': iid,
        'source': 'Isolate',
        'folder': 'Isolate',
        'path': path,
        'length': genome_length(path),
    }

print(f"  Total candidates: {len(candidates)}", file=sys.stderr)

# Skani edges
dedup_edges = []
strain_edges = []
def to_id(p):
    bn = os.path.basename(p)
    for ext in ('.fa.gz','.fna.gz','.fasta.gz','.fa','.fna','.fasta','.gz'):
        if bn.endswith(ext): bn = bn[:-len(ext)]
    return bn

with open(os.path.join(root, 'skani_all.tsv')) as f:
    rdr = csv.reader(f, delimiter='\t')
    hdr = next(rdr)
    i_r, i_q, i_a = hdr.index('Ref_file'), hdr.index('Query_file'), hdr.index('ANI')
    for row in rdr:
        if not row or len(row) <= max(i_r, i_q, i_a): continue
        a = float(row[i_a])
        ra, rb = to_id(row[i_r]), to_id(row[i_q])
        if ra == rb: continue
        if ra not in candidates or rb not in candidates: continue
        # Isolates are never collapsed: STAGE 0b called each of them present
        # on its own, and a >= 99% GlobDB relative left in the pool is there for
        # samples that do not carry the isolate (STAGE 3b).
        iso_pair = 'Isolate' in (candidates[ra]['source'], candidates[rb]['source'])
        if a >= ani_dedup and not iso_pair:
            dedup_edges.append((ra, rb))
        else:
            strain_edges.append((ra, rb, a))

# Union-find on dedup edges
parent = {r: r for r in candidates}
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb: parent[ra] = rb
for ra, rb in dedup_edges:
    union(ra, rb)

clusters = {}
for r in candidates:
    clusters.setdefault(find(r), []).append(r)

# Pick representative: an isolate first, then the longest genome; tiebreak by ref_id
def rank(r):
    c = candidates[r]
    # an isolate is the organism itself; never trade it for a longer relative
    return (c['source'] != 'Isolate', -c['length'], c['ref_id'])

keep = set()
dedup_log = []
for rep, members in clusters.items():
    if len(members) == 1:
        keep.add(members[0]); continue
    members_sorted = sorted(members, key=rank)
    keep.add(members_sorted[0])
    for m in members_sorted[1:]:
        kept = members_sorted[0]
        dedup_log.append((kept, m, candidates[kept]['length'], candidates[m]['length']))

with open(os.path.join(root, 'unified_dedup_dropped.tsv'), 'w') as out:
    out.write("kept\tdropped\tkept_length\tdropped_length\n")
    for k, d, kl, dl in dedup_log:
        out.write(f"{k}\t{d}\t{kl}\t{dl}\n")

with open(os.path.join(root, 'unified_strain_band.tsv'), 'w') as out:
    out.write("ref_a\tref_b\tani\n")
    for ra, rb, a in strain_edges:
        if ra in keep and rb in keep:
            out.write(f"{ra}\t{rb}\t{a:.2f}\n")

print("ref_id\tsource\tfolder\tpath\tlength")
for r in sorted(keep, key=lambda x: (candidates[x]['source'], -candidates[x]['length'], x)):
    c = candidates[r]
    print(f"{r}\t{c['source']}\t{c['folder']}\t{c['path']}\t{c['length']}")
PY

# The fungus joins after dedup: skani against bacteria would find nothing, and
# it is not in globdb_pool. Only when STAGE 0b found it somewhere - a genome no
# sample carries would still get JGI gene models injected by 5a.
if [[ -n "$FUNGAL_FASTA" ]]; then
    if [[ "$N_FUN_PRESENT" -gt 0 ]]; then
        if [[ "$FUNGAL_FASTA" == *.gz ]]; then FCAT=zcat; else FCAT=cat; fi
        flen=$($FCAT "$FUNGAL_FASTA" | awk '!/^>/ {n+=length($0)} END{print n+0}')
        printf '%s\tFungus\tFungus\t%s\t%s\n' "$FUNGAL_REF_ID" "$FUNGAL_FASTA" "$flen" >> unified_db.tsv
    else
        echo "  WARNING: $FUNGAL_REF_ID is in no sample; left out of unified_db.tsv."
        echo "           Clear FUNGAL_REF_ID before 5a, or it injects gene models for an absent genome."
    fi
fi

N_UNIFIED=$(tail -n +2 unified_db.tsv | wc -l)
N_DEDUP_DROP=$(tail -n +2 unified_dedup_dropped.tsv | wc -l)
N_MAG_KEPT=$(awk -F'\t' 'NR>1 && $2=="MAG"' unified_db.tsv | wc -l)
N_GDB_KEPT=$(awk -F'\t' 'NR>1 && $2=="GlobDB"' unified_db.tsv | wc -l)
N_FUN_KEPT=$(awk -F'\t' 'NR>1 && $2=="Fungus"' unified_db.tsv | wc -l)
N_ISO_KEPT=$(awk -F'\t' 'NR>1 && $2=="Isolate"' unified_db.tsv | wc -l)
echo "  Unified DB: $N_UNIFIED refs ($N_MAG_KEPT MAG + $N_ISO_KEPT Isolate + $N_GDB_KEPT GlobDB + $N_FUN_KEPT Fungus), dedup-dropped: $N_DEDUP_DROP"

# ============================================================================
# STAGE 6: Sketch unified database for second-round sylph
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 6: Sketching unified database for sylph..."

# The fungus stays out: STAGE 0b already profiled it, and STAGE 7 takes its rows
# from there. In round 2 it would come back as a second, conflicting row.
awk -F'\t' 'NR>1 && $2!="Fungus" {print $4}' unified_db.tsv > unified_paths.txt.new
if [[ -f unified_db.syldb ]] && cmp -s unified_paths.txt.new unified_paths.txt; then
    rm -f unified_paths.txt.new
    echo "  Reusing unified_db.syldb (same genome list)"
else
    mv unified_paths.txt.new unified_paths.txt
    sylph sketch -t "$THREADS" -c "$SYLPH_C" \
        -o unified_db \
        -l unified_paths.txt 2>&1 | tee sylph_sketch.log | tail -5
    echo "  Sketch done: unified_db.syldb"
fi

# ============================================================================
# STAGE 7: Per-sample reference lists
# ============================================================================
# WITH MAGs, this has to be a second sylph round. The cohort's own MAGs are not
# in GlobDB, so round 1 cannot say which sample contains which MAG; only
# profiling against the merged MAG+GlobDB pool can.
#
# WITHOUT MAGs there is nothing new to profile against. The pool IS the round-1
# hits, so each sample's sub-DB is just that sample's own round-1 rows. Measured
# on a 48-sample GlobDB profile: 1163 genomes detected, of which exactly 2 formed
# a >=99% ANI pair, so dedup took the pool from 1163 to 1161. Re-profiling every
# sample against a database that differs from round 1 by two genomes is work for
# nothing - `sylph profile` is already non-redundant within a sample (no genome
# appears twice), which is the property a second round would have supplied.
#
# The one thing dedup still buys is consistency: if two samples each picked a
# different member of a >=99% cluster, they would carry two ref_ids - and two
# disjoint blocks of genes - for one organism, and gene-level comparison across
# them would break. That hit 1 cluster in 1163 genomes, so it is not a reason to
# re-profile, but it is cheap to repair by rewriting dropped ids to the
# representative STAGE 5 kept. That is what the no-MAG branch does.
#
# Either branch writes ${REF_DIR}/<sample>.tsv in round-1's own column layout,
# which is what STAGE 8 reads.
REF_DIR="${OUT_ROOT}/sample_references"
mkdir -p "$REF_DIR"

mapfile -t SAMPLES < <(
    ls "$METAG_DIR"/*"$METAG_R1_SUFFIX" 2>/dev/null \
        | xargs -n1 basename \
        | sed "s/${METAG_R1_SUFFIX}\$//"
)
echo "  Samples found in metaG dir: ${#SAMPLES[@]}"

if [[ -n "$REFINEMENT_ROOT" ]]; then
    # ---------------- MAGs in the pool: a real second round ----------------
    echo "[$(date '+%F %T')] STAGE 7: second-round sylph (pool contains MAGs)..."
    N_TODO=0
    for sample in "${SAMPLES[@]}"; do
        [[ -s "${REF_DIR}/${sample}.tsv" ]] || N_TODO=$((N_TODO+1))
    done

    if [[ "$N_TODO" -eq 0 ]]; then
        echo "  all ${#SAMPLES[@]} samples already done, skipping"
    elif [[ -d "$SKDIR" ]] && [[ $(ls "$SKDIR"/*.sylsp 2>/dev/null | wc -l) -gt 0 ]]; then
        # STAGE 0 already sketched the reads; profiling them again from fastq
        # would repeat work that is on disk.
        echo "  reusing STAGE 0 sketches: one profile call over $(ls "$SKDIR"/*.sylsp | wc -l) samples"
        sylph profile unified_db.syldb "$SKDIR"/*.sylsp -t "$THREADS" -s 4 \
            > "${REF_DIR}/all_samples.tsv" 2>> sylph_round2.log \
            || { echo "ERROR: round-2 sylph failed, see sylph_round2.log"; exit 1; }
    else
        echo "  no STAGE 0 sketches, profiling each sample from its fastq"
        : > "${REF_DIR}/all_samples.tsv"
        for sample in "${SAMPLES[@]}"; do
            r1="${METAG_DIR}/${sample}${METAG_R1_SUFFIX}"
            r2="${METAG_DIR}/${sample}${METAG_R2_SUFFIX}"
            [[ -f "$r1" && -f "$r2" ]] || { echo "  [$sample] metaG fastq missing, skip"; continue; }
            echo "  [$sample] running sylph profile..."
            sylph profile unified_db.syldb -1 "$r1" -2 "$r2" -t "$THREADS" \
                >> "${REF_DIR}/all_samples.tsv" 2>> sylph_round2.log
        done
        # strip the repeated headers the appends produced
        awk 'NR==1 || $0 !~ /^Sample_file\t/' "${REF_DIR}/all_samples.tsv" > "${REF_DIR}/.t" \
            && mv "${REF_DIR}/.t" "${REF_DIR}/all_samples.tsv"
    fi
    SPLIT_SRC="${REF_DIR}/all_samples.tsv"
    DEDUP_MAP=""          # round 2 profiled the deduplicated pool directly
else
    # ---------------- no MAGs: reuse round 1 ----------------
    echo "[$(date '+%F %T')] STAGE 7: per-sample references from round 1 (no MAGs in pool)..."
    SPLIT_SRC="$SYLPH_TSV_R1"
    DEDUP_MAP="${OUT_ROOT}/unified_dedup_dropped.tsv"
fi

python3 - "$SPLIT_SRC" "$REF_DIR" "$MIN_EFF_COV_R2" "$DEDUP_MAP" \
    "$SYLPH_TIER_MAX" "$EXCLUDE_REFS" "$KNOWN_PROFILE" "$KNOWN_PRESENCE" "$SAME_SP" \
    "${SAMPLES[@]}" <<'PYREF'
import sys, os, csv
src, outdir, mincov, dedup_map = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]
tier_max, exclude = sys.argv[5], set(sys.argv[6].split())
known_profile, known_presence, same_sp_tsv = sys.argv[7], sys.argv[8], sys.argv[9]
samples = sys.argv[10:]

def gid(p):
    b = os.path.basename(p)
    for e in ('.fa.gz', '.fna.gz', '.fa', '.fna', '.fasta', '.gz'):
        if b.endswith(e):
            return b[:-len(e)]
    return b

# Same rule as STAGE 2&3; see SYLPH_TIER_MAX in the config.
def tier(row, ci):
    if row[ci['Eff_lambda']] == 'HIGH':
        # profile's --minimum-ani already keeps these >= 95; a query table does not
        return 1 if float(row[ci['Adjusted_ANI']]) >= 95 else 4
    cov = float(row[ci['Eff_cov']])
    try:
        lo = float(row[ci['ANI_5-95_percentile']].split('-')[0])
    except ValueError:
        lo = float('nan')
    if cov < 0.1 or not lo >= 95:
        return 4
    return 2 if cov >= 1 else 3

def num(x):
    try:
        return float(x)
    except ValueError:
        return float('-inf')

# Sample_file is whatever was profiled - a .sylsp sketch here, but a fastq in a
# table handed over from another run, and read suffixes vary too much to strip
# by guessing (_R1.fq.gz, .fq.1.gz, _1.fastq.gz, ...). Match against the sample
# names we already know instead, longest first so a prefix cannot shadow a
# longer name.
by_len = sorted(samples, key=len, reverse=True)
def sample_of(path):
    b = os.path.basename(path)
    for e in ('.paired.sylsp', '.sylsp'):
        if b.endswith(e):
            b = b[:-len(e)]
            break
    if b in samples:
        return b
    for s in by_len:
        if b.startswith(s):
            return s
    return b

# dropped -> kept, so a sample that hit a collapsed genome points at the
# representative instead of vanishing when STAGE 8 joins against unified_db.tsv
rep = {}
if dedup_map and os.path.isfile(dedup_map):
    with open(dedup_map) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            rep[r['dropped']] = r['kept']

# Where each genome's FASTA actually lives, taken from the table itself: a
# representative can sit in a different chunk directory from the genome it
# replaces, so its path cannot be synthesised from the dropped one's.
path_of = {}
with open(src) as f:
    rdr = csv.reader(f, delimiter='\t')
    h = next(rdr)
    ig = h.index('Genome_file')
    for r in rdr:
        if r:
            path_of.setdefault(gid(r[ig]), r[ig])

# The project's own genomes (STAGE 0b): which sample carries which. Their rows
# never come from the table being split; they are added below from STAGE 0b's
# profile. A GlobDB genome of the same species as an isolate the sample carries
# is a stand-in for it and is dropped there (STAGE 3b).
known, present = set(), set()
if os.path.isfile(known_presence):
    with open(known_presence) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            known.add(r['ref_id'])
            if r['present'] == 'yes':
                present.add((r['sample'], r['ref_id']))
same_sp = {}
if os.path.isfile(same_sp_tsv):
    with open(same_sp_tsv) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            same_sp[r['globdb']] = r['isolate']

with open(src) as f:
    rdr = csv.reader(f, delimiter='\t')
    hdr = next(rdr)
    i_s, i_g = hdr.index('Sample_file'), hdr.index('Genome_file')
    i_c, i_a = hdr.index('Eff_cov'), hdr.index('Taxonomic_abundance')
    ci = {n: hdr.index(n) for n in ['Eff_cov', 'Eff_lambda', 'ANI_5-95_percentile', 'Adjusted_ANI']}
    rows = {}
    n_in = n_lowcov = n_excl = n_shadow = n_remapped = n_merged = 0
    for r in rdr:
        if not r:
            continue
        n_in += 1
        try:
            sel = (tier(r, ci) <= int(tier_max[1:])) if tier_max else (float(r[i_c]) >= mincov)
        except (ValueError, IndexError):
            continue
        if not sel:
            n_lowcov += 1
            continue
        s = sample_of(r[i_s])
        g = gid(r[i_g])
        if g in exclude or g in known:
            n_excl += 1
            continue
        if g in same_sp and (s, same_sp[g]) in present:
            n_shadow += 1
            continue
        if g in rep:                       # rewrite to the kept representative
            r = list(r)
            g = rep[g]
            r[i_g] = path_of.get(g, r[i_g])
            n_remapped += 1
        per = rows.setdefault(s, {})
        if g in per:                       # both members of one cluster in this
            n_merged += 1                  # sample: keep the better-covered row
            old = per[g]                   # and carry the abundance across
            keep, drop = (r, old) if float(r[i_c]) > float(old[i_c]) else (old, r)
            keep = list(keep)
            keep[i_a] = f"{float(keep[i_a]) + float(drop[i_a]):.4f}"
            per[g] = keep
        else:
            per[g] = list(r)

# The project's own genomes, in every sample STAGE 0b called them present. They
# were profiled without the background, so their Taxonomic_abundance means
# nothing next to the GlobDB rows' and is written as NA.
n_known = 0
if present and os.path.isfile(known_profile):
    with open(known_profile) as f:
        frdr = csv.reader(f, delimiter='\t')
        fh = next(frdr)
        for fr in frdr:
            if not fr:
                continue
            s, g = sample_of(fr[fh.index('Sample_file')]), gid(fr[fh.index('Genome_file')])
            if (s, g) not in present:
                continue
            row = [fr[fh.index(n)] if n in fh else 'NA' for n in hdr]
            row[i_a] = 'NA'
            rows.setdefault(s, {})[g] = row
            n_known += 1

written = 0
for s in samples:
    with open(os.path.join(outdir, f'{s}.tsv'), 'w') as o:
        o.write('\t'.join(hdr) + '\n')
        for r in sorted(rows.get(s, {}).values(),
                        key=lambda x: -num(x[i_a])):
            o.write('\t'.join(r) + '\n')
    written += 1

sel_desc = f"tier > {tier_max}" if tier_max else f"below Eff_cov {mincov}"
sys.stderr.write(
    f"  {n_in} rows in, {n_lowcov} not selected ({sel_desc}), {n_excl} excluded, "
    f"{n_shadow} given way to an isolate the sample carries, "
    f"{n_remapped} remapped to a dedup representative, {n_merged} merged\n")
if known:
    sys.stderr.write(f"  {n_known} sample-genome pairs added from the project's own genomes\n")
sys.stderr.write(f"  wrote {written} per-sample reference lists\n")
found = set(rows)
missing = [s for s in samples if s not in found]
if missing:
    sys.stderr.write(f"  WARNING: no rows for {len(missing)} sample(s): {', '.join(missing)}\n")
extra = sorted(found - set(samples))
if extra:
    sys.stderr.write(f"  WARNING: {len(extra)} sample(s) in the profile are not in the metaG dir: {', '.join(extra)}\n")
PYREF

# ============================================================================
# STAGE 8: Per-sample sub-DB from second-round sylph hits
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 8: Per-sample sub-DB construction..."

build_one_sample() {
    local SAMPLE=$1
    local SAMPLE_DIR="${OUT_ROOT}/samples/${SAMPLE}"
    mkdir -p "$SAMPLE_DIR"
    cd "$SAMPLE_DIR"

    local SYLPH_R2="${OUT_ROOT}/sample_references/${SAMPLE}.tsv"
    [[ -s "$SYLPH_R2" ]] || { echo "  [$SAMPLE] no reference list, skip"; return; }

    {
        echo "[$(date '+%F %T')] === $SAMPLE ==="

        python3 - "$SAMPLE" "$OUT_ROOT" > references.tsv <<'PY'
import sys, os, csv

sample, root = sys.argv[1], sys.argv[2]

# Load unified db lookup: ref_id -> {source, folder, path}
unified = {}
with open(os.path.join(root, 'unified_db.tsv')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        unified[r['ref_id']] = r

# Parse second-round sylph result
sylph_path = os.path.join(root, 'sample_references', f'{sample}.tsv')
hits = []
with open(sylph_path) as f:
    rdr = csv.reader(f, delimiter='\t')
    hdr = next(rdr)
    ci = {n: hdr.index(n) for n in ['Genome_file','Taxonomic_abundance','Adjusted_ANI','Eff_cov']}
    # STAGE 7 already applied the selection (Eff_cov or tier); filtering again
    # on Eff_cov here would undo a tier choice.
    n_unmatched = 0
    for row in rdr:
        if len(row) < len(hdr):
            continue
        # Genome_file gives path; recover ref_id from basename
        gp = row[ci['Genome_file']]
        gid = os.path.basename(gp)
        for ext in ('.fa.gz','.fna.gz','.fasta.gz','.fa','.fna','.fasta','.gz'):
            if gid.endswith(ext): gid = gid[:-len(ext)]
        if gid in unified:
            ta = row[ci['Taxonomic_abundance']]
            hits.append({
                'ref_id': gid,
                'tax_abund': float(ta) if ta != 'NA' else None,
                'ani': float(row[ci['Adjusted_ANI']]),
                'eff_cov': float(row[ci['Eff_cov']]),
            })
        else:
            n_unmatched += 1
if n_unmatched:
    sys.stderr.write(f"  WARNING: {n_unmatched} selected genome(s) not in unified_db.tsv, left out\n")

# Output references for this sample
print("ref_id\tsource\tfolder\tpath\ttax_abund\tani\teff_cov")
for h in sorted(hits, key=lambda x: -(x['tax_abund'] if x['tax_abund'] is not None else -1)):
    u = unified[h['ref_id']]
    ta = f"{h['tax_abund']:.4f}" if h['tax_abund'] is not None else 'NA'
    print(f"{h['ref_id']}\t{u['source']}\t{u['folder']}\t{u['path']}\t{ta}\t{h['ani']:.2f}\t{h['eff_cov']:.3f}")
PY

        local N=$(tail -n +2 references.tsv | wc -l)
        local N_MAG=$(awk -F'\t' 'NR>1 && $2=="MAG"' references.tsv | wc -l)
        local N_GDB=$(awk -F'\t' 'NR>1 && $2=="GlobDB"' references.tsv | wc -l)
        local N_FUN=$(awk -F'\t' 'NR>1 && $2=="Fungus"' references.tsv | wc -l)
        local N_ISO=$(awk -F'\t' 'NR>1 && $2=="Isolate"' references.tsv | wc -l)
        echo "  Refs: $N (MAG: $N_MAG, Isolate: $N_ISO, GlobDB: $N_GDB, Fungus: $N_FUN)"

        # Concatenate FASTA + rename contigs
        : > ref.fna
        echo -e "new_contig\tref_id\tsource\toriginal_contig" > ref.contig_map.tsv

        while IFS=$'\t' read -r rid src folder rpath _; do
            [[ -f "$rpath" ]] || { echo "  MISSING: $rpath" >&2; continue; }
            local CAT="cat"; [[ "$rpath" == *.gz ]] && CAT="zcat"
            $CAT "$rpath" | awk -v rid="$rid" -v src="$src" -v mapfile="ref.contig_map.tsv" '
                /^>/ { orig=substr($1,2); new=rid"__"orig;
                       print ">"new; print new"\t"rid"\t"src"\t"orig >> mapfile; next }
                { print }
            ' >> ref.fna
        done < <(tail -n +2 references.tsv)

        local nc=$(grep -c "^>" ref.fna || true)
        local bp=$(awk '!/^>/ {n+=length($0)} END{print n}' ref.fna)
        echo "  ref.fna: $nc contigs, $(awk -v b=$bp 'BEGIN{printf "%.1f", b/1e6}') MB"

        bowtie2-build --threads "$THREADS" ref.fna ref > bt2_build.log 2>&1
        samtools faidx ref.fna 2>/dev/null

        echo "[$(date '+%F %T')] DONE: $SAMPLE"
    } | tee build.log
}

export -f build_one_sample
export OUT_ROOT THREADS MIN_EFF_COV MIN_EFF_COV_R2

# Subshell: build_one_sample cd's into the sample directory, and STAGE 9's
# relative paths (and its summary.tsv) would otherwise land in the last one.
for s in "${SAMPLES[@]}"; do
    ( build_one_sample "$s" )
done

# ============================================================================
# STAGE 9: Summary
# ============================================================================
echo
echo "[$(date '+%F %T')] STAGE 9: Summary"
{
    echo -e "sample\tn_refs\tn_MAG\tn_Isolate\tn_GlobDB\tn_Fungus\tref_Mb"
    for s in "${SAMPLES[@]}"; do
        d="samples/${s}"
        [[ -d "$d" ]] || continue
        n=$(tail -n +2 "$d/references.tsv" 2>/dev/null | wc -l || echo 0)
        nm=$(awk -F'\t' 'NR>1 && $2=="MAG"' "$d/references.tsv" 2>/dev/null | wc -l || echo 0)
        ng=$(awk -F'\t' 'NR>1 && $2=="GlobDB"' "$d/references.tsv" 2>/dev/null | wc -l || echo 0)
        nf=$(awk -F'\t' 'NR>1 && $2=="Fungus"' "$d/references.tsv" 2>/dev/null | wc -l || echo 0)
        ni=$(awk -F'\t' 'NR>1 && $2=="Isolate"' "$d/references.tsv" 2>/dev/null | wc -l || echo 0)
        if [[ -f "$d/ref.fna" ]]; then
            bp=$(awk '!/^>/ {n+=length($0)} END{printf "%.1f", n/1e6}' "$d/ref.fna")
        else
            bp=0
        fi
        echo -e "${s}\t${n}\t${nm}\t${ni}\t${ng}\t${nf}\t${bp}"
    done
} | column -t -s $'\t' | tee summary.tsv

echo
echo "[$(date '+%F %T')] PIPELINE COMPLETE."