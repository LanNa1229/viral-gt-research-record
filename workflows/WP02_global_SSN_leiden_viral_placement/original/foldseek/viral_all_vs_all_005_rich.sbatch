#!/bin/bash
#SBATCH --job-name=GTA005_allvall
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=168:00:00
#SBATCH --output=GTA005_allvall_%j.out
#SBATCH --error=GTA005_allvall_%j.err

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  sbatch viral_all_vs_all_005_rich.sbatch STRUCTURE_DIR RESULT_DIR

Example:
  sbatch viral_all_vs_all_005_rich.sbatch \
    /scratch/ln72030/005_GTA_extracted/005_GTA_virus_and_others_curated_20260723_manual/structures \
    /home/ln72030/005_GTA_virus_and_others_curated_20260723_manual/foldseek_run_20260724
EOF
}

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "ERROR: Submit this script with sbatch; do not run it on a login node." >&2
    usage >&2
    exit 2
fi

if [[ $# -ne 2 ]]; then
    usage >&2
    exit 2
fi

STRUCTURE_DIR="$1"
RESULT_DIR="$2"

if [[ ! -d "$STRUCTURE_DIR" ]]; then
    echo "ERROR: Structure directory does not exist: $STRUCTURE_DIR" >&2
    exit 1
fi

STRUCTURE_DIR="$(cd "$STRUCTURE_DIR" && pwd -P)"
mkdir -p "$RESULT_DIR"
RESULT_DIR="$(cd "$RESULT_DIR" && pwd -P)"

if ! command -v foldseek >/dev/null 2>&1; then
    echo "ERROR: foldseek is not available in PATH." >&2
    echo "Load the same Foldseek module or conda environment used for earlier runs, then resubmit." >&2
    exit 1
fi

if ! foldseek convertalis -h 2>&1 | grep "qtmscore" >/dev/null; then
    echo "ERROR: This Foldseek version does not advertise qtmscore/ttmscore output." >&2
    echo "Use a current Foldseek release before running the 005 workflow." >&2
    exit 1
fi

THREADS="${SLURM_CPUS_PER_TASK}"
SCRATCH_BASE="${SCRATCH:-/scratch/${USER}}"
RUN_DIR="${SCRATCH_BASE}/005_GTA_foldseek_${SLURM_JOB_ID}"

if [[ -e "$RUN_DIR" ]]; then
    echo "ERROR: Scratch run directory already exists: $RUN_DIR" >&2
    exit 1
fi
mkdir -p "$RUN_DIR"

DB_NAME="${RUN_DIR}/005_gta_db"
PREF_DB="${RUN_DIR}/005_gta_pref"
ALIGN_DB="${RUN_DIR}/005_gta_aln"
STRUCTURE_FILES="${RUN_DIR}/structure_files.txt"
STRUCTURE_STEMS="${RUN_DIR}/structure_stems.txt"
DB_KEYS="${RUN_DIR}/db_keys.txt"
PAIR_KEYS="${RUN_DIR}/all_pair_keys.tsv"

OUTPUT_M8="${RESULT_DIR}/005_GTA_all_vs_all_TM065_cov050_rich.m8"
OUTPUT_COLUMNS="${RESULT_DIR}/005_GTA_all_vs_all_TM065_cov050_rich.columns.txt"
OUTPUT_SUMMARY="${RESULT_DIR}/005_GTA_all_vs_all_TM065_cov050_summary.txt"

for output in "$OUTPUT_M8" "$OUTPUT_COLUMNS" "$OUTPUT_SUMMARY"; do
    if [[ -e "$output" ]]; then
        echo "ERROR: Refusing to overwrite existing output: $output" >&2
        exit 1
    fi
done

echo "Foldseek executable: $(command -v foldseek)"
echo "Foldseek version: $(foldseek version 2>&1 | head -n 1)"
echo "Structure directory: $STRUCTURE_DIR"
echo "Result directory: $RESULT_DIR"
echo "Scratch run directory: $RUN_DIR"
echo "Threads: $THREADS"

find "$STRUCTURE_DIR" -maxdepth 1 -type f \
    \( -iname '*.pdb' -o -iname '*.cif' \) \
    -printf '%f\n' |
    LC_ALL=C sort > "$STRUCTURE_FILES"

sed -E 's/\.[Pp][Dd][Bb]$//; s/\.[Cc][Ii][Ff]$//' \
    "$STRUCTURE_FILES" > "$STRUCTURE_STEMS"

N_STRUCTURES="$(wc -l < "$STRUCTURE_FILES")"
N_UNIQUE_STEMS="$(LC_ALL=C sort -u "$STRUCTURE_STEMS" | wc -l)"

if [[ "$N_STRUCTURES" -eq 0 ]]; then
    echo "ERROR: No top-level .pdb or .cif files were found." >&2
    exit 1
fi

if [[ "$N_STRUCTURES" -ne "$N_UNIQUE_STEMS" ]]; then
    echo "ERROR: Duplicate structure stems were found." >&2
    LC_ALL=C sort "$STRUCTURE_STEMS" | uniq -d >&2
    exit 1
fi

if [[ "$N_STRUCTURES" -ne 4854 ]]; then
    echo "ERROR: Expected 4,854 structures for dataset 005, found $N_STRUCTURES." >&2
    exit 1
fi

echo "Creating Foldseek database from $N_STRUCTURES structures..."
foldseek createdb "$STRUCTURE_DIR" "$DB_NAME" --threads "$THREADS"

N_DB_ENTRIES="$(wc -l < "${DB_NAME}.lookup")"
if [[ "$N_DB_ENTRIES" -ne "$N_STRUCTURES" ]]; then
    echo "ERROR: Foldseek DB has $N_DB_ENTRIES entries for $N_STRUCTURES files." >&2
    echo "This usually indicates an unexpected chain/file parsing issue." >&2
    exit 1
fi

cut -f1 "${DB_NAME}.lookup" | LC_ALL=C sort -n > "$DB_KEYS"
N_UNIQUE_DB_KEYS="$(LC_ALL=C sort -nu "$DB_KEYS" | wc -l)"
if [[ "$N_UNIQUE_DB_KEYS" -ne "$N_DB_ENTRIES" ]]; then
    echo "ERROR: Duplicate internal Foldseek DB keys were detected." >&2
    exit 1
fi

echo "Generating all ordered pairs..."
awk '
    BEGIN { OFS="\t" }
    NR==FNR { key[++n]=$1; next }
    { for (i=1; i<=n; i++) print $1, key[i] }
' "$DB_KEYS" "$DB_KEYS" > "$PAIR_KEYS"

EXPECTED_PAIRS=$((N_DB_ENTRIES * N_DB_ENTRIES))
OBSERVED_PAIRS="$(wc -l < "$PAIR_KEYS")"
if [[ "$OBSERVED_PAIRS" -ne "$EXPECTED_PAIRS" ]]; then
    echo "ERROR: Expected $EXPECTED_PAIRS pair keys, generated $OBSERVED_PAIRS." >&2
    exit 1
fi

echo "Creating explicit all-vs-all prefilter DB..."
foldseek tsv2db "$PAIR_KEYS" "$PREF_DB" --output-dbtype 7

echo "Running Foldseek TM-align structure alignment..."
foldseek structurealign \
    "$DB_NAME" "$DB_NAME" "$PREF_DB" "$ALIGN_DB" \
    --alignment-type 1 \
        -a 1 \
    -e inf \
    --tmscore-threshold 0.65 \
    --tmscore-threshold-mode 0 \
    -c 0.5 \
    --cov-mode 0 \
    --threads "$THREADS" \
    -v 3

FORMAT_COLUMNS="query,target,fident,alnlen,evalue,bits,alntmscore,qtmscore,ttmscore,qcov,tcov,qlen,tlen,qstart,qend,tstart,tend"

echo "Converting alignments to the 17-column rich m8 format..."
foldseek convertalis \
    "$DB_NAME" "$DB_NAME" "$ALIGN_DB" "$OUTPUT_M8" \
    --format-mode 4 \
    --format-output "$FORMAT_COLUMNS" \
    --threads "$THREADS"

printf '%s\n' "$FORMAT_COLUMNS" | tr ',' '\t' > "$OUTPUT_COLUMNS"



N_ALIGNMENTS="$(awk 'NR > 1 {n++} END {print n+0}' "$OUTPUT_M8")"
BAD_COLUMN_ROWS="$(awk -F '\t' 'NR > 1 && NF != 17 {n++} END {print n+0}' "$OUTPUT_M8")"
SELF_HITS="$(awk -F '\t' 'NR > 1 && $1==$2 {n++} END {print n+0}' "$OUTPUT_M8")"
UNIQUE_QUERIES="$(tail -n +2 "$OUTPUT_M8" | cut -f1 | LC_ALL=C sort -u | wc -l)"
UNIQUE_TARGETS="$(tail -n +2 "$OUTPUT_M8" | cut -f2 | LC_ALL=C sort -u | wc -l)"
if [[ "$N_ALIGNMENTS" -eq 0 ]]; then
    echo "ERROR: Foldseek produced an empty m8 file." >&2
    exit 1
fi

if [[ "$BAD_COLUMN_ROWS" -ne 0 ]]; then
    echo "ERROR: $BAD_COLUMN_ROWS m8 rows do not contain exactly 17 columns." >&2
    exit 1
fi

{
    echo "005 GTA FOLDSEEK ALL-VS-ALL SUMMARY"
    echo
    printf 'SLURM_JOB_ID\t%s\n' "$SLURM_JOB_ID"
    printf 'FOLDSEEK_VERSION\t%s\n' "$(foldseek version 2>&1 | head -n 1)"
    printf 'STRUCTURE_DIR\t%s\n' "$STRUCTURE_DIR"
    printf 'STRUCTURE_FILES\t%s\n' "$N_STRUCTURES"
    printf 'UNIQUE_STRUCTURE_STEMS\t%s\n' "$N_UNIQUE_STEMS"
    printf 'FOLDSEEK_DB_ENTRIES\t%s\n' "$N_DB_ENTRIES"
    printf 'ORDERED_PAIR_KEYS\t%s\n' "$OBSERVED_PAIRS"
    printf 'ALIGNMENT_ROWS\t%s\n' "$N_ALIGNMENTS"
    printf 'SELF_HITS\t%s\n' "$SELF_HITS"
    printf 'UNIQUE_QUERIES_IN_OUTPUT\t%s\n' "$UNIQUE_QUERIES"
    printf 'UNIQUE_TARGETS_IN_OUTPUT\t%s\n' "$UNIQUE_TARGETS"
    printf 'BAD_COLUMN_ROWS\t%s\n' "$BAD_COLUMN_ROWS"
    printf 'ALIGNMENT_TYPE\tTMalign_global\n'
    printf 'TM_THRESHOLD\t0.65\n'
    printf 'TM_THRESHOLD_MODE\t0\n'
    printf 'COVERAGE_THRESHOLD\t0.5\n'
    printf 'COVERAGE_MODE\t0\n'
    printf 'EVALUE_PREFILTER\tinf\n'
    printf 'OUTPUT_COLUMNS\t%s\n' "$FORMAT_COLUMNS"
    printf 'OUTPUT_M8\t%s\n' "$OUTPUT_M8"
    printf 'SCRATCH_RUN_DIR\t%s\n' "$RUN_DIR"
    printf 'RESULT\tPASS\n'
} | tee "$OUTPUT_SUMMARY"

echo
echo "First five alignment rows:"
head -n 5 "$OUTPUT_M8"
echo
echo "Done."
echo "Keep the scratch run directory until the m8 file has passed downstream validation:"
echo "  $RUN_DIR"
