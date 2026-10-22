#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash RUN_R699_DUAL_CONTEXT_MAPPING.sh [RUN_DIR] [R699_STRUCTURE] [MAPPING_SCRIPT]
#
# Defaults:
#   RUN_DIR        current directory
#   R699_STRUCTURE /home/rahnn/260801_Temp/9DYT_A.pdb
#   MAPPING_SCRIPT RUN_DIR/map_dual_context_cores_to_r699.py

RUN_DIR="${1:-$PWD}"
R699="${2:-/home/rahnn/260801_Temp/9DYT_A.pdb}"
MAP_SCRIPT="${3:-$RUN_DIR/map_dual_context_cores_to_r699.py}"

RUN_DIR="$(realpath "$RUN_DIR")"
R699="$(realpath "$R699")"
MAP_SCRIPT="$(realpath "$MAP_SCRIPT")"

C005_DIR="$RUN_DIR/C005"
C009_DIR="$RUN_DIR/C009"
OUT="$RUN_DIR/work/05_query_mapping"
MAP_OUT="$RUN_DIR/work/06_R699_dual_context_mapping"
LOG_DIR="$RUN_DIR/logs"

for required in \
  "$C005_DIR" \
  "$C009_DIR" \
  "$R699" \
  "$MAP_SCRIPT" \
  "$RUN_DIR/work/03_metrics/all_residue_metrics.tsv" \
  "$RUN_DIR/work/01_medoids/medoids.tsv" \
  "$RUN_DIR/work/02_alignments/alignments/within/C005/C005_aa.fa" \
  "$RUN_DIR/work/02_alignments/alignments/within/C009/C009_aa.fa"
do
  [[ -e "$required" ]] || { echo "ERROR: missing $required" >&2; exit 2; }
done

command -v foldmason >/dev/null || { echo "ERROR: foldmason not found in PATH" >&2; exit 2; }

mkdir -p "$OUT/tmp_C005" "$OUT/tmp_C009" "$MAP_OUT" "$LOG_DIR"

mapfile -t C005_INPUTS < <(
  find "$C005_DIR" -maxdepth 1 -type f \
    \( -iname '*.pdb' -o -iname '*.cif' -o -iname '*.mmcif' \) \
    -print | sort
)
mapfile -t C009_INPUTS < <(
  find "$C009_DIR" -maxdepth 1 -type f \
    \( -iname '*.pdb' -o -iname '*.cif' -o -iname '*.mmcif' \) \
    -print | sort
)

[[ ${#C005_INPUTS[@]} -eq 13 ]] || {
  echo "ERROR: expected 13 C005 reference structures; found ${#C005_INPUTS[@]}" >&2
  exit 2
}
[[ ${#C009_INPUTS[@]} -eq 15 ]] || {
  echo "ERROR: expected 15 C009 reference structures; found ${#C009_INPUTS[@]}" >&2
  exit 2
}

printf 'C005 references: %d; with R699: %d\n' "${#C005_INPUTS[@]}" "$(( ${#C005_INPUTS[@]} + 1 ))"
printf 'C009 references: %d; with R699: %d\n' "${#C009_INPUTS[@]}" "$(( ${#C009_INPUTS[@]} + 1 ))"

rm -f "$OUT"/C005_R699* "$OUT"/C009_R699*
rm -rf "$OUT/tmp_C005" "$OUT/tmp_C009"
mkdir -p "$OUT/tmp_C005" "$OUT/tmp_C009"

foldmason easy-msa \
  "${C005_INPUTS[@]}" \
  "$R699" \
  "$OUT/C005_R699" \
  "$OUT/tmp_C005" \
  --report-mode 1 \
  --refine-iters 0 \
  --pair-threshold 0.0 \
  --threads 8 \
  2>&1 | tee "$LOG_DIR/08_C005_R699_foldmason.log"

foldmason easy-msa \
  "${C009_INPUTS[@]}" \
  "$R699" \
  "$OUT/C009_R699" \
  "$OUT/tmp_C009" \
  --report-mode 1 \
  --refine-iters 0 \
  --pair-threshold 0.0 \
  --threads 8 \
  2>&1 | tee "$LOG_DIR/09_C009_R699_foldmason.log"

C005_N=$(grep -c '^>' "$OUT/C005_R699_aa.fa")
C009_N=$(grep -c '^>' "$OUT/C009_R699_aa.fa")
[[ "$C005_N" -eq 14 ]] || { echo "ERROR: C005 query MSA has $C005_N records, expected 14" >&2; exit 2; }
[[ "$C009_N" -eq 16 ]] || { echo "ERROR: C009 query MSA has $C009_N records, expected 16" >&2; exit 2; }

echo "C005+R699 AA MSA records: $C005_N"
echo "C009+R699 AA MSA records: $C009_N"

grep '^>' "$OUT/C005_R699_aa.fa" | grep -i '9DYT' || {
  echo "ERROR: R699 header not found in C005 query MSA" >&2; exit 2;
}
grep '^>' "$OUT/C009_R699_aa.fa" | grep -i '9DYT' || {
  echo "ERROR: R699 header not found in C009 query MSA" >&2; exit 2;
}

python "$MAP_SCRIPT" \
  --metrics "$RUN_DIR/work/03_metrics/all_residue_metrics.tsv" \
  --medoids "$RUN_DIR/work/01_medoids/medoids.tsv" \
  --c005-reference-aa "$RUN_DIR/work/02_alignments/alignments/within/C005/C005_aa.fa" \
  --c005-query-aa "$OUT/C005_R699_aa.fa" \
  --c009-reference-aa "$RUN_DIR/work/02_alignments/alignments/within/C009/C009_aa.fa" \
  --c009-query-aa "$OUT/C009_R699_aa.fa" \
  --query-id 9DYT_A \
  --query-structure "$R699" \
  --query-chain A \
  --core-flag passes_within_quality \
  --min-within-lddt 0.75 \
  --min-within-presence 0.70 \
  --outdir "$MAP_OUT"

echo
echo "Completed. Main outputs:"
printf '  %s\n' \
  "$MAP_OUT/R699_dual_context_classification.tsv" \
  "$MAP_OUT/R699_dual_context_summary.tsv" \
  "$MAP_OUT/R699_dual_context_mapping.pml"
