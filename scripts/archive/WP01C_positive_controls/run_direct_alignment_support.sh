#!/usr/bin/env bash
set -euo pipefail

# Direct structural support for corrected viral-GT spike controls.
#
# Usage:
#   ./run_direct_alignment_support.sh CONTROL_STRUCTURE_DIR [BASELINE_DB] \
#       [BASELINE_MEMBERSHIP_TSV] [OUTPUT_DIR]

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 CONTROL_STRUCTURE_DIR [BASELINE_DB] [BASELINE_MEMBERSHIP_TSV] [OUTPUT_DIR]" >&2
    exit 2
fi

CONTROL_STRUCTURE_DIR="$1"
BASELINE_DB="${2:-/home/ln72030/foldseek_virus_cluster_clean/database/protein_db}"
BASELINE_MEMBERSHIP="${3:-/home/ln72030/foldseek_virus_cluster_clean/accounting/cluster_membership.tsv}"
OUTPUT_DIR="${4:-/home/ln72030/foldseek_spike_runs/POS001_result/direct_alignment_support}"
THREADS="${THREADS:-64}"
GIANT_REP="${GIANT_REP:-A0A3Q8Q3U2}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTROL_IDS="$SCRIPT_DIR/control_ids.txt"
RAW_HITS="$OUTPUT_DIR/control_vs_baseline_all_hits.tsv"
TMP_DIR="$OUTPUT_DIR/tmp"

command -v foldseek >/dev/null 2>&1 || {
    echo "ERROR: foldseek is not available in PATH." >&2
    exit 1
}

[[ -d "$CONTROL_STRUCTURE_DIR" ]] || {
    echo "ERROR: control structure directory not found: $CONTROL_STRUCTURE_DIR" >&2
    exit 1
}
[[ -f "${BASELINE_DB}.dbtype" ]] || {
    echo "ERROR: baseline Foldseek DB not found: ${BASELINE_DB}.dbtype" >&2
    exit 1
}
[[ -f "$BASELINE_MEMBERSHIP" ]] || {
    echo "ERROR: baseline membership TSV not found: $BASELINE_MEMBERSHIP" >&2
    exit 1
}

mkdir -p "$OUTPUT_DIR" "$TMP_DIR"
foldseek version > "$OUTPUT_DIR/foldseek_version.txt"

FORMAT="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov,alntmscore,qtmscore,ttmscore,lddt,prob"

echo "[1/4] Exhaustive direct search: corrected controls versus baseline DB"
foldseek easy-search \
    "$CONTROL_STRUCTURE_DIR" \
    "$BASELINE_DB" \
    "$RAW_HITS" \
    "$TMP_DIR" \
    --exhaustive-search 1 \
    --alignment-type 2 \
    --max-seqs 100000 \
    -e 1000000 \
    -c 0.0 \
    --cov-mode 0 \
    --format-output "$FORMAT" \
    --threads "$THREADS"

echo "[2/4] Selecting representative and best non-representative giant-member hits"
python "$SCRIPT_DIR/summarize_direct_alignment_support.py" \
    --all-hits "$RAW_HITS" \
    --baseline-membership "$BASELINE_MEMBERSHIP" \
    --control-ids "$CONTROL_IDS" \
    --representative "$GIANT_REP" \
    --output-dir "$OUTPUT_DIR"

echo "[3/4] Rendering the primary alignment-coverage schematic"
MPLCONFIGDIR="$OUTPUT_DIR/mplconfig" \
python "$SCRIPT_DIR/make_alignment_coverage_figure.py" \
    --input "$OUTPUT_DIR/coverage_figure_input.tsv" \
    --output-dir "$OUTPUT_DIR"

echo "[4/4] Packaging results"
tar --no-same-owner -czf "$OUTPUT_DIR/direct_alignment_support_results.tar.gz" \
    -C "$OUTPUT_DIR" \
    control_vs_baseline_all_hits.tsv \
    control_representative_hits.tsv \
    control_best_giant_hits.tsv \
    control_best_nonrep_giant_hits.tsv \
    coverage_figure_input.tsv \
    control_search_qc.tsv \
    direct_alignment_coverage.png \
    direct_alignment_coverage.pdf \
    direct_alignment_coverage.svg \
    direct_alignment_coverage_legend.txt \
    foldseek_version.txt

echo "Completed: $OUTPUT_DIR"
echo "Upload this file for final figure review:"
echo "  $OUTPUT_DIR/direct_alignment_support_results.tar.gz"
