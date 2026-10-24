#!/usr/bin/env bash
set -euo pipefail

# Existing paths recovered from the supplied Foldseek scripts.
# Edit these variables or override them in the environment before execution.
PROJECT_ROOT="${PROJECT_ROOT:-/home/ln72030/foldseek_virus_cluster_clean}"
STRUCTURES_DIR="${STRUCTURES_DIR:-/home/ln72030/virus_pdb/virus_initial_hits}"
NETWORK_DIR="${NETWORK_DIR:-$PROJECT_ROOT/complete_cytoscape_network}"
OUTPUT_DIR="${OUTPUT_DIR:-$PROJECT_ROOT/singleton_length_plddt_analysis_v3}"
ORIGINAL_FASTA="${ORIGINAL_FASTA:-}"
THREAD_NOTE="This workflow parses structures serially; Foldseek pairwise recomputation is not required."

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NODES_TABLE="${NODES_TABLE:-$NETWORK_DIR/cytoscape_complete_nodes.tsv}"
NETWORK_ACCOUNTING="${NETWORK_ACCOUNTING:-$NETWORK_DIR/network_accounting.tsv}"
CLUSTER_TSV="${CLUSTER_TSV:-$PROJECT_ROOT/accounting/cluster_membership.tsv}"
FOLDSEEK_FASTA="${FOLDSEEK_FASTA:-$PROJECT_ROOT/accounting/database_entries.fasta}"

[[ -s "$NODES_TABLE" ]] || { echo "ERROR: missing nodes table: $NODES_TABLE" >&2; exit 2; }
[[ -d "$STRUCTURES_DIR" ]] || { echo "ERROR: missing structure directory: $STRUCTURES_DIR" >&2; exit 2; }

ARGS=(
  --nodes-table "$NODES_TABLE"
  --structures-dir "$STRUCTURES_DIR"
  --outdir "$OUTPUT_DIR"
  --bootstrap 500
)

[[ -s "$NETWORK_ACCOUNTING" ]] && ARGS+=(--network-accounting "$NETWORK_ACCOUNTING")
[[ -s "$CLUSTER_TSV" ]] && ARGS+=(--cluster-tsv "$CLUSTER_TSV")
[[ -s "$FOLDSEEK_FASTA" ]] && ARGS+=(--foldseek-fasta "$FOLDSEEK_FASTA")
[[ -n "$ORIGINAL_FASTA" && -s "$ORIGINAL_FASTA" ]] && ARGS+=(--original-fasta "$ORIGINAL_FASTA")

printf 'Project root:       %s\n' "$PROJECT_ROOT"
printf 'Structure directory:%s\n' "$STRUCTURES_DIR"
printf 'Nodes table:        %s\n' "$NODES_TABLE"
printf 'Output directory:   %s\n' "$OUTPUT_DIR"
printf '%s\n' "$THREAD_NOTE"

python "$SCRIPT_DIR/scripts/singleton_length_plddt_existing_project.py" "${ARGS[@]}"
