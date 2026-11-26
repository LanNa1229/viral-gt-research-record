#!/usr/bin/env bash
set -euo pipefail

# Robust small-cluster pairwise analysis for Foldseek.
#
# This version deliberately avoids createsubdb. It runs exhaustive all-vs-all
# alignment on the original Foldseek database (1,727 entries in this project),
# then filters the output to the selected 63 clusters / 185 proteins.
#
# Usage:
#   bash run_small_cluster_pairwise_v3_full_db.sh \
#     cluster_membership.tsv \
#     absent_reps.txt \
#     protein_db \
#     output_dir \
#     [threads]

MEMBERSHIP=${1:?cluster_membership.tsv required}
REPS=${2:?representative list required}
DB=${3:?Foldseek database prefix required}
OUT=${4:?output directory required}
THREADS=${5:-64}

mkdir -p "$OUT/work/tmp"

SELECTED_IDS="$OUT/selected_ids.txt"
SELECTED_MEMBERSHIP="$OUT/selected_cluster_membership.tsv"
FULL_FASTA="$OUT/full_database_structural_sequences.fasta"
SELECTED_FASTA="$OUT/selected_structural_sequences.fasta"

ALNDB="$OUT/work/full_all_vs_all_aln"
FULL_PAIRWISE="$OUT/full_all_vs_all_pairwise.tsv"
SELECTED_PAIRWISE="$OUT/pairwise_directed.tsv"
LOG="$OUT/run_small_cluster_pairwise_v3.log"

FORMAT_OUTPUT='query,target,qlen,tlen,alnlen,qcov,tcov,qstart,qend,tstart,tend,fident,evalue,bits,qtmscore,ttmscore,alntmscore,lddt'

: > "$LOG"

echo "============================================================"
echo "1. Selecting requested clusters and proteins"
echo "============================================================"

python - "$MEMBERSHIP" "$REPS" "$SELECTED_IDS" "$SELECTED_MEMBERSHIP" <<'PY'
import sys
from collections import defaultdict

membership_path, reps_path, ids_out, membership_out = sys.argv[1:]

requested = set()
with open(reps_path) as handle:
    for raw in handle:
        line = raw.strip()
        if line and not line.startswith("#"):
            requested.add(line.split()[0])

clusters = defaultdict(list)
member_to_rep = {}

with open(membership_path) as handle:
    for line_no, raw in enumerate(handle, 1):
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 2:
            raise SystemExit(
                f"{membership_path}:{line_no}: expected at least 2 columns"
            )
        rep, member = fields[0], fields[1]
        if member in member_to_rep and member_to_rep[member] != rep:
            raise SystemExit(
                f"{member} is assigned to multiple representatives"
            )
        if member not in member_to_rep:
            member_to_rep[member] = rep
            clusters[rep].append(member)

missing_reps = sorted(requested - set(clusters))
if missing_reps:
    raise SystemExit(
        f"{len(missing_reps)} requested representatives are absent from "
        f"the membership table. Examples: {', '.join(missing_reps[:10])}"
    )

selected_ids = []
with open(membership_out, "w") as out:
    for rep in sorted(requested):
        for member in clusters[rep]:
            out.write(f"{rep}\t{member}\n")
            selected_ids.append(member)

selected_unique = sorted(set(selected_ids))
with open(ids_out, "w") as out:
    for protein_id in selected_unique:
        out.write(protein_id + "\n")

print(f"Selected clusters: {len(requested)}")
print(f"Selected proteins: {len(selected_unique)}")
PY

echo
echo "============================================================"
echo "2. Exporting structure-derived sequences from the full DB"
echo "============================================================"

if [[ ! -s "$FULL_FASTA" ]]; then
    foldseek convert2fasta \
        "$DB" \
        "$FULL_FASTA" \
        2>&1 | tee -a "$LOG"
else
    echo "Reusing existing file: $FULL_FASTA"
fi

python - "$SELECTED_IDS" "$FULL_FASTA" "$SELECTED_FASTA" <<'PY'
import sys

ids_path, fasta_path, output_path = sys.argv[1:]

with open(ids_path) as handle:
    wanted = {line.strip() for line in handle if line.strip()}

found = set()
current_id = None
current_lines = []

def save(out):
    global current_id, current_lines
    if current_id in wanted:
        out.write(f">{current_id}\n")
        out.write("".join(current_lines))
        if current_lines and not current_lines[-1].endswith("\n"):
            out.write("\n")
        found.add(current_id)

with open(fasta_path) as inp, open(output_path, "w") as out:
    for raw in inp:
        if raw.startswith(">"):
            if current_id is not None:
                save(out)
            current_id = raw[1:].strip().split()[0]
            current_lines = []
        else:
            current_lines.append(raw)
    if current_id is not None:
        save(out)

missing = sorted(wanted - found)
if missing:
    raise SystemExit(
        f"{len(missing)} selected IDs were absent from the full FASTA. "
        f"Examples: {', '.join(missing[:10])}"
    )

print(f"Selected FASTA entries: {len(found)}")
PY

DB_ENTRIES=$(grep -c '^>' "$FULL_FASTA")
MAX_SEQS=$((DB_ENTRIES + 100))
if (( MAX_SEQS < 2000 )); then
    MAX_SEQS=2000
fi

echo
echo "Full database entries: $DB_ENTRIES"
echo "Foldseek --max-seqs:   $MAX_SEQS"

echo
echo "============================================================"
echo "3. Exhaustive all-vs-all global structural alignment"
echo "============================================================"

if [[ ! -f "${ALNDB}.dbtype" ]]; then
    rm -rf "$OUT/work/tmp"/*
    foldseek search \
        "$DB" \
        "$DB" \
        "$ALNDB" \
        "$OUT/work/tmp" \
        --exhaustive-search 1 \
        --alignment-type 1 \
        -a \
        -e 1e9 \
        -c 0.0 \
        --max-seqs "$MAX_SEQS" \
        --threads "$THREADS" \
        2>&1 | tee -a "$LOG"
else
    echo "Reusing existing alignment DB: $ALNDB"
fi

echo
echo "============================================================"
echo "4. Converting full all-vs-all results"
echo "============================================================"

if [[ ! -s "$FULL_PAIRWISE" ]]; then
    foldseek convertalis \
        "$DB" \
        "$DB" \
        "$ALNDB" \
        "$FULL_PAIRWISE" \
        --format-output "$FORMAT_OUTPUT" \
        2>&1 | tee -a "$LOG"
else
    echo "Reusing existing file: $FULL_PAIRWISE"
fi

echo
echo "============================================================"
echo "5. Retaining only within-cluster pairs for selected clusters"
echo "============================================================"

python - "$SELECTED_MEMBERSHIP" "$FULL_PAIRWISE" "$SELECTED_PAIRWISE" <<'PY'
import sys

membership_path, pairwise_path, output_path = sys.argv[1:]

member_to_rep = {}
with open(membership_path) as handle:
    for line_no, raw in enumerate(handle, 1):
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 2:
            raise SystemExit(
                f"{membership_path}:{line_no}: expected at least 2 columns"
            )
        rep, member = fields[0], fields[1]
        member_to_rep[member] = rep

kept = 0
self_rows = 0
cross_cluster_rows = 0
unselected_rows = 0

with open(pairwise_path) as inp, open(output_path, "w") as out:
    for line_no, raw in enumerate(inp, 1):
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 2:
            continue
        query, target = fields[0], fields[1]

        if query not in member_to_rep or target not in member_to_rep:
            unselected_rows += 1
            continue

        if member_to_rep[query] != member_to_rep[target]:
            cross_cluster_rows += 1
            continue

        out.write(raw)
        kept += 1
        if query == target:
            self_rows += 1

print(f"Selected within-cluster directed rows: {kept}")
print(f"  including self rows:                {self_rows}")
print(f"Filtered selected cross-cluster rows: {cross_cluster_rows}")
print(f"Filtered unselected rows:             {unselected_rows}")
PY

echo
echo "============================================================"
echo "6. Sanity checks"
echo "============================================================"

python - "$SELECTED_MEMBERSHIP" "$SELECTED_PAIRWISE" <<'PY'
import sys
from collections import defaultdict

membership_path, pairwise_path = sys.argv[1:]

clusters = defaultdict(list)
with open(membership_path) as handle:
    for raw in handle:
        rep, member = raw.rstrip("\n").split("\t")[:2]
        clusters[rep].append(member)

expected_undirected = sum(
    len(members) * (len(members) - 1) // 2
    for members in clusters.values()
)

directed_nonself = 0
undirected = set()
with open(pairwise_path) as handle:
    for raw in handle:
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 2:
            continue
        q, t = fields[0], fields[1]
        if q == t:
            continue
        directed_nonself += 1
        undirected.add(tuple(sorted((q, t))))

print(f"Expected undirected within-cluster pairs: {expected_undirected}")
print(f"Observed undirected within-cluster pairs: {len(undirected)}")
print(f"Observed directed non-self rows:           {directed_nonself}")

if len(undirected) < expected_undirected:
    print(
        "WARNING: Some expected within-cluster pairs were not reported. "
        "The summarization script will list them explicitly."
    )
PY

echo
echo "============================================================"
echo "Completed successfully"
echo "============================================================"
echo "Selected membership: $SELECTED_MEMBERSHIP"
echo "Selected FASTA:      $SELECTED_FASTA"
echo "Selected pairwise:   $SELECTED_PAIRWISE"
echo "Full pairwise file:  $FULL_PAIRWISE"
echo "Run log:             $LOG"
