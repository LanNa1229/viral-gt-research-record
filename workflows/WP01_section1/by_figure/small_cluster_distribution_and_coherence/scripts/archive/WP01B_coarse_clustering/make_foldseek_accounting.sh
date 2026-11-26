#!/bin/bash
set -euo pipefail

# ============================================================
# USER SETTINGS
# ============================================================

INPUT_PDB_DIR="/home/ln72030/foldseek_spike_runs/POS001"
OUTPUT_DIR="/home/ln72030/foldseek_spike_runs/POS001_result"

DB_NAME="$OUTPUT_DIR/database/protein_db"
CLUSTER_DB="$OUTPUT_DIR/cluster/cluster"

ACCOUNT_DIR="$OUTPUT_DIR/accounting"

# ============================================================

mkdir -p "$ACCOUNT_DIR"

CLUSTER_TSV="$ACCOUNT_DIR/cluster_membership.tsv"
DB_FASTA="$ACCOUNT_DIR/database_entries.fasta"

DB_IDS="$ACCOUNT_DIR/database_ids.txt"
CLUSTER_MEMBER_IDS="$ACCOUNT_DIR/cluster_member_ids.txt"
CLUSTER_SIZES="$ACCOUNT_DIR/cluster_sizes.tsv"

SINGLETON_IDS="$ACCOUNT_DIR/singleton_ids.txt"
MULTI_MEMBER_IDS="$ACCOUNT_DIR/multi_member_ids.txt"

MISSING_FROM_CLUSTER="$ACCOUNT_DIR/database_ids_missing_from_cluster.txt"
EXTRA_IN_CLUSTER="$ACCOUNT_DIR/cluster_ids_missing_from_database.txt"
DUPLICATE_MEMBERS="$ACCOUNT_DIR/members_assigned_to_multiple_clusters.txt"

ACCOUNT_TABLE="$ACCOUNT_DIR/accounting_table.tsv"

echo "============================================================"
echo "Foldseek cluster accounting"
echo "============================================================"

echo
echo "[1] Foldseek version"
foldseek version | tee "$ACCOUNT_DIR/foldseek_version.txt"

echo
echo "[2] Exporting cluster representative-member mapping"

foldseek createtsv \
    "$DB_NAME" \
    "$DB_NAME" \
    "$CLUSTER_DB" \
    "$CLUSTER_TSV"

echo
echo "[3] Exporting all sequences stored in Foldseek database"

foldseek convert2fasta \
    "$DB_NAME" \
    "$DB_FASTA"

echo
echo "[4] Extracting database identifiers"

grep '^>' "$DB_FASTA" \
    | sed 's/^>//' \
    | sed 's/[[:space:]].*$//' \
    | sort -u \
    > "$DB_IDS"

echo
echo "[5] Extracting cluster member identifiers"

cut -f2 "$CLUSTER_TSV" \
    | sort -u \
    > "$CLUSTER_MEMBER_IDS"

echo
echo "[6] Calculating cluster sizes"

awk -F'\t' '
{
    size[$1]++
}
END {
    for (rep in size) {
        print rep "\t" size[rep]
    }
}
' "$CLUSTER_TSV" \
    | sort -k2,2nr -k1,1 \
    > "$CLUSTER_SIZES"

echo
echo "[7] Extracting singleton identifiers"

awk -F'\t' '
NR == FNR {
    size[$1] = $2
    next
}
size[$1] == 1 {
    print $2
}
' "$CLUSTER_SIZES" "$CLUSTER_TSV" \
    | sort -u \
    > "$SINGLETON_IDS"

echo
echo "[8] Extracting members of clusters with size >= 2"

awk -F'\t' '
NR == FNR {
    size[$1] = $2
    next
}
size[$1] >= 2 {
    print $2
}
' "$CLUSTER_SIZES" "$CLUSTER_TSV" \
    | sort -u \
    > "$MULTI_MEMBER_IDS"

echo
echo "[9] Checking database IDs absent from cluster assignments"

comm -23 "$DB_IDS" "$CLUSTER_MEMBER_IDS" \
    > "$MISSING_FROM_CLUSTER"

echo
echo "[10] Checking cluster IDs absent from database"

comm -13 "$DB_IDS" "$CLUSTER_MEMBER_IDS" \
    > "$EXTRA_IN_CLUSTER"

echo
echo "[11] Checking members assigned to multiple representatives"

cut -f2 "$CLUSTER_TSV" \
    | sort \
    | uniq -d \
    > "$DUPLICATE_MEMBERS"

# ============================================================
# Counts
# ============================================================

INPUT_FILES=$(
    find "$INPUT_PDB_DIR" -maxdepth 1 -type f \
        \( -iname "*.pdb" \
        -o -iname "*.pdb.gz" \
        -o -iname "*.cif" \
        -o -iname "*.cif.gz" \
        -o -iname "*.mmcif" \
        -o -iname "*.mmcif.gz" \) \
        | wc -l
)

DB_ENTRIES=$(grep -c '^>' "$DB_FASTA" || true)

ASSIGNMENT_ROWS=$(wc -l < "$CLUSTER_TSV")
UNIQUE_CLUSTER_MEMBERS=$(wc -l < "$CLUSTER_MEMBER_IDS")
TOTAL_CLUSTERS=$(wc -l < "$CLUSTER_SIZES")

SINGLETON_CLUSTERS=$(
    awk -F'\t' '$2 == 1 {n++} END {print n+0}' "$CLUSTER_SIZES"
)

SINGLETON_MEMBERS=$(wc -l < "$SINGLETON_IDS")

MULTI_MEMBER_CLUSTERS=$(
    awk -F'\t' '$2 >= 2 {n++} END {print n+0}' "$CLUSTER_SIZES"
)

MULTI_MEMBER_MEMBERS=$(wc -l < "$MULTI_MEMBER_IDS")

MISSING_COUNT=$(wc -l < "$MISSING_FROM_CLUSTER")
EXTRA_COUNT=$(wc -l < "$EXTRA_IN_CLUSTER")
DUPLICATE_COUNT=$(wc -l < "$DUPLICATE_MEMBERS")

SELF_ASSIGNMENTS=$(
    awk -F'\t' '$1 == $2 {n++} END {print n+0}' "$CLUSTER_TSV"
)

PCT_SINGLETON=$(
    awk -v n="$SINGLETON_MEMBERS" -v d="$DB_ENTRIES" '
    BEGIN {
        if (d > 0) printf "%.2f", 100*n/d
        else print "NA"
    }'
)

PCT_MULTI=$(
    awk -v n="$MULTI_MEMBER_MEMBERS" -v d="$DB_ENTRIES" '
    BEGIN {
        if (d > 0) printf "%.2f", 100*n/d
        else print "NA"
    }'
)

# ============================================================
# Accounting table
# ============================================================

{
    printf "metric\tvalue\n"
    printf "input_structure_files\t%s\n" "$INPUT_FILES"
    printf "foldseek_database_entries\t%s\n" "$DB_ENTRIES"
    printf "cluster_assignment_rows\t%s\n" "$ASSIGNMENT_ROWS"
    printf "unique_cluster_members\t%s\n" "$UNIQUE_CLUSTER_MEMBERS"
    printf "total_clusters\t%s\n" "$TOTAL_CLUSTERS"
    printf "singleton_clusters\t%s\n" "$SINGLETON_CLUSTERS"
    printf "singleton_members\t%s\n" "$SINGLETON_MEMBERS"
    printf "singleton_members_percent\t%s\n" "$PCT_SINGLETON"
    printf "multi_member_clusters\t%s\n" "$MULTI_MEMBER_CLUSTERS"
    printf "multi_member_members\t%s\n" "$MULTI_MEMBER_MEMBERS"
    printf "multi_member_members_percent\t%s\n" "$PCT_MULTI"
    printf "representative_self_assignments\t%s\n" "$SELF_ASSIGNMENTS"
    printf "database_ids_missing_from_cluster\t%s\n" "$MISSING_COUNT"
    printf "cluster_ids_missing_from_database\t%s\n" "$EXTRA_COUNT"
    printf "members_assigned_multiple_times\t%s\n" "$DUPLICATE_COUNT"
} > "$ACCOUNT_TABLE"

echo
echo "============================================================"
echo "ACCOUNTING TABLE"
echo "============================================================"

column -t -s $'\t' "$ACCOUNT_TABLE" || cat "$ACCOUNT_TABLE"

echo
echo "Output directory:"
echo "$ACCOUNT_DIR"

echo
echo "Important files:"
echo "  $ACCOUNT_TABLE"
echo "  $CLUSTER_SIZES"
echo "  $SINGLETON_IDS"
echo "  $MULTI_MEMBER_IDS"
echo "  $MISSING_FROM_CLUSTER"

echo
echo "✅ Accounting completed."
