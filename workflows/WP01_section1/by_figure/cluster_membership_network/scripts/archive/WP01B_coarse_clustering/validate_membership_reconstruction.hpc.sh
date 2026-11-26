#!/usr/bin/env bash
set -euo pipefail

BASE=/home/ln72030/foldseek_virus_cluster_clean

RAW="$BASE/cluster/cluster_cluster.tsv"
SOURCE="$BASE/database/protein_db.source"
CANON="$BASE/accounting/cluster_membership.tsv"

REBUILT=$(mktemp /tmp/WP01B_rebuilt_membership.XXXXXX.tsv)
MISMATCH=$(mktemp /tmp/WP01B_membership_mismatch.XXXXXX.tsv)

cleanup() {
    rm -f "$REBUILT" "$MISMATCH"
}
trap cleanup EXIT

for f in "$RAW" "$SOURCE" "$CANON"; do
    if [[ ! -s "$f" ]]; then
        echo -e "FINAL_STATUS\tFAIL"
        echo -e "ERROR\tmissing_or_empty_file\t$f"
        exit 1
    fi
done

echo "WP01B MEMBERSHIP RECONSTRUCTION VALIDATION"
echo -e "validation_date\t$(date --iso-8601=seconds)"
echo

echo "### INPUT FILE IDENTITIES"
for f in "$RAW" "$SOURCE" "$CANON"; do
    stat -c '%n	size=%s	mtime=%y' "$f"
    sha256sum "$f"
done

echo
echo "### SOURCE-MAPPING INTEGRITY"
awk -F '\t' '
{
    rows++

    if (!($1 in numeric_id)) {
        numeric_id[$1] = 1
        unique_numeric_ids++
    } else {
        duplicate_numeric_ids++
    }

    if (!($2 in accession)) {
        accession[$2] = 1
        unique_accessions++
    } else {
        duplicate_accessions++
    }
}
END {
    print "source_rows\t" rows
    print "unique_numeric_ids\t" unique_numeric_ids
    print "unique_accessions\t" unique_accessions
    print "duplicate_numeric_ids\t" duplicate_numeric_ids + 0
    print "duplicate_accessions\t" duplicate_accessions + 0

    if (rows != 1727 || unique_numeric_ids != 1727 || unique_accessions != 1727 || duplicate_numeric_ids != 0 || duplicate_accessions != 0) exit 1
}
' "$SOURCE"

echo
echo "### REBUILD BIOLOGICAL MEMBERSHIP"
awk -F '\t' -v OFS='\t' '
NR == FNR {
    if (NF != 2) {
        malformed_source++
        next
    }

    id_to_accession[$1] = $2
    next
}
{
    if (NF != 2) {
        malformed_raw++
        next
    }

    if (!($2 in id_to_accession)) {
        print "UNMAPPED_MEMBER_ID", FNR, $1, $2 > "/dev/stderr"
        unmapped++
        next
    }

    print $1, id_to_accession[$2]
}
END {
    if (malformed_source != 0 || malformed_raw != 0 || unmapped != 0) {
        print "malformed_source_rows\t" malformed_source + 0 > "/dev/stderr"
        print "malformed_raw_rows\t" malformed_raw + 0 > "/dev/stderr"
        print "unmapped_member_ids\t" unmapped + 0 > "/dev/stderr"
        exit 1
    }
}
' "$SOURCE" "$RAW" > "$REBUILT"

echo -e "raw_membership_rows\t$(wc -l < "$RAW")"
echo -e "rebuilt_membership_rows\t$(wc -l < "$REBUILT")"
echo -e "canonical_membership_rows\t$(wc -l < "$CANON")"

echo
echo "### REBUILT AND CANONICAL CHECKSUMS"
sha256sum "$REBUILT" "$CANON"

echo
echo "### EXACT FILE COMPARISON"
if cmp -s "$REBUILT" "$CANON"; then
    echo -e "byte_for_byte_comparison\tPASS"
else
    echo -e "byte_for_byte_comparison\tFAIL"
    diff -u "$CANON" "$REBUILT" | sed -n '1,80p'
    exit 1
fi

paste "$CANON" "$REBUILT" |
awk -F '\t' -v OFS='\t' '
($1 != $3 || $2 != $4) {
    print NR, $1, $2, $3, $4
}
' > "$MISMATCH"

echo -e "mismatched_rows\t$(wc -l < "$MISMATCH")"

echo
echo "### CANONICAL MEMBERSHIP INTEGRITY"
awk -F '\t' '
{
    rows++

    if (!($1 in representative)) {
        representative[$1] = 1
        unique_representatives++
    }

    if (!($2 in member)) {
        member[$2] = 1
        unique_members++
    }

    assignment[$2]++

    if ($1 == $2)
        self_assignment[$1]++
}
END {
    for (m in assignment)
        if (assignment[m] > 1)
            multiple_assignments++

    for (r in representative) {
        if (!(r in member))
            representatives_missing_as_members++

        if (self_assignment[r] != 1)
            invalid_self_assignments++
    }

    print "membership_rows\t" rows
    print "unique_representatives\t" unique_representatives
    print "unique_members\t" unique_members
    print "multiple_assignments\t" multiple_assignments + 0
    print "representatives_missing_as_members\t" representatives_missing_as_members + 0
    print "invalid_representative_self_assignments\t" invalid_self_assignments + 0

    if (rows != 1727 || unique_representatives != 761 || unique_members != 1727 || multiple_assignments != 0 || representatives_missing_as_members != 0 || invalid_self_assignments != 0) exit 1
}
' "$CANON"

echo
echo -e "FINAL_STATUS\tPASS"
echo -e "CONCLUSION\traw membership plus protein_db.source reconstructs canonical biological membership exactly"
