#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <annotated-tag> [output-parent-directory]"
  exit 1
fi

TAG="$1"
ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

if ! git rev-parse --verify "$TAG^{tag}" >/dev/null 2>&1; then
  echo "ERROR: '$TAG' is not an annotated tag."
  echo "Create one with: git tag -a '$TAG' -m 'Milestone description'"
  exit 2
fi

STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
COMMIT="$(git rev-list -n 1 "$TAG")"
SHORT="$(git rev-parse --short=12 "$COMMIT")"
OUT_PARENT="${2:-${ROOT}/../viral_gt_snapshots}"
SNAPSHOT_ID="snapshot-${STAMP}-${SHORT}"
OUT="${OUT_PARENT}/${SNAPSHOT_ID}"

mkdir -p "$OUT"

git status --short > "$OUT/git_status.txt"
git log -1 --format=fuller "$COMMIT" > "$OUT/git_log.txt"

# Checksums for files tracked at the tagged commit.
TMP_TREE="$(mktemp -d)"
trap 'rm -rf "$TMP_TREE"' EXIT
git archive "$TAG" | tar -x -C "$TMP_TREE"
(
  cd "$TMP_TREE"
  find . -type f -print0 | LC_ALL=C sort -z | xargs -0 sha256sum
) > "$OUT/tracked_files.sha256"

git bundle create "$OUT/repository.bundle" --all
git archive --format=tar.gz --prefix="${TAG}/" -o "$OUT/repository_${TAG}.tar.gz" "$TAG"
git archive --format=zip --prefix="${TAG}/" -o "$OUT/repository_${TAG}.zip" "$TAG"

{
  printf "snapshot_id	%s
" "$SNAPSHOT_ID"
  printf "created_utc	%s
" "$STAMP"
  printf "created_by	%s
" "$(git config user.name)"
  printf "repository_name	%s
" "$(basename "$ROOT")"
  printf "commit	%s
" "$COMMIT"
  printf "annotated_tag	%s
" "$TAG"
  printf "working_tree_status	%s
" "$(test -s "$OUT/git_status.txt" && echo dirty || echo clean)"
  printf "institutional_storage_location	%s
" "$OUT_PARENT"
  printf "copy_verified_date	
"
  printf "verified_by	
"
  printf "related_manuscript_version	
"
  printf "notes	
"
} > "$OUT/SNAPSHOT_METADATA.tsv"

cat > "$OUT/external_datasets.tsv" <<'EOF'
dataset_id	data_tier	storage_location	file_manifest	manifest_sha256	size_bytes	access_classification	responsible_steward	related_notebook	related_claim
EOF

(
  cd "$OUT"
  sha256sum     repository.bundle     "repository_${TAG}.tar.gz"     "repository_${TAG}.zip"     tracked_files.sha256     SNAPSHOT_METADATA.tsv     external_datasets.tsv > SHA256SUMS
  {
    echo "Created UTC: $STAMP"
    echo "Commit: $COMMIT"
    echo "Tag: $TAG"
    echo
    sha256sum -c SHA256SUMS
    echo
    git bundle verify repository.bundle
  } > verification.txt 2>&1
)

echo "Snapshot created: $OUT"
echo "Review external_datasets.tsv, then copy the entire snapshot to UGA-managed storage."
