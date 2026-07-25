#!/usr/bin/env bash
set -euo pipefail

REPO_NAME="${1:-viral-gt-research-record}"
BASE_DIR="${2:-$HOME/projects}"
REPO_DIR="${BASE_DIR}/${REPO_NAME}"

mkdir -p "${REPO_DIR}"
cd "${REPO_DIR}"

if [ ! -d .git ]; then
  git init -b main
fi

mkdir -p   notebooks/{template,active,completed,retrospective}   src workflows configs   metadata/{schemas,curated,manifests}   data   results/{tables,networks,structures,logs}   figures/{source,panels,final,figure_manifests}   manuscript_claims reports environment archive_manifests   retrospective_records transfer_manifests docs scripts

touch   notebooks/active/.gitkeep   notebooks/completed/.gitkeep   notebooks/retrospective/.gitkeep   src/.gitkeep workflows/.gitkeep configs/.gitkeep   metadata/schemas/.gitkeep metadata/curated/.gitkeep metadata/manifests/.gitkeep   results/tables/.gitkeep results/networks/.gitkeep results/structures/.gitkeep results/logs/.gitkeep   figures/source/.gitkeep figures/panels/.gitkeep figures/final/.gitkeep figures/figure_manifests/.gitkeep   reports/.gitkeep transfer_manifests/.gitkeep

git config pull.ff only

echo "Repository scaffold created at: ${REPO_DIR}"
echo "Copy the starter template files into this directory, review them, then run:"
echo "  git add ."
echo '  git commit -m "Initialize reproducible Viral GT research record"'
echo ""
echo "Create a private remote through an authorized UGA-controlled organization or approved private service."
echo "Do not use a personal account as the sole institutional record."
