#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# C7 gated structural phylogeny
# 121 local C7-gated structures
# + 15 MGnify top hits
# + 8 LPS/O-antigen anchors
# - 1 duplicate Q5JBG6
# = 143 unique structures
# ------------------------------------------------------------

SOURCE_DIR="/home/REDACTED/C7_gated_constraint/C7_121_MGnify_LPS_merged"
SOURCE_STRUCTURES="${SOURCE_DIR}/structures"

FOLDTREE_BASE="/mnt/d/FoldTree"
RUN_NAME="C7_gated143_foldtree_run"
RUN_DIR="${FOLDTREE_BASE}/${RUN_NAME}"

CORES="${CORES:-16}"
RESUME="${RESUME:-0}"

echo "SOURCE_STRUCTURES=${SOURCE_STRUCTURES}"
echo "RUN_DIR=${RUN_DIR}"
echo "CORES=${CORES}"

# 1. Validate source
[[ -d "${SOURCE_STRUCTURES}" ]] || {
    echo "ERROR: source structure directory not found:"
    echo "  ${SOURCE_STRUCTURES}"
    exit 1
}

SOURCE_COUNT=$(
    find "${SOURCE_STRUCTURES}" -maxdepth 1 -type f \
        \( -iname '*.pdb' -o -iname '*.cif' -o -iname '*.mmcif' \
           -o -iname '*.pdb.gz' -o -iname '*.cif.gz' -o -iname '*.mmcif.gz' \) \
        | wc -l
)

if [[ "${SOURCE_COUNT}" -ne 143 ]]; then
    echo "ERROR: expected 143 source structures, observed ${SOURCE_COUNT}."
    echo "Inspect ${SOURCE_DIR}/merge_summary.json before running FoldTree."
    exit 1
fi

# 2. Find the same FoldTree workflow installation used under /mnt/d/FoldTree
if [[ -n "${FOLDTREE_WORKFLOW:-}" ]]; then
    WORKFLOW="${FOLDTREE_WORKFLOW}"
else
    mapfile -t WORKFLOW_CANDIDATES < <(
        find "${FOLDTREE_BASE}" -maxdepth 5 -type f \
            -path '*/workflow/fold_tree' 2>/dev/null | sort
    )

    if [[ "${#WORKFLOW_CANDIDATES[@]}" -eq 0 ]]; then
        echo "ERROR: could not find workflow/fold_tree under ${FOLDTREE_BASE}."
        echo "Set it explicitly, for example:"
        echo "  export FOLDTREE_WORKFLOW=/mnt/d/FoldTree/fold_tree/workflow/fold_tree"
        exit 1
    elif [[ "${#WORKFLOW_CANDIDATES[@]}" -gt 1 ]]; then
        echo "ERROR: multiple FoldTree workflows were found:"
        printf '  %s\n' "${WORKFLOW_CANDIDATES[@]}"
        echo "Set FOLDTREE_WORKFLOW explicitly to the Tree0 installation."
        exit 1
    fi

    WORKFLOW="${WORKFLOW_CANDIDATES[0]}"
fi

[[ -f "${WORKFLOW}" ]] || {
    echo "ERROR: FoldTree workflow does not exist: ${WORKFLOW}"
    exit 1
}

command -v snakemake >/dev/null 2>&1 || {
    echo "ERROR: snakemake is not available in the active environment."
    echo "Activate the same conda environment used for Tree0."
    exit 1
}

echo "FOLDTREE_WORKFLOW=${WORKFLOW}"
echo "SNAKEMAKE_VERSION=$(snakemake --version)"

# 3. Create a clean project folder in the same /mnt/d/FoldTree style
if [[ -e "${RUN_DIR}" && "${RESUME}" != "1" ]]; then
    echo "ERROR: run directory already exists:"
    echo "  ${RUN_DIR}"
    echo "To resume an existing run:"
    echo "  RESUME=1 bash $0"
    exit 1
fi

mkdir -p \
    "${RUN_DIR}/structures" \
    "${RUN_DIR}/00_metadata" \
    "${RUN_DIR}/99_logs"

# For a new run, synchronize exactly the 143 structures.
if [[ "${RESUME}" != "1" ]]; then
    rsync -a --delete \
        "${SOURCE_STRUCTURES}/" \
        "${RUN_DIR}/structures/"
fi

# FoldTree custom-structure mode requires a blank identifier file.
: > "${RUN_DIR}/identifiers.txt"

# Preserve provenance and exact input metadata.
for metadata_file in \
    merge_summary.json \
    foldtree_input_manifest.tsv \
    foldtree_input_ids.txt \
    foldtree_input_paths.txt
do
    if [[ -f "${SOURCE_DIR}/${metadata_file}" ]]; then
        cp -f \
            "${SOURCE_DIR}/${metadata_file}" \
            "${RUN_DIR}/00_metadata/${metadata_file}"
    fi
done

RUN_INPUT_COUNT=$(
    find "${RUN_DIR}/structures" -maxdepth 1 -type f \
        \( -iname '*.pdb' -o -iname '*.cif' -o -iname '*.mmcif' \
           -o -iname '*.pdb.gz' -o -iname '*.cif.gz' -o -iname '*.mmcif.gz' \) \
        | wc -l
)

if [[ "${RUN_INPUT_COUNT}" -ne 143 ]]; then
    echo "ERROR: FoldTree run directory contains ${RUN_INPUT_COUNT}, expected 143."
    exit 1
fi

find "${RUN_DIR}/structures" -maxdepth 1 -type f \
    -printf '%f\n' | sort \
    > "${RUN_DIR}/00_metadata/input_structure_filenames.txt"

sha256sum "${RUN_DIR}"/structures/* \
    > "${RUN_DIR}/00_metadata/input_structures.sha256"

cat > "${RUN_DIR}/00_metadata/run_design.txt" <<EOF
Analysis label: C7 gated 143-member FoldTree
Run name: ${RUN_NAME}
Input:
  121 C7-gated network structures
    - 23 Cluster 7 focal candidates
    - 98 neighbor/context proteins
  15 MGnify top C7-like candidates
  8 LPS/O-antigen functional anchors
  1 known overlap: Q5JBG6
Expected unique structures: 143
Source: ${SOURCE_STRUCTURES}
Workflow: ${WORKFLOW}
Cores: ${CORES}
customstructs=True
filter=False
EOF

# 4. Dry-run first
echo
echo "===== FoldTree dry run ====="
snakemake \
    --dry-run \
    --cores "${CORES}" \
    --use-conda \
    -s "${WORKFLOW}" \
    --config \
        folder="${RUN_DIR}" \
        filter=False \
        customstructs=True

# 5. Full FoldTree run
echo
echo "===== FoldTree full run ====="
set -o pipefail
snakemake \
    --cores "${CORES}" \
    --use-conda \
    --rerun-incomplete \
    --latency-wait 60 \
    --printshellcmds \
    -s "${WORKFLOW}" \
    --config \
        folder="${RUN_DIR}" \
        filter=False \
        customstructs=True \
    2>&1 | tee "${RUN_DIR}/99_logs/foldtree_run.log"

echo
echo "FoldTree completed."
echo "Run directory:"
echo "  ${RUN_DIR}"
echo
echo "Final-tree candidates:"
find "${RUN_DIR}" -type f \
    \( -name '*.nwk' -o -name '*.treefile' -o -name '*.contree' \) \
    | sort \
    | tee "${RUN_DIR}/00_metadata/tree_output_paths.txt"
