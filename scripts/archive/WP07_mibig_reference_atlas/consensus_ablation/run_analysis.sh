#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
for stage in construct_consensus compare_knn geometry validate_results write_report; do
  python3 -B "scripts/$stage.py" 2>&1 | tee "logs/$stage.log"
done
