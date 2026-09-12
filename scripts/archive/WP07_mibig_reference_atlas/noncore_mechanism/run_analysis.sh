#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
for stage in reconstruct mechanism compare_consensus plot_results validate write_report; do
  python3 -B "scripts/$stage.py" 2>&1 | tee "logs/$stage.log"
done
