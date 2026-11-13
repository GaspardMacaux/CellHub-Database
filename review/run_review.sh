#!/usr/bin/env bash
# Run the CellHubDB revision scripts (NARGAB-2026-310) from the repository root, one log per step
# Usage: bash review/run_review.sh            -> all steps in order
#        bash review/run_review.sh 01 02 03   -> only the listed steps
set -uo pipefail
cd "$(dirname "$0")/.."
mkdir -p review/logs
steps=("$@")
[ ${#steps[@]} -eq 0 ] && steps=(00 01 02 03 05 06 04 07)
for s in "${steps[@]}"; do
  script=$(ls review/${s}_*.R 2>/dev/null | head -n 1)
  if [ -z "$script" ]; then echo "No script for step $s"; exit 1; fi
  echo "=== Step $s: $script ($(date '+%H:%M:%S'))"
  if ! Rscript "$script" > "review/logs/${s}.log" 2>&1; then
    echo "Step $s failed, see review/logs/${s}.log"; tail -n 20 "review/logs/${s}.log"; exit 1
  fi
  echo "    done, log in review/logs/${s}.log"
done
