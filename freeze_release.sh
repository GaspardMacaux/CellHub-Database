#!/usr/bin/env bash
# Freeze the last pipeline run as release/<version>/ (never overwrites an existing release).
# Usage: bash freeze_release.sh <version>      e.g. bash freeze_release.sh 2.1
# Content: CellHubDB.rds, summary CSVs, app_databases/ (per-species merged + per-source files for Cell-Hub), reports/, gate/, report/, MD5SUMS
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
VERSION="${1:?usage: bash freeze_release.sh <version>}"
DEST="release/$VERSION"
if [[ -e "$DEST" ]]; then echo "ERROR: $DEST already exists; move it away before freezing again"; exit 1; fi
if [[ ! -f output/CellHubDB.rds ]]; then echo "ERROR: output/CellHubDB.rds not found; run run_pipeline.sh first"; exit 1; fi
built=$(Rscript -e 'cat(attr(readRDS("output/CellHubDB.rds"), "version"))')
if [[ "$built" != "$VERSION" ]]; then echo "ERROR: output/CellHubDB.rds has version '$built', not '$VERSION'"; exit 1; fi
if [[ ! -f "output/gate/v$VERSION/gate_report.csv" ]]; then echo "ERROR: no release gate report for v$VERSION"; exit 1; fi
if grep -q '^"FAIL"' "output/gate/v$VERSION/gate_report.csv"; then echo "ERROR: release gate has FAIL lines"; exit 1; fi
mkdir -p "$DEST"
cp output/CellHubDB.rds output/CellHubDB_species_summary.csv output/CellHubDB_source_contribution.csv "$DEST/"
cp -r data/databases "$DEST/app_databases"
cp -r output/reports "$DEST/reports"
cp -r "output/gate/v$VERSION" "$DEST/gate"
cp -r "output/report/$VERSION" "$DEST/report"
(cd "$DEST" && find . -type f ! -name MD5SUMS -exec md5sum {} + | sort -k2 > MD5SUMS)
echo "Frozen: $DEST"
md5sum "$DEST/CellHubDB.rds"
