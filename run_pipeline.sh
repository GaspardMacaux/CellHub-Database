#!/usr/bin/env bash
# Rebuild CellHubDB end to end from the CellHub_DB root.
# Usage: bash run_pipeline.sh              -> geneInfo, curation, QC, finalize, build, release gate, report (uses existing data/converted/)
#        bash run_pipeline.sh --convert    -> also re-runs every converter on data/raw/
# Inputs: data/raw/ (source files), config/ (manual tables), data/geneinfo/gtf/ (cached Ensembl GTFs)
# Outputs: output/CellHubDB.rds (release candidate), output/reports/ (one console log per step), output/gate/, output/report/
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
mkdir -p output/reports
P=pipeline
VERSION=2.1

# Run one R script with its arguments and keep its full console output in output/reports/<name>.txt
step() {
  local log="output/reports/$1.txt"; shift
  echo "== $(date +%H:%M:%S) $1"
  Rscript "$@" 2>&1 | tee "$log"
}

step provenance "$P/01_collect_source_provenance.R"
if [[ "${1:-}" == "--convert" ]]; then
  for conv in "$P"/02_convert/convert_*.R; do step "convert_$(basename "$conv" .R)" "$conv"; done
  step qc_analyze_converted "$P/03_analyze_converted.R"
fi
# Genome-wide geneInfo from Ensembl GTF; always run, because convert_cellchatdb.R resets the human, mouse and zebrafish tables
step geneinfo "$P/04_build_geneinfo.R"
step curation_log_all_sources "$P/05_curate.R"
step qc_check_before_freeze "$P/06_check_before_freeze.R"
step qc_verify_integrity "$P/07_verify_integrity.R"
step qc_analyze_databases "$P/08_analyze_databases.R"
# CellChat schema, release repairs, interaction types, pathways and pair-level confidence
step finalize "$P/09_finalize.R"
step build "$P/10_build.R"
step stats "$P/11_stats.R"
step upset "$P/12_upset_plot.R"
# Release gate (stops the pipeline on FAIL) and numbers report, both computed from the built file
step release_gate "$P/13_release_gate.R" output/CellHubDB.rds data/databases "v$VERSION"
step report "$P/14_database_report.R" output/CellHubDB.rds "$VERSION"
echo "== Done. Release candidate: output/CellHubDB.rds (freeze it with: bash freeze_release.sh $VERSION)"
