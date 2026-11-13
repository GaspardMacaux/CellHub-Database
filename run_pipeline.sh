#!/usr/bin/env bash
# Rebuild CellHubDB end to end from the CellHub_DB root. Logs and reports are written to reports/.
# Usage: bash run_pipeline.sh              -> geneInfo, curation, QC, finalize, build, release gate, report (uses existing converted/)
#        bash run_pipeline.sh --convert    -> also re-runs every converter on raw_databases/
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
mkdir -p reports
S=database_curation
VERSION=2.1

# Run one R script with its arguments and keep its full console output in reports/<name>.txt
step() {
  local log="reports/$1.txt"; shift
  echo "== $(date +%H:%M:%S) $1"
  Rscript "$@" 2>&1 | tee "$log"
}

step provenance "$S/collect_source_provenance.R"
if [[ "${1:-}" == "--convert" ]]; then
  for conv in "$S"/convert_databases/convert_*.R; do step "convert_$(basename "$conv" .R)" "$conv"; done
  step qc_analyze_converted "$S/analyze_converted.R"
fi
# Genome-wide geneInfo from Ensembl GTF; always run, because convert_cellchatdb.R resets the human, mouse and zebrafish tables
step geneinfo "$S/build_geneinfo.R"
step curation_log_all_sources "$S/curate_CellHub_DB.R"
step qc_check_before_freeze "$S/check_before_freeze.R"
step qc_verify_integrity "$S/verify_integrity.R"
step qc_analyze_databases "$S/analyze_databases.R"
# CellChat schema, release repairs, interaction types, pathways and pair-level confidence
step finalize "$S/finalize_cellhubdb.R"
step build "$S/build_cellhubdb.R"
step stats "$S/cellhubdb_stats.R"
step upset "$S/upset_plot.R"
# Release gate (stops the pipeline on FAIL) and numbers report, both computed from the released file
step release_gate review/09_release_gate.R CellHubDB.rds databases "v$VERSION"
step report review/11_database_report.R CellHubDB.rds "$VERSION"

# Move generated tables and reports from the root into reports/ (CellHubDB.rds and summary CSVs stay at the root)
for f in CellHubDB_Statistics_Report.md source_provenance.csv package_versions.csv; do [[ -e "$f" ]] && mv -f "$f" reports/; done
echo "== Done. Release files: CellHubDB.rds, CellHubDB_species_summary.csv, CellHubDB_source_contribution.csv"
