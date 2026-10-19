#!/usr/bin/env bash
# Rebuild CellHubDB end to end from the CellHub_DB root. Logs and reports are written to reports/.
# Usage: bash run_pipeline.sh              -> curation, QC, build, statistics (uses existing converted/)
#        bash run_pipeline.sh --convert    -> also re-runs every converter on raw_databases/
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
mkdir -p reports
S=database_curation

# Run one R script and keep its full console output in reports/<name>.txt
step() {
  local script="$1" log="reports/$2.txt"
  echo "== $(date +%H:%M:%S) $script"
  Rscript "$script" 2>&1 | tee "$log"
}

step "$S/collect_source_provenance.R" provenance
if [[ "${1:-}" == "--convert" ]]; then
  for conv in "$S"/convert_databases/convert_*.R; do step "$conv" "convert_$(basename "$conv" .R)"; done
  step "$S/analyze_converted.R" qc_analyze_converted
fi
step "$S/curate_CellHub_DB.R" curation_log_all_sources
step "$S/check_before_freeze.R" qc_check_before_freeze
step "$S/verify_integrity.R" qc_verify_integrity
step "$S/analyze_databases.R" qc_analyze_databases
step "$S/build_cellhubdb.R" build
step "$S/cellhubdb_stats.R" stats
step "$S/upset_plot.R" upset

# Move generated tables and reports from the root into reports/ (CellHubDB.rds and summary CSVs stay at the root)
for f in CellHubDB_Statistics_Report.md source_provenance.csv package_versions.csv; do [[ -e "$f" ]] && mv -f "$f" reports/; done
echo "== Done. Release files: CellHubDB.rds, CellHubDB_species_summary.csv, CellHubDB_source_contribution.csv"
