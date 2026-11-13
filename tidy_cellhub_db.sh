#!/usr/bin/env bash
# Tidy the CellHub_DB folder after the v2.0 build. Nothing is deleted: development files move to Archives/.
# Usage: bash tidy_cellhub_db.sh            -> dry run, prints what would move
#        bash tidy_cellhub_db.sh --apply    -> performs the moves
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
APPLY=0; [[ "${1:-}" == "--apply" ]] && APPLY=1
ARCH="Archives/dev_2026-10-19"
cd "$ROOT"

# Move one file (if it exists) into a destination folder, never overwriting an existing file
move() {
  local src="$1" dest="$2"
  [[ -e "$src" ]] || return 0
  if [[ $APPLY -eq 1 ]]; then mkdir -p "$dest" && mv -n "$src" "$dest/" && echo "moved   $src -> $dest/"
  else echo "would move  $src -> $dest/"; fi
}

echo "== Patch scripts and one-off diagnostics -> $ARCH/patches"
for f in database_curation/patch_*.py database_curation/curate_step5_fixed.R database_curation/check_complexes.R \
         database_curation/*.bak database_curation/*~ database_curation/*.orig; do move "$f" "$ARCH/patches"; done

echo "== Superseded GaspouDB documents -> $ARCH/docs_v1"
for f in integration_plan.md CellHub_DB_manuscript_section.md GaspouDB_README.md \
         database_curation/integration_plan.md; do move "$f" "$ARCH/docs_v1"; done

echo "== QC logs, statistics and provenance -> reports/"
for f in qc_*.txt curation_log_all_sources.txt CellHubDB_Statistics_Report.md \
         source_provenance.csv package_versions.csv Table_source_databases.csv; do move "$f" reports; done

echo "== Files left at the root (check manually, nothing else is moved)"
for f in *; do
  case "$f" in
    CellHubDB.rds|CellHubDB_species_summary.csv|CellHubDB_source_contribution.csv|README.md|CellHub_DB_README.md|run_pipeline.sh|tidy_cellhub_db.sh) ;;
    databases|raw_databases|converted|geneinfo|database_curation|reports|figures|Archives|Script) ;;
    *) echo "  unclassified: $f" ;;
  esac
done
[[ -d Script ]] && echo "  note: Script/ still exists; list it (ls Script) before deciding whether to archive it"
[[ $APPLY -eq 0 ]] && echo "Dry run only. Re-run with --apply to perform the moves."
exit 0
