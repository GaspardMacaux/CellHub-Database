#!/usr/bin/env bash
# Reorganise the CellHub_DB folder (2026-11-13). Run from the CellHub_DB root.
# Usage: bash migrate_cellhub_db.sh <new_code.tgz>            -> dry run: prints every operation, changes nothing
#        bash migrate_cellhub_db.sh <new_code.tgz> --apply    -> performs the moves, writes undo_migration.sh
# Nothing is deleted except macOS metadata files (._*, .DS_Store) and R history files; every move is reversed by undo_migration.sh.
set -euo pipefail
CODE_TGZ="${1:?usage: bash migrate_cellhub_db.sh <new_code.tgz> [--apply]}"
APPLY=0; [[ "${2:-}" == "--apply" ]] && APPLY=1
CODE_TGZ="$(cd "$(dirname "$CODE_TGZ")" && pwd)/$(basename "$CODE_TGZ")"
UNDO=undo_migration.sh
MOVED=()

# Preconditions: old layout present, new layout absent, code archive readable, git state committed
[[ -f run_pipeline.sh && -d database_curation && -d review ]] || { echo "ERROR: run this from the CellHub_DB root (old layout not found)"; exit 1; }
for d in data output pipeline R config revision archive; do [[ -e "$d" ]] && { echo "ERROR: $d/ already exists: migration already done or partial"; exit 1; }; done
tar tzf "$CODE_TGZ" >/dev/null || { echo "ERROR: cannot read $CODE_TGZ"; exit 1; }
if [[ -d .git && -n "$(git status --porcelain --untracked-files=no)" ]]; then echo "ERROR: uncommitted changes in git; commit them first"; git status --short | head; exit 1; fi
[[ -d .git ]] || echo "WARNING: not a git repository; only undo_migration.sh can revert the code changes"

# Move one path (creating the destination folder); missing sources are reported and skipped
mv_op() {
  local src="$1" dst="$2"
  if [[ ! -e "$src" ]]; then echo "  skip (absent)  $src"; return 0; fi
  echo "  move  $src  ->  $dst"
  if (( APPLY )); then mkdir -p "$(dirname "$dst")"; mv -n "$src" "$dst"; if [[ -e "$src" || ! -e "$dst" ]]; then echo "ERROR: move failed: $src"; exit 1; fi; MOVED+=("$src|$dst"); fi
  return 0
}
# Delete macOS metadata and history files only
rm_op() {
  echo "  delete  $1"
  if (( APPLY )); then rm -f -- "$1"; fi
  return 0
}

if (( APPLY )); then echo "=== APPLYING migration"; else echo "=== DRY RUN (nothing is changed; add --apply to perform)"; fi

echo "--- 1. Data -> data/"
mv_op raw_databases/extract_all_databases.R archive/code/raw_databases/extract_all_databases.R
mv_op raw_databases data/raw
mv_op converted data/converted
mv_op geneinfo data/geneinfo
mv_op databases data/databases

echo "--- 2. Last pipeline run -> output/"
for f in CellHubDB.rds CellHubDB_species_summary.csv CellHubDB_source_contribution.csv CellHubDB_Statistics_Report.md source_provenance.csv package_versions.csv; do mv_op "$f" "output/$f"; done
mv_op figures output/figures
mv_op reports output/reports
mv_op review/09_gate output/gate
mv_op review/11_report output/report

echo "--- 3. Revision results -> revision/results/"
for d in 04_benchmark_v2.1 06_symbols 12_evidence 13_conservation_v2.1 logs; do mv_op "review/$d" "revision/results/$d"; done

echo "--- 4. Superseded releases, old code and documents -> archive/"
mv_op release/2.1_superseded_338ba8a1 archive/release/2.1_superseded_338ba8a1
mv_op release_v2.1 archive/release/release_v2.1_dir_2.0.2
mv_op app_databases_v2.0.1 archive/release/app_databases_v2.0.1
mv_op Archives archive/Archives_before_2026-11-13
mv_op database_curation archive/code/database_curation
mv_op review archive/code/review
mv_op run_pipeline.sh archive/code/run_pipeline_old.sh
mv_op tidy_cellhub_db.sh archive/code/tidy_cellhub_db.sh
mv_op README.md archive/docs/README_v2.0.md
mv_op CellHub_DB_README.md archive/docs/GaspouDB_README.md
mv_op inventory.txt archive/docs/inventory_2026-11-13.txt
mv_op hgnc_complete_set.txt archive/data/hgnc_complete_set.txt
for f in cellhub_db_scripts_*.tgz; do if [[ -e "$f" && "$PWD/$f" != "$CODE_TGZ" ]]; then mv_op "$f" "archive/docs/$f"; fi; done

echo "--- 5. macOS metadata and R history (deleted)"
while IFS= read -r -d '' f; do rm_op "$f"; done < <(find . \( -name '._*' -o -name '.DS_Store' \) -not -path './.git/*' -print0)
for f in .Rhistory .Rapp.history; do if [[ -e "$f" ]]; then rm_op "$f"; fi; done

echo "--- 6. New code tree (from $(basename "$CODE_TGZ"))"
tar tzf "$CODE_TGZ" | grep -v '/$' | sed 's/^/  add  /'
if (( APPLY )); then
  if [[ -e .gitignore ]]; then mkdir -p archive/docs; mv .gitignore archive/docs/gitignore_old; fi
  tar xzf "$CODE_TGZ"
  # Undo script: remove the new code files, then reverse every move (last first)
  { echo "#!/usr/bin/env bash"; echo "# Reverts migrate_cellhub_db.sh (deleted macOS metadata files are not restored)"; echo "set -euo pipefail"
    tar tzf "$CODE_TGZ" | grep -v '/$' | sed 's/^/rm -f -- "/; s/$/"/'
    for (( i=${#MOVED[@]}-1; i>=0; i-- )); do IFS='|' read -r s d <<< "${MOVED[$i]}"; echo "mkdir -p \"$(dirname "$s")\" && mv -n \"$d\" \"$s\""; done
    echo 'if [[ -e archive/docs/gitignore_old ]]; then mv archive/docs/gitignore_old .gitignore; fi'
    echo 'find pipeline R config revision output data archive -type d -empty -delete 2>/dev/null || true'
    echo 'echo "Migration reverted"'; } > "$UNDO"
  echo "=== Done. Undo with: bash $UNDO"
  echo "=== Next: bash run_pipeline.sh, then compare output/CellHubDB.rds with release/2.1/CellHubDB.rds"
else
  echo "=== Dry run finished: check the list above, then re-run with --apply"
fi
