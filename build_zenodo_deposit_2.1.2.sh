#!/usr/bin/env bash
# Build the files of the Zenodo deposit of CellHubDB 2.1.2 (new version of doi:10.5281/zenodo.24197571) in output/zenodo_2.1.2/.
# Same layout as the 2.1 deposit. Run from the CellHub_DB root after the v2.1.2 tag exists. Reads release/2.1.2/, the 2.1 raw-source
# archive (unchanged sources) and the git tag.
set -euo pipefail
V=2.1.2; OUT=output/zenodo_$V; REL=release/$V
EXPECTED_MD5=c08462e1beeb9ae9013a87a86fe106d1
[[ -f "$REL/CellHubDB.rds" ]] || { echo "ERROR: $REL/CellHubDB.rds not found"; exit 1; }
[[ "$(md5sum < "$REL/CellHubDB.rds" | cut -d' ' -f1)" == "$EXPECTED_MD5" ]] || { echo "ERROR: $REL/CellHubDB.rds is not the frozen 2.1.2 (MD5 differs)"; exit 1; }
(cd "$REL" && md5sum -c MD5SUMS --quiet) || { echo "ERROR: $REL/MD5SUMS does not match the files"; exit 1; }
git rev-parse -q --verify "refs/tags/v$V" >/dev/null || { echo "ERROR: tag v$V missing (commit the README and patches, then: git tag v$V)"; exit 1; }
command -v zip >/dev/null || { echo "ERROR: zip not installed (sudo apt install zip)"; exit 1; }
[[ -e "$OUT" ]] && { echo "ERROR: $OUT already exists; remove it to rebuild"; exit 1; }
RAW_21=output/zenodo_2.1/CellHubDB_2.1_raw_archived_sources.zip
[[ -f "$RAW_21" ]] || { echo "ERROR: $RAW_21 not found (raw sources are unchanged since 2.1: take the file from the 2.1 deposit)"; exit 1; }
mkdir -p "$OUT"
ROOT=$PWD
# 1. Database bundle
cp "$REL/CellHubDB.rds" "$OUT/CellHubDB_$V.rds"
# 2. Per-species files loaded by Cell-Hub (merged + per source, audit CSVs)
(cd "$REL" && zip -qr -X "$ROOT/$OUT/CellHubDB_${V}_app_databases.zip" app_databases)
# 3. Release reports, including the two patch reports and change notes
STAGE=$(mktemp -d); mkdir -p "$STAGE/CellHubDB_${V}_reports"
cp -r "$REL/reports" "$STAGE/CellHubDB_${V}_reports/pipeline_reports"
for d in report_11 report; do [[ -d "$REL/$d" ]] && cp -r "$REL/$d" "$STAGE/CellHubDB_${V}_reports/numbers_report_2.1"; done
[[ -d "$REL/gate" ]] && cp -r "$REL/gate" "$STAGE/CellHubDB_${V}_reports/release_gate_2.1"
cp "$REL"/CHANGES_*.md "$STAGE/CellHubDB_${V}_reports/"
cp "$REL/MD5SUMS" "$STAGE/CellHubDB_${V}_reports/release_MD5SUMS.txt"
(cd "$STAGE" && zip -qr -X "$ROOT/$OUT/CellHubDB_${V}_reports.zip" "CellHubDB_${V}_reports")
rm -rf "$STAGE"
# 4. Raw source files: unchanged since 2.1, same archive under its 2.1 name
cp "$RAW_21" "$OUT/"
# 5. Code at the release tag (pipeline + patch scripts)
git archive --format=zip --prefix="CellHub-Database-$V/" -o "$OUT/CellHub-Database-v$V-code.zip" "v$V"
# 6. Checksums of the deposit
(cd "$OUT" && md5sum CellHubDB_$V.rds *.zip > MD5SUMS.txt && cat MD5SUMS.txt)
ls -lh "$OUT"
