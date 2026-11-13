#!/usr/bin/env bash
# Build the files of the Zenodo deposit of CellHubDB 2.1 (doi:10.5281/zenodo.24197571) in output/zenodo_2.1/.
# Run from the CellHub_DB root after the v2.1 tag exists. Reads release/2.1/, output/gate/v2.1/, data/raw/ and the git tag.
set -euo pipefail
V=2.1; OUT=output/zenodo_$V; REL=release/$V
EXPECTED_MD5=84419abc3ee00cd516da00b99b0ba425
[[ -f "$REL/CellHubDB.rds" ]] || { echo "ERROR: $REL/CellHubDB.rds not found"; exit 1; }
[[ "$(md5sum < "$REL/CellHubDB.rds" | cut -d' ' -f1)" == "$EXPECTED_MD5" ]] || { echo "ERROR: $REL/CellHubDB.rds is not the frozen 2.1 (MD5 differs)"; exit 1; }
git rev-parse -q --verify "refs/tags/v$V" >/dev/null || { echo "ERROR: tag v$V missing (create it before building the code archive)"; exit 1; }
command -v zip >/dev/null || { echo "ERROR: zip not installed (sudo apt install zip)"; exit 1; }
[[ -e "$OUT" ]] && { echo "ERROR: $OUT already exists; remove it to rebuild"; exit 1; }
mkdir -p "$OUT"
ROOT=$PWD

# 1. Database bundle
cp "$REL/CellHubDB.rds" "$OUT/CellHubDB_$V.rds"
# 2. Per-species files loaded by Cell-Hub (merged + per source, audit CSVs)
(cd "$REL" && zip -qr -X "$ROOT/$OUT/CellHubDB_${V}_app_databases.zip" app_databases)
# 3. Release reports: numbers report, release gate, step logs, finalize repairs and removal audits
STAGE=$(mktemp -d); mkdir -p "$STAGE/CellHubDB_${V}_reports"
cp -r "$REL/reports" "$STAGE/CellHubDB_${V}_reports/pipeline_reports"
[[ -d "$REL/report_11" ]] && cp -r "$REL/report_11" "$STAGE/CellHubDB_${V}_reports/numbers_report"
[[ -d "$REL/report" ]] && cp -r "$REL/report" "$STAGE/CellHubDB_${V}_reports/numbers_report"
if [[ -d "$REL/gate" ]]; then cp -r "$REL/gate" "$STAGE/CellHubDB_${V}_reports/release_gate"; else cp -r "output/gate/v$V" "$STAGE/CellHubDB_${V}_reports/release_gate"; fi
cp "$REL/MD5SUMS" "$STAGE/CellHubDB_${V}_reports/release_MD5SUMS.txt"
(cd "$STAGE" && zip -qr -X "$ROOT/$OUT/CellHubDB_${V}_reports.zip" "CellHubDB_${V}_reports")
# 4. Raw source files without a versioned public download whose licence allows redistribution (CellTalkDB GPL-3.0, connectomeDB2025 MIT)
RAW="$STAGE/CellHubDB_${V}_raw_archived_sources"; mkdir -p "$RAW/CellTalkDB" "$RAW/ConnectomeDB2025" "$RAW/LICENSES"
cp data/raw/CellTalkDB/human_lr_pair.txt data/raw/CellTalkDB/mouse_lr_pair.txt "$RAW/CellTalkDB/"
cp data/raw/ConnectomeDB2025/all_species.csv "$RAW/ConnectomeDB2025/"
cp config/source_retrieval.tsv "$RAW/"
curl -fsSL https://raw.githubusercontent.com/ZJUFanLab/CellTalkDB/7276cd4e4dcdd3c5d3a719b962a4fd6a8bfc4b29/LICENSE -o "$RAW/LICENSES/CellTalkDB_LICENSE_GPL-3.0.txt"
cat > "$RAW/LICENSES/connectomeDB2025_NOTICE.txt" <<'NOTICE'
connectomeDB2025 (https://connectomedb.org), file all_species.csv downloaded on 2026-08-20.
Licence statement of the connectomeDB download page (github.com/bioinfo-YCU/ConnectomeDB, downloads.qmd):
"All ConnectomeDB ligand-receptor interaction data are freely available under the MIT License and are provided AS IS, without any warranty."
Citation: Liu P, Maezono S, et al. connectomeDB2025: a rigorously curated, multi-species resource of experimentally supported
ligand-receptor interactions. Nucleic Acids Research 54(D1):D546-D554 (2026). doi:10.1093/nar/gkaf1108

MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files
(the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge,
publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH
THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
NOTICE
{ echo "PlantCellChatDB raw files are NOT included: no licence is stated by their authors."
  echo "Download them from the PlantCellChat web site (https://bioinformatics.fafu.edu.cn/PlantCellChat/) and check them against these MD5:"
  awk -F'\t' 'NR > 1 && $1 ~ /^PlantCellChatDB\/(ath|gmx|osa|sly|zma)\.csv$/ {print $2 "  " $1}' config/source_retrieval.tsv
  echo "Then run: Rscript pipeline/00_download_sources.R"; } > "$RAW/PlantCellChatDB_NOT_INCLUDED.txt"
(cd "$RAW" && md5sum CellTalkDB/* ConnectomeDB2025/* > MD5SUMS.txt)
(cd "$STAGE" && zip -qr -X "$ROOT/$OUT/CellHubDB_${V}_raw_archived_sources.zip" "CellHubDB_${V}_raw_archived_sources")
# 5. Code at the release tag
git archive --format=zip --prefix="CellHub-Database-$V/" -o "$OUT/CellHub-Database-v$V-code.zip" "v$V"
rm -rf "$STAGE"
# 6. Checksums of the deposit
(cd "$OUT" && md5sum CellHubDB_$V.rds *.zip > MD5SUMS.txt && cat MD5SUMS.txt)
ls -lh "$OUT"
grep -c . "$OUT/MD5SUMS.txt" | xargs -I{} echo "{} files checksummed in $OUT/MD5SUMS.txt"
