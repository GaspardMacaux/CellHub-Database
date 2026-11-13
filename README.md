# Cell-Hub Database (CellHubDB) — Multi-Species Curated Ligand–Receptor Interaction Database

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24197571.svg)](https://doi.org/10.5281/zenodo.24197571)

**Version:** 2.1
**Build date:** 2026-11-13
**Release file:** `CellHubDB_2.1.rds` (MD5 `84419abc3ee00cd516da00b99b0ba425`), on Zenodo: [doi:10.5281/zenodo.24197571](https://doi.org/10.5281/zenodo.24197571) and attached to the GitHub release v2.1
**Content:** 76,582 interactions · 21 species · 10 source databases · pair-level confidence annotation
**Format:** CellChat-native RDS (ready to use in Cell-Hub or in CellChat)
**Preprint:** doi:10.64898/2026.09.23.753685

---

## Overview

CellHubDB consolidates 10 independently curated ligand–receptor (LR) resources into one standardized, CellChat-compatible format. It provides:

- **CellChat-native files.** One `.rds` file per species and source, a merged file per species, and a single `CellHubDB.rds` bundle. Every file passes a full CellChat run in the release gate.
- **21 species in one schema.** 11 mammals, 3 other vertebrates, *Drosophila*, *C. elegans* and 5 plants.
- **Per-interaction provenance.** A `sources` column lists every database that contains the pair, and `curated_sources` lists those whose evidence for this pair is curated.
- **Pair-level confidence.** Each merged interaction carries a confidence tier (`high`, `medium`, `low`) computed from the evidence of each supporting source, plus flags for prediction-only pairs, complexes and corroborated ligand/receptor roles.
- **Standardized curation.** Species-specific nomenclature, genome-wide Ensembl gene tables, protected complex names, a definition for every heteromeric partner, and audit files for every removed entry.
- **Reproducible build.** Every raw file is re-obtained from a pinned origin and checked by checksum; one command rebuilds the release bit for bit.
- **No-code access.** Species, source and confidence selection through Cell-Hub's graphical interface.

---

## Species coverage

| Species | Organism | Sources | Merged interactions | Multi-source | Complexes | Tiers high / medium / low |
|---|---|---|---|---|---|---|
| `human` | *Homo sapiens* | CellChatDB v2, CellPhoneDB, CellTalkDB, connectomeDB2025, NicheNet v2, NeuronChatDB | 11,503 | 3,834 (33.3%) | 534 | 2,117 / 7,159 / 2,227 |
| `mouse` | *Mus musculus* | CellChatDB v2, CellTalkDB, connectomeDB2025, NicheNet v2, NeuronChatDB | 9,800 | 2,955 (30.2%) | 337 | 906 / 4,766 / 4,128 |
| `arabidopsis` | *Arabidopsis thaliana* | PlantCellChatDB, PlantPhoneDB | 7,191 | 39 (0.5%) | 0 | 38 / 3,940 / 3,213 |
| `zebrafish` | *Danio rerio* | CellChatDB v2, connectomeDB2025 | 4,947 | 284 (5.7%) | 175 | 0 / 17 / 4,930 |
| `chimp` | *Pan troglodytes* | connectomeDB2025 | 3,444 | – | 0 | 0 / 0 / 3,444 |
| `macaque` | *Macaca mulatta* | connectomeDB2025 | 3,417 | – | 0 | 0 / 0 / 3,417 |
| `rat` | *Rattus norvegicus* | connectomeDB2025 | 3,315 | – | 0 | 0 / 158 / 3,157 |
| `cow` | *Bos taurus* | connectomeDB2025 | 3,295 | – | 0 | 0 / 3 / 3,292 |
| `horse` | *Equus caballus* | connectomeDB2025 | 3,285 | – | 0 | 0 / 0 / 3,285 |
| `dog` | *Canis lupus familiaris* | connectomeDB2025 | 3,249 | – | 0 | 0 / 0 / 3,249 |
| `pig` | *Sus scrofa* | connectomeDB2025 | 3,244 | – | 0 | 0 / 0 / 3,244 |
| `marmoset` | *Callithrix jacchus* | connectomeDB2025 | 3,183 | – | 0 | 0 / 0 / 3,183 |
| `sheep` | *Ovis aries* | connectomeDB2025 | 3,121 | – | 0 | 0 / 0 / 3,121 |
| `frog` | *Xenopus tropicalis* | connectomeDB2025 | 2,776 | – | 0 | 0 / 21 / 2,755 |
| `chicken` | *Gallus gallus* | connectomeDB2025 | 2,593 | – | 0 | 0 / 28 / 2,565 |
| `maize` | *Zea mays* | PlantCellChatDB | 1,913 | – | 0 | 0 / 855 / 1,058 |
| `rice` | *Oryza sativa* | PlantCellChatDB | 1,810 | – | 0 | 0 / 189 / 1,621 |
| `drosophila` | *Drosophila melanogaster* | FlyPhoneDB2 | 1,798 | – | 0 | 0 / 310 / 1,488 |
| `soybean` | *Glycine max* | PlantCellChatDB | 1,412 | – | 0 | 0 / 245 / 1,167 |
| `tomato` | *Solanum lycopersicum* | PlantCellChatDB | 1,041 | – | 0 | 0 / 139 / 902 |
| `celegans` | *Caenorhabditis elegans* | Armingol et al. 2022 | 245 | – | 0 | 0 / 181 / 64 |
| **Total** | | | **76,582** | **7,112** | | **3,061 / 18,011 / 55,510** |

Most vertebrate entries outside human and mouse are orthology projections (connectomeDB2025 "Inferred"; the CellChatDB zebrafish database is projected from mouse) and are therefore in the `low` tier. See *Confidence annotation* and *Known limitations*.

---

## Source databases

| Source | Version / origin | Species | Interactions after curation | Reference |
|---|---|---|---|---|
| CellChatDB v2 | CellChat 2.2.0.9001, commit `75253cd` | human, mouse, zebrafish | 3,233 / 3,379 / 2,463 | Jin et al., 2024 |
| CellPhoneDB | `ventolab/cellphonedb-data` tag v5.0.0 | human | 2,909 | Troulé et al., 2025 |
| CellTalkDB | text release of tcm.zju.edu.cn/celltalkdb (archived) | human, mouse | 3,398 / 2,033 | Shao et al., 2021 |
| connectomeDB2025 | `all_species.csv` export of connectomedb.org (archived) | 14 vertebrates | 44,468 in total | Liu et al., 2025 |
| NicheNet v2 LR network | Zenodo 10229222, `lr_network_*_allInfo_30112033.rds` | human, mouse | 5,086 / 5,596 | Browaeys et al., 2023 |
| NeuronChatDB | NeuronChat 1.0.0, commit `36e06b2` | human, mouse | 471 / 452 | Zhao et al., 2023 |
| FlyPhoneDB2 | FlyPhone 0.1, commit `49649b6` | *Drosophila* | 1,798 | FlyPhoneDB2 (2025) |
| *C. elegans* LR pairs | `LewisLabUCSD/Ligand-Receptor-Pairs` commit `c0ce816` | *C. elegans* | 245 | Armingol et al., 2022 |
| PlantPhoneDB | same repository and commit (Xu et al. LR tables) | *Arabidopsis* | 3,514 | Xu et al., 2022 |
| PlantCellChatDB | web release of bioinformatics.fafu.edu.cn/PlantCellChat (archived) | *Arabidopsis*, maize, rice, soybean, tomato | 3,716 / 1,913 / 1,810 / 1,412 / 1,041 | PlantCellChatDB |

Maize identifiers are converted from B73 v4 to v5 with the MaizeGDB cross-reference `B73v4_to_B73v5.tsv` (one-to-one correspondences only). Every raw file, its origin and its MD5 are listed in `config/source_retrieval.tsv` and `config/sources_manifest.tsv` (see *Reproducibility*).

---

## Confidence annotation

Merged files (`CellHubDB_<species>_merged.rds` and `CellHubDB.rds`) carry pair-level annotation columns computed in `pipeline/09_finalize.R`:

| Column | Description |
|---|---|
| `curated_sources`, `n_curated_sources` | Sources whose evidence for this pair is curated |
| `predicted_only` | TRUE when no source gives curated evidence (orthology-inferred or computationally predicted only) |
| `confidence_tier` | `high` = at least two curated sources; `medium` = one curated source; `low` = no curated source |
| `is_complex` | Ligand or receptor is a multi-subunit complex |
| `role_corroborated` | Ligand used as a ligand and receptor used as a receptor by at least two curated sources (NA for single-source species) |

A source's evidence for a pair is curated unless the source marks it otherwise:

| Source | Counted as curated |
|---|---|
| connectomeDB2025 | "Direct" evidence only ("Inferred" = orthology projection) |
| CellChatDB v2 | human and mouse; the zebrafish database is an orthology projection of the mouse database and is not counted |
| CellPhoneDB, CellTalkDB, NeuronChatDB | all entries |
| FlyPhoneDB2 | high-confidence pairs |
| *C. elegans* (Armingol) | "Literature supported" pairs |
| PlantPhoneDB | literature, IntAct, BioGRID, interactome2.0 and plant.MAP entries |
| PlantCellChatDB | "experimentally validated" entries |
| NicheNet v2 | never (computationally integrated network) |

**Recommended set (human, mouse):** `high`, plus `medium` with `role_corroborated == TRUE`. The tier grades the evidence that a pair exists; it does not measure cell-type specificity, and a `high` pair is not guaranteed to be active in a given dataset. For species where nearly all pairs are `low`, use all tiers and interpret results as orthology-based hypotheses.

---

## File structure

```
README.md, run_pipeline.sh, freeze_release.sh
R/cellhubdb_functions.R        # shared functions: CellChat schema, repairs, pathways, pair-level confidence
pipeline/                      # numbered in run order (see "Curation pipeline")
config/                        # source_retrieval.tsv, sources_manifest.tsv, complex_reconstructions.tsv
data/                          # not versioned on GitHub (archived on Zenodo)
├── raw/                       #   source files as downloaded or extracted from R packages
├── converted/                 #   one standardized table per source and species
├── geneinfo/                  #   Ensembl-based gene tables (gtf/ = cached Ensembl GTFs)
└── databases/<species>/       #   working copy: CellHubDB_<species>_<source>.rds, CellHubDB_<species>_merged.rds, audit CSVs
output/                        # last pipeline run: CellHubDB.rds, summary CSVs, figures/, reports/ (logs), gate/, report/
release/<version>/             # frozen releases: CellHubDB.rds, app_databases/, reports, gate, report, MD5SUMS
revision/                      # analyses for the manuscript revisions; results/ holds their outputs
archive/                       # superseded scripts, intermediate releases and outputs (kept, not used)
```

All scripts are run from the repository root. `release/<version>/app_databases/` holds the per-species folders loaded by Cell-Hub.

---

## Data format

Each database is a CellChat-compatible list:

```r
db <- readRDS("release/2.1/app_databases/mouse/CellHubDB_mouse_merged.rds")
db$interaction  # LR pairs with annotations, provenance and confidence
db$complex      # Multi-subunit complex definitions (subunit_1, subunit_2, ...)
db$cofactor     # Agonists, antagonists and co-receptors (CellChatDB-derived)
db$geneInfo     # Genome-wide gene table (Ensembl 112; Ensembl Plants 59)
```

### Interaction table columns

| Column | Description |
|---|---|
| `interaction_name` | Primary key `LIGAND_RECEPTOR` |
| `interaction_name_2` | Readable form, e.g. `Tgfb1 - (Tgfbr1+Tgfbr2)` |
| `ligand`, `receptor` | Gene symbol or complex name (species nomenclature) |
| `pathway_name` | Pathway from the source, inherited from the same human or mouse pair, or `Unassigned` |
| `interaction_type` | `Secreted Signaling`, `Cell-Cell Contact`, `ECM-Receptor`, `Non-protein Signaling` or NA |
| `annotation` | Functional annotation |
| `evidence` | PMIDs, KEGG identifiers or DOIs |
| `evidence_type` | connectomeDB2025 `Direct` or `Inferred` (NA for other sources) |
| `source_db` | Source whose record was retained |
| `sources` | All sources containing the pair (comma-separated) |
| `curated_sources` … `role_corroborated` | Confidence annotation (merged files only, see above) |

### Gene identifier conventions

| Species | Convention | Example |
|---|---|---|
| human, chimp, macaque, marmoset, cow, dog, horse, pig, sheep, chicken | HGNC-style uppercase | `TGFB1` |
| mouse, rat | MGI/RGD (first letter uppercase) | `Tgfb1` |
| zebrafish | ZFIN lowercase; genes without a ZFIN symbol keep their upper-case Ensembl name | `tgfb1`, `FAT4` |
| frog | Xenbase, as provided (lowercase) | `tgfb1` |
| drosophila | FlyBase, as provided | `dpp`, `Toll` |
| celegans | WormBase, as provided | `daf-7` |
| arabidopsis | TAIR AGI locus | `AT1G01900` |
| rice | RAP-DB locus | `Os01g0958100` |
| maize | B73 v5 gene model | `Zm00001eb375620` |
| soybean | Phytozome Wm82 (converted from `GLYMA_`) | `Glyma.07G056900` |
| tomato | ITAG gene model | `Solyc10g076410.1` |

Plant identifiers must match the annotation used to build your count matrix. Gene models from another annotation release will not match.

---

## Curation pipeline

```
bash run_pipeline.sh [--convert]        full rebuild into output/; bash freeze_release.sh <version> copies it to release/<version>/

pipeline/00_download_sources.R          re-obtain every raw file from its pinned origin and check it   → data/raw/
pipeline/01_collect_source_provenance.R → output/source_provenance.csv, output/package_versions.csv
pipeline/02_convert/convert_*.R         one converter per source                     → data/converted/
pipeline/03_analyze_converted.R         QC of converted tables
pipeline/04_build_geneinfo.R            Ensembl gene tables                          → data/geneinfo/
pipeline/05_curate.R                    cleaning, nomenclature, validation, merge     → data/databases/<species>/
pipeline/06–08_*.R                      QC (pre-freeze checks, integrity, annotation coverage)
pipeline/09_finalize.R                  CellChat schema, repairs, pathways, confidence (in place)
pipeline/10_build.R                     → output/CellHubDB.rds + summary CSVs
pipeline/11_stats.R, 12_upset_plot.R    statistics report and figures
pipeline/13_release_gate.R              release gate (stops on FAIL)                  → output/gate/
pipeline/14_database_report.R           numbers report                                → output/report/
```

### Curation steps

1. **Whitespace and empty values.** Whitespace is trimmed and empty strings become NA. Entries without a ligand or receptor are removed.
2. **Pathways and annotations.** Source names used as pathways, and generic annotations, become NA.
3. **Interaction types.** Variant spellings are mapped to the four canonical types.
4. **Unresolved UniProt accessions.** Tokens that are neither known symbols nor complex names are removed and audited.
5. **Nomenclature.** Symbols are normalized per species; complex names keep their exact spelling. Maize B73 v4 identifiers are converted to v5; interactions keeping a v4 identifier are removed and audited.
6. **Symbol validation.** Symbols absent from the genome-wide Ensembl gene table of the species are removed and audited.
7. **Identifier rebuild.** Interaction names are rebuilt and case-induced duplicates removed.
8. **Evidence.** References are normalized.
9. **Merge.** Pairs are merged on a canonical key that ignores case and, for non-oriented sources, partner order, using the source priority below. Metadata from lower-priority sources fill the gaps in the retained record, and a `sources` provenance column is added.
10. **Complex tables.** Curated definitions are kept. Heteromeric partners without a definition are defined from their underscore-separated subunits. Zebrafish complexes used by CellChatDB but undefined in its release are reconstructed from paralog combinations (`config/complex_reconstructions.tsv`).
11. **Finalization** (`09_finalize.R`). CellChat schema completed on every file; undefined cofactors, cofactor genes in another case, unresolvable partners and complexes with a numeric subunit are repaired or removed and audited; `interaction_type` is restored from the CellChat annotation; pathways are filled; pair-level confidence is computed.

### Source priority (merge)

```
CellChatDB_v2 > connectomeDB2025 > CellPhoneDB > NeuronChatDB >
FlyPhoneDB2 = cell2cell_Armingol = PlantPhoneDB = PlantCellChatDB > CellTalkDB > NicheNet v2
```

### Quality control

| Script | Checks |
|---|---|
| `pipeline/07_verify_integrity.R` | Traceability from converted to per-source to merged files (exhaustive), multi-source consistency, row names |
| `pipeline/08_analyze_databases.R` | Per-species annotation coverage, casing, complex/cofactor/geneInfo tables |
| `pipeline/11_stats.R` | Statistics report (`output/CellHubDB_Statistics_Report.md`) |
| `pipeline/12_upset_plot.R` | Source overlap figures |
| `pipeline/13_release_gate.R` | Every reviewer check plus a CellChat smoke test on every file; the release is not frozen on FAIL |

Release 2.1 gate: 0 FAIL, 15 WARN (documented below), 1,034 PASS.

---

## Reproducibility

```bash
Rscript pipeline/00_download_sources.R --install   # fetch data/raw/ from pinned origins, check MD5 or content
bash run_pipeline.sh --convert                      # rebuild output/CellHubDB.rds
bash freeze_release.sh 2.1                          # copy to release/2.1/ with MD5SUMS
```

- **Raw sources.** 25 of the 33 checked raw files are re-obtained from pinned origins (GitHub commits or tags, Zenodo records, R packages at fixed commits, MaizeGDB) and verified by MD5 or, for re-serialised R objects, by content. CellTalkDB, connectomeDB2025 and PlantCellChatDB have no versioned download (8 raw files, checked by MD5). The CellTalkDB (GPL-3.0) and connectomeDB2025 (MIT) files are distributed in the Zenodo deposit ([doi:10.5281/zenodo.24197571](https://doi.org/10.5281/zenodo.24197571)) with their licences; the PlantCellChatDB files carry no licence and are not redistributed: obtain them from the PlantCellChat web site and check them against the MD5 in `config/source_retrieval.tsv`.
- **Gene tables.** Ensembl 112 (vertebrates, *Drosophila*, *C. elegans*) and Ensembl Plants 59 GTFs; files and MD5 in `data/geneinfo/gtf_manifest.tsv`.
- **Build.** Running the pipeline on these inputs reproduces `release/2.1/CellHubDB.rds` bit for bit (same MD5 when built on the same day; `attr(, "build_date")` is the only date-dependent field).

---

## Known limitations

- **Orthology-inferred content.** In 12 of the 14 vertebrates, fewer than 5% of pairs have curated support (connectomeDB2025 projects human and mouse pairs to the other species; CellChatDB zebrafish is projected from mouse). Overlap between these species reflects the projection, not independent evidence of conservation.
- **Single-source pairs.** About two thirds of human and mouse pairs come from one database. Use the confidence tier and the recommended set to restrict analyses.
- **Gene coverage.** 80–85% of interactions are scorable (all genes present in the Ensembl gene table) for dog, frog, horse, soybean and zebrafish; symbols absent from Ensembl names cannot be matched to count matrices built on Ensembl annotation.
- **Interaction types in plants.** 82–97% of plant interactions have no canonical interaction type, because the plant sources do not provide one.
- **Pathways.** Pairs without a source pathway inherit the pathway of the same human or mouse pair when it exists, otherwise `Unassigned` (CellChat requires a pathway for every interaction).
- **Zebrafish complexes.** `ITGA6_ITGB1_1/2` cannot be reconstructed unambiguously and their interactions are removed (audited in the finalize reports).
- **Reciprocal pairs.** The merge key ignores partner order for non-oriented sources, so a pair listed in both directions within such a source is kept once.

---

## Usage

### In Cell-Hub

1. Select the species.
2. Load the merged database, or tick individual sources; filter by confidence tier if needed.
3. Run CellChat. When several sources are selected, they are merged with the same rules on the fly and the tiers are recomputed on the selected sources.

### In R

```r
library(CellChat)
cellhubdb <- readRDS("CellHubDB.rds")
attr(cellhubdb, "version")              # "2.1"
db <- cellhubdb[["mouse"]]              # merged mouse database
cellchat <- createCellChat(seurat_obj, group.by = "cell_type")
# Recommended set: high tier, plus medium tier with corroborated ligand and receptor roles
keep <- db$interaction$confidence_tier == "high" | (db$interaction$confidence_tier == "medium" & db$interaction$role_corroborated %in% TRUE)
db_rec <- db
db_rec$interaction <- db$interaction[keep, ]
cellchat@DB <- db_rec
```

---

## Changes from v2.0

- **Confidence annotation.** New pair-level columns (`curated_sources`, `confidence_tier`, `predicted_only`, `is_complex`, `role_corroborated`) from the evidence of each source; CellChatDB zebrafish counted as orthology-inferred.
- **Genome-wide gene tables.** `geneInfo` built from Ensembl 112 / Ensembl Plants 59 for all 21 species; symbol validation no longer depends on the CellChatDB gene table (DEFB4A, S100A8/A9, P2RY12 and other valid symbols are kept).
- **CellChat readiness of every file.** Complete complex and cofactor tables, repaired cofactors and partners, interaction types restored from the CellChat annotation (for example 429 human ECM-Receptor interactions), pathways for every interaction. All 21 merged and all per-source files pass a CellChat run in the release gate.
- **Maize.** Identifiers converted from B73 v4 to v5 (MaizeGDB); interactions keeping a v4 identifier removed (maize 2,268 → 1,913; scorable 57.6% → 91.9%).
- **Zebrafish complexes.** 10 complexes used but undefined in CellChatDB zebrafish reconstructed from paralog combinations.
- **Reproducibility.** Pinned retrieval of every raw source, single build script, release gate and frozen releases with checksums.

v2.0 had rebuilt the database from versioned raw sources and fixed defects of v1.0 (complex names altered by case normalization, an over-broad UniProt filter, zebrafish validation against a non-genome-wide gene table, uppercased *Xenopus* and crop identifiers).

Previous versions are kept as tags of this repository: [v2.0](https://github.com/GaspardMacaux/CellHub-Database/tree/v2.0) (version analysed in the first review round) and [v1.0](https://github.com/GaspardMacaux/CellHub-Database/tree/v1.0) (version described in the preprint). They are superseded by 2.1.

---

## Source licences

Each source keeps its own licence and terms of use; cite the source articles when you use their content.

| Source | Licence |
|---|---|
| CellChatDB v2 | GPL-3 |
| CellPhoneDB | MIT (CellPhoneDB software) |
| CellTalkDB | GPL-3.0 |
| connectomeDB2025 | MIT ("freely available under the MIT License and provided AS IS", connectomedb.org) |
| NicheNet v2 (nichenetr, multinichenetr) | GPL-3 |
| NeuronChatDB | GPL-3 |
| FlyPhoneDB2 | MIT |
| *C. elegans* LR pairs (Armingol et al.), PlantPhoneDB, PlantCellChatDB | no licence stated; content redistributed in derived form with citation |

The curation pipeline and CellHubDB release files are distributed under GPL-3.0.

---

## Citation

> Macaux G. et al. (2026). *bioRxiv* doi:10.64898/2026.09.23.753685.

Please also cite the original source databases used in your analysis.

## License

The curation pipeline, conversion scripts and quality-control code are released under the GNU General Public License v3.0 (GPL-3.0), consistent with Cell-Hub. The interaction data are derived from the source databases listed above, each of which retains its original license and terms of use; users redistributing derived data should comply with the terms of every source they use.
