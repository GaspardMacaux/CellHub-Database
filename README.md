# CellHub Database (CellHubDB) — Multi-Species Curated Ligand–Receptor Interaction Database

**Version:** 2.0
**Build date:** 2026-10-19
**Content:** 73,589 curated interactions · 21 species · 10 source databases · 34,280 distinct LR pairs
**Format:** CellChat-native RDS (ready to use in Cell-Hub or in CellChat)
**Preprint:** doi:10.64898/2026.09.23.753685

---

## Overview

CellHubDB consolidates 10 independently curated ligand–receptor (LR) resources into one standardized, CellChat-compatible format. It provides:

- **CellChat-native files.** One `.rds` file per species and source, plus a merged file per species and a single `CellHubDB.rds` bundle.
- **21 species in one schema.** 11 mammals, 3 other vertebrates, *Drosophila*, *C. elegans* and 5 plants.
- **Per-interaction provenance.** A `sources` column lists every database that contains the pair.
- **Standardized curation.** Species-specific nomenclature, protected complex names, a complex definition for every heteromeric partner (human and mouse), and audit files for every removed entry.
- **No-code access.** Species and source selection through Cell-Hub's graphical interface.

---

## Species coverage

| Species | Organism | Sources | Merged interactions | Complexes |
|---|---|---|---|---|
| `human` | *Homo sapiens* | CellChatDB v2, CellPhoneDB, CellTalkDB, connectomeDB2025, NicheNet v2, NeuronChatDB | 11,089 | 534 |
| `mouse` | *Mus musculus* | CellChatDB v2, CellTalkDB, connectomeDB2025, NicheNet v2, NeuronChatDB | 9,400 | 337 |
| `arabidopsis` | *Arabidopsis thaliana* | PlantCellChatDB, PlantPhoneDB | 7,071 | 0 |
| `zebrafish` | *Danio rerio* | CellChatDB v2, connectomeDB2025 | 4,756 | 179 |
| `chimp` | *Pan troglodytes* | connectomeDB2025 | 3,222 | 0 |
| `macaque` | *Macaca mulatta* | connectomeDB2025 | 3,209 | 0 |
| `rat` | *Rattus norvegicus* | connectomeDB2025 | 3,103 | 0 |
| `cow` | *Bos taurus* | connectomeDB2025 | 3,082 | 0 |
| `horse` | *Equus caballus* | connectomeDB2025 | 3,069 | 0 |
| `pig` | *Sus scrofa* | connectomeDB2025 | 3,036 | 0 |
| `dog` | *Canis lupus familiaris* | connectomeDB2025 | 3,035 | 0 |
| `marmoset` | *Callithrix jacchus* | connectomeDB2025 | 2,992 | 1 |
| `sheep` | *Ovis aries* | connectomeDB2025 | 2,910 | 0 |
| `frog` | *Xenopus tropicalis* | connectomeDB2025 | 2,600 | 0 |
| `chicken` | *Gallus gallus* | connectomeDB2025 | 2,449 | 0 |
| `maize` | *Zea mays* | PlantCellChatDB | 2,268 | see note |
| `rice` | *Oryza sativa* | PlantCellChatDB | 1,810 | 0 |
| `drosophila` | *Drosophila melanogaster* | FlyPhoneDB2 | 1,790 | 0 |
| `soybean` | *Glycine max* | PlantCellChatDB | 1,412 | 3 |
| `tomato` | *Solanum lycopersicum* | PlantCellChatDB | 1,041 | 0 |
| `celegans` | *Caenorhabditis elegans* | Armingol et al. 2022 | 245 | 0 |
| **Total** | | | **73,589** | |

In human and mouse, 33.8% and 31.0% of interactions are supported by at least two sources. In zebrafish the figure is 5.9%, and in *Arabidopsis* 0.6%.

---

## Source databases

| Source | Version / access | Species | Interactions after curation | Reference |
|---|---|---|---|---|
| CellChatDB v2 | CellChat 2.2.0.9001 (75253cd), 2026-08-21 | human, mouse, zebrafish | 3,233 / 3,379 / 2,514 | Jin et al., 2024 |
| CellPhoneDB | v5.0.0 data release (2023-10-31) | human | 2,909 | Troulé et al., 2025 |
| CellTalkDB | web download, 2026-08-20 | human, mouse | 3,392 / 2,026 | Shao et al., 2021 |
| connectomeDB2025 | `all_species.csv`, 2026-08-20 | 14 vertebrates | 44,433 in total | Liu et al., 2025 |
| NicheNet v2 LR network | nichenetr 2.2.1.1 / multinichenetr 2.1.0, 2026-08-21 | human, mouse | 5,074 / 5,575 | Browaeys et al., 2023 |
| NeuronChatDB | NeuronChat 1.0.0 (36e06b2), 2026-08-21 | human, mouse | 468 / 452 | Zhao et al., 2023 |
| FlyPhoneDB2 | v2.5, 2026-08-21 | *Drosophila* | 1,798 | FlyPhoneDB2 (2025) |
| *C. elegans* LR pairs | Armingol et al. supplementary, 2026-08-21 | *C. elegans* | 245 | Armingol et al., 2022 |
| PlantPhoneDB | Xu et al. LR tables (2022-07-26) | *Arabidopsis* | 3,514 | Xu et al., 2022 |
| PlantCellChatDB | web download, 2026-08-24 | *Arabidopsis*, maize, rice, soybean, tomato | 5,481 / 2,316 / 1,864 / 1,457 / 1,069 | PlantCellChatDB |

The MD5 checksum, size and date of every raw file are listed in `source_provenance.csv`. Pinned package versions are in `package_versions.csv`.

---

## File structure

```
CellHubDB.rds                          # Named list: one merged CellChat DB per species
CellHubDB_species_summary.csv          # Interactions, sources, multi-source %, complexes per species
CellHubDB_source_contribution.csv      # Supported / exclusive interactions per source and species
databases/
├── <species>/
│   ├── CellHubDB_<species>_<source>.rds   # One file per source
│   └── CellHubDB_<species>_merged.rds     # All sources merged and deduplicated
├── removed_uniprot_entries.csv            # Audit: unresolved UniProt entries (4)
└── removed_invalid_genesymbols.csv        # Audit: symbols absent from the reference gene table (185)
```

---

## Data format

Each database is a CellChat-compatible list:

```r
db <- readRDS("databases/mouse/CellHubDB_mouse_merged.rds")
db$interaction  # LR pairs with annotations and provenance
db$complex      # Multi-subunit complex definitions (subunit_1, subunit_2, ...)
db$cofactor     # Agonists, antagonists and co-receptors (CellChatDB-derived)
db$geneInfo     # Reference gene table (human, mouse, zebrafish)
```

### Interaction table columns

| Column | Description |
|---|---|
| `interaction_name` | Primary key `LIGAND_RECEPTOR` |
| `interaction_name_2` | Readable form, e.g. `Tgfb1 - (Tgfbr1+Tgfbr2)` |
| `ligand`, `receptor` | Gene symbol or complex name (species nomenclature) |
| `pathway_name` | Signaling pathway when provided by the source, otherwise NA |
| `interaction_type` | `Secreted Signaling`, `Cell-Cell Contact`, `ECM-Receptor`, `Non-protein Signaling` or NA |
| `annotation` | Functional annotation |
| `evidence` | PMIDs, KEGG identifiers or DOIs |
| `source_db` | Source whose record was retained |
| `sources` | All sources containing the pair (comma-separated) |

### Gene identifier conventions

| Species | Convention | Example |
|---|---|---|
| human, chimp, macaque, marmoset, cow, dog, horse, pig, sheep, chicken | HGNC-style uppercase | `TGFB1` |
| mouse, rat | MGI/RGD (first letter uppercase) | `Tgfb1` |
| zebrafish | ZFIN lowercase | `tgfb1` |
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
raw_databases/                          raw files (RDS, CSV, TSV, XLSX)
      │  collect_source_provenance.R    → source_provenance.csv, package_versions.csv
      ▼
database_curation/convert_databases/    one converter per source
      ▼
converted/                              92,032 source-level entries
      ▼
database_curation/curate_CellHub_DB.R   cleaning, nomenclature, validation, merge
      ▼
databases/<species>/                    91,199 curated entries → 73,589 merged
      ▼
database_curation/build_cellhubdb.R     → CellHubDB.rds + summary CSVs
```

### Curation steps

1. **Whitespace and empty values.** Whitespace is trimmed and empty strings become NA. Entries without a ligand or receptor are removed (1).
2. **Pathways and annotations.** Source names used as pathways, and generic annotations, become NA.
3. **Interaction types.** Variant spellings are mapped to the four canonical types.
4. **Unresolved UniProt accessions.** Tokens that are neither known symbols nor complex names are removed and audited (4).
5. **Nomenclature.** Symbols are normalized per species; complex names keep their exact spelling.
6. **Symbol validation (human and mouse only).** Symbols absent from the genome-wide reference gene table are removed and audited (185).
7. **Identifier rebuild.** Interaction names are rebuilt and case-induced duplicates removed (643).
8. **Evidence.** References are normalized.
9. **Merge.** Pairs are merged on a canonical key that ignores case and partner order, using the source priority below. Metadata from lower-priority sources fill the gaps in the retained record, and a `sources` provenance column is added.
10. **Complex tables.** Curated definitions are kept. Heteromeric partners without a definition are defined from their underscore-separated subunits.

### Source priority (merge)

```
CellChatDB_v2 > connectomeDB2025 > CellPhoneDB > NeuronChatDB >
FlyPhoneDB2 = cell2cell_Armingol = PlantPhoneDB = PlantCellChatDB > CellTalkDB > NicheNet v2
```

### Quality control

| Script | Checks |
|---|---|
| `verify_integrity.R` | Traceability from converted to per-source to merged files (exhaustive), multi-source consistency, row names |
| `analyze_databases.R` | Per-species annotation coverage, casing, complex/cofactor/geneInfo tables |
| `cellhubdb_stats.R` | Statistics report (`CellHubDB_Statistics_Report.md`) |
| `upset_plot.R` | Source overlap figures |

---

## Known limitations

- **Unresolved zebrafish complexes.** 9 named zebrafish complexes (`inhibin a`, `activin ab`, `il23 complex`, `il27 complex`, `npnt complex 1/2`, `il11r complex 1/2`, `gp complex`) lack subunit definitions in the CellChatDB zebrafish release. CellChat ignores them.
- **Pathway annotations** are available only for interactions from CellChatDB, PlantCellChatDB, PlantPhoneDB, FlyPhoneDB2 and the *C. elegans* list. connectomeDB2025-only species have no pathway labels.
- **Orthology-inferred entries.** Several resources (connectomeDB2025 non-human species, PlantPhoneDB, part of PlantCellChatDB) include interactions inferred by orthology. See each source for its evidence levels.
- **Reference gene table.** Symbol validation uses the CellChatDB gene table. Symbols that are newer than that table (for example `CRIPTO`, `ACKR5`) or retired aliases (`Il1f6`) are excluded and listed in the audit file.
- **Reciprocal pairs.** The merge key ignores partner order, so a pair listed in both directions within a source is kept once.

---

## Usage

### In Cell-Hub

1. Select the species.
2. Load the merged database, or tick individual sources.
3. Run CellChat. When several sources are selected, they are merged with the same rules on the fly.

### In R

```r
library(CellChat)
cellhubdb <- readRDS("CellHubDB.rds")
attr(cellhubdb, "version")              # "2.0"
db <- cellhubdb[["mouse"]]              # merged mouse database
cellchat <- createCellChat(seurat_obj, group.by = "cell_type")
cellchat@DB <- db
# Restrict to interactions supported by at least two sources
db_hc <- db
db_hc$interaction <- db$interaction[grepl(",", db$interaction$sources), ]
cellchat@DB <- db_hc
```

---

## Changes from v1.0

v2.0 rebuilds the database from versioned raw sources and fixes defects present in v1.0:

- **Complex names.** Case normalization altered complex names (for example `TGFbR1_R2` became `Tgfbr1_R2`), and heteromeric partners from non-CellChat sources had no complex definition. As a result, CellChat silently ignored 1,758 of 9,302 mouse and 1,964 of 10,950 human v1.0 interactions, including canonical TGF-β receptor signaling.
- **UniProt filter.** An overly broad UniProt filter removed 515 valid interactions (for example GABA-A receptor complexes, S100A8/A9 and P2RY12).
- **Zebrafish validation.** Zebrafish entries were validated against a gene table that is not genome-wide, which removed 1,425 valid interactions.
- **Identifier casing.** *Xenopus* and crop identifiers were uppercased, so they did not match count matrices.
- **New content.** CellChatDB zebrafish and PlantCellChatDB are now included, adding maize, rice, soybean and tomato.

---

## Citation

> Macaux G. et al. (2026). *bioRxiv* doi:10.64898/2026.09.23.753685.

Please also cite the original source databases used in your analysis.

## License

The curation pipeline, conversion scripts and quality-control code are released under the GNU General Public License v3.0 (GPL-3.0), consistent with Cell-Hub. The interaction data are derived from the source databases listed above, each of which retains its original license and terms of use; users redistributing derived data should comply with the terms of every source they use.
