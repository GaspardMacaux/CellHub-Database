# GaspouDB — Multi-Species Curated Ligand-Receptor Interaction Database

**Version:** 2.0  
**Release date:** August 2026  
**Total interactions:** 60,456 (deduplicated, across 17 species)  
**Sources:** 10 curated databases  
**Format:** CellChat-native RDS (ready-to-use in Cell-Hub and CellChat)

---

## Overview

GaspouDB is a consolidated, multi-species ligand-receptor (LR) interaction database designed for cell-cell communication analysis. It integrates 10 independently curated source databases into a unified, standardized format compatible with CellChat and Cell-Hub.

Unlike existing aggregators (OmniPath, LIANA+), GaspouDB provides:

- **CellChat-native format** — load a single `.rds` file, ready for analysis
- **17 species** in a unified schema — from human and mouse to drosophila, C. elegans, and Arabidopsis
- **Per-source traceability** — every interaction tracks its origin database(s)
- **No-code access** — users select species and databases through Cell-Hub's graphical interface
- **Standardized curation** — consistent gene symbol conventions, deduplicated cross-source interactions, normalized evidence references

---

## Species Coverage

| Species | Common name | Sources | Interactions (merged) |
|---|---|---|---|
| `human` | Homo sapiens | CellChatDB, CellPhoneDB, CellTalkDB, connectomeDB2025, MultiNicheNet, NeuronChatDB | 10,775 |
| `mouse` | Mus musculus | CellChatDB, CellTalkDB, connectomeDB2025, MultiNicheNet, NeuronChatDB | 9,194 |
| `arabidopsis` | Arabidopsis thaliana | PlantPhoneDB | 3,394 |
| `chimp` | Pan troglodytes | connectomeDB2025 | 3,190 |
| `macaque` | Macaca mulatta | connectomeDB2025 | 3,185 |
| `rat` | Rattus norvegicus | connectomeDB2025 | 3,103 |
| `cow` | Bos taurus | connectomeDB2025 | 3,050 |
| `horse` | Equus caballus | connectomeDB2025 | 3,037 |
| `dog` | Canis lupus familiaris | connectomeDB2025 | 3,004 |
| `pig` | Sus scrofa | connectomeDB2025 | 3,004 |
| `marmoset` | Callithrix jacchus | connectomeDB2025 | 2,961 |
| `sheep` | Ovis aries | connectomeDB2025 | 2,891 |
| `frog` | Xenopus tropicalis | connectomeDB2025 | 2,601 |
| `zebrafish` | Danio rerio | CellChatDB, connectomeDB2025 | 2,591 |
| `chicken` | Gallus gallus | connectomeDB2025 | 2,441 |
| `drosophila` | Drosophila melanogaster | FlyPhoneDB2 | 1,790 |
| `celegans` | Caenorhabditis elegans | cell2cell (Armingol et al.) | 245 |

---

## Source Databases

| Source | Reference | Species | Interactions | Evidence type |
|---|---|---|---|---|
| **CellChatDB v2** | Jin et al. (2024) | human, mouse, zebrafish | ~3,300/species | KEGG refs, PMIDs |
| **CellPhoneDB v5** | Troulé et al. (2023) | human | 2,911 | PMIDs, UniProt |
| **CellTalkDB** | Shao et al. (2021) | human, mouse | 3,398 H / 2,033 M | PMIDs |
| **connectomeDB2025** | Liu et al. (2025) | 14 vertebrates | ~3,200/species | Direct / Inferred (orthology) |
| **MultiNicheNet** | Browaeys et al. (2023) | human, mouse | 5,098 H / 5,681 M | PMIDs, multi-resource scores |
| **NeuronChatDB** | Zhao et al. (2023) | human, mouse | 190 H / 183 M | Curated (neuro-specific) |
| **FlyPhoneDB2** | (2025) | drosophila | 1,798 | Annotation + literature |
| **cell2cell** | Armingol et al. (2022) | C. elegans | 245 | DOIs, orthology + STRING |
| **PlantPhoneDB** | Xu et al. (2022) | Arabidopsis | 3,514 | STRING, IntAct, BioGRID, literature |
| **MRCLinkdb** | (reserved) | human, mouse | — | Metabolite-receptor (future) |

---

## File Structure

```
databases/
├── <species>/
│   ├── GaspouDB_<species>_<source>.rds    # Individual source databases
│   └── GaspouDB_<species>_merged.rds      # All sources merged and deduplicated
│
├── human/
│   ├── GaspouDB_human_CellChatDB_v2.rds
│   ├── GaspouDB_human_CellPhoneDB.rds
│   ├── GaspouDB_human_CellTalkDB.rds
│   ├── GaspouDB_human_connectomeDB2025.rds
│   ├── GaspouDB_human_MultiNicheNet.rds
│   ├── GaspouDB_human_NeuronChatDB.rds
│   └── GaspouDB_human_merged.rds           # 10,775 interactions
│
├── mouse/
│   ├── GaspouDB_mouse_CellChatDB_v2.rds
│   ├── GaspouDB_mouse_CellTalkDB.rds
│   ├── GaspouDB_mouse_connectomeDB2025.rds
│   ├── GaspouDB_mouse_MultiNicheNet.rds
│   ├── GaspouDB_mouse_NeuronChatDB.rds
│   └── GaspouDB_mouse_merged.rds           # 9,194 interactions
│
├── drosophila/
│   ├── GaspouDB_drosophila_FlyPhoneDB2.rds
│   └── GaspouDB_drosophila_merged.rds      # 1,790 interactions
│
├── celegans/
│   ├── GaspouDB_celegans_cell2cell_Armingol.rds
│   └── GaspouDB_celegans_merged.rds        # 245 interactions
│
├── arabidopsis/
│   ├── GaspouDB_arabidopsis_PlantPhoneDB.rds
│   └── GaspouDB_arabidopsis_merged.rds     # 3,394 interactions
│
└── ... (12 additional vertebrate species)
```

---

## Data Format

Each `.rds` file is a CellChat-compatible R list with four elements:

```r
db <- readRDS("databases/mouse/GaspouDB_mouse_merged.rds")
db$interaction  # data.frame — LR pairs with annotations
db$complex      # data.frame — multi-subunit complex definitions
db$cofactor     # data.frame — cofactor/co-receptor annotations
db$geneInfo     # data.frame — gene symbol reference table
```

### Interaction Table Columns

| Column | Description |
|---|---|
| `interaction_name` | Primary key: `LIGAND_RECEPTOR` (underscore-joined) |
| `interaction_name_2` | Human-readable: `Ligand - Receptor` or `Ligand - (Sub1+Sub2)` |
| `ligand` | Gene symbol of the ligand (species-appropriate casing) |
| `receptor` | Gene symbol of the receptor (or complex: `Sub1_Sub2`) |
| `pathway_name` | Signaling pathway (e.g., TGFb, WNT, CXCL) — NA if unknown |
| `interaction_type` | `Secreted Signaling`, `Cell-Cell Contact`, `ECM-Receptor`, `Non-protein Signaling` |
| `annotation` | Detailed functional annotation |
| `evidence` | Bibliographic references (PMIDs, KEGG, DOIs) |
| `source_db` | Winning source database for this interaction |
| `sources` | All databases containing this interaction (comma-separated) |

### Gene Symbol Conventions

| Species | Convention | Example |
|---|---|---|
| Mouse | MGI — first uppercase | `Tgfb1`, `Bmpr1a` |
| Human | HGNC — all uppercase | `TGFB1`, `BMPR1A` |
| Zebrafish | ZFIN — all lowercase | `tgfb1`, `bmpr1a` |
| Rat | RGD — first uppercase | `Tgfb1`, `Bmpr1a` |
| Drosophila | FlyBase — as published | `dpp`, `Toll`, `sNPF` |
| C. elegans | WormBase — lowercase + hyphen | `daf-7`, `sma-6` |
| Arabidopsis | TAIR — AGI locus IDs | `AT1G01900` |

---

## Curation Pipeline

### Architecture

```
raw_databases/           → Source-specific raw files (RDS, CSV, TSV, Excel)
       │
       ▼
convert_<source>.R       → Per-source conversion to standardized format
       │
       ▼
converted/               → Intermediate files (one RDS per source × species)
       │
       ▼
curate_gaspoudb.R        → Cleaning, normalization, deduplication
       │
       ▼
databases/<species>/     → Final output (individual + merged per species)
```

### Curation Steps

1. **Empty string normalization** — `""` → `NA` across all columns
2. **Pathway name cleaning** — source names used as pathways → `NA`
3. **Annotation cleaning** — generic values like `"Signaling"` → `NA`
4. **Interaction type standardization** — variant spellings normalized to 4 canonical types
5. **UniProt ID removal** — unresolvable accessions flagged and removed (saved to audit file)
6. **Gene symbol casing** — normalized per species convention using geneInfo lookup
7. **Interaction name rebuild** — regenerated from normalized ligand/receptor
8. **Evidence normalization** — PMID formatting standardized, aberrant entries cleaned
9. **Cross-source deduplication** — canonical L-R pairs merged with source priority

### Source Priority (for deduplication)

When the same L-R pair exists in multiple databases, the highest-priority source determines directionality. Metadata (evidence, annotation, pathway) from lower-priority sources is merged into the winning row.

```
CellChatDB_v2 (10) > CellChatDB_v1 (9) > connectomeDB2025 (8) > CellPhoneDB (7) >
NeuronChatDB (6) > SingleCellSignalR (5) > ICELLNET (5) > MRCLinkdb (4) >
cell2cell (3) > PlantPhoneDB (3) > FlyPhoneDB2 (3) > CellTalkDB (2) > MultiNicheNet (1)
```

---

## Usage in Cell-Hub

### Graphical Interface

1. **Select species** from the dropdown menu
2. **Choose databases**: load the pre-merged database (all sources, deduplicated) or select individual sources via checkboxes
3. **Load** — the database is assigned to the CellChat analysis pipeline

### Programmatic Usage (R)

```r
# Load merged database for mouse
db <- readRDS("databases/mouse/GaspouDB_mouse_merged.rds")

# Use with CellChat
cellchat <- createCellChat(seurat_obj, group.by = "cell_type")
cellchat@DB <- db

# Or load a single source
db_cellchat <- readRDS("databases/mouse/GaspouDB_mouse_CellChatDB_v2.rds")
```

---

## Reproducibility

### Pipeline Scripts

| Script | Purpose |
|---|---|
| `extract_all_databases.R` | Extract raw data from R packages + document manual downloads |
| `convert_cellchatdb.R` | CellChatDB v2 → GaspouDB format |
| `convert_cellphonedb.R` | CellPhoneDB v5 → GaspouDB (resolves UniProt IDs and complexes) |
| `convert_celltalkdb.R` | CellTalkDB → GaspouDB |
| `convert_connectomedb2025.R` | connectomeDB2025 → GaspouDB (14 species from single CSV) |
| `convert_multinichenet.R` | NicheNet v2 LR network → GaspouDB |
| `convert_neuronchatdb.R` | NeuronChatDB → GaspouDB |
| `convert_flyphonedb2.R` | FlyPhoneDB2 → GaspouDB |
| `convert_celegans.R` | C. elegans Armingol et al. → GaspouDB |
| `convert_plantphonedb.R` | PlantPhoneDB (Arabidopsis) → GaspouDB |
| `curate_gaspoudb.R` | Main curation pipeline (cleaning + merging) |
| `analyze_converted.R` | Quality control: inspect all converted files |

### Running the Full Pipeline

```bash
cd /path/to/GaspouDB/

# 1. Extract data from R packages
Rscript database_curation/extract_all_databases.R

# 2. Convert each source to standardized format
Rscript database_curation/convert_cellchatdb.R
Rscript database_curation/convert_connectomedb2025.R
Rscript database_curation/convert_cellphonedb.R
Rscript database_curation/convert_celltalkdb.R
Rscript database_curation/convert_multinichenet.R
Rscript database_curation/convert_neuronchatdb.R
Rscript database_curation/convert_flyphonedb2.R
Rscript database_curation/convert_celegans.R
Rscript database_curation/convert_plantphonedb.R

# 3. Quality control
Rscript database_curation/analyze_converted.R

# 4. Run curation pipeline (cleaning + merging)
Rscript database_curation/curate_gaspoudb.R
```

---

## Citation

If you use GaspouDB, please cite:

> Macaux G. et al. (2026). Cell-Hub: a no-code platform for single-cell RNA sequencing analysis with integrated multi-species ligand-receptor databases. *NAR Genomics and Bioinformatics* (submitted).

Additionally, please cite the original source databases used in your analysis (see Source Databases table above).

---

## License

GaspouDB curation pipeline and integration code: MIT License.  
Individual source databases retain their original licenses. See each source's documentation for terms of use.
