# CellHub Database

A consolidated multi-species ligand-receptor interaction resource for cell-cell communication analysis.

CellHub Database integrates **10 independently curated ligand-receptor databases** into **67,584 curated interactions** across **21 species** in CellChat-native format.

## Quick start

```r
# Load the database
db <- readRDS("CellHubDB.rds")

# See available species
names(db)
# [1] "arabidopsis" "celegans"  "chicken"    "chimp"     "cow"       "dog"
# [7] "drosophila"  "frog"      "horse"      "human"     "macaque"   "maize"
# [13] "marmoset"   "mouse"     "pig"        "rat"       "rice"      "sheep"
# [19] "soybean"    "tomato"    "zebrafish"

# Assign to a CellChat object
cellchat@DB <- db[["mouse"]]
```

No additional package installation is required beyond CellChat itself.

## Source databases

| Source | Species | Reference |
|---|---|---|
| CellChatDB v2 | Human, Mouse, Zebrafish | Jin et al., *Nat. Protoc.* 2025 |
| CellPhoneDB v5 | Human | Troulé et al., *Nat. Protoc.* 2025 |
| NicheNet v2 | Human, Mouse | Browaeys et al., *Nat. Methods* 2020 |
| CellTalkDB | Human, Mouse | Shao et al., *Brief. Bioinform.* 2021 |
| connectomeDB2025 | 14 vertebrates | Liu et al., *Nucleic Acids Res.* 2026 |
| NeuronChatDB | Human, Mouse | Zhao et al., *Nat. Commun.* 2023 |
| FlyPhoneDB2 | Drosophila | Qadiri et al., *Comput. Struct. Biotechnol. J.* 2025 |
| cell2cell | C. elegans | Armingol et al., *PLoS Comput. Biol.* 2022 |
| PlantPhoneDB | Arabidopsis | Xu et al., *Plant Biotechnol. J.* 2022 |
| PlantCellChatDB | 5 plant species | Liu et al., *Plant J.* 2026 |

## Species coverage

**Mammals:** Human (6 sources, 10,950 interactions), Mouse (5 sources, 9,302), Rat, Chimp, Macaque, Marmoset, Cow, Dog, Horse, Pig, Sheep (connectomeDB2025)

**Other vertebrates:** Chicken, Frog, Zebrafish (1–2 sources)

**Invertebrates:** Drosophila melanogaster (FlyPhoneDB2), Caenorhabditis elegans (cell2cell)

**Plants:** Arabidopsis thaliana (2 sources, 3,716), Maize, Rice, Soybean, Tomato (PlantCellChatDB)

## Cross-source validation

For human and mouse, 34% and 31% of interactions respectively are confirmed by at least two independent databases.

## RDS file structure

Each species entry is a list with four elements, directly compatible with `cellchat@DB`:

- **interaction** — Ligand-receptor pairs with provenance tracking
- **complex** — Multi-subunit complex definitions
- **cofactor** — Cofactors and co-receptors
- **geneInfo** — Gene symbol reference table

The interaction table includes: `interaction_name`, `ligand`, `receptor`, `pathway_name`, `interaction_type`, `evidence` (PMIDs), `source_db` (winning source), and `provenance` (all contributing databases).

## Standalone usage

```r
db <- readRDS("CellHubDB.rds")
cellchat@DB <- db[["human"]]
```

## Usage in Cell-Hub

CellHub Database is integrated into [Cell-Hub](https://github.com/GaspardMacaux/Cell-Hub), an open-source R/Shiny application for no-code scRNA-seq analysis. Within Cell-Hub, users select a species and choose between loading the merged database or individual sources via a graphical interface.

## Repository structure

```
CellHub-Database/
├── CellHubDB.rds                    # Main database (all 21 species, merged)
├── databases/                       # Per-species RDS files (merged + individual sources)
├── scripts/                         # Complete curation pipeline (10 conversion + QC)
├── geneinfo/                        # Gene symbol reference tables
├── figures/                         # Article figures
├── LICENSE                          # MIT
└── README.md
```

## Reproducing the pipeline

```bash
# 1. Place raw source files in raw_databases/
# 2. Run conversion scripts (one per source)
Rscript scripts/convert_cellchatdb.R
Rscript scripts/convert_cellphonedb.R
# ... (all 10 sources)

# 3. Run curation pipeline
Rscript scripts/curate_CellHub_DB.R

# 4. Verify integrity
Rscript scripts/verify_integrity.R
```

## Citation

If you use CellHub Database, please cite:

> Macaux, G. & Maire, P. CellHub Database: a consolidated multi-species ligand-receptor interaction resource for cell-cell communication analysis. (2026).

## License

MIT License. See [LICENSE](LICENSE).

## Contact

Gaspard Macaux — Institut Cochin, INSERM U1016, Université Paris Cité
