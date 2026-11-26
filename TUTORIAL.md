# CellHub Database tutorial

This tutorial shows how to use CellHub Database (CellHubDB) release 2.1.2 in four ways:
- in R with CellChat;
- in Python with pandas or LIANA+;
- through the Cell-Hub graphical interface;
- by rebuilding the release from the original sources.

Release 2.1.2 is archived on Zenodo ([doi:10.5281/zenodo.24283592](https://doi.org/10.5281/zenodo.24283592)) and attached to the GitHub release [v2.1.2](https://github.com/GaspardMacaux/CellHub-Database/releases/tag/v2.1.2). It contains 76,582 ligand–receptor interactions for 21 species, merged from 10 sources.

| File | Content | MD5 |
|---|---|---|
| `CellHubDB_2.1.2.rds` | All 21 species, one CellChat-format database per species | c08462e1beeb9ae9013a87a86fe106d1 |
| `CellHubDB_2.1.2_tsv.zip` | The same content as tab-separated tables, plus a LIANA+ resource per species | f9aef2d43d43549d25e376d553c7d8aa |
| `CellHubDB_2.1.2_app_databases.zip` | Per-species and per-source files used by the Cell-Hub application | |

---

## 1. What is in a species database

`CellHubDB_2.1.2.rds` is a list with one element per species (`human`, `mouse`, `zebrafish`, `drosophila`, `arabidopsis`, ...). Each element holds the four tables of a CellChat database:

- `interaction`: one row per ligand–receptor interaction. It has the CellChat columns (`interaction_name`, `pathway_name`, `ligand`, `receptor`, `agonist`, `antagonist`, `co_A_receptor`, `co_I_receptor`, `evidence`, `annotation`, `interaction_name_2`) and the CellHubDB columns below.
- `complex`: multi-subunit ligands and receptors with their subunits. It is filled for human, mouse and zebrafish only; the sources of the other species list receptors subunit by subunit.
- `cofactor`: agonists, antagonists and co-receptors.
- `geneInfo`: the gene reference table, genome-wide (Ensembl 112 or Ensembl Plants 59).

CellHubDB columns of the `interaction` table:

| Column | Meaning |
|---|---|
| `sources` | Sources supporting the interaction, separated by ", " |
| `n_sources` | Number of supporting sources |
| `n_curated_sources` | Number of sources supporting the interaction with curated evidence |
| `confidence_tier` | `high` (2 or more curated sources), `medium` (1 curated source), `low` (no curated source: orthology-inferred or computationally predicted only) |
| `predicted_only` | TRUE when no source gives curated evidence (same as tier `low`) |
| `role_corroborated` | TRUE when the ligand gene is used as a ligand, and the receptor gene as a receptor, by 2 or more curated sources (human and mouse; NA for the other species) |
| `is_complex` | TRUE when the ligand or the receptor is a multi-subunit complex |

A source counts as curated for a pair unless it flags the pair as inferred or predicted: connectomeDB2025 "Inferred" pairs, CellChatDB zebrafish pairs projected from mouse, FlyPhoneDB2 pairs below high confidence, PlantPhoneDB orthology or STRING pairs, PlantCellChatDB pairs without experimental validation, and all NicheNet v2 pairs.

---

## 2. In R with CellChat

### 2.1 Load a species

```r
library(CellChat)
db_all <- readRDS("CellHubDB_2.1.2.rds")
names(db_all)                                   # the 21 species
db <- db_all[["human"]]
table(db$interaction$confidence_tier)           # high / medium / low
```

### 2.2 Choose the interactions to keep

The tiers grade the evidence that the sources give for a pair. They do not say whether the pair is active in your tissue.

```r
# Subset a CellChat database to some interactions; complexes and cofactors are kept whole, as CellChat looks them up by name
subset_db <- function(db, keep) { db$interaction <- db$interaction[keep, , drop = FALSE]; db }
# Recommended set for human and mouse: high tier + medium tier with corroborated roles
db_rec <- subset_db(db, db$interaction$confidence_tier == "high" | (db$interaction$confidence_tier == "medium" & db$interaction$role_corroborated %in% TRUE))
# Curated evidence only (no predicted or orthology-inferred pair)
db_cur <- subset_db(db, !db$interaction$predicted_only)
# Direct (non-complex) interactions only
db_direct <- subset_db(db_rec, !db_rec$interaction$is_complex)
nrow(db_rec$interaction); nrow(db_cur$interaction); nrow(db_direct$interaction)
```

Recommendations:
- **Human and mouse**: the recommended set (high + medium with corroborated roles) favours specificity. High + medium (`db_cur`) is broader.
- **Other species**: nearly all pairs are orthology projections, so most of them are in the low tier. Keep all tiers and read the results as hypotheses.
- **Low tier**: it includes intracellular, protease–substrate and same-cell pairs, which CellChat can score because it uses expression alone. Use it only to explore.

### 2.3 Find the interactions of a list of genes

```r
# Interactions whose ligand or receptor, or one of their complex subunits, is in a gene list (case-insensitive)
genes <- c("EGFR", "NOTCH1", "DLL1", "VCAM1")
partner_genes <- function(p) if (p %in% rownames(db$complex)) unique(c(p, unlist(db$complex[p, ]))) else c(p, strsplit(p, "_")[[1]])
hit <- vapply(seq_len(nrow(db$interaction)), function(i) any(toupper(genes) %in% toupper(c(partner_genes(db$interaction$ligand[i]), partner_genes(db$interaction$receptor[i])))), logical(1))
db$interaction[hit, c("interaction_name", "pathway_name", "confidence_tier", "role_corroborated", "sources")]
```

### 2.4 Run CellChat with a species database

The example below uses the human blood cells of the demonstration data (`pbmc_1k_v3_filtered_feature_bc_matrix.h5`, [doi:10.5281/zenodo.24198709](https://doi.org/10.5281/zenodo.24198709)), annotated in Seurat. Any Seurat object with a cell-type column works the same way.

```r
library(Seurat)
seu <- readRDS("pbmc_annotated.rds")              # Seurat object with a "cell_type" column
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])          # Seurat v5: join the layers before CellChat
seu$cell_type <- droplevels(factor(seu$cell_type))
cellchat <- createCellChat(object = seu, group.by = "cell_type", assay = "RNA")
cellchat@DB <- db_rec                             # CellHubDB instead of CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
calls <- subsetCommunication(cellchat, thresh = 0.05)
# Evidence of each significant call, from the database
calls <- merge(calls, db_rec$interaction[, c("interaction_name", "confidence_tier", "role_corroborated", "sources")], by = "interaction_name")
head(calls[order(calls$prob, decreasing = TRUE), ])
```

The same code runs for every species of the release, including *Drosophila*, *C. elegans* and the five plants: each file passes a complete CellChat run before release. For species without a complex table, CellChat scores each subunit of a heteromeric receptor on its own. It cannot require all subunits to be co-expressed, so interpret receptors such as integrins with caution.

In your methods, report the release (2.1.2), the species and the tiers you used.

---

## 3. In Python (pandas, LIANA+)

Unzip `CellHubDB_2.1.2_tsv.zip`. Each species folder holds the four CellChat tables (`interaction.tsv`, `complex.tsv`, `cofactor.tsv`, `geneInfo.tsv`) and `resource_liana.tsv`. The LIANA+ file has one row per ligand–receptor gene pair, with complex subunits joined by "_" (the LIANA+ notation) and the CellHubDB annotation. `MD5SUMS` checks every file.

```python
import pandas as pd
res = pd.read_csv("tsv/human/resource_liana.tsv", sep="\t")
# Recommended set for human and mouse: high tier + medium tier with corroborated roles
rec = res[(res.confidence_tier == "high") | ((res.confidence_tier == "medium") & (res.role_corroborated == True))]
print(len(res), len(rec))

import liana as li
# adata: AnnData with log-normalized expression and a "cell_type" column
li.mt.rank_aggregate(adata, groupby="cell_type", resource=rec[["ligand", "receptor"]], use_raw=False)
adata.uns["liana_res"].head()
```

When the same gene pair appears in several CellChat interactions (a complex written differently, agonist or cofactor variants), `resource_liana.tsv` keeps one row for it: the best-supported one. This is why it has fewer rows than `interaction.tsv` for human and mouse.

---

## 4. In Cell-Hub (no code)

CellHubDB 2.1.2 is built into the Cell-Hub Docker image (`gaspardmacaux/cell-hub:2.0.0`); see the [Cell-Hub tutorial](https://github.com/GaspardMacaux/Cell-Hub/blob/main/TUTORIAL.md), section 3.

1. **CellChat → Load dataset**: click **Use current object** (an object analysed in Cell-Hub) or load a Seurat object (`.rds`).
2. **Select Species**, then **All databases (merged)** or **Select individual databases**. With selected sources, the tiers are recomputed on the sources kept.
3. **Confidence tiers to keep**. The default is set per species: high + medium with confirmed roles for human and mouse, all tiers for the other species. Click **Load Database**. **Database Information** gives the number of interactions per tier.
4. **Genes of interest**: paste a list of genes, click **Add list**, then **Show interactions in the database**. This lists their interactions without running CellChat, like section 2.3.
5. **Ligand-Receptor → Create & Analyze**, then the bubble and chord plots. Both can be restricted to a gene list.

---

## 5. Rebuild the release from the sources

The raw sources are not stored in the repository. A download script gets them from pinned origins and checks each file against its MD5 (`config/source_retrieval.tsv`). From an empty clone:

```bash
git clone https://github.com/GaspardMacaux/CellHub-Database.git && cd CellHub-Database
Rscript pipeline/00_download_sources.R --install     # installs the R packages at pinned versions and downloads the sources
bash run_pipeline.sh --convert                       # conversion, harmonization, curation, consolidation, CellChat checks
md5sum release/*/CellHubDB*.rds
```

Three sources have no versioned download: CellTalkDB, connectomeDB2025 and PlantCellChatDB. The CellTalkDB and connectomeDB2025 files are deposited on Zenodo. PlantCellChatDB states no licence and is not redistributed: download its files as explained by `00_download_sources.R`, which then checks their checksums.

The revision analyses are in `revision/`:
- `04_benchmark_cellchat.R`: CellChat on the muscle dataset with each database, and permuted labels;
- `04b_musc_shortlist.R`: new interactions of muscle stem cells, for the literature check;
- `04c_merged_calls_by_source.R`: calls of the merged database by evidence;
- `13_conservation_direct_inferred.R`: conservation between species, curated and inferred pairs separated;
- `15_placeholders_2.1.2.R`: the 2.1.2 patch.

`pipeline/15_export_tsv.R` writes the tab-separated tables of section 3.

---

## 6. Licences and citation

The pipeline and the release files are GPL-3.0. Each source keeps its own licence (Table 3 of the article and `config/source_retrieval.tsv`); cite the sources you rely on, in addition to CellHubDB.

Questions and problems: https://github.com/GaspardMacaux/CellHub-Database/issues
