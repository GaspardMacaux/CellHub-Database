# Step 04 (Reviewer 2, major 2c + 3): CellChat benchmark on the MuSC/myofiber dataset:
# default CellChatDB, each source, merged, consensus sets, label-permutation negative control
source("revision/revision_helpers.R")
suppressPackageStartupMessages({ library(Seurat); library(CellChat) })
seurat_path <- "/media/gaspard/Gasp/Article Cell-Hub/WT_2M.rds"; label_col <- "cluster_names"; focus_group <- "MuSC"; species <- "mouse"
# DoubletFinder calls: nuclei labelled "Doublet" are excluded before the benchmark (set doublet_col <- NULL to keep them)
doublet_col <- "DF.classifications_0.25_0.005_1023"
db_dir <- file.path("release/2.1/app_databases", species); annotated_path <- "release/2.1/CellHubDB.rds"; out_dir <- "revision/results/04_benchmark_v2.1"
max_cells_per_group <- 500; n_perm <- 5; nboot <- 100; n_workers <- 6
cellchat_seed <- 1  # same seed.use as the Cell-Hub CellChat module that produced Figure 3
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("CellChat version: ", as.character(packageVersion("CellChat")))
future::plan("multisession", workers = n_workers); options(future.globals.maxSize = 8 * 1024^3)

# Expression data: joined layers, clean factor levels, stratified subsampling with seed 42
obj <- readRDS(seurat_path)
DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
if (!is.null(doublet_col)) {
  n_before <- ncol(obj); obj <- subset(obj, cells = colnames(obj)[obj@meta.data[[doublet_col]] != "Doublet"])
  message("Doublets removed: ", n_before - ncol(obj), " of ", n_before, " nuclei")
}
obj@meta.data[] <- lapply(obj@meta.data, function(x) if (is.factor(x)) droplevels(x) else x)
print(table(obj@meta.data[[label_col]]))
set.seed(42)
cells <- unlist(lapply(split(colnames(obj), as.character(obj@meta.data[[label_col]])), function(cl) if (length(cl) > max_cells_per_group) sample(cl, max_cells_per_group) else cl), use.names = FALSE)
data_norm <- LayerData(obj, assay = "RNA", layer = "data")[, cells]
meta <- data.frame(labels = droplevels(factor(obj@meta.data[cells, label_col])), row.names = cells)
rm(obj); gc(verbose = FALSE)

# Database variants
schema <- get_cellchat_schema()
merged <- readRDS(annotated_path)[[species]]
db_variants <- list(CellChatDB_default = CellChat::CellChatDB.mouse)
src_files <- list.files(db_dir, pattern = paste0("^CellHubDB_", species, "_.*\\.rds$"), full.names = TRUE)
src_files <- src_files[!grepl("_merged\\.rds$", src_files)]
# Per-source files: schema fixed and, when their own geneInfo is empty, the corrected species geneInfo (same reference for every variant)
for (f in src_files) {
  x <- standardize_cellchat_schema(readRDS(f), schema)
  if (nrow(x$geneInfo) == 0) x$geneInfo <- merged$geneInfo
  db_variants[[paste0("source_", sub(paste0("^CellHubDB_", species, "_(.*)\\.rds$"), "\\1", basename(f)))]] <- x
}
# Merged variants: full set, without NicheNet-only pairs, Reviewer 2 consensus (>= 2 sources, NicheNet included),
# >= 2 curated sources (tier high), and curated pairs whose roles are corroborated (tier high + medium with role_corroborated)
db_variants$merged_all <- merged
db_variants$merged_no_predicted_only <- subset_cellchat_db(merged, !merged$interaction$predicted_only)
db_variants$consensus_2plus <- subset_cellchat_db(merged, merged$interaction$n_sources >= 2)
db_variants$tier_high <- subset_cellchat_db(merged, merged$interaction$confidence_tier == "high")
db_variants$curated_role_corroborated <- subset_cellchat_db(merged, merged$interaction$confidence_tier == "high" |
                                                              (merged$interaction$confidence_tier == "medium" & merged$interaction$role_corroborated %in% TRUE))

# Real labels and permuted labels (group sizes preserved) for each variant
results <- list(); summary_rows <- list()
for (v in names(db_variants)) {
  message("Variant: ", v, " (", nrow(db_variants[[v]]$interaction), " interactions)")
  real <- tryCatch(run_cellchat_db(data_norm, meta, db_variants[[v]], seed = cellchat_seed, nboot = nboot), error = function(e) { message("  failed: ", conditionMessage(e)[1]); NULL })
  perm <- lapply(seq_len(n_perm), function(i) {
    set.seed(42 + i); meta_p <- meta; meta_p$labels <- sample(meta$labels)
    tryCatch(run_cellchat_db(data_norm, meta_p, db_variants[[v]], seed = cellchat_seed, nboot = nboot), error = function(e) NULL)
  })
  results[[v]] <- list(real = real, perm = perm)
  s_real <- summarise_calls(real); s_perm <- vapply(perm, function(p) unname(summarise_calls(p)["n_calls"]), numeric(1))
  summary_rows[[v]] <- data.frame(variant = v, n_db_interactions = nrow(db_variants[[v]]$interaction), n_calls = s_real[["n_calls"]], n_lr = s_real[["n_lr"]],
                                  n_pathways = s_real[["n_pathways"]], perm_calls_mean = mean(s_perm), perm_calls_max = max(s_perm),
                                  perm_to_real_ratio = round(mean(s_perm) / max(s_real[["n_calls"]], 1), 3))
  saveRDS(results, file.path(out_dir, "benchmark_raw_results.rds"))
}
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(out_dir, "benchmark_summary.csv"), row.names = FALSE)
print(summary_df)

# Merged-database calls annotated with provenance; flags for predicted-only and for pairs absent from the default CellChatDB
info_cols <- intersect(c("interaction_name", "sources", "n_sources", "n_curated_sources", "predicted_only", "is_complex", "role_corroborated", "confidence_tier", "pathway_origin"), colnames(merged$interaction))
lr_key <- function(df) tolower(paste(df$ligand, df$receptor, sep = "|"))
annotate_calls <- function(res) {
  if (is.null(res) || nrow(res) == 0) return(NULL)
  out <- merge(res, merged$interaction[, info_cols], by = "interaction_name", all.x = TRUE)
  out$in_default_db <- lr_key(out) %in% lr_key(CellChat::CellChatDB.mouse$interaction)
  out$in_default_result <- if (is.null(results$CellChatDB_default$real)) FALSE else lr_key(out) %in% lr_key(results$CellChatDB_default$real)
  return(out)
}
merged_calls <- annotate_calls(results$merged_all$real)
consensus_calls <- annotate_calls(results$consensus_2plus$real)
if (is.null(merged_calls)) stop("No significant call with the merged database: check the CellChat run before going further")
write.csv(merged_calls, file.path(out_dir, "calls_merged_annotated.csv"), row.names = FALSE)
if (!is.null(consensus_calls)) write.csv(consensus_calls, file.path(out_dir, "calls_consensus2_annotated.csv"), row.names = FALSE)
message("Merged calls: ", nrow(merged_calls), " | predicted-only (NicheNet): ", sum(merged_calls$predicted_only, na.rm = TRUE),
        " | absent from default CellChatDB: ", sum(!merged_calls$in_default_db))
print(table(tier = merged_calls$confidence_tier, in_default_db = merged_calls$in_default_db, useNA = "ifany"))

# Interactions sent or received by the focus population, for manual checking against the literature (Figure 3e)
focus <- merged_calls[as.character(merged_calls$source) == focus_group | as.character(merged_calls$target) == focus_group, ]
focus <- focus[order(focus$in_default_db, -focus$n_sources, focus$pval), ]
write.csv(focus, file.path(out_dir, paste0("calls_", focus_group, "_for_manual_check.csv")), row.names = FALSE)
