# Step 04c (Reviewer 2, major 3): significant calls of the merged database classified by the evidence of their pair
# (curated support / NicheNet v2 only / other non-curated sources). Same data, subsampling, seed and CellChat settings
# as step 04 (merged_all, real labels only, no permutation). Run from the CellHub_DB folder.
source("revision/revision_helpers.R")
suppressPackageStartupMessages({ library(Seurat); library(CellChat) })
seurat_path <- "/media/gaspard/Gasp/Article Cell-Hub/WT_2M.rds"; label_col <- "cluster_names"; species <- "mouse"
doublet_col <- "DF.classifications_0.25_0.005_1023"
annotated_path <- "release/2.1/CellHubDB.rds"; out_dir <- "revision/results/04c_merged_by_source"
max_cells_per_group <- 500; nboot <- 100; n_workers <- 6; cellchat_seed <- 1
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("CellChat version: ", as.character(packageVersion("CellChat")), " | database: ", annotated_path, " (MD5 ", unname(tools::md5sum(annotated_path)), ")")
future::plan("multisession", workers = n_workers); options(future.globals.maxSize = 8 * 1024^3)
# Expression data: identical to step 04 (joined layers, doublets removed, clean factor levels, stratified subsampling with seed 42)
obj <- readRDS(seurat_path)
DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
n_before <- ncol(obj); obj <- subset(obj, cells = colnames(obj)[obj@meta.data[[doublet_col]] != "Doublet"])
message("Doublets removed: ", n_before - ncol(obj), " of ", n_before, " nuclei")
obj@meta.data[] <- lapply(obj@meta.data, function(x) if (is.factor(x)) droplevels(x) else x)
set.seed(42)
cells <- unlist(lapply(split(colnames(obj), as.character(obj@meta.data[[label_col]])), function(cl) if (length(cl) > max_cells_per_group) sample(cl, max_cells_per_group) else cl), use.names = FALSE)
data_norm <- LayerData(obj, assay = "RNA", layer = "data")[, cells]
meta <- data.frame(labels = droplevels(factor(obj@meta.data[cells, label_col])), row.names = cells)
rm(obj); gc(verbose = FALSE)
# Merged database, real labels only
merged <- readRDS(annotated_path)[[species]]
calls <- run_cellchat_db(data_norm, meta, merged, seed = cellchat_seed, nboot = nboot)
message("Merged calls: ", nrow(calls), " (Table 2: 5,397)")
# Evidence class of each call, from the pair annotation of the database
inter <- merged$interaction
key <- if ("interaction_name" %in% colnames(inter)) as.character(inter$interaction_name) else rownames(inter)
idx <- match(as.character(calls$interaction_name), key)
if (anyNA(idx)) stop(sum(is.na(idx)), " calls not found in the merged database")
tier <- tolower(as.character(inter$confidence_tier[idx]))
srcs <- strsplit(as.character(inter$sources[idx]), "\\s*[;,|]\\s*")
nichenet_only <- vapply(srcs, function(s) length(s) == 1 && grepl("NicheNet", s), logical(1))
calls$sources <- vapply(srcs, paste, character(1), collapse = " + ")
calls$confidence_tier <- tier
calls$evidence_class <- ifelse(tier %in% c("high", "medium"), "curated support", ifelse(nichenet_only, "NicheNet v2 only", "no curated support, other sources"))
write.csv(calls, file.path(out_dir, "calls_merged_by_source.csv"), row.names = FALSE)
tab <- table(calls$evidence_class)
print(data.frame(evidence = names(tab), calls = as.integer(tab), percent = round(100 * as.integer(tab) / nrow(calls), 1)))
message("\nSources of the non-curated calls that are not NicheNet-only:")
print(sort(table(calls$sources[calls$evidence_class == "no curated support, other sources"]), decreasing = TRUE))
# Same classification restricted to the calls absent from the default CellChatDB result (the 4,392 added calls)
lr_key <- function(df) tolower(paste(df$ligand, df$receptor, sep = "|"))
default_calls <- run_cellchat_db(data_norm, meta, CellChat::CellChatDB.mouse, seed = cellchat_seed, nboot = nboot)
added <- !paste(lr_key(calls), calls$source, calls$target) %in% paste(lr_key(default_calls), default_calls$source, default_calls$target)
message("\nDefault CellChatDB calls: ", nrow(default_calls), " (Table 2: 1,005) | merged calls absent from them: ", sum(added), " (expected 4,392)")
print(table(added_call_evidence = calls$evidence_class[added]))
