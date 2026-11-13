# Step 12: inventory of the per-interaction evidence kept from each source (to classify every pair as curated, orthology-inferred or predicted)
# Reads the per-source files of the app databases; for every evidence-like column: top values and keyword counts
# Usage: Rscript revision/12_evidence_inventory.R [app_databases_dir]
args <- commandArgs(trailingOnly = TRUE)
db_dir <- if (length(args) >= 1) args[1] else "data/databases"
out_dir <- "revision/results/12_evidence"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# Columns that are never evidence
skip_cols <- c("interaction_name", "interaction_name_2", "ligand", "receptor", "agonist", "antagonist", "co_A_receptor", "co_I_receptor", "sources", "n_sources",
               "n_curated_sources", "predicted_only", "is_complex", "ligand_role_sources", "receptor_role_sources", "role_corroborated", "confidence_tier", "pathway_origin",
               "pathway_name", "interaction_type", "ligand.symbol", "receptor.symbol", "ligand.keyword", "receptor.keyword", "version")
# Keywords that separate literature / experiment from prediction
keywords <- c(PMID = "PMID|pubmed|^[0-9]{6,9}$", literature = "literature|curated|manual", experimental = "experiment|IntAct|BioGRID|MINT|HPRD|DIP|interactome|plant\\.MAP|MIST",
              STRING = "STRING", ortholog = "ortholog|orthology|DIOPT|inferred|homolog", predicted = "predict|AlphaFold|ipTM|LIS|putative|computational",
              rank_high = "^high$|high rank|high confidence", rank_moderate = "^moderate$|moderate", rank_low = "^low$|low rank|low confidence",
              Direct = "^Direct$|:Direct", Inferred = "^Inferred$|:Inferred", KEGG = "KEGG", DOI = "doi|10\\.[0-9]{4}")

files <- list.files(db_dir, pattern = "^CellHubDB_.*\\.rds$", recursive = TRUE, full.names = TRUE)
files <- files[!grepl("_merged\\.rds$", files)]
top_rows <- list(); kw_rows <- list()
for (f in files) {
  it <- readRDS(f)$interaction
  sp <- basename(dirname(f)); src <- sub(paste0("^CellHubDB_", sp, "_(.*)\\.rds$"), "\\1", basename(f))
  for (cl in setdiff(colnames(it), skip_cols)) {
    v <- as.character(it[[cl]]); v[is.na(v) | v == ""] <- "<empty>"
    # Multi-valued cells (";" or "," separated) counted token by token
    tokens <- trimws(unlist(strsplit(v, "[;,|]")))
    tab <- sort(table(tokens), decreasing = TRUE)
    top_rows[[length(top_rows) + 1]] <- data.frame(species = sp, source = src, column = cl, n_rows = nrow(it), n_distinct = length(tab),
                                                   value = substr(names(tab)[seq_len(min(15, length(tab)))], 1, 80), count = as.integer(tab[seq_len(min(15, length(tab)))]), stringsAsFactors = FALSE)
    kw <- vapply(keywords, function(p) sum(grepl(p, v, ignore.case = TRUE)), integer(1))
    kw_rows[[length(kw_rows) + 1]] <- data.frame(species = sp, source = src, column = cl, n_rows = nrow(it), n_empty = sum(v == "<empty>"), t(kw), stringsAsFactors = FALSE)
  }
}
top_df <- do.call(rbind, top_rows); kw_df <- do.call(rbind, kw_rows)
write.csv(top_df, file.path(out_dir, "evidence_top_values.csv"), row.names = FALSE)
write.csv(kw_df, file.path(out_dir, "evidence_keywords.csv"), row.names = FALSE)

# Console: keyword counts for the sources whose curation status is in question, then their top values
focus <- c("PlantPhoneDB", "PlantCellChatDB", "FlyPhoneDB2", "cell2cell_Armingol", "connectomeDB2025", "CellTalkDB")
kw_show <- kw_df[kw_df$source %in% focus & rowSums(kw_df[, names(keywords)]) > 0, ]
kw_show <- kw_show[!duplicated(paste(kw_show$source, kw_show$column)) | kw_show$source != "connectomeDB2025", ]
print(kw_show, row.names = FALSE)
for (s in setdiff(focus, "connectomeDB2025")) {
  x <- top_df[top_df$source == s & !duplicated(paste(top_df$source, top_df$column, top_df$value)), ]
  x <- x[x$species == x$species[1], ]
  message("\n== ", s, " (", x$species[1], ")"); print(x[, c("column", "n_distinct", "value", "count")], row.names = FALSE)
}
message("\nFull tables: ", out_dir)
