# Step 10: diagnosis of the release gate warnings (zebrafish cofactors / complexes / mouse genes, arabidopsis cofactor, malformed mouse symbol)
# Tells for each problem whether it comes from the original source database (CellChat's own DB) or from the CellHubDB build
# Usage: Rscript review/10_diagnose_gate_warnings.R [app_databases_dir]
suppressPackageStartupMessages(library(CellChat))
args <- commandArgs(trailingOnly = TRUE)
db_dir <- if (length(args) >= 1) args[1] else "app_databases_v2.0.1"
out_dir <- "review/10_gate_warnings"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cof_cols_all <- c("agonist", "antagonist", "co_A_receptor", "co_I_receptor")

# Cofactor identifiers referenced by interactions but missing from the cofactor table
undefined_cofactors <- function(db) {
  cols <- intersect(cof_cols_all, colnames(db$interaction))
  ids <- setdiff(unique(unlist(db$interaction[, cols, drop = FALSE])), c("", NA))
  return(setdiff(ids, rownames(db$cofactor)))
}

# Rows of a table where at least one cell matches a pattern
rows_matching <- function(tab, pattern) {
  if (is.null(tab) || nrow(tab) == 0) return(tab)
  keep <- apply(tab, 1, function(r) any(grepl(pattern, r)))
  return(tab[keep, , drop = FALSE])
}

# Interaction partners (ligand or receptor) of a database
partners <- function(db) unique(c(db$interaction$ligand, db$interaction$receptor))

# Interactions using a given partner, with provenance columns when present
interactions_with <- function(db, genes) {
  it <- db$interaction
  hit <- it$ligand %in% genes | it$receptor %in% genes
  cols <- intersect(c("interaction_name", "ligand", "receptor", "pathway_name", "source_db", "sources"), colnames(it))
  return(it[hit, cols, drop = FALSE])
}

summary_rows <- list()
add <- function(species, check, n_cellhub, n_source, detail) {
  summary_rows[[length(summary_rows) + 1]] <<- data.frame(species = species, check = check, n_cellhub = n_cellhub, n_in_cellchat_source = n_source,
                                                          detail = paste(head(detail, 15), collapse = ", "), stringsAsFactors = FALSE)
}

############## Zebrafish ##############
z <- readRDS(file.path(db_dir, "zebrafish", "CellHubDB_zebrafish_merged.rds"))
ref_z <- CellChat::CellChatDB.zebrafish

# 1. Undefined cofactor references: defined in CellChat's own zebrafish DB or not
undef_z <- undefined_cofactors(z)
add("zebrafish", "undefined cofactor references", length(undef_z), sum(undef_z %in% rownames(ref_z$cofactor)), undef_z)
message("Zebrafish undefined cofactors: ", length(undef_z), " | defined in CellChatDB.zebrafish: ", sum(undef_z %in% rownames(ref_z$cofactor)))

# 2. Mouse-case genes (Abcd1 style) in the cofactor table, here and in CellChat's own zebrafish DB
cof_mouse <- rows_matching(z$cofactor, "^[A-Z][a-z]"); cof_mouse_ref <- rows_matching(ref_z$cofactor, "^[A-Z][a-z]")
add("zebrafish", "cofactor rows with mouse-case genes", nrow(cof_mouse), nrow(cof_mouse_ref), rownames(cof_mouse))
write.csv(cbind(cofactor_id = rownames(cof_mouse), cof_mouse), file.path(out_dir, "zebrafish_cofactor_mouse_case.csv"), row.names = FALSE)

# 3. Purely numeric "subunits" in the complex table
cpx_num <- rows_matching(z$complex, "^[0-9]+$"); cpx_num_ref <- rows_matching(ref_z$complex, "^[0-9]+$")
add("zebrafish", "complex rows with numeric subunits", nrow(cpx_num), nrow(cpx_num_ref), rownames(cpx_num))
write.csv(cbind(complex_id = rownames(cpx_num), cpx_num), file.path(out_dir, "zebrafish_complex_numeric_subunits.csv"), row.names = FALSE)

# 4. Mouse MHC genes (h2-*) used as partners
p_z <- partners(z); p_ref <- partners(ref_z)
h2 <- grep("^h2-", p_z, value = TRUE, ignore.case = TRUE); h2_ref <- grep("^h2-", p_ref, value = TRUE, ignore.case = TRUE)
add("zebrafish", "mouse MHC (h2-*) partners", length(h2), length(h2_ref), h2)

# 5. Named complexes (with spaces) used as partners but not defined in the complex table
undef_cpx <- p_z[grepl(" ", p_z) & !p_z %in% rownames(z$complex)]
undef_cpx_ref <- p_ref[grepl(" ", p_ref) & !p_ref %in% rownames(ref_z$complex)]
add("zebrafish", "undefined named complexes (with spaces)", length(undef_cpx), length(undef_cpx_ref), undef_cpx)

# Interactions affected by 4 and 5, with their source, to decide filter vs fix
write.csv(interactions_with(z, c(h2, undef_cpx)), file.path(out_dir, "zebrafish_interactions_h2_undefined_complexes.csv"), row.names = FALSE)

############## Arabidopsis ##############
a <- readRDS(file.path(db_dir, "arabidopsis", "CellHubDB_arabidopsis_merged.rds"))
undef_a <- undefined_cofactors(a)
add("arabidopsis", "undefined cofactor references", length(undef_a), NA, undef_a)
cols_a <- intersect(cof_cols_all, colnames(a$interaction))
hit_a <- apply(a$interaction[, cols_a, drop = FALSE], 1, function(r) any(r %in% undef_a))
write.csv(a$interaction[hit_a, intersect(c("interaction_name", "ligand", "receptor", cols_a, "source_db", "sources"), colnames(a$interaction)), drop = FALSE],
          file.path(out_dir, "arabidopsis_undefined_cofactor_interactions.csv"), row.names = FALSE)

############## Mouse ##############
m <- readRDS(file.path(db_dir, "mouse", "CellHubDB_mouse_merged.rds"))
# Partners with characters that are never part of a mouse gene symbol or complex name
bad_m <- partners(m)[grepl("[^A-Za-z0-9_.]", partners(m)) & !partners(m) %in% rownames(m$complex)]
add("mouse", "malformed partner symbols", length(bad_m), NA, bad_m)
write.csv(interactions_with(m, bad_m), file.path(out_dir, "mouse_malformed_symbols.csv"), row.names = FALSE)

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(out_dir, "gate_warnings_summary.csv"), row.names = FALSE)
print(summary_df, row.names = FALSE)
message("n_in_cellchat_source > 0 means the problem already exists in CellChat's own database (source defect); 0 means it comes from the CellHubDB build")
message("Detailed tables written to ", out_dir)
