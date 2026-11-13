# Step 01b (diagnostic only, nothing here is released): which missing piece breaks CellChat on each species
# A = schema fixed only; B = A + minimal geneInfo built from the database genes; C = B + NA pathway_name replaced by "Unassigned"
source("review/revision_helpers.R")
suppressPackageStartupMessages(library(CellChat))
future::plan("sequential")
db_path <- "review/01_audit/CellHubDB_schema_fixed.rds"; out_dir <- "review/01_audit"
db <- readRDS(db_path)
message("CellChat version: ", as.character(packageVersion("CellChat")))

# Minimal geneInfo listing only the genes used by the interactions (diagnostic, not a validated reference)
add_minimal_geneinfo <- function(sp_db) {
  if (nrow(sp_db$geneInfo) == 0) sp_db$geneInfo <- data.frame(Symbol = extract_partner_genes(sp_db), stringsAsFactors = FALSE)
  return(sp_db)
}

# Missing pathway labels replaced by a placeholder so that computeCommunProbPathway can aggregate them
fill_missing_pathway <- function(sp_db, placeholder = "Unassigned") {
  pw <- sp_db$interaction$pathway_name
  sp_db$interaction$pathway_name[is.na(pw) | pw == ""] <- placeholder
  return(sp_db)
}

# Pathway-level diagnostics: NA count, pathways with a single interaction, pathway names that are not plain strings
pathway_diag <- function(sp_db, species) {
  pw <- sp_db$interaction$pathway_name
  ok <- pw[!is.na(pw) & pw != ""]
  return(data.frame(species = species, n_pathway_na = sum(is.na(pw) | pw == ""), n_pathways = length(unique(ok)),
                    n_singleton_pathways = sum(table(ok) == 1), pathway_is_factor = is.factor(pw), stringsAsFactors = FALSE))
}

variants <- list(A_schema = db, B_geneinfo = lapply(db, add_minimal_geneinfo), C_geneinfo_pathway = lapply(lapply(db, add_minimal_geneinfo), fill_missing_pathway))
diag_rows <- list()
for (v in names(variants)) {
  for (sp in names(variants[[v]])) {
    message("Diagnostic ", v, ": ", sp)
    r <- run_cellchat_smoke(variants[[v]][[sp]], sp)
    r$variant <- v
    diag_rows[[paste(v, sp)]] <- r
  }
}
diag_df <- do.call(rbind, diag_rows)
diag_df <- diag_df[, c("variant", setdiff(colnames(diag_df), "variant"))]
write.csv(diag_df, file.path(out_dir, "diagnostic_variants.csv"), row.names = FALSE)

pw_df <- do.call(rbind, lapply(names(db), function(sp) pathway_diag(db[[sp]], sp)))
write.csv(pw_df, file.path(out_dir, "diagnostic_pathways.csv"), row.names = FALSE)

# One line per species: last step reached in each variant
last_step <- function(r) { st <- c("subsetData", "identifyOverExpressedGenes", "identifyOverExpressedInteractions", "computeCommunProb", "filterCommunication", "computeCommunProbPathway"); s <- unlist(r[st]); if (all(s == "ok")) "ALL_OK" else names(s)[which(s != "ok")[1]] }
wide <- reshape(data.frame(variant = diag_df$variant, species = diag_df$species, stop_at = vapply(seq_len(nrow(diag_df)), function(i) last_step(diag_df[i, ]), character(1))),
                idvar = "species", timevar = "variant", direction = "wide")
print(wide, row.names = FALSE)
print(pw_df, row.names = FALSE)
