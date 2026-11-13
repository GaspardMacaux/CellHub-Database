# Step 06 (Reviewer 2, other 4): v1 vs v2 presence of the genes cited by Reviewer 2, and outdated / rescuable human symbols (HGNC)
source("revision/revision_helpers.R")
# v1.0 = the databases/ folder on the desktop copy (checked: arabidopsis 3,716 interactions); audit file of removed entries from the pipeline
v1_human_path <- "/home/gaspard/Bureau/CellHub_DB/databases/human/CellHubDB_human_merged.rds"; v2_path <- "release/2.1/CellHubDB.rds"
audit_removed_path <- Sys.getenv("CELLHUB_AUDIT", "audit/removed_entries.csv")
if (!file.exists(v1_human_path)) stop("v1 human database not found: ", v1_human_path)
out_dir <- "revision/results/06_symbols"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
genes_check <- c("S100A8", "S100A9", "P2RY12", "C10orf99", "GPR15LG", "DEFB4A", "DEFB4B")

# Number of interactions involving each gene (as partner or complex subunit, case-insensitive)
gene_presence <- function(sp_db, genes) {
  it <- sp_db$interaction
  partner_genes <- lapply(seq_len(nrow(it)), function(i) toupper(c(partner_subunits(it$ligand[i], sp_db$complex), partner_subunits(it$receptor[i], sp_db$complex))))
  return(vapply(toupper(genes), function(g) sum(vapply(partner_genes, function(p) g %in% p, logical(1))), integer(1)))
}
v1 <- readRDS(v1_human_path); v2 <- readRDS(v2_path)$human
message("v1 human interactions: ", nrow(v1$interaction), " (expected 10,950) | v2: ", nrow(v2$interaction), " (expected 11,479)")
presence <- data.frame(gene = genes_check, n_v1 = gene_presence(v1, genes_check), n_v2 = gene_presence(v2, genes_check))
write.csv(presence, file.path(out_dir, "reviewer_genes_v1_v2.csv"), row.names = FALSE)
print(presence)

# HGNC status of every human token in v2 (current map downloaded when possible, bundled map otherwise)
hgnc_map <- tryCatch(HGNChelper::getCurrentHumanMap(), error = function(e) NULL)
tokens <- extract_partner_genes(v2); tokens <- tokens[!grepl("-", tokens)]  # metabolite-enzyme composites are not HGNC symbols
chk <- if (is.null(hgnc_map)) HGNChelper::checkGeneSymbols(tokens, species = "human") else HGNChelper::checkGeneSymbols(tokens, map = hgnc_map)
outdated <- chk[!chk$Approved, ]
write.csv(outdated, file.path(out_dir, "v2_human_non_approved_symbols.csv"), row.names = FALSE)
message("Non-approved human tokens in v2: ", nrow(outdated), " (with a suggested symbol: ", sum(!is.na(outdated$Suggested.Symbol)), ")")

# Removed entries whose symbol could have been renamed instead of dropped
if (file.exists(audit_removed_path)) {
  removed <- read.csv(audit_removed_path, stringsAsFactors = FALSE)
  if (!all(c("ligand", "receptor") %in% colnames(removed))) stop("Audit file needs ligand and receptor columns; found: ", paste(colnames(removed), collapse = ", "))
  rem_genes <- unique(c(removed$ligand, removed$receptor)); rem_genes <- rem_genes[!is.na(rem_genes) & !grepl("[-_]", rem_genes)]
  rchk <- if (is.null(hgnc_map)) HGNChelper::checkGeneSymbols(rem_genes, species = "human") else HGNChelper::checkGeneSymbols(rem_genes, map = hgnc_map)
  rescuable <- rchk[!rchk$Approved & !is.na(rchk$Suggested.Symbol), ]
  write.csv(rescuable, file.path(out_dir, "removed_but_rescuable.csv"), row.names = FALSE)
  message("Removed symbols rescuable through HGNC: ", nrow(rescuable))
} else {
  message("Audit file not found: ", audit_removed_path)
}
