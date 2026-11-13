# Step 07 (Reviewer 1 Python interface; Reviewer 2 other 2): tool-agnostic TSV export per species
# (LIANA+-style ligand/receptor columns, complex subunits joined by "_") and MD5 manifest of the release files
source("review/revision_helpers.R")
db_path <- "review/03_confidence/CellHubDB_annotated.rds"; out_dir <- "release_v2.1"
dir.create(file.path(out_dir, "tsv"), recursive = TRUE, showWarnings = FALSE)
db <- readRDS(db_path)
canon <- c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor", "Non-protein Signaling")

# interaction_type rebuilt from the CellChat annotation column when it holds a canonical category (the original conversion lost
# most ECM-Receptor labels of CellChatDB entries: 1 instead of 424 in human); other entries keep their mapped type
harmonize_interaction_type <- function(sp_db) {
  it <- sp_db$interaction
  if (!"interaction_type" %in% colnames(it)) it$interaction_type <- NA_character_
  if ("annotation" %in% colnames(it)) {
    use_ann <- !is.na(it$annotation) & it$annotation %in% canon
    it$interaction_type[use_ann] <- it$annotation[use_ann]
  }
  sp_db$interaction <- it
  return(sp_db)
}
count_types <- function(db) do.call(rbind, lapply(names(db), function(sp) {
  t <- db[[sp]]$interaction$interaction_type; t[is.na(t) | t == ""] <- "Unlabelled"; t[!t %in% c(canon, "Unlabelled")] <- "Other"
  data.frame(species = sp, t(as.integer(table(factor(t, levels = c(canon, "Other", "Unlabelled"))))), check.names = FALSE) |> setNames(c("species", canon, "Other", "Unlabelled"))
}))
before <- count_types(db)
db <- lapply(db, harmonize_interaction_type)
after <- count_types(db)
types_report <- merge(before, after, by = "species", suffixes = c("_before", "_after"))
write.csv(types_report, file.path(out_dir, "interaction_type_harmonization.csv"), row.names = FALSE)
print(after, row.names = FALSE)
saveRDS(db, file.path(out_dir, "CellHubDB.rds"))

# Complex partners expanded to their subunits joined by "_"
expand_partner <- function(p, cpx) paste(partner_subunits(p, cpx), collapse = "_")
for (sp in names(db)) {
  it <- db[[sp]]$interaction; cpx <- db[[sp]]$complex
  keep_cols <- intersect(c("interaction_name", "pathway_name", "pathway_origin", "interaction_type", "annotation", "sources", "n_sources", "n_curated_sources", "predicted_only",
                           "is_complex", "role_corroborated", "confidence_tier", "evidence_type", "evidence"), colnames(it))
  out <- data.frame(ligand = vapply(it$ligand, expand_partner, character(1), cpx = cpx, USE.NAMES = FALSE),
                    receptor = vapply(it$receptor, expand_partner, character(1), cpx = cpx, USE.NAMES = FALSE),
                    ligand_cellchat = it$ligand, receptor_cellchat = it$receptor, it[, keep_cols, drop = FALSE], stringsAsFactors = FALSE)
  write.table(out, file.path(out_dir, "tsv", paste0("CellHubDB_", sp, ".tsv")), sep = "\t", row.names = FALSE, quote = TRUE)
}

# MD5 manifest of every release file (to archive with the Zenodo deposit)
files <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
files <- files[basename(files) != "MD5SUMS.tsv"]
manifest <- data.frame(file = substring(normalizePath(files), nchar(normalizePath(out_dir)) + 2), size_bytes = file.info(files)$size, md5 = unname(tools::md5sum(files)))
write.table(manifest, file.path(out_dir, "MD5SUMS.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(manifest)
