# Step 08 (interim Cell-Hub fix, v2.0.1): CellChat-ready database files for the Cell-Hub app, same folder layout and file names as
# databases/<species>/ (merged + per-source), built from the corrected and repaired database of steps 07 / 07b
source("review/revision_helpers.R")
source("review/repair_functions.R")
suppressPackageStartupMessages(library(CellChat))
src_root <- "databases"; corrected_path <- "release_v2.1/CellHubDB.rds"; out_root <- "app_databases_v2.0.1"; placeholder <- "Unassigned"
canon <- c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor", "Non-protein Signaling")
db <- readRDS(corrected_path)
schema <- get_cellchat_schema()

# Per-source file made CellChat-ready: schema, species geneInfo, pathway labels taken from the corrected merged pair,
# interaction_type from annotation, then the same repairs as the merged database (step 07b)
fix_source_file <- function(x, merged, schema, species) {
  x <- standardize_cellchat_schema(x, schema)
  x$geneInfo <- merged$geneInfo
  it <- x$interaction
  key <- tolower(paste(it$ligand, it$receptor, sep = "|")); mkey <- tolower(paste(merged$interaction$ligand, merged$interaction$receptor, sep = "|"))
  missing_pw <- is.na(it$pathway_name) | it$pathway_name == ""
  it$pathway_name[missing_pw] <- merged$interaction$pathway_name[match(key[missing_pw], mkey)]
  it$pathway_name[is.na(it$pathway_name) | it$pathway_name == ""] <- placeholder
  if ("annotation" %in% colnames(it)) {
    if (!"interaction_type" %in% colnames(it)) it$interaction_type <- NA_character_
    use_ann <- !is.na(it$annotation) & it$annotation %in% canon
    it$interaction_type[use_ann] <- it$annotation[use_ann]
  }
  x$interaction <- it
  return(repair_species_db(x, species))
}

# Checks that Cell-Hub relies on: non-empty geneInfo, at least two pathways, full cofactor schema
check_file <- function(x, sp, file) {
  pw <- unique(x$interaction$pathway_name)
  data.frame(species = sp, file = file, n_interactions = nrow(x$interaction), n_geneinfo = nrow(x$geneInfo), n_pathways = length(pw),
             cofactor_ncol = ncol(x$cofactor), ok = nrow(x$geneInfo) > 0 && length(pw) >= 2 && ncol(x$cofactor) > 1, stringsAsFactors = FALSE)
}

checks <- list(); removed_src <- list()
for (sp in names(db)) {
  out_dir <- file.path(out_root, sp); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  merged_file <- paste0("CellHubDB_", sp, "_merged.rds")
  saveRDS(db[[sp]], file.path(out_dir, merged_file))
  checks[[length(checks) + 1]] <- check_file(db[[sp]], sp, merged_file)
  src_files <- list.files(file.path(src_root, sp), pattern = paste0("^CellHubDB_", sp, "_.*\\.rds$"))
  for (f in src_files[!grepl("_merged\\.rds$", src_files)]) {
    res <- fix_source_file(readRDS(file.path(src_root, sp, f)), db[[sp]], schema, sp)
    saveRDS(res$db, file.path(out_dir, f))
    checks[[length(checks) + 1]] <- check_file(res$db, sp, f)
    if (!is.null(res$removed) && nrow(res$removed) > 0) removed_src[[length(removed_src) + 1]] <- cbind(file = f, res$removed)
  }
}
checks_df <- do.call(rbind, checks)
write.csv(checks_df, file.path(out_root, "app_database_checks.csv"), row.names = FALSE)
if (length(removed_src) > 0) write.csv(do.call(rbind, removed_src), file.path(out_root, "removed_interactions_per_source.csv"), row.names = FALSE)
print(checks_df, row.names = FALSE)
message("Files written: ", nrow(checks_df), " | failing checks: ", sum(!checks_df$ok))
