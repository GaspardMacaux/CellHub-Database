################################################################################
#          convert_cellchatdb.R — CellChatDB v2 → CellHub Database                     #
#                                                                              #
# Input:  raw_databases/CellChat/CellChatDB_v2_<species>.rds                   #
#         (mouse, human, zebrafish — as distributed with CellChat)             #
#         database_curation/complex_reconstructions.tsv                        #
# Output: converted/cellchatdb_v2_<species>.rds                                #
#         geneinfo/geneinfo_<species>.rds                                      #
#         geneinfo/complexes.rds, geneinfo/cofactors.rds                       #
#                                                                              #
# CellChatDB is already in CellChat-native format. This script extracts the    #
# interaction table and adds the mandatory species/source_db columns.          #
# It also exports geneInfo, complex, and cofactor tables for the pipeline.     #
# Complexes used by CellChatDB interactions but missing from its complex table #
# (zebrafish numbered paralog variants) are added from complex_reconstructions.tsv
################################################################################

INPUT_DIR   <- "raw_databases/CellChat"
OUTPUT_DIR  <- "converted"
GENEINFO_DIR <- "geneinfo"
SOURCE_NAME <- "CellChatDB_v2"
RECONSTRUCTION_FILE <- "database_curation/complex_reconstructions.tsv"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(GENEINFO_DIR)) dir.create(GENEINFO_DIR, recursive = TRUE)

complexes_all  <- list()
cofactors_all  <- list()
reconstructions <- if (file.exists(RECONSTRUCTION_FILE)) read.delim(RECONSTRUCTION_FILE, stringsAsFactors = FALSE, colClasses = "character") else NULL

# Add reconstructed complex definitions to a CellChat complex table (only complexes it does not define); returns the table
add_reconstructed_complexes <- function(complex_df, rec, species) {
  rec <- rec[rec$species == species & !rec$complex %in% rownames(complex_df), , drop = FALSE]
  if (nrow(rec) == 0) return(complex_df)
  complex_df[] <- lapply(complex_df, as.character)
  sub_cols <- grep("^subunit_", colnames(rec), value = TRUE)
  for (cl in setdiff(sub_cols, colnames(complex_df))) complex_df[[cl]] <- rep("", nrow(complex_df))
  added <- as.data.frame(matrix("", nrow = nrow(rec), ncol = ncol(complex_df), dimnames = list(rec$complex, colnames(complex_df))), stringsAsFactors = FALSE)
  for (cl in sub_cols) added[[cl]] <- ifelse(is.na(rec[[cl]]), "", rec[[cl]])
  cat(sprintf("  complexes: %d reconstructed definitions added (%s)\n", nrow(rec), paste(rec$complex, collapse = ", ")))
  return(rbind(complex_df, added))
}

# All species shipped with CellChat; zebrafish was previously extracted but never converted
for (species in c("mouse", "human", "zebrafish")) {
  input_file <- file.path(INPUT_DIR, sprintf("CellChatDB_v2_%s.rds", species))
  if (!file.exists(input_file)) {
    cat(sprintf("[CellChatDB] SKIPPED %s — file not found: %s\n", species, input_file))
    next
  }

  cat(sprintf("[CellChatDB] Loading %s...\n", species))
  db <- readRDS(input_file)
  inter <- db$interaction

  # Add mandatory columns
  inter$species   <- species
  inter$source_db <- SOURCE_NAME

  # Ensure standard columns exist
  if (!"pathway_name" %in% colnames(inter)) inter$pathway_name <- NA_character_
  if (!"interaction_type" %in% colnames(inter)) inter$interaction_type <- NA_character_
  if (!"annotation" %in% colnames(inter)) inter$annotation <- NA_character_

  # Save interaction table
  out_path <- file.path(OUTPUT_DIR, sprintf("cellchatdb_v2_%s.rds", species))
  saveRDS(inter, out_path)
  cat(sprintf("  %s: %d interactions -> %s\n", species, nrow(inter), out_path))

  # Save geneInfo (used by the curation pipeline for casing and symbol validation)
  gi_path <- file.path(GENEINFO_DIR, sprintf("geneinfo_%s.rds", species))
  saveRDS(db$geneInfo, gi_path)
  cat(sprintf("  geneInfo: %d genes -> %s\n", nrow(db$geneInfo), gi_path))

  # Collect complex/cofactor tables, completed with the reconstructed definitions
  complex_df <- db$complex
  if (!is.null(reconstructions)) complex_df <- add_reconstructed_complexes(complex_df, reconstructions, species)
  # Complexes still used by interactions without any definition (reported, handled by the finalize step)
  partners <- unique(c(inter$ligand, inter$receptor))
  undefined <- partners[grepl("_", partners) & !partners %in% rownames(complex_df)]
  if (length(undefined) > 0) cat(sprintf("  complexes used but undefined: %d (%s)\n", length(undefined), paste(head(undefined, 10), collapse = ", ")))
  complexes_all[[species]]  <- complex_df
  cofactors_all[[species]]  <- db$cofactor
}

# Save complex and cofactor tables (keyed by species)
saveRDS(complexes_all, file.path(GENEINFO_DIR, "complexes.rds"))
saveRDS(cofactors_all, file.path(GENEINFO_DIR, "cofactors.rds"))
cat("[CellChatDB] Complex/cofactor tables saved.\n")
cat("[CellChatDB] Done!\n")
