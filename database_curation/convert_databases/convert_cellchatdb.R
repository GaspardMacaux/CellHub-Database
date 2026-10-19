################################################################################
#          convert_cellchatdb.R — CellChatDB v2 → CellHub Database                     #
#                                                                              #
# Input:  raw_databases/CellChat/CellChatDB_v2_<species>.rds                   #
#         (mouse, human, zebrafish — as distributed with CellChat)             #
# Output: converted/cellchatdb_v2_<species>.rds                                #
#         geneinfo/geneinfo_<species>.rds                                      #
#         geneinfo/complexes.rds, geneinfo/cofactors.rds                       #
#                                                                              #
# CellChatDB is already in CellChat-native format. This script extracts the    #
# interaction table and adds the mandatory species/source_db columns.          #
# It also exports geneInfo, complex, and cofactor tables for the pipeline.     #
################################################################################

INPUT_DIR   <- "raw_databases/CellChat"
OUTPUT_DIR  <- "converted"
GENEINFO_DIR <- "geneinfo"
SOURCE_NAME <- "CellChatDB_v2"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(GENEINFO_DIR)) dir.create(GENEINFO_DIR, recursive = TRUE)

complexes_all  <- list()
cofactors_all  <- list()

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
  
  # Collect complex/cofactor tables
  complexes_all[[species]]  <- db$complex
  cofactors_all[[species]]  <- db$cofactor
}

# Save complex and cofactor tables (keyed by species)
saveRDS(complexes_all, file.path(GENEINFO_DIR, "complexes.rds"))
saveRDS(cofactors_all, file.path(GENEINFO_DIR, "cofactors.rds"))
cat("[CellChatDB] Complex/cofactor tables saved.\n")
cat("[CellChatDB] Done!\n")