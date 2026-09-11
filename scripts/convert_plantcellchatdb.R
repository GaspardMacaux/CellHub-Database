################################################################################
#     convert_plantcellchatdb.R — PlantCellChatDB → CellHub Database                   #
#                                                                              #
# Input:  raw_databases/PlantCellChat/ath.csv  (Arabidopsis thaliana)           #
#         raw_databases/PlantCellChat/osa.csv  (Oryza sativa / rice)           #
#         raw_databases/PlantCellChat/sly.csv  (Solanum lycopersicum / tomato) #
#         raw_databases/PlantCellChat/zma.csv  (Zea mays / maize)             #
#         raw_databases/PlantCellChat/gmx.csv  (Glycine max / soybean)        #
# Output: converted/plantcellchatdb_<species>.rds (one per species)           #
#                                                                              #
# PlantCellChatDB provides 12,717 LR pairs across 5 plant species in a        #
# CellChat-like format (with Agonist, Antagonist, CoAreceptor, etc.).         #
# Includes hormone/signaling molecule interactions (ABA, auxin, JA, SA, etc.) #
# with PMIDs and KEGG references as evidence.                                 #
# Reference: Liu et al. (2026) Plant J. PMID:42119054                         #
################################################################################
setwd("/media/gaspard/Gasp/CellHub Database")
INPUT_DIR   <- "raw_databases/PlantCellChatDB"
OUTPUT_DIR  <- "converted"
SOURCE_NAME <- "PlantCellChatDB"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Species file mapping: filename prefix -> CellHub Database species code
SPECIES_FILES <- c(
  "ath" = "arabidopsis",
  "osa" = "rice",
  "sly" = "tomato",
  "zma" = "maize",
  "gmx" = "soybean"
)

cat("[PlantCellChatDB] Loading CSVs...\n")

for (prefix in names(SPECIES_FILES)) {
  sp_code <- SPECIES_FILES[prefix]
  input_file <- file.path(INPUT_DIR, paste0(prefix, ".csv"))
  
  if (!file.exists(input_file)) {
    cat(sprintf("  %s (%s): SKIPPED — file not found: %s\n", sp_code, prefix, input_file))
    next
  }
  
  lr <- read.csv(input_file, stringsAsFactors = FALSE)
  cat(sprintf("  %s (%s): %d interactions loaded\n", sp_code, prefix, nrow(lr)))
  
  # Map Interaction_type to CellChat canonical types
  map_type <- function(itype) {
    if (is.na(itype) || itype == "") return(NA_character_)
    itype_lower <- tolower(itype)
    if (grepl("signaling molecule", itype_lower)) return("Non-protein Signaling")
    if (grepl("secreted", itype_lower)) return("Secreted Signaling")
    if (grepl("cell-cell|contact|adhesion", itype_lower)) return("Cell-Cell Contact")
    if (grepl("ecm|matrix", itype_lower)) return("ECM-Receptor")
    return(NA_character_)
  }
  
  converted <- data.frame(
    interaction_name   = lr$Interaction_name,
    interaction_name_2 = paste(lr$Ligand, "-", lr$Receptor),
    ligand             = lr$Ligand,
    receptor           = lr$Receptor,
    species            = sp_code,
    source_db          = SOURCE_NAME,
    pathway_name       = ifelse(lr$Signal == "" | is.na(lr$Signal), NA_character_, lr$Signal),
    interaction_type   = sapply(lr$Interaction_type, map_type),
    annotation         = ifelse(lr$Origin_of_interaction == "" | is.na(lr$Origin_of_interaction),
                                NA_character_, lr$Origin_of_interaction),
    evidence           = ifelse(lr$Evidence == "" | is.na(lr$Evidence),
                                NA_character_, lr$Evidence),
    # CellChat-compatible extra columns
    agonist            = ifelse(lr$Agonist == "" | is.na(lr$Agonist),
                                NA_character_, lr$Agonist),
    antagonist         = ifelse(lr$Antagonist == "" | is.na(lr$Antagonist),
                                NA_character_, lr$Antagonist),
    co_A_receptor      = ifelse(lr$CoAreceptor == "" | is.na(lr$CoAreceptor),
                                NA_character_, lr$CoAreceptor),
    co_I_receptor      = ifelse(lr$CoIreceptor == "" | is.na(lr$CoIreceptor),
                                NA_character_, lr$CoIreceptor),
    complex1           = ifelse(lr$Complex1 == "" | is.na(lr$Complex1),
                                NA_character_, lr$Complex1),
    complex2           = ifelse(lr$Complex2 == "" | is.na(lr$Complex2),
                                NA_character_, lr$Complex2),
    signal             = lr$Signal,
    stringsAsFactors   = FALSE
  )
  
  # Remove duplicates
  n_before <- nrow(converted)
  converted <- converted[!duplicated(converted$interaction_name), ]
  n_after <- nrow(converted)
  
  out_path <- file.path(OUTPUT_DIR, sprintf("plantcellchatdb_%s.rds", sp_code))
  saveRDS(converted, out_path)
  
  # Stats
  n_signaling_mol <- sum(converted$interaction_type == "Non-protein Signaling", na.rm = TRUE)
  n_secreted <- sum(converted$interaction_type == "Secreted Signaling", na.rm = TRUE)
  n_evidence <- sum(!is.na(converted$evidence))
  n_pathways <- length(unique(converted$pathway_name[!is.na(converted$pathway_name)]))
  
  cat(sprintf("    %d -> %d interactions (-%d dups)\n", n_before, n_after, n_before - n_after))
  cat(sprintf("    Signaling molecules: %d | Secreted: %d | Pathways: %d | Evidence: %d\n",
              n_signaling_mol, n_secreted, n_pathways, n_evidence))
  cat(sprintf("    -> %s\n", basename(out_path)))
}

cat("[PlantCellChatDB] Done!\n")