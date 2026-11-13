################################################################################
#          convert_connectomedb2025.R — connectomeDB2025 → CellHub Database            #
#                                                                              #
# Input:  data/raw/all_species.csv  (from https://connectomedb.org)       #
# Output: converted/connectomedb2025_<species>.rds  (one per species)          #
#                                                                              #
# connectomeDB2025 provides 3,579 curated LR pairs across 14 vertebrate       #
# species, each with experimental evidence (Direct) or orthology mapping       #
# (Inferred). Gene symbols are already in species-correct casing.              #
# Reference: Liu P et al. (2026) Nucleic Acids Res. PMID:41171146             #
################################################################################

INPUT_PATH  <- "data/raw/ConnectomeDB2025/all_species.csv"
OUTPUT_DIR  <- "data/converted"
SOURCE_NAME <- "connectomeDB2025"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

cat(sprintf("[connectomeDB2025] Loading %s...\n", INPUT_PATH))
raw <- read.csv(INPUT_PATH, stringsAsFactors = FALSE)
cat(sprintf("[connectomeDB2025] Loaded: %d rows | %d columns\n", nrow(raw), ncol(raw)))
cat(sprintf("[connectomeDB2025] Species: %s\n",
  paste(unique(raw$Species), collapse = ", ")))

# ==============================================================================
# Species name mapping: connectomeDB names -> CellHub Database species codes
# ==============================================================================

SPECIES_MAP <- c(
  "human"     = "human",
  "mouse"     = "mouse",
  "rat"       = "rat",
  "zebrafish" = "zebrafish",
  "chicken"   = "chicken",
  "frog"      = "frog",
  "chimp"     = "chimp",
  "macaque"   = "macaque",
  "marmoset"  = "marmoset",
  "horse"     = "horse",
  "pig"       = "pig",
  "cow"       = "cow",
  "dog"       = "dog",
  "sheep"     = "sheep"
)

# ==============================================================================
# Interaction type mapping: ligand/receptor location -> CellChat categories
# ==============================================================================
# CellChat uses 4 categories: Secreted Signaling, Cell-Cell Contact,
# ECM-Receptor, Non-protein Signaling.
# connectomeDB2025 uses subcellular localization (secreted, cell membrane, etc.)

map_interaction_type <- function(ligand_loc, receptor_loc) {
  lig <- tolower(ligand_loc)
  rec <- tolower(receptor_loc)
  # Both membrane-bound -> Cell-Cell Contact (juxtacrine)
  if (grepl("cell membrane", lig) && !grepl("secreted", lig) &&
      grepl("cell membrane", rec)) {
    return("Cell-Cell Contact")
  }
  # Ligand is secreted -> Secreted Signaling
  if (grepl("secreted", lig)) {
    return("Secreted Signaling")
  }
  # GPI-anchored ligand -> Cell-Cell Contact
  if (grepl("gpi", lig)) {
    return("Cell-Cell Contact")
  }
  # Fallback
  return(NA_character_)
}

# ==============================================================================
# Extract primary gene symbol from "GeneSymbol (alias1, alias2, ...)" format
# ==============================================================================

extract_primary_symbol <- function(symbol_str) {
  if (is.na(symbol_str) || symbol_str == "") return(NA_character_)
  # Clean HTML artifacts from connectomeDB2025 CSV export
  # e.g. ', Cd, ska26, skax26)"">Lrp6 (ska, Cd, ska26, skax26)'
  clean <- symbol_str
  clean <- gsub("<[^>]+>", "", clean)        # Strip HTML tags
  clean <- gsub('"">', "", clean)             # Strip escaped quote artifacts
  clean <- gsub("^[^A-Za-z0-9]+", "", clean) # Strip leading garbage
  # Take everything before the first "("
  primary <- trimws(strsplit(clean, "\\(")[[1]][1])
  primary <- sub("\\s+$", "", primary)
  primary
}

# ==============================================================================
# Convert all rows to CellHub Database intermediate format
# ==============================================================================

cat("[connectomeDB2025] Converting to CellHub Database format...\n")

converted <- data.frame(
  interaction_name = character(nrow(raw)),
  interaction_name_2 = character(nrow(raw)),
  ligand = character(nrow(raw)),
  receptor = character(nrow(raw)),
  species = character(nrow(raw)),
  source_db = rep(SOURCE_NAME, nrow(raw)),
  pathway_name = rep(NA_character_, nrow(raw)),
  interaction_type = character(nrow(raw)),
  annotation = rep(NA_character_, nrow(raw)),
  evidence = character(nrow(raw)),
  # Extra columns from connectomeDB2025
  connectome_id = raw$Interaction.ID,
  evidence_type = raw$Evidence,
  ligand_ensembl = raw$Ligand.ENSEMBL.ID,
  receptor_ensembl = raw$Receptor.ENSEMBL.ID,
  ligand_location = raw$Ligand.Location,
  receptor_location = raw$Receptor.Location,
  human_ligand = character(nrow(raw)),
  human_receptor = character(nrow(raw)),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(raw))) {
  # Extract primary gene symbols from LR Pair column (clean, space-separated)
  # This is more reliable than Ligand/Receptor Symbols which can contain
  # HTML artifacts from connectomeDB2025's export
  lr_pair <- raw$LR.Pair[i]
  lr_parts <- strsplit(lr_pair, " ")[[1]]
  
  if (length(lr_parts) >= 2) {
    lig <- lr_parts[1]
    rec <- lr_parts[2]
  } else {
    # Fallback to Ligand/Receptor Symbols columns
    lig <- extract_primary_symbol(raw$Ligand.Symbols[i])
    rec <- extract_primary_symbol(raw$Receptor.Symbols[i])
  }
  
  converted$ligand[i] <- lig
  converted$receptor[i] <- rec
  
  # Build interaction names
  converted$interaction_name[i] <- paste(lig, rec, sep = "_")
  converted$interaction_name_2[i] <- paste(lig, "-", rec)
  
  # Map species
  sp_raw <- raw$Species[i]
  converted$species[i] <- ifelse(sp_raw %in% names(SPECIES_MAP),
                                 SPECIES_MAP[sp_raw], sp_raw)
  
  # Map interaction type from location
  converted$interaction_type[i] <- map_interaction_type(
    raw$Ligand.Location[i], raw$Receptor.Location[i])
  
  # Evidence
  converted$evidence[i] <- paste0(SOURCE_NAME, ":", raw$Evidence[i])
  
  # Human orthologs
  converted$human_ligand[i] <- extract_primary_symbol(raw$Human.Ligand.Symbols[i])
  converted$human_receptor[i] <- extract_primary_symbol(raw$Human.Receptor.Symbols[i])
}
# ==============================================================================
# Split by species and save
# ==============================================================================

species_list <- unique(converted$species)
cat(sprintf("[connectomeDB2025] Splitting into %d species...\n", length(species_list)))

for (sp in species_list) {
  sp_data <- converted[converted$species == sp, ]

  # Remove within-species duplicates (same interaction_name)
  n_before <- nrow(sp_data)
  sp_data <- sp_data[!duplicated(sp_data$interaction_name), ]
  n_after <- nrow(sp_data)

  out_path <- file.path(OUTPUT_DIR, sprintf("connectomedb2025_%s.rds", sp))
  saveRDS(sp_data, out_path)

  n_direct <- sum(sp_data$evidence_type == "Direct")
  n_inferred <- sum(sp_data$evidence_type == "Inferred")

  cat(sprintf("  %s: %d interactions (%d direct, %d inferred",
    sp, n_after, n_direct, n_inferred))
  if (n_before != n_after) cat(sprintf(", %d dups removed", n_before - n_after))
  cat(sprintf(") -> %s\n", basename(out_path)))
}

# ==============================================================================
# Summary
# ==============================================================================

cat("\n[connectomeDB2025] Conversion complete!\n")
cat(sprintf("  Total input rows: %d\n", nrow(raw)))
cat(sprintf("  Total output interactions: %d\n", nrow(converted)))
cat(sprintf("  Species covered: %d\n", length(species_list)))
cat(sprintf("  Output files in: %s/\n", OUTPUT_DIR))
