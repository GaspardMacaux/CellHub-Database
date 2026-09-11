################################################################################
#          convert_flyphonedb2.R — FlyPhoneDB2 → CellHub Database                      #
#                                                                              #
# Input:  raw_databases/Flyphone/FlyPhoneDB2_all.rds                           #
# Output: converted/flyphonedb2_drosophila.rds                                 #
#                                                                              #
# FlyPhoneDB2 v2.5 (2025) contains 1,804 Drosophila melanogaster LR pairs     #
# with 3 confidence tiers: high (316), moderate (651), low (837).              #
# Gene symbols follow FlyBase convention (Burs, sNPF, dpp, Toll).             #
# Includes mammalian orthologs and pathway annotations.                        #
# Reference: FullStackGoogler/FlyPhoneDB2 (GitHub)                            #
################################################################################

INPUT_DIR   <- "raw_databases/Flyphone"
OUTPUT_DIR  <- "converted"
SOURCE_NAME <- "FlyPhoneDB2"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

cat("[FlyPhoneDB2] Loading...\n")
db <- readRDS(file.path(INPUT_DIR, "FlyPhoneDB2_all.rds"))

# Use V2A (all 1,804 pairs including all confidence levels)
# Keep the rank column so users can filter by confidence in Cell-Hub
v2a <- db$V2A

converted <- data.frame(
  interaction_name   = paste(v2a$Gene_secreted, v2a$Gene_receptor, sep = "_"),
  interaction_name_2 = paste(v2a$Gene_secreted, "-", v2a$Gene_receptor),
  ligand             = v2a$Gene_secreted,
  receptor           = v2a$Gene_receptor,
  species            = "drosophila",
  source_db          = SOURCE_NAME,
  pathway_name       = ifelse(v2a$Pathway_receptor == "other", NA_character_, v2a$Pathway_receptor),
  interaction_type   = NA_character_,
  annotation         = v2a$Ligand_annotation,
  evidence           = ifelse(is.na(v2a$More_information) | v2a$More_information == "",
                               NA_character_, v2a$More_information),
  # FlyPhoneDB2-specific columns
  flyphone_confidence    = v2a$rank,
  flyphone_source        = v2a$Source,
  flyphone_fbgn_ligand   = v2a$FBgn_secreted,
  flyphone_fbgn_receptor = v2a$FBgn_receptor,
  stringsAsFactors       = FALSE
)

# Derive interaction_type from Ligand_annotation where possible
secreted_keywords <- c("secreted", "neuropeptide", "cytokine", "growth factor",
                        "hormone", "extracellular matrix", "morphogen")
membrane_keywords <- c("membrane", "transmembrane", "GPI-anchored",
                        "cell surface", "adhesion")

for (i in seq_len(nrow(converted))) {
  ann <- tolower(converted$annotation[i])
  if (is.na(ann)) next
  if (any(sapply(secreted_keywords, function(k) grepl(k, ann)))) {
    converted$interaction_type[i] <- "Secreted Signaling"
  } else if (any(sapply(membrane_keywords, function(k) grepl(k, ann)))) {
    converted$interaction_type[i] <- "Cell-Cell Contact"
  }
}

# Remove duplicates (same L-R pair from different annotation sources)
n_before <- nrow(converted)
converted <- converted[!duplicated(converted$interaction_name), ]
n_after <- nrow(converted)

out_path <- file.path(OUTPUT_DIR, "flyphonedb2_drosophila.rds")
saveRDS(converted, out_path)

cat(sprintf("  Drosophila: %d -> %d interactions (-%d dups)\n",
  n_before, n_after, n_before - n_after))
cat(sprintf("  Confidence: high=%d, moderate=%d, low=%d\n",
  sum(converted$flyphone_confidence == "high"),
  sum(converted$flyphone_confidence == "moderate"),
  sum(converted$flyphone_confidence == "low")))
cat(sprintf("  -> %s\n", out_path))
cat("[FlyPhoneDB2] Done!\n")
