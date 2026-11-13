################################################################################
#          convert_plantphonedb.R — PlantPhoneDB → CellHub Database                    #
#                                                                              #
# Input:  data/raw/Plants-2022-Xu-LR-pairs/Athaliana-*.csv                #
# Output: converted/plantphonedb_arabidopsis.rds                               #
#                                                                              #
# PlantPhoneDB provides 3,657 LR pairs for Arabidopsis thaliana, curated from #
# 7 resources (plant.MAP, Interactome2.0, IntAct, BioGRID, Text-mining,       #
# STRING, Orthologs). Gene identifiers are AGI locus IDs (AT1G01900),         #
# which is the standard TAIR convention for Arabidopsis.                       #
# Only Arabidopsis is integrated (other species use UniProt/assembly IDs       #
# that cannot be resolved without biomaRt).                                    #
# Reference: Xu et al. (2022) Plant J.                                        #
################################################################################

INPUT_DIR   <- "data/raw/Plants-2022-Xu-LR-pairs"
OUTPUT_DIR  <- "data/converted"
SOURCE_NAME <- "PlantPhoneDB"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Find the Arabidopsis file
at_files <- list.files(INPUT_DIR, pattern = "Athaliana.*\\.csv$", full.names = TRUE)
if (length(at_files) == 0) {
  stop(sprintf("No Arabidopsis CSV found in %s", INPUT_DIR))
}
input_file <- at_files[1]

cat(sprintf("[PlantPhoneDB] Loading %s...\n", basename(input_file)))
lr <- read.csv(input_file, stringsAsFactors = FALSE)
cat(sprintf("  Loaded: %d rows | Columns: %s\n", nrow(lr), paste(colnames(lr), collapse = ", ")))

converted <- data.frame(
  interaction_name   = paste(lr$Ligands, lr$Receptors, sep = "_"),
  interaction_name_2 = paste(lr$Ligands, "-", lr$Receptors),
  ligand             = lr$Ligands,
  receptor           = lr$Receptors,
  species            = "arabidopsis",
  source_db          = SOURCE_NAME,
  pathway_name       = NA_character_,
  interaction_type   = NA_character_,
  annotation         = NA_character_,
  evidence           = lr$source,  # "STRING", "IntAct", "BioGRID", "orthologs", etc.
  plantphonedb_source = lr$source,
  stringsAsFactors   = FALSE
)

# Remove duplicates
n_before <- nrow(converted)
converted <- converted[!duplicated(converted$interaction_name), ]
n_after <- nrow(converted)

out_path <- file.path(OUTPUT_DIR, "plantphonedb_arabidopsis.rds")
saveRDS(converted, out_path)

cat(sprintf("  Arabidopsis: %d -> %d interactions (-%d dups)\n",
  n_before, n_after, n_before - n_after))
cat(sprintf("  Evidence sources: %s\n",
  paste(names(table(converted$plantphonedb_source)), table(converted$plantphonedb_source),
        sep = "=", collapse = ", ")))
cat(sprintf("  -> %s\n", out_path))
cat("[PlantPhoneDB] Done!\n")
