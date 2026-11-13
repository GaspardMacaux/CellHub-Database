################################################################################
#          convert_celegans.R — C. elegans Armingol LR pairs → CellHub Database        #
#                                                                              #
# Input:  data/raw/Celegans_LR/Celegans-Curated-LR-pairs.xlsx            #
# Output: converted/celegans_celegans.rds                                      #
#                                                                              #
# The only curated LR database for C. elegans. 245 manually curated pairs     #
# with WormBase gene IDs, functional categories, confidence levels, and DOIs. #
# Gene symbols follow WormBase convention (daf-7, sma-6, arf-1.2).            #
# Reference: Armingol et al. (2022) PLOS Comp Biol. PMID:36395331            #
################################################################################

INPUT_DIR   <- "data/raw/Celegans_LR"
OUTPUT_DIR  <- "data/converted"
SOURCE_NAME <- "cell2cell_Armingol"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

cat("[C. elegans] Loading Excel file...\n")

# Read Excel — requires readxl package
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' required. Install with: install.packages('readxl')")
}

input_file <- file.path(INPUT_DIR, "Celegans-Curated-LR-pairs.xlsx")
# Also try alternative filename from Armingol compendium
if (!file.exists(input_file)) {
  input_file <- file.path(INPUT_DIR, "Celegans-2022-Armingol-LR-pairs.xlsx")
}
if (!file.exists(input_file)) {
  # Try listing the directory to find the file
  files <- list.files(INPUT_DIR, pattern = "\\.xlsx$", full.names = TRUE)
  if (length(files) > 0) input_file <- files[1]
}

lr <- readxl::read_excel(input_file)
cat(sprintf("  Loaded: %d rows | Columns: %s\n", nrow(lr), paste(colnames(lr), collapse = ", ")))

# Map Location_type to interaction_type
map_location_type <- function(loc) {
  if (is.na(loc)) return(NA_character_)
  loc <- tolower(loc)
  if (grepl("secreted", loc)) return("Secreted Signaling")
  if (grepl("membrane", loc)) return("Cell-Cell Contact")
  return(NA_character_)
}

converted <- data.frame(
  interaction_name   = paste(lr$Ligand_symbol, lr$Receptor_symbol, sep = "_"),
  interaction_name_2 = paste(lr$Ligand_symbol, "-", lr$Receptor_symbol),
  ligand             = lr$Ligand_symbol,
  receptor           = lr$Receptor_symbol,
  species            = "celegans",
  source_db          = SOURCE_NAME,
  pathway_name       = lr$`LR Function`,
  interaction_type   = sapply(lr$Location_type, map_location_type),
  annotation         = paste(lr$`L Function`, "->", lr$`R Function`),
  evidence           = lr$Confidence_evidence,
  # C. elegans specific
  ligand_wormbase    = lr$Ligand_WB,
  receptor_wormbase  = lr$Receptor_WB,
  pair_confidence    = lr$Pair_confidence,
  location_type      = lr$Location_type,
  stringsAsFactors   = FALSE
)

# Remove duplicates
n_before <- nrow(converted)
converted <- converted[!duplicated(converted$interaction_name), ]
n_after <- nrow(converted)

out_path <- file.path(OUTPUT_DIR, "celegans_celegans.rds")
saveRDS(converted, out_path)

cat(sprintf("  C. elegans: %d -> %d interactions (-%d dups)\n",
  n_before, n_after, n_before - n_after))
cat(sprintf("  Pathways: %s\n",
  paste(names(table(converted$pathway_name)), collapse = ", ")))
cat(sprintf("  Confidence: %s\n",
  paste(names(table(converted$pair_confidence)), table(converted$pair_confidence),
        sep = "=", collapse = ", ")))
cat(sprintf("  -> %s\n", out_path))
cat("[C. elegans] Done!\n")
