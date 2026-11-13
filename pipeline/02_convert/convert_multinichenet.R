################################################################################
#          convert_multinichenet.R — NicheNet v2 LR network → CellHub Database         #
#                                                                              #
# Input:  data/raw/Multinichenet/MultiNicheNet_lr_network_human.rds       #
#         data/raw/Multinichenet/MultiNicheNet_lr_network_mouse.rds       #
# Output: converted/multinichenet_human.rds, converted/multinichenet_mouse.rds #
#                                                                              #
# NicheNet v2 LR network (from Zenodo) contains 5,281 human + 5,877 mouse     #
# pairs with source databases, reference PMIDs, and resource counts.           #
################################################################################

INPUT_DIR   <- "data/raw/Multinichenet"
OUTPUT_DIR  <- "data/converted"
SOURCE_NAME <- "MultiNicheNet"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

for (species in c("human", "mouse")) {
  input_file <- file.path(INPUT_DIR, sprintf("MultiNicheNet_lr_network_%s.rds", species))
  if (!file.exists(input_file)) {
    cat(sprintf("[MultiNicheNet] SKIPPED %s — file not found: %s\n", species, input_file))
    next
  }

  cat(sprintf("[MultiNicheNet] Loading %s...\n", species))
  lr <- readRDS(input_file)

  converted <- data.frame(
    interaction_name   = paste(lr$ligand, lr$receptor, sep = "_"),
    interaction_name_2 = paste(lr$ligand, "-", lr$receptor),
    ligand             = lr$ligand,
    receptor           = lr$receptor,
    species            = species,
    source_db          = SOURCE_NAME,
    pathway_name       = NA_character_,
    interaction_type   = NA_character_,
    annotation         = NA_character_,
    # NicheNet provides actual references (PMIDs) and source DB info
    evidence           = ifelse(is.na(lr$references), NA_character_, lr$references),
    nichenet_database   = lr$database,
    nichenet_sources    = lr$sources,
    n_resources         = lr$n_resources,
    n_references        = lr$n_references,
    is_directed         = lr$is_directed,
    is_stimulation      = lr$is_stimulation,
    is_inhibition       = lr$is_inhibition,
    stringsAsFactors    = FALSE
  )

  # Remove duplicates (same L-R pair from multiple NicheNet sub-databases)
  n_before <- nrow(converted)
  converted <- converted[!duplicated(converted$interaction_name), ]
  n_after <- nrow(converted)

  out_path <- file.path(OUTPUT_DIR, sprintf("multinichenet_%s.rds", species))
  saveRDS(converted, out_path)

  cat(sprintf("  %s: %d -> %d interactions (-%d dups) -> %s\n",
    species, n_before, n_after, n_before - n_after, out_path))
}

cat("[MultiNicheNet] Done!\n")
