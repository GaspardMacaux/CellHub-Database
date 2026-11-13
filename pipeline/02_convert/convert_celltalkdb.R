################################################################################
#          convert_celltalkdb.R — CellTalkDB → CellHub Database                        #
#                                                                              #
# Input:  data/raw/CellTalkDB/human_lr_pair.txt                           #
#         data/raw/CellTalkDB/mouse_lr_pair.txt                           #
# Output: converted/celltalkdb_human.rds, converted/celltalkdb_mouse.rds       #
#                                                                              #
# CellTalkDB is a manually curated database with PMID evidence for each pair. #
# TSV format with gene symbols, Ensembl IDs, and PubMed references.           #
# Reference: Shao et al. (2020) Nucleic Acids Res.                            #
################################################################################

INPUT_DIR   <- "data/raw/CellTalkDB"
OUTPUT_DIR  <- "data/converted"
SOURCE_NAME <- "CellTalkDB"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

for (species in c("human", "mouse")) {
  input_file <- file.path(INPUT_DIR, sprintf("%s_lr_pair.txt", species))
  if (!file.exists(input_file)) {
    cat(sprintf("[CellTalkDB] SKIPPED %s — file not found: %s\n", species, input_file))
    next
  }

  cat(sprintf("[CellTalkDB] Loading %s...\n", species))
  lr <- read.delim(input_file, stringsAsFactors = FALSE, fileEncoding = "latin1")

  # Normalize PMID evidence: bare numbers -> "PMID:XXXXXXXX"
  evidence <- lr$evidence
  evidence <- ifelse(
    grepl("^[0-9]+$", evidence),
    paste0("PMID:", evidence),
    evidence
  )

  converted <- data.frame(
    interaction_name   = lr$lr_pair,
    interaction_name_2 = paste(lr$ligand_gene_symbol, "-", lr$receptor_gene_symbol),
    ligand             = lr$ligand_gene_symbol,
    receptor           = lr$receptor_gene_symbol,
    species            = species,
    source_db          = SOURCE_NAME,
    pathway_name       = NA_character_,
    interaction_type   = NA_character_,
    annotation         = NA_character_,
    evidence           = evidence,
    ligand_ensembl     = lr$ligand_ensembl_gene_id,
    receptor_ensembl   = lr$receptor_ensembl_gene_id,
    stringsAsFactors   = FALSE
  )

  out_path <- file.path(OUTPUT_DIR, sprintf("celltalkdb_%s.rds", species))
  saveRDS(converted, out_path)
  cat(sprintf("  %s: %d interactions -> %s\n", species, nrow(converted), out_path))
}

cat("[CellTalkDB] Done!\n")
