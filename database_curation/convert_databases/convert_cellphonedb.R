################################################################################
#          convert_cellphonedb.R — CellPhoneDB v5 → CellHub Database                   #
#                                                                              #
# Input:  raw_databases/data_interaction_input.csv                             #
#         raw_databases/data_gene_input.csv                                    #
#         raw_databases/data_complex_input.csv                                 #
#         raw_databases/data_protein_input.csv                                 #
# Output: converted/cellphonedb_human.rds                                      #
#                                                                              #
# CellPhoneDB v5 uses UniProt accessions and complex names as partner IDs.     #
# This converter resolves them to HGNC gene symbols using gene_input.csv       #
# and complex_input.csv. The classification column provides rich pathway       #
# annotations (Signaling by Glutamate, Adhesion by Cadherin, etc.)             #
#                                                                              #
# CellPhoneDB is human-only (manually curated for Homo sapiens).              #
# Reference: Efremova M et al. (2020) Nat Protoc. + Troulé K et al. (2023)   #
################################################################################

INPUT_DIR   <- "raw_databases/cellphonedb-data-5.0.0/data"
OUTPUT_DIR  <- "converted"
SOURCE_NAME <- "CellPhoneDB"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ==============================================================================
# Load mapping tables
# ==============================================================================

cat("[CellPhoneDB] Loading mapping tables...\n")

# UniProt -> HGNC gene symbol
gene_input <- read.csv(file.path(INPUT_DIR, "gene_input.csv"), 
                        stringsAsFactors = FALSE)
uniprot_to_gene <- setNames(gene_input$hgnc_symbol, gene_input$uniprot)
cat(sprintf("  gene_input: %d UniProt -> HGNC mappings\n", nrow(gene_input)))

# Complex name -> list of subunit UniProt IDs
complex_input <- read.csv(file.path(INPUT_DIR, "complex_input.csv"),
                           stringsAsFactors = FALSE)
# Build complex lookup: name -> vector of gene symbols
complex_to_genes <- list()
for (i in seq_len(nrow(complex_input))) {
  cname <- complex_input$complex_name[i]
  subunits <- c()
  for (j in 1:5) {
    col <- paste0("uniprot_", j)
    uid <- complex_input[[col]][i]
    if (!is.na(uid) && uid != "") {
      # Resolve UniProt subunit to gene symbol
      gene <- uniprot_to_gene[uid]
      if (!is.na(gene)) {
        subunits <- c(subunits, gene)
      } else {
        subunits <- c(subunits, uid)  # Keep UniProt ID if unresolvable
      }
    }
  }
  complex_to_genes[[cname]] <- subunits
}
cat(sprintf("  complex_input: %d complexes resolved\n", length(complex_to_genes)))

# Protein annotations: for deriving interaction_type
protein_input <- read.csv(file.path(INPUT_DIR, "protein_input.csv"),
                           stringsAsFactors = FALSE)
# Build protein property lookup by UniProt ID
protein_props <- list()
for (i in seq_len(nrow(protein_input))) {
  uid <- protein_input$uniprot[i]
  protein_props[[uid]] <- list(
    secreted      = toupper(protein_input$secreted[i]) == "TRUE",
    transmembrane = toupper(protein_input$transmembrane[i]) == "TRUE",
    peripheral    = toupper(protein_input$peripheral[i]) == "TRUE",
    receptor      = toupper(protein_input$receptor[i]) == "TRUE",
    integrin      = toupper(protein_input$integrin[i]) == "TRUE"
  )
}
# Same for complexes
complex_props <- list()
for (i in seq_len(nrow(complex_input))) {
  cname <- complex_input$complex_name[i]
  complex_props[[cname]] <- list(
    secreted      = toupper(complex_input$secreted[i]) == "TRUE",
    transmembrane = toupper(complex_input$transmembrane[i]) == "TRUE",
    receptor      = toupper(complex_input$receptor[i]) == "TRUE",
    integrin      = toupper(complex_input$integrin[i]) == "TRUE"
  )
}
cat(sprintf("  protein_input: %d protein annotations\n", nrow(protein_input)))

# ==============================================================================
# Helper: resolve a partner ID to gene symbol(s)
# ==============================================================================

# A partner can be:
#   - A UniProt accession (P12830 -> CDH1)
#   - A complex name (integrin_a2b1_complex -> ITGA2_ITGB1)
#   - A non-standard ID (HLAA, HLAB, etc.)
resolve_partner <- function(partner_id) {
  if (is.na(partner_id) || partner_id == "") return(NA_character_)

  # Try as UniProt ID first
  gene <- uniprot_to_gene[partner_id]
  if (!is.na(gene)) return(gene)

  # Try as complex name
  if (partner_id %in% names(complex_to_genes)) {
    subunits <- complex_to_genes[[partner_id]]
    if (length(subunits) > 0) return(paste(subunits, collapse = "_"))
  }

  # Fallback: return as-is (non-standard IDs like HLAA, HLAB)
  return(partner_id)
}

# ==============================================================================
# Helper: derive interaction_type from partner properties
# ==============================================================================

get_interaction_type <- function(partner_a, partner_b, directionality) {
  # Use directionality as primary signal
  if (directionality == "Gap-Gap") return("Cell-Cell Contact")
  if (directionality == "Adhesion-Adhesion") return("Cell-Cell Contact")

  # Check protein/complex properties for the partners
  props_a <- protein_props[[partner_a]]
  if (is.null(props_a)) props_a <- complex_props[[partner_a]]
  props_b <- protein_props[[partner_b]]
  if (is.null(props_b)) props_b <- complex_props[[partner_b]]

  # If either partner is secreted -> Secreted Signaling
  if (!is.null(props_a) && isTRUE(props_a$secreted)) return("Secreted Signaling")
  if (!is.null(props_b) && isTRUE(props_b$secreted)) return("Secreted Signaling")

  # If both transmembrane -> Cell-Cell Contact
  if (!is.null(props_a) && !is.null(props_b)) {
    if (isTRUE(props_a$transmembrane) && isTRUE(props_b$transmembrane)) {
      return("Cell-Cell Contact")
    }
  }

  # If integrin involved -> ECM-Receptor
  if (!is.null(props_a) && isTRUE(props_a$integrin)) return("ECM-Receptor")
  if (!is.null(props_b) && isTRUE(props_b$integrin)) return("ECM-Receptor")

  # Default based on directionality
  if (directionality == "Ligand-Receptor") return("Secreted Signaling")

  return(NA_character_)
}

# ==============================================================================
# Load and convert interactions
# ==============================================================================

cat("[CellPhoneDB] Loading interactions...\n")

interactions <- read.csv(file.path(INPUT_DIR, "interaction_input.csv"),
                          stringsAsFactors = FALSE)
cat(sprintf("  %d raw interactions loaded\n", nrow(interactions)))

cat("[CellPhoneDB] Resolving partners to gene symbols...\n")

n_resolved <- 0
n_failed <- 0
results <- list()

for (i in seq_len(nrow(interactions))) {
  row <- interactions[i, ]

  # Resolve partners to gene symbols
  ligand_raw  <- row$partner_a
  receptor_raw <- row$partner_b
  ligand  <- resolve_partner(ligand_raw)
  receptor <- resolve_partner(receptor_raw)

  # Skip if either partner couldn't be resolved
  if (is.na(ligand) || is.na(receptor)) {
    n_failed <- n_failed + 1
    next
  }

  # For Adhesion-Adhesion, Ligand-Ligand, Gap-Gap: partner_a and partner_b
  # have symmetric roles. We keep partner_a as "ligand" and partner_b as
  # "receptor" for consistency, with the directionality column to distinguish.
  dir <- row$directionality

  # Derive interaction_type from protein properties + directionality
  int_type <- get_interaction_type(ligand_raw, receptor_raw, dir)

  # Use classification as annotation (rich pathway names)
  annotation <- row$classification
  if (is.na(annotation) || annotation == "") annotation <- NA_character_

  # Build evidence from source column
  evidence <- row$source
  if (is.na(evidence) || evidence == "") evidence <- NA_character_

  # CellChat uses underscore for complex subunit separators
  # CellPhoneDB uses + in interactors column — we already resolved via
  # complex_to_genes which joins with underscore

  # Build interaction name
  int_name <- paste(ligand, receptor, sep = "_")

  # Build human-readable name
  int_name_2 <- if (grepl("_", receptor)) {
    paste0(ligand, " - (", gsub("_", "+", receptor), ")")
  } else {
    paste(ligand, "-", receptor)
  }

  results[[length(results) + 1]] <- data.frame(
    interaction_name   = int_name,
    interaction_name_2 = int_name_2,
    ligand             = ligand,
    receptor           = receptor,
    species            = "human",
    source_db          = SOURCE_NAME,
    pathway_name       = NA_character_,  # CellPhoneDB has no pathway grouping
    interaction_type   = int_type,
    annotation         = annotation,
    evidence           = evidence,
    directionality     = dir,
    modulatory_effect  = ifelse(row$modulatory_effect == "", NA_character_, row$modulatory_effect),
    cpdb_version       = row$version,
    stringsAsFactors   = FALSE
  )

  n_resolved <- n_resolved + 1
}

cat(sprintf("  Resolved: %d | Failed: %d\n", n_resolved, n_failed))

# Combine all results
converted <- do.call(rbind, results)

# Remove duplicate interaction_names (can happen if two complexes resolve
# to the same gene combination)
n_before <- nrow(converted)
converted <- converted[!duplicated(converted$interaction_name), ]
n_after <- nrow(converted)
if (n_before != n_after) {
  cat(sprintf("  Removed %d duplicate interaction_names after resolution\n", n_before - n_after))
}

# ==============================================================================
# Save
# ==============================================================================

out_path <- file.path(OUTPUT_DIR, "cellphonedb_human.rds")
saveRDS(converted, out_path)

cat(sprintf("\n[CellPhoneDB] Conversion complete!\n"))
cat(sprintf("  Output: %s\n", out_path))
cat(sprintf("  Interactions: %d\n", nrow(converted)))
cat(sprintf("  Unique ligands: %d\n", length(unique(converted$ligand))))
cat(sprintf("  Unique receptors: %d\n", length(unique(converted$receptor))))

# Stats
cat("\n  Interaction type distribution:\n")
print(table(converted$interaction_type, useNA = "always"))

cat("\n  Directionality distribution:\n")
print(table(converted$directionality, useNA = "always"))

cat("\n  Top 15 annotations:\n")
ann_table <- sort(table(converted$annotation), decreasing = TRUE)
print(head(ann_table, 15))

cat("\n  Sample rows:\n")
print(head(converted[, c("interaction_name", "ligand", "receptor", "interaction_type", "annotation")], 10))
