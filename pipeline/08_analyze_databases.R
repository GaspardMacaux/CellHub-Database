################################################################################
#              analyze_databases.R — Final CellHub Database Quality Control            #
#                                                                              #
# Analyzes all files in databases/ (individual + merged) and reports:          #
#   - Per-species summary with source counts and interaction totals            #
#   - Merged vs sum-of-individuals consistency check                           #
#   - Cross-species overlap detection (same gene pair in multiple species)     #
#   - Gene symbol convention check per species                                 #
#   - Interaction type coverage                                                #
#   - Evidence/pathway/annotation coverage rates                               #
#   - Duplicate detection within and across files                              #
#   - Empty string / NA consistency                                            #
#   - Multi-source interaction statistics                                      #
################################################################################
DATABASES_DIR <- "data/databases"

cat("============================================================\n")
cat("  CellHub Database — Final Database Quality Control\n")
cat("============================================================\n\n")

# Expected gene casing conventions
CASING_RULES <- c(
  "mouse"       = "MGI",
  "human"       = "HGNC",
  "zebrafish"   = "ZFIN",
  "rat"         = "MGI",
  "drosophila"  = "FlyBase",
  "celegans"    = "WormBase",
  "arabidopsis" = "TAIR",
  "rice"        = "AsIs",
  "tomato"      = "AsIs",
  "maize"       = "AsIs",
  "soybean"     = "AsIs",
  "chimp"       = "HGNC",
  "macaque"     = "HGNC",
  "marmoset"    = "HGNC",
  "cow"         = "HGNC",
  "dog"         = "HGNC",
  "horse"       = "HGNC",
  "pig"         = "HGNC",
  "sheep"       = "HGNC",
  "chicken"     = "HGNC",
  "frog"        = "AsIs"
)

species_dirs <- list.dirs(DATABASES_DIR, recursive = FALSE)
species_list <- basename(species_dirs)

cat(sprintf("Found %d species: %s\n\n", length(species_list), paste(species_list, collapse = ", ")))

# ==============================================================================
# 1. PER-SPECIES ANALYSIS
# ==============================================================================

all_species_stats <- list()
all_issues <- list()
grand_total <- 0

for (sp in sort(species_list)) {
  sp_dir <- file.path(DATABASES_DIR, sp)
  files <- list.files(sp_dir, pattern = "\\.rds$", full.names = TRUE)
  merged_file <- files[grepl("_merged\\.rds$", files)]
  individual_files <- files[!grepl("_merged\\.rds$", files)]

  cat("============================================================\n")
  cat(sprintf("  SPECIES: %s\n", toupper(sp)))
  cat("============================================================\n")

  # --- Load merged ---
  if (length(merged_file) == 0) {
    cat("  !! NO MERGED FILE FOUND\n\n")
    all_issues[[sp]] <- c(all_issues[[sp]], "No merged file")
    next
  }

  merged_db <- readRDS(merged_file)
  merged <- merged_db$interaction
  n_merged <- nrow(merged)
  grand_total <- grand_total + n_merged

  cat(sprintf("  Merged: %d interactions\n", n_merged))

  # --- Load individuals ---
  individual_counts <- c()
  individual_total <- 0
  for (f in individual_files) {
    db <- readRDS(f)
    src_name <- sub(paste0("CellHubDB_", sp, "_"), "", sub("\\.rds$", "", basename(f)))
    n <- nrow(db$interaction)
    individual_counts[src_name] <- n
    individual_total <- individual_total + n
    cat(sprintf("    %s: %d interactions\n", src_name, n))
  }

  # --- Check: merged should be <= sum of individuals (dedup removes some) ---
  n_dedup <- individual_total - n_merged
  dedup_pct <- if (individual_total > 0) 100 * n_dedup / individual_total else 0
  cat(sprintf("  Sum of individuals: %d | Merged: %d | Deduplicated: %d (%.1f%%)\n",
    individual_total, n_merged, n_dedup, dedup_pct))
  if (n_merged > individual_total) {
    cat("  !! WARNING: merged has MORE interactions than sum of individuals\n")
    all_issues[[sp]] <- c(all_issues[[sp]], "Merged > sum of individuals")
  }

  # --- Check: empty strings ---
  n_empty <- sum(sapply(merged, function(col) {
    if (is.character(col)) sum(col == "", na.rm = TRUE) else 0
  }))
  if (n_empty > 0) {
    cat(sprintf("  !! EMPTY STRINGS: %d found in merged\n", n_empty))
    all_issues[[sp]] <- c(all_issues[[sp]], paste("Empty strings:", n_empty))
  }

  # --- Check: duplicate interaction_names ---
  n_dup <- sum(duplicated(merged$interaction_name))
  if (n_dup > 0) {
    cat(sprintf("  !! DUPLICATE NAMES: %d in merged\n", n_dup))
    all_issues[[sp]] <- c(all_issues[[sp]], paste("Duplicate names:", n_dup))
  }

  # --- Check: rownames consistency ---
  if (!all(rownames(merged) == merged$interaction_name)) {
    cat("  !! ROWNAMES MISMATCH\n")
    all_issues[[sp]] <- c(all_issues[[sp]], "Rownames mismatch")
  }

  # --- Check: NA in mandatory columns ---
  for (col in c("ligand", "receptor", "species", "source_db")) {
    n_na <- sum(is.na(merged[[col]]))
    if (n_na > 0) {
      cat(sprintf("  !! NA in mandatory column '%s': %d\n", col, n_na))
      all_issues[[sp]] <- c(all_issues[[sp]], paste("NA in", col))
    }
  }

  # --- Check: species column consistency ---
  species_vals <- unique(merged$species)
  if (length(species_vals) != 1 || species_vals != sp) {
    cat(sprintf("  !! SPECIES MISMATCH: expected '%s', found '%s'\n",
      sp, paste(species_vals, collapse = ", ")))
    all_issues[[sp]] <- c(all_issues[[sp]], "Species mismatch")
  }

  # --- Gene symbol casing check ---
  convention <- CASING_RULES[sp]
  if (!is.na(convention) && convention %in% c("MGI", "HGNC", "ZFIN")) {
    # Sample 100 ligand symbols for casing check
    sample_genes <- head(unique(merged$ligand), 100)
    # Remove complex names (contain _)
    simple_genes <- sample_genes[!grepl("_", sample_genes)]
    # Remove metabolite-enzyme composites (contain -)
    simple_genes <- simple_genes[!grepl("-", simple_genes)]

    if (length(simple_genes) > 0) {
      if (convention == "MGI") {
        # First uppercase, rest lowercase: Tgfb1
        bad_case <- simple_genes[!grepl("^[A-Z][a-z0-9]+$", simple_genes) &
                                  !grepl("^[A-Z][a-z]+[0-9]+[a-z]*$", simple_genes)]
      } else if (convention == "HGNC") {
        # All uppercase: TGFB1
        bad_case <- simple_genes[!grepl("^[A-Z0-9]+$", simple_genes)]
      } else if (convention == "ZFIN") {
        # All lowercase: tgfb1
        bad_case <- simple_genes[!grepl("^[a-z0-9]+$", simple_genes)]
      }

      if (length(bad_case) > 0) {
        cat(sprintf("  CASING: %d / %d simple gene symbols don't follow %s convention\n",
          length(bad_case), length(simple_genes), convention))
        cat(sprintf("    Examples: %s\n", paste(head(bad_case, 10), collapse = ", ")))
      } else {
        cat(sprintf("  CASING: OK (%s convention)\n", convention))
      }
    }
  }

  # --- Interaction type coverage ---
  type_table <- table(merged$interaction_type, useNA = "always")
  n_typed <- sum(!is.na(merged$interaction_type))
  cat(sprintf("  INTERACTION TYPE: %d / %d (%.0f%%) annotated\n",
    n_typed, n_merged, 100 * n_typed / n_merged))
  for (t in names(type_table)) {
    label <- if (is.na(t)) "NA" else t
    cat(sprintf("    %-25s %5d (%4.1f%%)\n", label, type_table[t],
      100 * type_table[t] / n_merged))
  }

  # --- Pathway coverage ---
  n_pathway <- sum(!is.na(merged$pathway_name))
  n_unique_pw <- length(unique(merged$pathway_name[!is.na(merged$pathway_name)]))
  cat(sprintf("  PATHWAYS: %d / %d (%.0f%%) annotated | %d unique pathways\n",
    n_pathway, n_merged, 100 * n_pathway / n_merged, n_unique_pw))

  # --- Evidence coverage ---
  if ("evidence" %in% colnames(merged)) {
    n_evidence <- sum(!is.na(merged$evidence))
    cat(sprintf("  EVIDENCE: %d / %d (%.0f%%) have references\n",
      n_evidence, n_merged, 100 * n_evidence / n_merged))
  }

  # --- Annotation coverage ---
  if ("annotation" %in% colnames(merged)) {
    n_annot <- sum(!is.na(merged$annotation))
    cat(sprintf("  ANNOTATION: %d / %d (%.0f%%) annotated\n",
      n_annot, n_merged, 100 * n_annot / n_merged))
  }

  # --- Multi-source interactions ---
  if ("sources" %in% colnames(merged)) {
    n_multi <- sum(grepl(",", merged$sources), na.rm = TRUE)
    cat(sprintf("  MULTI-SOURCE: %d / %d (%.0f%%) validated by 2+ databases\n",
      n_multi, n_merged, 100 * n_multi / n_merged))

    if (n_multi > 0) {
      # Count source combinations
      multi_entries <- merged[grepl(",", merged$sources), ]
      source_combos <- sort(table(multi_entries$sources), decreasing = TRUE)
      cat("    Top source combinations:\n")
      for (combo in head(names(source_combos), 5)) {
        cat(sprintf("      %-50s %5d\n", combo, source_combos[combo]))
      }
    }
  }

  # --- Complex / cofactor / geneInfo tables ---
  cat(sprintf("  COMPLEX TABLE: %d entries\n", nrow(merged_db$complex)))
  cat(sprintf("  COFACTOR TABLE: %d entries\n", nrow(merged_db$cofactor)))
  if (!is.null(merged_db$geneInfo)) {
    cat(sprintf("  GENEINFO TABLE: %d genes\n", nrow(merged_db$geneInfo)))

    # Gene coverage: how many ligand/receptor genes are in geneInfo?
    all_genes <- unique(c(
      unlist(strsplit(merged$ligand[!grepl("-", merged$ligand)], "_")),
      unlist(strsplit(merged$receptor[!grepl("-", merged$receptor)], "_"))
    ))
    in_gi <- sum(all_genes %in% merged_db$geneInfo$Symbol)
    cat(sprintf("  GENE COVERAGE: %d / %d (%.1f%%) of interaction genes in geneInfo\n",
      in_gi, length(all_genes), 100 * in_gi / length(all_genes)))
  } else {
    cat("  GENEINFO TABLE: MISSING\n")
    all_issues[[sp]] <- c(all_issues[[sp]], "No geneInfo")
  }

  # Store stats
  all_species_stats[[sp]] <- list(
    n_merged = n_merged,
    n_sources = length(individual_files),
    n_multi = if ("sources" %in% colnames(merged)) sum(grepl(",", merged$sources), na.rm = TRUE) else 0,
    n_typed = n_typed,
    n_pathway = n_pathway,
    n_evidence = if ("evidence" %in% colnames(merged)) sum(!is.na(merged$evidence)) else 0
  )

  cat("\n")
}

# ==============================================================================
# 2. GLOBAL SUMMARY TABLE
# ==============================================================================

cat("\n============================================================\n")
cat("  GLOBAL SUMMARY\n")
cat("============================================================\n\n")

cat(sprintf("%-15s %6s %5s %7s %7s %7s %7s\n",
  "Species", "Merged", "Srcs", "Multi%", "Type%", "Path%", "Evid%"))
cat(paste(rep("-", 70), collapse = ""), "\n")

for (sp in sort(names(all_species_stats))) {
  s <- all_species_stats[[sp]]
  cat(sprintf("%-15s %6d %5d %6.0f%% %6.0f%% %6.0f%% %6.0f%%\n",
    sp, s$n_merged, s$n_sources,
    100 * s$n_multi / s$n_merged,
    100 * s$n_typed / s$n_merged,
    100 * s$n_pathway / s$n_merged,
    100 * s$n_evidence / s$n_merged))
}
cat(paste(rep("-", 70), collapse = ""), "\n")
cat(sprintf("%-15s %6d %5s\n", "TOTAL", grand_total, ""))

# ==============================================================================
# 3. CROSS-SPECIES CONSISTENCY
# ==============================================================================

cat("\n============================================================\n")
cat("  CROSS-SPECIES CHECKS\n")
cat("============================================================\n\n")

# Check if human orthologs exist for mouse interactions (via connectomeDB2025)
if (all(c("human", "mouse") %in% names(all_species_stats))) {
  human_db <- readRDS(file.path(DATABASES_DIR, "human", "CellHubDB_human_merged.rds"))
  mouse_db <- readRDS(file.path(DATABASES_DIR, "mouse", "CellHubDB_mouse_merged.rds"))

  # Canonical pairs (case-insensitive) for overlap detection
  human_pairs <- unique(paste(toupper(human_db$interaction$ligand),
                               toupper(human_db$interaction$receptor), sep = "|||"))
  mouse_pairs <- unique(paste(toupper(mouse_db$interaction$ligand),
                               toupper(mouse_db$interaction$receptor), sep = "|||"))

  overlap <- sum(mouse_pairs %in% human_pairs)
  cat(sprintf("  Human-Mouse ortholog overlap: %d pairs (%.0f%% of mouse in human)\n",
    overlap, 100 * overlap / length(mouse_pairs)))
}

# ==============================================================================
# 4. SOURCE DATABASE CONTRIBUTION
# ==============================================================================

cat("\n============================================================\n")
cat("  SOURCE CONTRIBUTION (human + mouse)\n")
cat("============================================================\n\n")

for (sp in c("human", "mouse")) {
  if (!sp %in% species_list) next
  sp_dir <- file.path(DATABASES_DIR, sp)
  merged_db <- readRDS(file.path(sp_dir, sprintf("CellHubDB_%s_merged.rds", sp)))
  merged <- merged_db$interaction

  cat(sprintf("  %s (%d interactions):\n", toupper(sp), nrow(merged)))

  # Count by source_db (winning source)
  src_table <- sort(table(merged$source_db), decreasing = TRUE)
  for (src in names(src_table)) {
    cat(sprintf("    %-25s %5d (%4.1f%%) — winning source\n",
      src, src_table[src], 100 * src_table[src] / nrow(merged)))
  }

  # Count by presence in sources column (any contribution)
  if ("sources" %in% colnames(merged)) {
    all_sources <- unique(unlist(strsplit(merged$sources, ", ")))
    cat("    Contribution (present in any role):\n")
    for (src in sort(all_sources)) {
      n_present <- sum(grepl(src, merged$sources, fixed = TRUE))
      cat(sprintf("      %-25s %5d (%4.1f%%)\n",
        src, n_present, 100 * n_present / nrow(merged)))
    }
  }
  cat("\n")
}

# ==============================================================================
# 5. ISSUES SUMMARY
# ==============================================================================

cat("============================================================\n")
cat("  ISSUES SUMMARY\n")
cat("============================================================\n\n")

if (length(all_issues) == 0) {
  cat("  No issues found. All checks passed.\n")
} else {
  for (sp in names(all_issues)) {
    cat(sprintf("  %s:\n", toupper(sp)))
    for (issue in all_issues[[sp]]) {
      cat(sprintf("    !! %s\n", issue))
    }
  }
}

cat("\n============================================================\n")
cat("  Analysis complete.\n")
cat("============================================================\n")
