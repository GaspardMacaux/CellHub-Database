################################################################################
#         verify_integrity.R — Deep integrity check for CellHub Database               #
#                                                                              #
# Compares converted/ (pre-curation) vs databases/ (post-curation) to verify   #
# that cleaning and merging did not introduce data corruption:                 #
#   - Column shifts (ligand in receptor column, etc.)                          #
#   - L-R pair consistency (same pair, same direction after curation)           #
#   - Interaction name matches ligand_receptor                                 #
#   - Merged entries traceable back to individual source files                  #
#   - No phantom interactions (in merged but not in any individual)            #
#   - No lost interactions (in individual but not in merged, unexplained)       #
################################################################################

CONVERTED_DIR <- "data/converted"
DATABASES_DIR <- "data/databases"

cat("============================================================\n")
cat("  CellHub Database — Deep Integrity Verification\n")
cat("============================================================\n\n")

issues_found <- 0

# ==============================================================================
# CHECK 1: interaction_name == ligand_receptor in all files
# ==============================================================================

cat("=== CHECK 1: interaction_name consistency ===\n")

all_db_files <- list.files(DATABASES_DIR, pattern = "\\.rds$",
                           recursive = TRUE, full.names = TRUE)

for (f in all_db_files) {
  db <- readRDS(f)
  inter <- db$interaction
  # Rebuild expected interaction_name from ligand + receptor
  expected <- paste(inter$ligand, inter$receptor, sep = "_")
  mismatches <- which(inter$interaction_name != expected)
  if (length(mismatches) > 0) {
    cat(sprintf("  !! %s: %d name mismatches\n", basename(f), length(mismatches)))
    for (i in head(mismatches, 3)) {
      cat(sprintf("    Row %d: name='%s' but ligand='%s' receptor='%s' -> expected='%s'\n",
                  i, inter$interaction_name[i], inter$ligand[i], inter$receptor[i], expected[i]))
    }
    issues_found <- issues_found + length(mismatches)
  }
}
cat(sprintf("  Result: %d mismatches found\n\n",
            sum(sapply(all_db_files, function(f) {
              db <- readRDS(f); inter <- db$interaction
              sum(inter$interaction_name != paste(inter$ligand, inter$receptor, sep = "_"))
            }))))

# ==============================================================================
# CHECK 2: No column shifts — ligand and receptor should be gene symbols,
# not pathway names, evidence, or other metadata
# ==============================================================================

cat("=== CHECK 2: Column content sanity ===\n")

for (f in all_db_files) {
  db <- readRDS(f)
  inter <- db$interaction
  fname <- sub(paste0(DATABASES_DIR, "/"), "", f)
  
  # Ligand should not contain spaces (except complex names from CellChat)
  lig_with_spaces <- inter$ligand[grepl(" ", inter$ligand) & !is.na(inter$ligand)]
  # Filter out known CellChat complex names
  known_complex_patterns <- c("ACTIVIN", "INHIBIN", "COMPLEX", "complex")
  lig_suspicious <- lig_with_spaces[!grepl(paste(known_complex_patterns, collapse = "|"), lig_with_spaces)]
  if (length(lig_suspicious) > 0) {
    cat(sprintf("  !! %s: %d ligands with unexpected spaces\n", fname, length(lig_suspicious)))
    cat(sprintf("    Examples: %s\n", paste(head(unique(lig_suspicious), 5), collapse = ", ")))
    issues_found <- issues_found + length(lig_suspicious)
  }
  
  # Receptor should not contain PMID, KEGG, or URL patterns
  rec_metadata <- inter$receptor[grepl("PMID|KEGG|http|doi|10\\.", inter$receptor) & !is.na(inter$receptor)]
  if (length(rec_metadata) > 0) {
    cat(sprintf("  !! %s: %d receptors look like metadata\n", fname, length(rec_metadata)))
    cat(sprintf("    Examples: %s\n", paste(head(unique(rec_metadata), 3), collapse = ", ")))
    issues_found <- issues_found + length(rec_metadata)
  }
  
  # Species should be a single consistent value per file
  species_vals <- unique(inter$species)
  if (length(species_vals) > 1) {
    cat(sprintf("  !! %s: multiple species in one file: %s\n", fname, paste(species_vals, collapse = ", ")))
    issues_found <- issues_found + 1
  }
  
  # source_db should not be NA
  na_source <- sum(is.na(inter$source_db))
  if (na_source > 0) {
    cat(sprintf("  !! %s: %d rows with NA source_db\n", fname, na_source))
    issues_found <- issues_found + na_source
  }
}
cat(sprintf("  Result: column content check complete\n\n"))

# ==============================================================================
# CHECK 3: Individual source files — every interaction should exist in merged
# ==============================================================================

cat("=== CHECK 3: Individual -> Merged traceability ===\n")

species_dirs <- list.dirs(DATABASES_DIR, recursive = FALSE)

for (sp_dir in species_dirs) {
  sp <- basename(sp_dir)
  merged_file <- file.path(sp_dir, sprintf("CellHubDB_%s_merged.rds", sp))
  if (!file.exists(merged_file)) next
  
  merged_db <- readRDS(merged_file)
  merged <- merged_db$interaction
  
  # Build canonical pairs from merged (case-insensitive, order-insensitive)
  merged_canonical <- vapply(seq_len(nrow(merged)), function(i) {
    pair <- sort(c(tolower(merged$ligand[i]), tolower(merged$receptor[i])))
    paste(pair, collapse = "|||")
  }, character(1))
  
  individual_files <- list.files(sp_dir, pattern = "\\.rds$", full.names = TRUE)
  individual_files <- individual_files[!grepl("_merged\\.rds$", individual_files)]
  
  for (f in individual_files) {
    db <- readRDS(f)
    ind <- db$interaction
    src_name <- sub(paste0("CellHubDB_", sp, "_"), "", sub("\\.rds$", "", basename(f)))
    
    # Build canonical pairs from individual
    ind_canonical <- vapply(seq_len(nrow(ind)), function(i) {
      pair <- sort(c(tolower(ind$ligand[i]), tolower(ind$receptor[i])))
      paste(pair, collapse = "|||")
    }, character(1))
    
    # Every individual pair should be in merged (canonical match)
    missing_in_merged <- ind_canonical[!ind_canonical %in% merged_canonical]
    if (length(missing_in_merged) > 0) {
      cat(sprintf("  !! %s / %s: %d pairs from individual NOT in merged\n",
                  sp, src_name, length(missing_in_merged)))
      # Show examples
      missing_idx <- which(!ind_canonical %in% merged_canonical)
      for (i in head(missing_idx, 3)) {
        cat(sprintf("    %s -> %s\n", ind$ligand[i], ind$receptor[i]))
      }
      issues_found <- issues_found + length(missing_in_merged)
    } else {
      cat(sprintf("  OK %s / %s: all %d pairs found in merged\n", sp, src_name, nrow(ind)))
    }
  }
}
cat("\n")

# ==============================================================================
# CHECK 4: Merged — no phantom interactions (every pair traceable to a source)
# ==============================================================================

cat("=== CHECK 4: Merged -> Individual traceability ===\n")

for (sp_dir in species_dirs) {
  sp <- basename(sp_dir)
  merged_file <- file.path(sp_dir, sprintf("CellHubDB_%s_merged.rds", sp))
  if (!file.exists(merged_file)) next
  
  merged_db <- readRDS(merged_file)
  merged <- merged_db$interaction
  
  # Collect all canonical pairs from all individual source files
  individual_files <- list.files(sp_dir, pattern = "\\.rds$", full.names = TRUE)
  individual_files <- individual_files[!grepl("_merged\\.rds$", individual_files)]
  
  all_individual_canonical <- c()
  for (f in individual_files) {
    db <- readRDS(f)
    ind <- db$interaction
    ind_canonical <- vapply(seq_len(nrow(ind)), function(i) {
      pair <- sort(c(tolower(ind$ligand[i]), tolower(ind$receptor[i])))
      paste(pair, collapse = "|||")
    }, character(1))
    all_individual_canonical <- c(all_individual_canonical, ind_canonical)
  }
  all_individual_canonical <- unique(all_individual_canonical)
  
  # Check merged against all individuals
  merged_canonical <- vapply(seq_len(nrow(merged)), function(i) {
    pair <- sort(c(tolower(merged$ligand[i]), tolower(merged$receptor[i])))
    paste(pair, collapse = "|||")
  }, character(1))
  
  phantom <- merged_canonical[!merged_canonical %in% all_individual_canonical]
  if (length(phantom) > 0) {
    cat(sprintf("  !! %s: %d phantom pairs in merged (not in any individual)\n", sp, length(phantom)))
    phantom_idx <- which(!merged_canonical %in% all_individual_canonical)
    for (i in head(phantom_idx, 3)) {
      cat(sprintf("    %s -> %s (source: %s)\n",
                  merged$ligand[i], merged$receptor[i], merged$source_db[i]))
    }
    issues_found <- issues_found + length(phantom)
  } else {
    cat(sprintf("  OK %s: all %d merged pairs traceable to source files\n", sp, nrow(merged)))
  }
}
cat("\n")

# ==============================================================================
# CHECK 5: Converted vs Database — spot check random samples
# ==============================================================================

cat("=== CHECK 5: Converted vs Database spot check ===\n")

converted_files <- list.files(CONVERTED_DIR, pattern = "\\.rds$", full.names = TRUE)

for (f in converted_files) {
  conv <- readRDS(f)
  if (!is.data.frame(conv)) next
  if (nrow(conv) == 0) next
  
  src_name <- unique(conv$source_db)[1]
  species_in_file <- unique(conv$species)
  
  for (sp in species_in_file) {
    sp_conv <- conv[conv$species == sp, ]
    if (nrow(sp_conv) == 0) next
    
    # Find matching database file
    db_file <- file.path(DATABASES_DIR, sp, sprintf("CellHubDB_%s_%s.rds",
                                                    sp, gsub("[^A-Za-z0-9_]", "_", src_name)))
    if (!file.exists(db_file)) next
    
    db <- readRDS(db_file)
    db_inter <- db$interaction
    
    # Sample 10 random rows from converted and verify they exist in database
    sample_idx <- sample(seq_len(nrow(sp_conv)), min(10, nrow(sp_conv)))
    n_found <- 0
    n_mismatch <- 0
    
    for (i in sample_idx) {
      lig <- sp_conv$ligand[i]
      rec <- sp_conv$receptor[i]
      
      # Find in database (case-insensitive)
      match_idx <- which(tolower(db_inter$ligand) == tolower(lig) &
                           tolower(db_inter$receptor) == tolower(rec))
      
      if (length(match_idx) > 0) {
        n_found <- n_found + 1
        # Verify key fields match
        db_row <- db_inter[match_idx[1], ]
        if (!is.na(sp_conv$source_db[i]) && !is.na(db_row$source_db)) {
          if (sp_conv$source_db[i] != db_row$source_db) {
            n_mismatch <- n_mismatch + 1
          }
        }
      }
    }
    
    if (n_found < length(sample_idx)) {
      missing <- length(sample_idx) - n_found
      cat(sprintf("  !! %s / %s: %d / %d sampled pairs NOT found in database file\n",
                  sp, src_name, missing, length(sample_idx)))
      issues_found <- issues_found + missing
    } else {
      cat(sprintf("  OK %s / %s: %d / %d sampled pairs verified\n",
                  sp, src_name, n_found, length(sample_idx)))
    }
  }
}
cat("\n")

# ==============================================================================
# CHECK 6: Multi-source merged entries — sources column consistency
# ==============================================================================

cat("=== CHECK 6: Multi-source entries validation ===\n")

for (sp in c("human", "mouse")) {
  merged_file <- file.path(DATABASES_DIR, sp, sprintf("CellHubDB_%s_merged.rds", sp))
  if (!file.exists(merged_file)) next
  
  merged <- readRDS(merged_file)$interaction
  
  if (!"sources" %in% colnames(merged)) next
  
  multi <- merged[grepl(",", merged$sources), ]
  if (nrow(multi) == 0) next
  
  # For each multi-source entry, verify:
  #   1. source_db (winning source) is listed in sources
  #   2. All listed sources have the species in their individual files
  
  n_winner_missing <- 0
  n_source_invalid <- 0
  
  for (i in seq_len(min(nrow(multi), 200))) {
    winner <- multi$source_db[i]
    all_sources <- trimws(strsplit(multi$sources[i], ",")[[1]])
    
    # Winner should be in sources list
    if (!winner %in% all_sources) {
      n_winner_missing <- n_winner_missing + 1
    }
    
    # Each source should have an individual file for this species
    for (src in all_sources) {
      safe_src <- gsub("[^A-Za-z0-9_]", "_", src)
      src_file <- file.path(DATABASES_DIR, sp, sprintf("CellHubDB_%s_%s.rds", sp, safe_src))
      if (!file.exists(src_file)) {
        n_source_invalid <- n_source_invalid + 1
      }
    }
  }
  
  cat(sprintf("  %s: checked %d multi-source entries\n", sp, min(nrow(multi), 200)))
  if (n_winner_missing > 0) {
    cat(sprintf("    !! %d entries where winning source not in sources list\n", n_winner_missing))
    issues_found <- issues_found + n_winner_missing
  }
  if (n_source_invalid > 0) {
    cat(sprintf("    !! %d source references with no matching individual file\n", n_source_invalid))
    issues_found <- issues_found + n_source_invalid
  }
  if (n_winner_missing == 0 && n_source_invalid == 0) {
    cat(sprintf("    OK — all winners in sources, all sources have files\n"))
  }
}
cat("\n")

# ==============================================================================
# CHECK 7: Rowname integrity across all files
# ==============================================================================

cat("=== CHECK 7: Rowname integrity ===\n")

for (f in all_db_files) {
  db <- readRDS(f)
  inter <- db$interaction
  fname <- sub(paste0(DATABASES_DIR, "/"), "", f)
  
  # Rownames should equal interaction_name
  if (!all(rownames(inter) == inter$interaction_name)) {
    cat(sprintf("  !! %s: rownames don't match interaction_name\n", fname))
    mismatches <- which(rownames(inter) != inter$interaction_name)
    cat(sprintf("    %d mismatches. First: rowname='%s', name='%s'\n",
                length(mismatches), rownames(inter)[mismatches[1]], inter$interaction_name[mismatches[1]]))
    issues_found <- issues_found + 1
  }
  
  # No duplicate rownames
  dup_rn <- sum(duplicated(rownames(inter)))
  if (dup_rn > 0) {
    cat(sprintf("  !! %s: %d duplicate rownames\n", fname, dup_rn))
    issues_found <- issues_found + dup_rn
  }
}
cat(sprintf("  Result: rowname check complete\n\n"))


# ==============================================================================
# CHECK 8: Exhaustive converted vs database pair comparison
# ==============================================================================
# For each converted file, verify that EVERY L-R pair exists in the
# corresponding database file with the same ligand and receptor
# (after case normalization). Catches: swapped ligand/receptor,
# corrupted gene names, lost interactions.

cat("=== CHECK 8: Exhaustive converted -> database pair comparison ===\n")

converted_files <- list.files(CONVERTED_DIR, pattern = "\\.rds$", full.names = TRUE)

for (f in converted_files) {
  conv <- readRDS(f)
  if (!is.data.frame(conv) || nrow(conv) == 0) next
  
  species_in_file <- unique(conv$species)
  src_name <- unique(conv$source_db)[1]
  
  for (sp in species_in_file) {
    sp_conv <- conv[conv$species == sp, ]
    if (nrow(sp_conv) == 0) next
    
    safe_src <- gsub("[^A-Za-z0-9_]", "_", src_name)
    db_file <- file.path(DATABASES_DIR, sp, sprintf("CellHubDB_%s_%s.rds", sp, safe_src))
    if (!file.exists(db_file)) {
      cat(sprintf("  SKIP %s / %s: no database file found\n", sp, src_name))
      next
    }
    
    db <- readRDS(db_file)
    db_inter <- db$interaction
    
    # Build canonical pairs for both (case-insensitive)
    conv_pairs <- paste(tolower(sp_conv$ligand), tolower(sp_conv$receptor), sep = "|||")
    db_pairs <- paste(tolower(db_inter$ligand), tolower(db_inter$receptor), sep = "|||")
    
    # Check: every converted pair should be in database
    missing_in_db <- conv_pairs[!conv_pairs %in% db_pairs]
    n_missing <- length(unique(missing_in_db))
    
    # Check: every database pair should be in converted
    extra_in_db <- db_pairs[!db_pairs %in% conv_pairs]
    n_extra <- length(unique(extra_in_db))
    
    # Check: same ligand-receptor direction (not swapped)
    # For pairs present in both, verify ligand stayed as ligand
    common_pairs <- intersect(conv_pairs, db_pairs)
    n_swapped <- 0
    for (pair in head(unique(common_pairs), 500)) {
      parts <- strsplit(pair, "\\|\\|\\|")[[1]]
      conv_idx <- which(conv_pairs == pair)[1]
      db_idx <- which(db_pairs == pair)[1]
      # Compare original case ligand/receptor
      if (tolower(sp_conv$ligand[conv_idx]) != tolower(db_inter$ligand[db_idx])) {
        n_swapped <- n_swapped + 1
      }
    }
    
    # Report
    if (n_missing == 0 && n_extra == 0 && n_swapped == 0) {
      cat(sprintf("  OK %s / %s: %d pairs — all match, none swapped\n",
                  sp, src_name, nrow(sp_conv)))
    } else {
      if (n_missing > 0) {
        cat(sprintf("  !! %s / %s: %d pairs in converted NOT in database\n",
                    sp, src_name, n_missing))
        missing_idx <- which(!conv_pairs %in% db_pairs)
        for (i in head(missing_idx, 3)) {
          cat(sprintf("    MISSING: %s -> %s\n", sp_conv$ligand[i], sp_conv$receptor[i]))
        }
        issues_found <- issues_found + n_missing
      }
      if (n_extra > 0) {
        cat(sprintf("  !! %s / %s: %d pairs in database NOT in converted (phantom)\n",
                    sp, src_name, n_extra))
        issues_found <- issues_found + n_extra
      }
      if (n_swapped > 0) {
        cat(sprintf("  !! %s / %s: %d pairs with SWAPPED ligand/receptor\n",
                    sp, src_name, n_swapped))
        issues_found <- issues_found + n_swapped
      }
    }
  }
}
cat("\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("============================================================\n")
cat("  VERIFICATION SUMMARY\n")
cat("============================================================\n\n")

if (issues_found == 0) {
  cat("  ALL CHECKS PASSED — no data corruption detected.\n")
  cat("  The cleaning and merging pipeline preserved data integrity.\n")
} else {
  cat(sprintf("  !! %d ISSUES FOUND — review the details above.\n", issues_found))
}

cat("\n============================================================\n")
cat("  Verification complete.\n")
cat("============================================================\n")