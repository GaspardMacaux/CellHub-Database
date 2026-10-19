################################################################################
#              analyze_converted.R — Inspect all converted source files        #
#                                                                              #
# Reads every RDS in converted/ and reports:                                   #
#   - File size, number of rows/columns                                        #
#   - Column names and types                                                   #
#   - Mandatory columns present/missing                                        #
#   - NA counts per column                                                     #
#   - Unique values for key columns (species, source_db, interaction_type)     #
#   - Sample rows                                                              #
#   - Cross-file consistency checks                                            #
################################################################################

CONVERTED_DIR <- "converted"

cat("============================================================\n")
cat("  CellHub Database — Converted Files Analyzer\n")
cat("============================================================\n\n")

# Mandatory columns that every converted file MUST have
MANDATORY_COLS <- c("interaction_name", "ligand", "receptor", "species", "source_db")

# Expected columns (not mandatory but should be present for a complete DB)
EXPECTED_COLS <- c("interaction_name", "interaction_name_2", "ligand", "receptor",
                   "species", "source_db", "pathway_name", "interaction_type",
                   "annotation", "evidence")

# Valid interaction types
VALID_TYPES <- c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor",
                 "Non-protein Signaling")

files <- sort(list.files(CONVERTED_DIR, pattern = "\\.rds$", full.names = TRUE))
cat(sprintf("Found %d files in %s/\n\n", length(files), CONVERTED_DIR))

# Collect global stats for summary table
global_stats <- data.frame(
  file = character(),
  species = character(),
  source = character(),
  rows = integer(),
  cols = integer(),
  size_kb = numeric(),
  mandatory_ok = logical(),
  na_ligand = integer(),
  na_receptor = integer(),
  na_pathway = integer(),
  na_type = integer(),
  na_evidence = integer(),
  stringsAsFactors = FALSE
)

all_columns_seen <- c()
issues <- list()

for (f in files) {
  fname <- basename(f)
  fsize <- file.info(f)$size

  cat("============================================================\n")
  cat(sprintf("  FILE: %s (%.1f KB)\n", fname, fsize / 1024))
  cat("============================================================\n")

  df <- readRDS(f)

  if (!is.data.frame(df)) {
    cat(sprintf("  !! NOT A DATA.FRAME — class: %s. Skipping.\n\n", paste(class(df), collapse = ", ")))
    issues[[fname]] <- "Not a data.frame"
    next
  }

  # --- Dimensions ---
  cat(sprintf("  Rows: %d | Columns: %d\n", nrow(df), ncol(df)))

  # --- Species and source ---
  species_vals <- if ("species" %in% colnames(df)) unique(df$species) else "MISSING"
  source_vals  <- if ("source_db" %in% colnames(df)) unique(df$source_db) else "MISSING"
  cat(sprintf("  Species: %s\n", paste(species_vals, collapse = ", ")))
  cat(sprintf("  Source:  %s\n", paste(source_vals, collapse = ", ")))

  # --- Mandatory columns check ---
  missing_mandatory <- setdiff(MANDATORY_COLS, colnames(df))
  mandatory_ok <- length(missing_mandatory) == 0
  if (!mandatory_ok) {
    cat(sprintf("  !! MISSING MANDATORY COLUMNS: %s\n", paste(missing_mandatory, collapse = ", ")))
    issues[[fname]] <- c(issues[[fname]], paste("Missing:", paste(missing_mandatory, collapse = ", ")))
  } else {
    cat("  Mandatory columns: OK\n")
  }

  # --- Expected columns check ---
  missing_expected <- setdiff(EXPECTED_COLS, colnames(df))
  if (length(missing_expected) > 0) {
    cat(sprintf("  Missing expected columns: %s\n", paste(missing_expected, collapse = ", ")))
  }

  # --- Extra columns (source-specific) ---
  extra_cols <- setdiff(colnames(df), EXPECTED_COLS)
  if (length(extra_cols) > 0) {
    cat(sprintf("  Extra columns: %s\n", paste(extra_cols, collapse = ", ")))
  }

  # --- Column types ---
  cat("  Column types:\n")
  for (col in colnames(df)) {
    cat(sprintf("    %-30s %s\n", col, paste(class(df[[col]]), collapse = "/")))
  }

  # --- NA counts ---
  cat("  NA/missing counts:\n")
  na_counts <- sapply(colnames(df), function(col) sum(is.na(df[[col]])))
  for (col in names(na_counts)) {
    pct <- 100 * na_counts[col] / nrow(df)
    flag <- if (col %in% MANDATORY_COLS && na_counts[col] > 0) " !! MANDATORY" else ""
    if (na_counts[col] > 0 || col %in% MANDATORY_COLS) {
      cat(sprintf("    %-30s %5d / %d (%5.1f%%)%s\n",
        col, na_counts[col], nrow(df), pct, flag))
    }
  }

  # --- Empty strings (should be 0 after curation) ---
  empty_counts <- sapply(colnames(df), function(col) {
    if (is.character(df[[col]])) sum(df[[col]] == "", na.rm = TRUE) else 0
  })
  total_empty <- sum(empty_counts)
  if (total_empty > 0) {
    cat(sprintf("  !! EMPTY STRINGS FOUND: %d total\n", total_empty))
    for (col in names(empty_counts[empty_counts > 0])) {
      cat(sprintf("    %-30s %d empty strings\n", col, empty_counts[col]))
    }
    issues[[fname]] <- c(issues[[fname]], paste("Empty strings:", total_empty))
  }

  # --- Duplicate interaction_names ---
  if ("interaction_name" %in% colnames(df)) {
    n_dup <- sum(duplicated(df$interaction_name))
    if (n_dup > 0) {
      cat(sprintf("  !! DUPLICATE interaction_name: %d\n", n_dup))
      issues[[fname]] <- c(issues[[fname]], paste("Duplicates:", n_dup))
    }
  }

  # --- interaction_type distribution ---
  if ("interaction_type" %in% colnames(df)) {
    cat("  Interaction type distribution:\n")
    type_table <- table(df$interaction_type, useNA = "always")
    for (t in names(type_table)) {
      label <- if (is.na(t)) "NA" else t
      valid <- if (is.na(t) || t %in% VALID_TYPES) "" else " !! NON-STANDARD"
      cat(sprintf("    %-30s %5d%s\n", label, type_table[t], valid))
    }
  }

  # --- Pathway coverage ---
  if ("pathway_name" %in% colnames(df)) {
    n_pw <- sum(!is.na(df$pathway_name))
    n_unique_pw <- length(unique(df$pathway_name[!is.na(df$pathway_name)]))
    cat(sprintf("  Pathway coverage: %d / %d (%.0f%%) | %d unique pathways\n",
      n_pw, nrow(df), 100 * n_pw / nrow(df), n_unique_pw))
  }

  # --- Evidence coverage ---
  if ("evidence" %in% colnames(df)) {
    n_ev <- sum(!is.na(df$evidence))
    cat(sprintf("  Evidence coverage: %d / %d (%.0f%%)\n",
      n_ev, nrow(df), 100 * n_ev / nrow(df)))
  }

  # --- Sample rows ---
  cat("  Sample rows (first 3):\n")
  display_cols <- intersect(c("interaction_name", "ligand", "receptor",
                               "interaction_type", "pathway_name", "evidence"), colnames(df))
  for (i in seq_len(min(3, nrow(df)))) {
    vals <- sapply(display_cols, function(col) {
      v <- df[[col]][i]
      if (is.na(v)) return("NA")
      if (nchar(v) > 30) return(paste0(substr(v, 1, 27), "..."))
      return(v)
    })
    cat(sprintf("    %s\n", paste(display_cols, "=", vals, collapse = " | ")))
  }

  # --- Collect global stats ---
  global_stats <- rbind(global_stats, data.frame(
    file = fname,
    species = paste(species_vals, collapse = ","),
    source = paste(source_vals, collapse = ","),
    rows = nrow(df),
    cols = ncol(df),
    size_kb = round(fsize / 1024, 1),
    mandatory_ok = mandatory_ok,
    na_ligand = na_counts["ligand"],
    na_receptor = na_counts["receptor"],
    na_pathway = if ("pathway_name" %in% names(na_counts)) na_counts["pathway_name"] else NA,
    na_type = if ("interaction_type" %in% names(na_counts)) na_counts["interaction_type"] else NA,
    na_evidence = if ("evidence" %in% names(na_counts)) na_counts["evidence"] else NA,
    stringsAsFactors = FALSE
  ))

  all_columns_seen <- unique(c(all_columns_seen, colnames(df)))
  cat("\n")
}

# ==============================================================================
# GLOBAL SUMMARY
# ==============================================================================

cat("\n============================================================\n")
cat("  GLOBAL SUMMARY\n")
cat("============================================================\n\n")

# Summary table
cat(sprintf("%-45s %6s %8s %4s %10s %10s %10s\n",
  "File", "Rows", "Species", "Cols", "NA_pathway", "NA_type", "NA_evid"))
cat(paste(rep("-", 120), collapse = ""), "\n")
for (i in seq_len(nrow(global_stats))) {
  r <- global_stats[i, ]
  cat(sprintf("%-45s %6d %8s %4d %9d %10d %10d\n",
    r$file, r$rows, r$species, r$cols,
    ifelse(is.na(r$na_pathway), 0, r$na_pathway),
    ifelse(is.na(r$na_type), 0, r$na_type),
    ifelse(is.na(r$na_evidence), 0, r$na_evidence)))
}

# Totals
cat(paste(rep("-", 120), collapse = ""), "\n")
cat(sprintf("%-45s %6d %8s\n", "TOTAL",
  sum(global_stats$rows), paste0(length(unique(global_stats$species)), " species")))

# All columns seen across all files
cat(sprintf("\nAll columns seen across all files (%d):\n", length(all_columns_seen)))
cat(paste("  ", sort(all_columns_seen), collapse = "\n"), "\n")

# Issues
if (length(issues) > 0) {
  cat("\n!! ISSUES FOUND:\n")
  for (fname in names(issues)) {
    cat(sprintf("  %s: %s\n", fname, paste(issues[[fname]], collapse = "; ")))
  }
} else {
  cat("\nNo issues found.\n")
}

cat("\nAnalysis complete.\n")
