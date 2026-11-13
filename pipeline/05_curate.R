################################################################################
#                     CellHub Database Curation Pipeline v3.0                          #
#                                                                              #
# Multi-species, multi-source architecture.                                    #
#                                                                              #
# Workflow:                                                                    #
#   1. Load all converted source files (one RDS per source, standardized)      #
#   2. Apply common cleaning to all entries                                    #
#   3. Split by species                                                        #
#   4. Within each species: save individual source files + merged              #
#                                                                              #
# Input:  converted/*.rds  (output of convert_*.R scripts)                     #
# Output: databases/<species>/CellHubDB_<species>_<source>.rds                  #
#         databases/<species>/CellHubDB_<species>_merged.rds                    #
################################################################################

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Directory containing converted source files (one RDS per source DB).
# Each file must be a data.frame with at least the mandatory columns:
#   interaction_name, ligand, receptor, species, source_db,
#   pathway_name, interaction_type, annotation, evidence
CONVERTED_DIR <- "data/converted"

# Root output directory. Subdirectories per species will be created.
OUTPUT_DIR <- "data/databases"

# Path to geneInfo tables (one per species). Each is a data.frame with at
# least: Symbol, EntryID.uniprot. Used for gene symbol normalization.
GENEINFO_DIR <- "data/geneinfo"

# Source priority for dedup during merge: higher number = wins on conflicts.
# When the same L-R pair exists in multiple sources, the highest-priority
# source determines directionality. Metadata (evidence, annotation) from
# lower-priority sources is still merged into the winning row.
SOURCE_PRIORITY <- c(
  "MultiNicheNet"       = 1,
  "CellTalkDB"          = 2,
  "cell2cell_Armingol"  = 3,
  "PlantPhoneDB"        = 3,
  "FlyPhoneDB2"         = 3,
  "MRCLinkdb"           = 4,
  "PlantCellChatDB"     = 4,
  "SingleCellSignalR"   = 5,
  "ICELLNET"            = 5,
  "NeuronChatDB"        = 6,
  "CellPhoneDB"         = 7,
  "connectomeDB2025"    = 8,
  "CellChatDB_v1"       = 9,
  "CellChatDB_v2"       = 10
)

# Gene symbol casing convention per species.
# "MGI"   = first letter uppercase, rest lowercase (Tgfb1)
# "HGNC"  = all uppercase (TGFB1)
# "ZFIN"  = all lowercase (tgfb1)
# "FlyBase" = kept as-is (complex rules: dpp vs Toll)
# "WormBase" = kept as-is (daf-2, lin-12)
# "TAIR"  = all uppercase (CLV3)
SPECIES_CASING <- c(
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

# ==============================================================================
# HELPERS
# ==============================================================================

# Print a timestamped log message tagged with the current step name.
log_step <- function(step_name, msg) {
  timestamp <- format(Sys.time(), "%H:%M:%S")
  cat(sprintf("[%s] STEP %s | %s\n", timestamp, step_name, msg))
}

# Print total interaction count, unique L-R pairs, and per-source breakdown.
log_db_state <- function(df, step_name) {
  log_step(step_name, sprintf("State: %d interactions | %d unique L-R pairs | %d species | %d sources",
    nrow(df),
    length(unique(paste(df$ligand, df$receptor, sep = "->"))),
    length(unique(df$species)),
    length(unique(df$source_db))
  ))
  # Per-species breakdown
  sp_counts <- table(df$species)
  for (sp in names(sp_counts)) {
    src_in_sp <- unique(df$source_db[df$species == sp])
    log_step(step_name, sprintf("  %s: %d interactions (%d sources: %s)",
      sp, sp_counts[sp], length(src_in_sp), paste(src_in_sp, collapse = ", ")))
  }
}

# Print N randomly sampled rows for visual spot-checking.
log_sample <- function(df, step_name, n = 3) {
  cols <- c("interaction_name", "ligand", "receptor", "species", "source_db", "interaction_type")
  cols <- cols[cols %in% colnames(df)]
  idx <- sample(seq_len(nrow(df)), min(n, nrow(df)))
  log_step(step_name, "Sample rows:")
  for (i in idx) {
    cat(sprintf("    %s\n", paste(cols, "=", df[i, cols], collapse = " | ")))
  }
}

# Convert a gene symbol to the correct casing for a given species.
# For MGI (mouse): first letter uppercase, rest lowercase per subunit.
# For HGNC (human) / TAIR (plant): all uppercase.
# For ZFIN (zebrafish): all lowercase.
# For FlyBase / WormBase: keep as-is (complex organism-specific rules).
# Only normalizes subunits found in the geneInfo lookup table; unknown
# names are kept unchanged to avoid corrupting aliases or shorthands.
normalize_gene_case <- function(gene_name, lookup, convention) {
  if (is.na(gene_name)) return(NA_character_)
  if (convention %in% c("FlyBase", "WormBase", "AsIs")) return(gene_name)  # Keep as-is
  parts <- strsplit(gene_name, "_")[[1]]
  normalized <- vapply(parts, function(p) {
    # Try to match in lookup (case-insensitive)
    upper_p <- toupper(p)
    if (upper_p %in% names(lookup)) {
      lookup[[upper_p]]  # Use canonical symbol from geneInfo
    } else {
      # Not in lookup — apply convention blindly only for HGNC/ZFIN/TAIR
      # where the rule is simple and universal
      switch(convention,
        "HGNC" = toupper(p),
        "ZFIN" = tolower(p),
        "TAIR" = toupper(p),
        p  # Default: keep as-is
      )
    }
  }, character(1))
  paste(normalized, collapse = "_")
}

# Build a CellChat-compatible DB list from an interaction data.frame.
# CellChat expects: db$interaction (rownames = interaction_name),
# db$complex, db$cofactor, db$geneInfo.
build_cellchat_db <- function(interaction_df, complex_df = NULL,
                               cofactor_df = NULL, geneinfo_df = NULL) {
  rownames(interaction_df) <- interaction_df$interaction_name
  all_names <- unique(c(interaction_df$ligand, interaction_df$receptor))
  if (is.null(complex_df) || nrow(complex_df) == 0) {
    complex_df <- data.frame(subunit_1 = character(0), subunit_2 = character(0),
                             stringsAsFactors = FALSE)
  }
  # Keep only the complexes used by this subset
  complex_df <- complex_df[rownames(complex_df) %in% all_names, , drop = FALSE]
  # Multi-subunit partners without a complex definition (heteromers from CellPhoneDB,
  # connectomeDB2025...) are silently dropped by CellChat: define them from their
  # underscore-separated subunits. Metabolite composites ("-") and names with spaces
  # need a curated definition and are not generated.
  undefined <- all_names[!is.na(all_names) & grepl("_", all_names) & !grepl("[- ]", all_names) &
                           !(all_names %in% rownames(complex_df))]
  if (length(undefined) > 0) {
    parts <- strsplit(undefined, "_")
    n_cols <- max(ncol(complex_df), max(lengths(parts)))
    for (k in seq_len(n_cols)) {
      col <- paste0("subunit_", k)
      if (!col %in% colnames(complex_df)) complex_df[[col]] <- rep("", nrow(complex_df))
    }
    # CellChat expects empty subunit slots as "" (not NA)
    generated <- as.data.frame(t(vapply(parts, function(p) c(p, rep("", n_cols - length(p))),
                                        character(n_cols))), stringsAsFactors = FALSE)
    colnames(generated) <- paste0("subunit_", seq_len(n_cols))
    rownames(generated) <- undefined
    for (col in setdiff(colnames(complex_df), colnames(generated))) generated[[col]] <- ""
    complex_df <- rbind(complex_df, generated[, colnames(complex_df), drop = FALSE])
    cat(sprintf("    complex table: %d curated + %d generated from subunits\n",
                nrow(complex_df) - length(undefined), length(undefined)))
  }
  if (is.null(cofactor_df) || nrow(cofactor_df) == 0) {
    cofactor_df <- data.frame(cofactor1 = character(0), stringsAsFactors = FALSE)
  }
  if (is.null(geneinfo_df)) {
    geneinfo_df <- data.frame(Symbol = character(0), stringsAsFactors = FALSE)
  }
  return(list(
    interaction = interaction_df,
    complex     = complex_df,
    cofactor    = cofactor_df,
    geneInfo    = geneinfo_df
  ))
}

# ==============================================================================
# STEP 0: Load all converted source files
# ==============================================================================

log_step("0-LOAD", sprintf("Scanning %s/ for converted source files...", CONVERTED_DIR))

if (!dir.exists(CONVERTED_DIR)) {
  stop(sprintf("Directory '%s' not found. Run convert_*.R scripts first.", CONVERTED_DIR))
}

source_files <- list.files(CONVERTED_DIR, pattern = "\\.rds$", full.names = TRUE)
log_step("0-LOAD", sprintf("Found %d source files", length(source_files)))

if (length(source_files) == 0) {
  stop("No .rds files found in converted/. Nothing to process.")
}

# Load all sources into a list, then merge with column alignment.
# Different sources may have different columns (connectomeDB2025 has
# connectome_id, evidence_type, etc. that CellChatDB doesn't). We keep
# all columns — missing values become NA.
source_list <- lapply(source_files, function(f) {
  log_step("0-LOAD", sprintf("  Loading %s...", basename(f)))
  df <- readRDS(f)
  # Validate mandatory columns
  mandatory <- c("interaction_name", "ligand", "receptor", "species", "source_db")
  missing <- setdiff(mandatory, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf("File %s missing mandatory columns: %s", basename(f), paste(missing, collapse = ", ")))
  }
  log_step("0-LOAD", sprintf("    -> %d rows | species: %s | source: %s",
    nrow(df), paste(unique(df$species), collapse = ","), paste(unique(df$source_db), collapse = ",")))
  df
})

# Align columns: collect all column names, add missing ones as NA per source
all_colnames <- unique(unlist(lapply(source_list, colnames)))
source_list <- lapply(source_list, function(df) {
  for (col in setdiff(all_colnames, colnames(df))) {
    df[[col]] <- NA
  }
  df[, all_colnames, drop = FALSE]
})
all_interactions <- do.call(rbind, source_list)

log_step("0-LOAD", sprintf("Total loaded: %d interactions", nrow(all_interactions)))
log_db_state(all_interactions, "0-LOAD")

# Load geneInfo tables (one per species, if available)
geneinfo_tables <- list()
if (dir.exists(GENEINFO_DIR)) {
  gi_files <- list.files(GENEINFO_DIR, pattern = "\\.rds$", full.names = TRUE)
  for (f in gi_files) {
    # Filename convention: geneinfo_mouse.rds, geneinfo_human.rds, etc.
    sp <- sub("geneinfo_", "", sub("\\.rds$", "", basename(f)))
    geneinfo_tables[[sp]] <- readRDS(f)
    log_step("0-LOAD", sprintf("  GeneInfo loaded for %s: %d genes", sp, nrow(geneinfo_tables[[sp]])))
  }
}

# Load complex/cofactor tables if available
complex_tables <- list()
cofactor_tables <- list()
if (file.exists(file.path(GENEINFO_DIR, "complexes.rds"))) {
  complex_tables <- readRDS(file.path(GENEINFO_DIR, "complexes.rds"))
  log_step("0-LOAD", sprintf("  Complex tables loaded for: %s", paste(names(complex_tables), collapse = ", ")))
}
if (file.exists(file.path(GENEINFO_DIR, "cofactors.rds"))) {
  cofactor_tables <- readRDS(file.path(GENEINFO_DIR, "cofactors.rds"))
  log_step("0-LOAD", sprintf("  Cofactor tables loaded for: %s", paste(names(cofactor_tables), collapse = ", ")))
}

# ==============================================================================
# STEP 1: Normalize empty strings and whitespace-only to NA
# ==============================================================================

log_step("1-NA", "Normalizing empty strings and whitespace-only to NA...")

count_before <- sum(sapply(all_interactions, function(col) {
  if (is.character(col)) sum(col == "" | grepl("^\\s+$", col), na.rm = TRUE) else 0
}))

all_interactions[] <- lapply(all_interactions, function(col) {
  if (is.character(col)) {
    col[col == "" | grepl("^\\s+$", col)] <- NA_character_
  }
  col
})

log_step("1-NA", sprintf("Cleaned %d empty/whitespace entries", count_before))
# Leading/trailing spaces in identifiers prevent exact matching in CellChat (e.g. "Os07g0173700 ")
all_interactions$ligand   <- trimws(all_interactions$ligand)
all_interactions$receptor <- trimws(all_interactions$receptor)
# Remove rows where ligand or receptor is NA (broken source entries)
n_na_lr <- sum(is.na(all_interactions$ligand) | is.na(all_interactions$receptor))
if (n_na_lr > 0) {
  log_step("1-NA", sprintf("Removing %d interactions with NA ligand or receptor", n_na_lr))
  all_interactions <- all_interactions[!is.na(all_interactions$ligand) & 
                                         !is.na(all_interactions$receptor), ]
}
# ==============================================================================
# STEP 1-ID: Repair species-specific identifier formats
# ==============================================================================
# CellChat reads "_" as a complex separator, so identifiers carrying one must be rewritten
# before complexes are generated from subunits.
log_step("1-ID", "Repairing species-specific identifier formats...")
# Optional two-column file (B73 v4 gene ID, B73 v5 gene ID, no header) to convert legacy maize IDs
MAIZE_V4_TO_V5 <- file.path("data/raw", "PlantCellChatDB", "maize_v4_to_v5.tsv")
maize_map <- character(0)
if (file.exists(MAIZE_V4_TO_V5)) {
  mm <- read.delim(MAIZE_V4_TO_V5, header = FALSE, stringsAsFactors = FALSE)[, 1:2]
  mm <- unique(mm[mm[[1]] != "" & mm[[2]] != "", ])
  # Keep 1:1 correspondences only (split or merged gene models are ambiguous)
  mm <- mm[!(mm[[1]] %in% mm[[1]][duplicated(mm[[1]])]) & !(mm[[2]] %in% mm[[2]][duplicated(mm[[2]])]), ]
  maize_map <- setNames(mm[[2]], mm[[1]])
  log_step("1-ID", sprintf("  maize: %d one-to-one v4 -> v5 correspondences loaded", length(maize_map)))
}
repair_ids <- function(x, sp) {
  if (sp == "maize") {
    x <- sub("^ZEAMMB73_", "", x)  # legacy prefix: ZEAMMB73_Zm00001d048611 -> Zm00001d048611
    hit <- x %in% names(maize_map)
    x[hit] <- maize_map[x[hit]]
  }
  if (sp == "soybean") x <- sub("^GLYMA_", "Glyma.", x)  # Phytozome Wm82 form: Glyma.07G056900
  x
}
for (sp in intersect(c("maize", "soybean"), unique(all_interactions$species))) {
  idx <- which(all_interactions$species == sp)
  for (col in c("ligand", "receptor")) {
    old_ids <- all_interactions[[col]][idx]
    new_ids <- repair_ids(old_ids, sp)
    log_step("1-ID", sprintf("  %s %s: %d identifiers rewritten", sp, col, sum(old_ids != new_ids)))
    all_interactions[[col]][idx] <- new_ids
  }
}
if ("maize" %in% all_interactions$species) {
  mz <- all_interactions$species == "maize"
  v4 <- unique(grep("^Zm00001d", c(all_interactions$ligand[mz], all_interactions$receptor[mz]), value = TRUE))
  log_step("1-ID", sprintf("  maize: %d distinct B73 v4 identifiers have no one-to-one B73 v5 gene model", length(v4)))
  # Interactions keeping a v4 identifier can never be scored against a v5 count matrix or the v5 geneInfo: removed and logged
  drop_v4 <- mz & (all_interactions$ligand %in% v4 | all_interactions$receptor %in% v4)
  if (any(drop_v4)) {
    if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
    removed_v4 <- all_interactions[drop_v4, ]
    removed_v4$removal_reason <- "B73 v4 identifier without a one-to-one B73 v5 gene model (MaizeGDB B73v4_to_B73v5.tsv)"
    write.csv(removed_v4, file.path(OUTPUT_DIR, "removed_maize_v4_identifiers.csv"), row.names = FALSE)
    log_step("1-ID", sprintf("  maize: %d interactions removed (logged in removed_maize_v4_identifiers.csv)", sum(drop_v4)))
    all_interactions <- all_interactions[!drop_v4, ]
  }
}
# Underscore identifiers in species without curated complex tables will become generated complexes
no_cx <- !(all_interactions$species %in% names(complex_tables)) &
  (grepl("_", all_interactions$ligand) | grepl("_", all_interactions$receptor))
if (any(no_cx)) {
  # No curated complexes exist for these species: an underscore identifier is malformed, not a heteromer
  log_step("1-ID", sprintf("  Removing %d entries with malformed '_' identifiers in species without curated complex tables:", sum(no_cx)))
  print(head(all_interactions[no_cx, c("species", "source_db", "ligand", "receptor")], 20), row.names = FALSE)
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  write.csv(all_interactions[no_cx, ], file.path(OUTPUT_DIR, "removed_malformed_identifiers.csv"), row.names = FALSE)
  all_interactions <- all_interactions[!no_cx, ]
}

# ==============================================================================
# STEP 2: Clean pathway_name — remove source names used as pathway
# ==============================================================================

log_step("2-PATHWAY", "Cleaning non-informative pathway_name values...")

# These are database names that were erroneously stored as pathway names
fake_pathways <- c("MultiNiche", "CellTalk", "CellPhoneDB", "MultiNicheNet",
                   "NeuronChat", "connectomeDB", "ICELLNET", "SingleCellSignalR",
                   "MRCLinkdb", "FlyPhoneDB", "PlantPhoneDB", "cell2cell")
n_fake <- sum(all_interactions$pathway_name %in% fake_pathways, na.rm = TRUE)
all_interactions$pathway_name[all_interactions$pathway_name %in% fake_pathways] <- NA_character_

log_step("2-PATHWAY", sprintf("Replaced %d source-as-pathway entries with NA", n_fake))

# ==============================================================================
# STEP 3: Clean annotation — remove generic/uninformative values
# ==============================================================================

log_step("3-ANNOT", "Cleaning generic annotation values...")

generic_annotations <- c("Signaling", "signaling", "interaction", "Interaction")
n_generic <- sum(all_interactions$annotation %in% generic_annotations, na.rm = TRUE)
all_interactions$annotation[all_interactions$annotation %in% generic_annotations] <- NA_character_

log_step("3-ANNOT", sprintf("Replaced %d generic annotations with NA", n_generic))

# ==============================================================================
# STEP 4: Clean interaction_type — standardize and set missing to NA
# ==============================================================================

log_step("4-TYPE", "Standardizing interaction_type values...")

# Map variant spellings to canonical values
type_mapping <- c(
  "Unknown"               = NA_character_,
  "unknown"               = NA_character_,
  "Secreted Signaling"    = "Secreted Signaling",
  "secreted"              = "Secreted Signaling",
  "Cell-Cell Contact"     = "Cell-Cell Contact",
  "cell-cell"             = "Cell-Cell Contact",
  "ECM-Receptor"          = "ECM-Receptor",
  "ecm"                   = "ECM-Receptor",
  "Non-protein Signaling" = "Non-protein Signaling",
  "non-protein"           = "Non-protein Signaling",
  "Metabolite"            = "Non-protein Signaling"
)

for (old_val in names(type_mapping)) {
  match_idx <- which(all_interactions$interaction_type == old_val)
  if (length(match_idx) > 0) {
    all_interactions$interaction_type[match_idx] <- type_mapping[old_val]
  }
}

# Anything not in the canonical set -> NA
valid_types <- c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor", "Non-protein Signaling")
invalid <- !is.na(all_interactions$interaction_type) & !(all_interactions$interaction_type %in% valid_types)
if (sum(invalid) > 0) {
  log_step("4-TYPE", sprintf("WARNING: %d entries with non-standard interaction_type set to NA", sum(invalid)))
  log_step("4-TYPE", sprintf("  Values: %s", paste(unique(all_interactions$interaction_type[invalid]), collapse = ", ")))
  all_interactions$interaction_type[invalid] <- NA_character_
}

log_step("4-TYPE", "Interaction type distribution:")
print(table(all_interactions$interaction_type, useNA = "always"))

# ==============================================================================
# STEP 5: Remove entries with unresolvable UniProt IDs
# ==============================================================================

log_step("5-UNIPROT", "Identifying entries with unresolvable UniProt IDs...")

uniprot_pattern <- "^[A-Z][0-9][A-Z0-9]{3}[0-9]$"

# Many valid identifiers match the UniProt accession pattern (S100A8, S100A9, P2RY12,
# GABA-A subunit codes such as A1B3G2 inside complex names), so a token is flagged
# only if it is neither a known gene symbol (any species with geneInfo) nor part of a complex name
known_symbols <- unique(toupper(unlist(lapply(geneinfo_tables, function(gi) as.character(gi$Symbol)))))
known_complexes <- unique(unlist(lapply(complex_tables, rownames)))
log_step("5-UNIPROT", sprintf("Reference: %d known gene symbols, %d complex names",
                              length(known_symbols), length(known_complexes)))

# Tokens of one ligand or receptor that look like UniProt accessions and are not known identifiers
unresolved_tokens <- function(partner) {
  if (is.na(partner) || partner %in% known_complexes) return(character(0))
  parts <- strsplit(partner, "_")[[1]]
  return(parts[grepl(uniprot_pattern, parts) & !(toupper(parts) %in% known_symbols)])
}

uniprot_tokens <- lapply(seq_len(nrow(all_interactions)), function(i) {
  c(unresolved_tokens(all_interactions$ligand[i]), unresolved_tokens(all_interactions$receptor[i]))
})
has_uniprot <- lengths(uniprot_tokens) > 0

n_uniprot <- sum(has_uniprot)
log_step("5-UNIPROT", sprintf("Entries with unresolvable UniProt IDs: %d", n_uniprot))

if (n_uniprot > 0) {
  removed <- all_interactions[has_uniprot, ]
  # Keep the offending tokens so every removal can be audited
  removed$unresolved_tokens <- vapply(uniprot_tokens[has_uniprot], paste, character(1), collapse = ", ")
  log_step("5-UNIPROT", "By source:")
  print(table(removed$source_db))
  
  removed_path <- file.path(OUTPUT_DIR, "removed_uniprot_entries.csv")
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  write.csv(removed, removed_path, row.names = FALSE)
  log_step("5-UNIPROT", sprintf("Removed entries saved to: %s", removed_path))
  
  all_interactions <- all_interactions[!has_uniprot, ]
}

log_step("5-UNIPROT", sprintf("After removal: %d interactions", nrow(all_interactions)))
# ==============================================================================
# STEP 6: Normalize gene symbol casing per species
# ==============================================================================

log_step("6-CASE", "Normalizing gene symbol casing per species...")

for (sp in unique(all_interactions$species)) {
  convention <- SPECIES_CASING[sp]
  if (is.na(convention)) {
    log_step("6-CASE", sprintf("  %s: no casing convention defined, skipping", sp))
    next
  }
  if (convention %in% c("FlyBase", "WormBase", "AsIs")) {
    log_step("6-CASE", sprintf("  %s: %s convention — keeping as-is", sp, convention))
    next
  }

  sp_idx <- which(all_interactions$species == sp)
  log_step("6-CASE", sprintf("  %s: applying %s convention to %d interactions...", sp, convention, length(sp_idx)))

  # Build lookup from geneInfo if available
  lookup <- character(0)
  if (sp %in% names(geneinfo_tables)) {
    gi <- geneinfo_tables[[sp]]
    lookup <- setNames(gi$Symbol, toupper(gi$Symbol))
  }
  # Complex names (TGFbR1_R2, Activin AB, IL22RA1_IL20RB) must keep the exact spelling of the complex table
  sp_complex_names <- if (sp %in% names(complex_tables)) rownames(complex_tables[[sp]]) else character(0)
  complex_lookup <- setNames(sp_complex_names, toupper(sp_complex_names))

  n_changed <- 0
  for (col in c("ligand", "receptor")) {
    old_vals <- all_interactions[[col]][sp_idx]
    new_vals <- vapply(old_vals, normalize_gene_case, character(1),
                       lookup = lookup, convention = convention)
    is_cx <- !is.na(old_vals) & toupper(old_vals) %in% names(complex_lookup)
    new_vals[is_cx] <- complex_lookup[toupper(old_vals[is_cx])]
    n_changed <- n_changed + sum(old_vals != new_vals, na.rm = TRUE)
    all_interactions[[col]][sp_idx] <- new_vals
  }

  log_step("6-CASE", sprintf("    Corrected: %d gene names", n_changed))
}



# ==============================================================================
# STEP 6B: Validate gene symbols against geneInfo (with alias rescue)
# ==============================================================================
# A ligand/receptor name is valid if it is a known complex name, if all its underscore-separated
# subunits are reference symbols, or if it is a metabolite-enzyme composite (GABA-GAD1) whose
# enzyme is a reference symbol. Failing names are first rescued through an alias table (retired
# aliases such as Il1f6 -> Il36a, symbols newer than the reference such as CRIPTO -> TDGF1).
# Only names that cannot be rescued are removed. Renamed and removed entries are audited.
ALIAS_DIR <- file.path("data/raw", "CellTalkDB")

# Alias (uppercase) -> reference symbol, from a Synonyms column of geneInfo (if any) and the
# NCBI gene_info table shipped with CellTalkDB. Ambiguous aliases are dropped.
build_alias_map <- function(gi, gene_info_path) {
  valid_upper <- setNames(gi$Symbol, toupper(gi$Symbol))
  rows <- list()
  syn_col <- grep("^synonyms?$", colnames(gi), ignore.case = TRUE, value = TRUE)
  if (length(syn_col) > 0) {
    rows[["geneInfo"]] <- data.frame(id = paste0("gi", seq_len(nrow(gi))),
                                     names = paste(gi$Symbol, gi[[syn_col[1]]], sep = "|"), stringsAsFactors = FALSE)
  }
  if (file.exists(gene_info_path)) {
    ni <- read.delim(gene_info_path, quote = "", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
    colnames(ni) <- sub("^#", "", colnames(ni))
    s_col <- grep("^symbol$", colnames(ni), ignore.case = TRUE, value = TRUE)
    y_col <- grep("^synonyms?$", colnames(ni), ignore.case = TRUE, value = TRUE)
    if (length(s_col) > 0 && length(y_col) > 0) {
      rows[["ncbi"]] <- data.frame(id = paste0("ncbi", seq_len(nrow(ni))),
                                   names = paste(ni[[s_col[1]]], ni[[y_col[1]]], sep = "|"), stringsAsFactors = FALSE)
    } else {
      log_step("6B-GENESYM", sprintf("    %s: no Symbol/Synonyms columns (found: %s)", basename(gene_info_path),
                                     paste(head(colnames(ni), 12), collapse = ", ")))
    }
  } else {
    log_step("6B-GENESYM", sprintf("    alias table not found: %s", gene_info_path))
  }
  if (length(rows) == 0) return(character(0))
  tab <- do.call(rbind, rows)
  split_names <- strsplit(tab$names, "[|,;]")
  long <- data.frame(id = rep(tab$id, lengths(split_names)), name = trimws(unlist(split_names)), stringsAsFactors = FALSE)
  long <- long[!long$name %in% c("", "-", "NA"), ]
  long$upper <- toupper(long$name)
  long$official <- unname(valid_upper[long$upper])  # NA unless the name is itself a reference symbol
  # Each record must point to exactly one reference symbol
  target_per_id <- tapply(long$official, long$id, function(v) {
    v <- unique(v[!is.na(v)])
    if (length(v) == 1) v else NA_character_
  })
  long$target <- unname(target_per_id[long$id])
  long <- long[!is.na(long$target) & is.na(long$official), ]
  # Drop aliases that point to several different reference symbols
  n_targets <- tapply(long$target, long$upper, function(v) length(unique(v)))
  long <- long[unname(n_targets[long$upper]) == 1, ]
  long <- long[!duplicated(long$upper), ]
  setNames(long$target, long$upper)
}

log_step("6B-GENESYM", "Validating gene symbols against geneInfo...")

n_before_genesym <- nrow(all_interactions)
flagged <- rep(FALSE, nrow(all_interactions))
unmatched_names <- rep("", nrow(all_interactions))
unmatched_log <- list()
renamed_log <- list()

for (sp in unique(all_interactions$species)) {
  sp_idx <- which(all_interactions$species == sp)
  if (!sp %in% names(geneinfo_tables) || is.null(geneinfo_tables[[sp]])) {
    log_step("6B-GENESYM", sprintf("  %s: no geneInfo available — skipping validation", sp))
    next
  }
  # Only genome-wide references are used for validation: the zebrafish geneInfo shipped with
  # CellChatDB covers part of the genome and would discard valid symbols (itgb1a, fgfr4, epha7...)
  if (!sp %in% c("human", "mouse")) {
    log_step("6B-GENESYM", sprintf("  %s: geneInfo is not genome-wide — skipping validation", sp))
    next
  }
  valid_symbols <- geneinfo_tables[[sp]]$Symbol
  valid_complexes <- character(0)
  if (sp %in% names(complex_tables) && !is.null(complex_tables[[sp]])) valid_complexes <- rownames(complex_tables[[sp]])
  valid_complexes_upper <- toupper(valid_complexes)
  alias_map <- build_alias_map(geneinfo_tables[[sp]], file.path(ALIAS_DIR, sprintf("%s_gene_info.txt", sp)))
  log_step("6B-GENESYM", sprintf("  %s: %d aliases available for rescue", sp, length(alias_map)))

  # Returns the valid (possibly renamed) name, or NA when it cannot be resolved
  resolve_name <- function(name) {
    if (toupper(name) %in% valid_complexes_upper) return(name)
    parts <- strsplit(name, "_")[[1]]
    fixed <- ifelse(parts %in% valid_symbols, parts, unname(alias_map[toupper(parts)]))
    if (!anyNA(fixed)) return(paste(fixed, collapse = "_"))
    # Metabolite-enzyme composite: the enzyme after the last hyphen is the expression proxy
    if (grepl("-", name)) {
      h <- strsplit(name, "-")[[1]]
      enzyme <- h[length(h)]
      if (enzyme %in% valid_symbols) return(name)
      if (!is.na(alias_map[toupper(enzyme)])) return(paste(c(h[-length(h)], alias_map[[toupper(enzyme)]]), collapse = "-"))
    }
    NA_character_
  }

  for (i in sp_idx) {
    lig <- all_interactions$ligand[i]
    rec <- all_interactions$receptor[i]
    new_lig <- resolve_name(lig)
    new_rec <- resolve_name(rec)
    if (is.na(new_lig) || is.na(new_rec)) {
      flagged[i] <- TRUE
      bad <- c(if (is.na(new_lig)) lig, if (is.na(new_rec)) rec)
      unmatched_names[i] <- paste(unique(bad), collapse = ", ")
      for (name in bad) unmatched_log[[name]] <- c(unmatched_log[[name]], sp)
    } else if (new_lig != lig || new_rec != rec) {
      renamed_log[[length(renamed_log) + 1]] <- data.frame(species = sp, source_db = all_interactions$source_db[i],
        old_ligand = lig, old_receptor = rec, ligand = new_lig, receptor = new_rec, stringsAsFactors = FALSE)
      all_interactions$ligand[i] <- new_lig
      all_interactions$receptor[i] <- new_rec
    }
  }
  log_step("6B-GENESYM", sprintf("  %s: %d / %d interactions flagged (%d gene symbols + %d complex names in reference)",
                                 sp, sum(flagged[sp_idx]), length(sp_idx), length(valid_symbols), length(valid_complexes)))
}

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Audit renamed entries
if (length(renamed_log) > 0) {
  renamed_df <- do.call(rbind, renamed_log)
  renamed_path <- file.path(OUTPUT_DIR, "renamed_alias_genesymbols.csv")
  write.csv(renamed_df, renamed_path, row.names = FALSE)
  log_step("6B-GENESYM", sprintf("  Rescued through aliases: %d interactions -> %s", nrow(renamed_df), renamed_path))
  print(table(renamed_df$source_db))
}

# Report most common unmatched symbols
if (length(unmatched_log) > 0) {
  unmatched_counts <- sort(sapply(unmatched_log, length), decreasing = TRUE)
  log_step("6B-GENESYM", sprintf("  Unique unmatched symbols: %d", length(unmatched_counts)))
  log_step("6B-GENESYM", "  Top 30 unmatched symbols:")
  for (sym in head(names(unmatched_counts), 30)) {
    cat(sprintf("    %-25s (%d occurrences, species: %s)\n", sym, unmatched_counts[sym],
                paste(unique(unmatched_log[[sym]]), collapse = ", ")))
  }
}

# Remove unresolved interactions and save them to the audit file
if (sum(flagged) > 0) {
  removed_genesym <- all_interactions[flagged, ]
  removed_genesym$unmatched_symbols <- unmatched_names[flagged]
  removed_path <- file.path(OUTPUT_DIR, "removed_invalid_genesymbols.csv")
  write.csv(removed_genesym, removed_path, row.names = FALSE)
  log_step("6B-GENESYM", sprintf("  Removed entries saved to: %s", removed_path))
  log_step("6B-GENESYM", "  By source:")
  print(table(removed_genesym$source_db))
  all_interactions <- all_interactions[!flagged, ]
}

n_after_genesym <- nrow(all_interactions)
log_step("6B-GENESYM", sprintf("  Gene symbol validation: %d -> %d interactions (-%d removed)",
                               n_before_genesym, n_after_genesym, n_before_genesym - n_after_genesym))

# ==============================================================================
# STEP 7: Rebuild interaction_name and interaction_name_2
# ==============================================================================

log_step("7-NAMES", "Rebuilding interaction_name and interaction_name_2...")

all_interactions$interaction_name <- paste(all_interactions$ligand,
                                           all_interactions$receptor, sep = "_")

all_interactions$interaction_name_2 <- ifelse(
  grepl("_", all_interactions$receptor),
  paste0(all_interactions$ligand, " - (", gsub("_", "+", all_interactions$receptor), ")"),
  paste(all_interactions$ligand, "-", all_interactions$receptor)
)

# Check for duplicate interaction_names WITHIN the same species+source
# (can happen after case normalization: "TGFB1_TGFBR1" and "Tgfb1_Tgfbr1" collide)
dup_key <- paste(all_interactions$species, all_interactions$source_db,
                 all_interactions$interaction_name, sep = "|||")
n_intra_dup <- sum(duplicated(dup_key))
if (n_intra_dup > 0) {
  log_step("7-NAMES", sprintf("Removing %d intra-source duplicates (from case normalization)", n_intra_dup))
  all_interactions <- all_interactions[!duplicated(dup_key), ]
}

log_step("7-NAMES", sprintf("After rebuild: %d interactions", nrow(all_interactions)))
log_sample(all_interactions, "7-NAMES")

# ==============================================================================
# STEP 8: Normalize evidence column
# ==============================================================================

log_step("8-EVIDENCE", "Normalizing evidence column formatting...")

n_before <- sum(!is.na(all_interactions$evidence))

# Clean aberrant suffixes and free-text noise
all_interactions$evidence <- gsub(";NO$", "", all_interactions$evidence)
all_interactions$evidence <- gsub(";also_called_[^;]+", "", all_interactions$evidence)
all_interactions$evidence <- gsub("Induces_.*_in_", "", all_interactions$evidence)

# Normalize PMID prefix whitespace
all_interactions$evidence <- gsub("PMID:\\s+", "PMID:", all_interactions$evidence)

# Source-name-as-evidence -> NA
source_evidence <- c("MultiNiche database", "MultiNicheNet", "NeuronChat",
                     "CellTalkDB", "CellPhoneDB", "PlantPhoneDB")
all_interactions$evidence[all_interactions$evidence %in% source_evidence] <- NA_character_

# Normalize bare PMIDs to "PMID:XXXXXXXX"
is_bare <- grepl("^[0-9]{7,8}(,[0-9]{7,8})*$", all_interactions$evidence) &
           !is.na(all_interactions$evidence)
if (sum(is_bare) > 0) {
  all_interactions$evidence[is_bare] <- sapply(all_interactions$evidence[is_bare], function(x) {
    paste(paste0("PMID:", trimws(strsplit(x, ",")[[1]])), collapse = "; ")
  })
}

n_after <- sum(!is.na(all_interactions$evidence))
log_step("8-EVIDENCE", sprintf("Evidence: %d -> %d non-NA entries", n_before, n_after))

# ==============================================================================
# STEP 9: Add 'sources' column
# ==============================================================================

log_step("9-SOURCES", "Adding 'sources' column (initialized from source_db)...")
all_interactions$sources <- all_interactions$source_db

# ==============================================================================
# STEP 10: Split by species, save individual + merged per species
# ==============================================================================

log_step("10-OUTPUT", "Splitting by species and saving...")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

species_list <- unique(all_interactions$species)

for (sp in species_list) {
  log_step("10-OUTPUT", sprintf("=== Processing %s ===", toupper(sp)))

  sp_dir <- file.path(OUTPUT_DIR, sp)
  if (!dir.exists(sp_dir)) dir.create(sp_dir, recursive = TRUE)

  sp_data <- all_interactions[all_interactions$species == sp, ]
  sp_sources <- unique(sp_data$source_db)

  # Get species-specific geneInfo, complex, cofactor tables
  sp_geneinfo <- geneinfo_tables[[sp]]  # Can be NULL
  sp_complex  <- complex_tables[[sp]]   # Can be NULL
  sp_cofactor <- cofactor_tables[[sp]]  # Can be NULL

  # --- Save individual source files ---
  for (src in sp_sources) {
    src_data <- sp_data[sp_data$source_db == src, ]
    rownames(src_data) <- src_data$interaction_name

    src_db <- build_cellchat_db(src_data, sp_complex, sp_cofactor, sp_geneinfo)

    safe_src <- gsub("[^A-Za-z0-9_]", "_", src)
    out_path <- file.path(sp_dir, sprintf("CellHubDB_%s_%s.rds", sp, safe_src))
    saveRDS(src_db, out_path)

    log_step("10-OUTPUT", sprintf("  %s: %d interactions -> %s", src, nrow(src_data), basename(out_path)))
  }

  # --- Merge + dedup within species ---
  log_step("10-OUTPUT", sprintf("  Merging %d sources for %s...", length(sp_sources), sp))

  n_before_merge <- nrow(sp_data)

  # Merge key: lowercase ligand|||receptor in the source orientation, so that reciprocal pairs
  # (EFNB2 -> EPHB4 and EPHB4 -> EFNB2) stay distinct. CellPhoneDB partners are not oriented:
  # a CellPhoneDB pair is flipped when only its reverse exists in an oriented source.
  # CellPhoneDB (partner_a/partner_b) and PlantCellChatDB (undirected PPI) do not encode direction
  NON_ORIENTED_SOURCES <- c("CellPhoneDB", "PlantCellChatDB")
  key_fwd <- paste(tolower(sp_data$ligand), tolower(sp_data$receptor), sep = "|||")
  key_rev <- paste(tolower(sp_data$receptor), tolower(sp_data$ligand), sep = "|||")
  is_unoriented <- sp_data$source_db %in% NON_ORIENTED_SOURCES
  oriented_keys <- key_fwd[!is_unoriented]
  to_flip <- is_unoriented & !(key_fwd %in% oriented_keys) & (key_rev %in% oriented_keys)
  if (any(to_flip)) {
    old_lig <- sp_data$ligand[to_flip]
    sp_data$ligand[to_flip] <- sp_data$receptor[to_flip]
    sp_data$receptor[to_flip] <- old_lig
    sp_data$interaction_name[to_flip] <- paste(sp_data$ligand[to_flip], sp_data$receptor[to_flip], sep = "_")
    sp_data$interaction_name_2[to_flip] <- ifelse(grepl("_", sp_data$receptor[to_flip]),
      paste0(sp_data$ligand[to_flip], " - (", gsub("_", "+", sp_data$receptor[to_flip]), ")"),
      paste(sp_data$ligand[to_flip], "-", sp_data$receptor[to_flip]))
    key_fwd[to_flip] <- key_rev[to_flip]
    log_step("10-OUTPUT", sprintf("  Oriented %d CellPhoneDB pairs on oriented sources", sum(to_flip)))
  }
  # Non-oriented pairs with no oriented counterpart use an order-free key, so that A-B and B-A
  # listed by an undirected source collapse into one entry
  key_sorted <- ifelse(tolower(sp_data$ligand) < tolower(sp_data$receptor), key_fwd, key_rev)
  free <- is_unoriented & !(key_fwd %in% oriented_keys)
  key_fwd[free] <- key_sorted[free]
  canonical <- key_fwd

  sp_data$canonical_key <- canonical
  sp_data$priority <- ifelse(sp_data$source_db %in% names(SOURCE_PRIORITY),
                              SOURCE_PRIORITY[sp_data$source_db], 0)

  dup_keys <- unique(canonical[duplicated(canonical)])
  is_dup <- canonical %in% dup_keys

  unique_part <- sp_data[!is_dup, ]
  dup_part <- sp_data[is_dup, ]

  if (nrow(dup_part) > 0) {
    deduped_list <- lapply(split(dup_part, dup_part$canonical_key), function(group) {
      group <- group[order(-group$priority), ]
      winner <- group[1, , drop = FALSE]
      # Merge sources
      winner$sources <- paste(unique(group$source_db), collapse = ", ")
      # Merge evidence
      all_ev <- unique(na.omit(group$evidence))
      if (length(all_ev) > 0) winner$evidence <- paste(all_ev, collapse = "; ")
      # Fill NA gaps from lower-priority sources
      for (col in c("pathway_name", "annotation", "interaction_type")) {
        if (is.na(winner[[col]])) {
          fill <- na.omit(group[[col]])
          if (length(fill) > 0) winner[[col]] <- fill[1]
        }
      }
      winner
    })
    deduped_part <- do.call(rbind, deduped_list)
  } else {
    deduped_part <- dup_part[0, ]
  }

  merged_data <- rbind(unique_part, deduped_part)
  merged_data$canonical_key <- NULL
  merged_data$priority <- NULL

  # Remove any remaining name collisions
  if (any(duplicated(merged_data$interaction_name))) {
    merged_data <- merged_data[!duplicated(merged_data$interaction_name), ]
  }
  rownames(merged_data) <- merged_data$interaction_name

  n_after_merge <- nrow(merged_data)
  n_multi <- sum(grepl(",", merged_data$sources))

  log_step("10-OUTPUT", sprintf("  Merge: %d -> %d interactions (-%d dups, %d multi-source)",
    n_before_merge, n_after_merge, n_before_merge - n_after_merge, n_multi))

  # Save merged DB
  merged_db <- build_cellchat_db(merged_data, sp_complex, sp_cofactor, sp_geneinfo)
  merged_path <- file.path(sp_dir, sprintf("CellHubDB_%s_merged.rds", sp))
  saveRDS(merged_db, merged_path)

  log_step("10-OUTPUT", sprintf("  Merged DB: %d interactions -> %s", n_after_merge, basename(merged_path)))
}

# ==============================================================================
# STEP 11: Final integrity checks (all species)
# ==============================================================================

log_step("11-CHECK", "Running final integrity checks...")

for (sp in species_list) {
  sp_dir <- file.path(OUTPUT_DIR, sp)
  merged_path <- file.path(sp_dir, sprintf("CellHubDB_%s_merged.rds", sp))
  if (!file.exists(merged_path)) next
  
  merged_db <- readRDS(merged_path)
  inter <- merged_db$interaction
  
  n_empty <- sum(sapply(inter, function(col) {
    if (is.character(col)) sum(col == "", na.rm = TRUE) else 0
  }))
  n_dup <- sum(duplicated(inter$interaction_name))
  # Same rule as STEP 5: known gene symbols and complex names are not UniProt accessions
  n_up <- sum(vapply(seq_len(nrow(inter)), function(i) {
    length(c(unresolved_tokens(inter$ligand[i]), unresolved_tokens(inter$receptor[i]))) > 0
  }, logical(1)))
  rn_ok <- all(rownames(inter) == inter$interaction_name)
  
  log_step("11-CHECK", sprintf("  %s: %d interactions | empty=%d | dups=%d | uniprot=%d | rownames=%s",
                               toupper(sp), nrow(inter), n_empty, n_dup, n_up, if (rn_ok) "OK" else "FAIL"))
}

# ==============================================================================
# STEP 12: Summary
# ==============================================================================

log_step("12-SUMMARY", "")
log_step("12-SUMMARY", "==================== PIPELINE COMPLETE ====================")
log_step("12-SUMMARY", "")

# List all output files
all_outputs <- list.files(OUTPUT_DIR, pattern = "CellHubDB_.*\\.rds$",
                          recursive = TRUE, full.names = TRUE)
for (f in sort(all_outputs)) {
  fsize <- file.info(f)$size
  db <- readRDS(f)
  n <- nrow(db$interaction)
  rel_path <- sub(paste0(OUTPUT_DIR, "/"), "", f)
  log_step("12-SUMMARY", sprintf("  %-50s %6d interactions  (%5.0f KB)", rel_path, n, fsize / 1024))
}

log_step("12-SUMMARY", "")
log_step("12-SUMMARY", sprintf("Total output files: %d", length(all_outputs)))
log_step("12-SUMMARY", "Done!")
