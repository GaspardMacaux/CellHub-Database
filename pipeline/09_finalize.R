################################################################################
#   finalize_cellhubdb.R — make every curated database CellChat-ready and      #
#   annotate it, after curate_CellHub_DB.R and before build_cellhubdb.R.       #
#                                                                              #
# Applied to databases/<species>/CellHubDB_<species>_<source|merged>.rds:      #
#   1. CellChat schema (complete complex / cofactor tables, "" for empty refs) #
#   2. release repairs (undefined cofactors, cofactor genes in another case,   #
#      unresolvable partners, complexes with a numeric subunit)                #
#   3. interaction_type from the CellChat annotation column                    #
#   4. pathway labels: source, inherited from the human / mouse pair, or       #
#      "Unassigned" (CellChat needs a pathway for every interaction)           #
#   5. pair-level confidence (merged files): curated_sources, n_curated_sources,#
#      predicted_only, is_complex, role_corroborated, confidence_tier          #
# Input / output: databases/ (files overwritten); logs in reports/finalize_*   #
# Usage: Rscript pipeline/09_finalize.R (from the CellHub_DB root)
################################################################################

source("R/cellhubdb_functions.R")
suppressPackageStartupMessages(library(CellChat))
DATABASE_DIR <- "data/databases"
REPORT_DIR <- "output/reports"
PLACEHOLDER <- "Unassigned"
# Sources computationally inferred for every pair
PREDICTED_SOURCES <- c("MultiNicheNet")
dir.create(REPORT_DIR, showWarnings = FALSE)
schema <- get_cellchat_schema()

# Print a timestamped log message
log_msg <- function(msg) {
  cat(sprintf("[%s] FINALIZE | %s\n", format(Sys.time(), "%H:%M:%S"), msg))
  return(invisible(NULL))
}

# Per-source pathway labels taken from the merged pair of the same species, placeholder for the rest
fill_source_pathways <- function(x, merged) {
  it <- x$interaction
  key <- pair_key(it); mkey <- pair_key(merged$interaction)
  missing_pw <- is.na(it$pathway_name) | it$pathway_name == ""
  it$pathway_name[missing_pw] <- merged$interaction$pathway_name[match(key[missing_pw], mkey)]
  it$pathway_name[is.na(it$pathway_name) | it$pathway_name == ""] <- PLACEHOLDER
  x$interaction <- it
  return(x)
}

species <- sort(list.dirs(DATABASE_DIR, full.names = FALSE, recursive = FALSE))
files_of <- function(sp) list.files(file.path(DATABASE_DIR, sp), pattern = paste0("^CellHubDB_", sp, "_.*\\.rds$"))
log_msg(sprintf("%d species: %s", length(species), paste(species, collapse = ", ")))

# 1-3. Schema, repairs and interaction types on every file; the merged geneInfo is shared by the per-source files
dbs <- list(); repairs <- list(); removed <- list()
for (sp in species) {
  dbs[[sp]] <- list()
  for (f in files_of(sp)) {
    x <- standardize_cellchat_schema(readRDS(file.path(DATABASE_DIR, sp, f)), schema)
    r <- repair_species_db(x, sp)
    if (!is.null(r$repairs) && nrow(r$repairs) > 0) repairs[[length(repairs) + 1]] <- cbind(file = f, r$repairs)
    if (!is.null(r$removed) && nrow(r$removed) > 0) removed[[length(removed) + 1]] <- cbind(file = f, r$removed)
    dbs[[sp]][[f]] <- harmonize_interaction_type(r$db)
  }
  merged_file <- paste0("CellHubDB_", sp, "_merged.rds")
  if (!merged_file %in% names(dbs[[sp]])) stop(sp, ": merged file missing")
  for (f in setdiff(names(dbs[[sp]]), merged_file)) dbs[[sp]][[f]]$geneInfo <- dbs[[sp]][[merged_file]]$geneInfo
}

# 4. Pathways: lookups from the human and mouse merged databases, applied to every merged database, then to the per-source files
merged_name <- function(sp) paste0("CellHubDB_", sp, "_merged.rds")
lookups <- list(human = pathway_lookup(dbs$human[[merged_name("human")]]$interaction), mouse = pathway_lookup(dbs$mouse[[merged_name("mouse")]]$interaction))
for (sp in species) {
  dbs[[sp]][[merged_name(sp)]] <- propagate_pathways(dbs[[sp]][[merged_name(sp)]], lookups, PLACEHOLDER)
  for (f in setdiff(names(dbs[[sp]]), merged_name(sp))) dbs[[sp]][[f]] <- fill_source_pathways(dbs[[sp]][[f]], dbs[[sp]][[merged_name(sp)]])
}

# 5. Pair-level confidence on the merged databases, from the evidence of each per-source file
unmatched <- list(); support_tab <- list()
for (sp in species) {
  src_files <- setdiff(names(dbs[[sp]]), merged_name(sp))
  srcs <- sub(paste0("^CellHubDB_", sp, "_(.*)\\.rds$"), "\\1", src_files)
  support <- setNames(lapply(seq_along(src_files), function(k) source_support_status(dbs[[sp]][[src_files[k]]], srcs[k])), srcs)
  for (s in names(support)) support_tab[[length(support_tab) + 1]] <- data.frame(species = sp, source = s, pairs = length(support[[s]]), curated = sum(support[[s]]))
  m <- annotate_confidence(dbs[[sp]][[merged_name(sp)]], PREDICTED_SOURCES, support)
  u <- attr(m, "unmatched"); attr(m, "unmatched") <- NULL
  if (length(u) > 0) unmatched[[sp]] <- data.frame(species = sp, source = names(u), pairs_not_found = as.integer(u))
  dbs[[sp]][[merged_name(sp)]] <- m
}

# Checks CellChat relies on, then write the files back
checks <- list()
for (sp in species) for (f in names(dbs[[sp]])) {
  x <- dbs[[sp]][[f]]; pw <- unique(x$interaction$pathway_name)
  ok <- nrow(x$geneInfo) > 0 && length(pw) >= 2 && !any(duplicated(x$interaction$interaction_name))
  checks[[length(checks) + 1]] <- data.frame(species = sp, file = f, interactions = nrow(x$interaction), geneinfo = nrow(x$geneInfo), pathways = length(pw), ok = ok)
  rownames(x$interaction) <- x$interaction$interaction_name
  saveRDS(x, file.path(DATABASE_DIR, sp, f))
}
checks_df <- do.call(rbind, checks)
write.csv(checks_df, file.path(REPORT_DIR, "finalize_checks.csv"), row.names = FALSE)
if (length(repairs) > 0) write.csv(do.call(rbind, repairs), file.path(REPORT_DIR, "finalize_repairs.csv"), row.names = FALSE)
if (length(removed) > 0) write.csv(do.call(rbind, removed), file.path(REPORT_DIR, "finalize_removed_interactions.csv"), row.names = FALSE)
write.csv(do.call(rbind, support_tab), file.path(REPORT_DIR, "finalize_pair_support.csv"), row.names = FALSE)
if (length(unmatched) > 0) write.csv(do.call(rbind, unmatched), file.path(REPORT_DIR, "finalize_pairs_not_found.csv"), row.names = FALSE)

removed_df <- if (length(removed) > 0) do.call(rbind, removed) else NULL
if (!is.null(removed_df)) { log_msg("Interactions removed (merged files):"); print(table(removed_df$species[grepl("_merged", removed_df$file)], removed_df$reason[grepl("_merged", removed_df$file)])) }
if (length(unmatched) > 0) { log_msg("Merged pairs not found in their source file (source default used):"); print(do.call(rbind, unmatched), row.names = FALSE) }
print(checks_df[grepl("_merged", checks_df$file), ], row.names = FALSE)
if (!all(checks_df$ok)) stop("Finalize checks failed: ", paste(checks_df$file[!checks_df$ok], collapse = ", "))
log_msg(sprintf("Done: %d files finalized", nrow(checks_df)))
