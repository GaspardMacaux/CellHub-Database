################################################################################
#   build_cellhubdb.R — assemble the per-species merged databases into the     #
#   single CellHubDB.rds distributed with Cell-Hub and the GitHub repository,  #
#   and compute the statistics reported in the manuscript.                     #
#                                                                              #
# Input:  databases/<species>/CellHubDB_<species>_merged.rds                   #
# Output: CellHubDB.rds                      (named list, one CellChat DB per species)
#         CellHubDB_species_summary.csv      (interactions, sources, multi-source support)
#         CellHubDB_source_contribution.csv  (interactions contributed by each source)
# Usage:  Rscript database_curation/build_cellhubdb.R  (from the CellHub_DB root)
################################################################################

DATABASE_DIR <- "databases"
OUTPUT_RDS   <- "CellHubDB.rds"
DB_VERSION   <- "2.1"

# Sources supporting each interaction, from the comma-separated 'sources' column written by the merge
split_sources <- function(sources) {
  return(lapply(strsplit(as.character(sources), ",\\s*"), function(x) unique(x[nzchar(x)])))
}

# One summary row per species: size, sources present, and support by >= 2 sources
summarize_species <- function(sp, db) {
  inter <- db$interaction
  src_list <- split_sources(inter$sources)
  n_src <- lengths(src_list)
  return(data.frame(
    species            = sp,
    interactions       = nrow(inter),
    sources            = paste(sort(unique(unlist(src_list))), collapse = ", "),
    n_sources          = length(unique(unlist(src_list))),
    multi_source       = sum(n_src >= 2),
    pct_multi_source   = round(100 * mean(n_src >= 2), 1),
    complexes          = if (is.null(db$complex)) 0L else nrow(db$complex),
    stringsAsFactors   = FALSE
  ))
}

# Interactions supported by each source in each species (an interaction counts for every source supporting it)
source_contribution <- function(sp, db) {
  src_list <- split_sources(db$interaction$sources)
  all_src <- sort(unique(unlist(src_list)))
  rows <- lapply(all_src, function(s) {
    has_s <- vapply(src_list, function(x) s %in% x, logical(1))
    data.frame(species = sp, source = s, supported = sum(has_s),
               exclusive = sum(has_s & lengths(src_list) == 1), stringsAsFactors = FALSE)
  })
  return(do.call(rbind, rows))
}

merged_files <- list.files(DATABASE_DIR, pattern = "^CellHubDB_.+_merged\\.rds$", recursive = TRUE, full.names = TRUE)
if (length(merged_files) == 0) stop("No merged databases found under ", DATABASE_DIR, "/ - run curate_CellHub_DB.R first")
species <- sub("^CellHubDB_(.+)_merged\\.rds$", "\\1", basename(merged_files))
cat(sprintf("[Build] %d species found: %s\n", length(species), paste(sort(species), collapse = ", ")))

cellhubdb <- setNames(lapply(merged_files, readRDS), species)
cellhubdb <- cellhubdb[sort(names(cellhubdb))]

# Each entry must be directly assignable to cellchat@DB
required <- c("interaction", "complex", "cofactor", "geneInfo")
for (sp in names(cellhubdb)) {
  missing <- setdiff(required, names(cellhubdb[[sp]]))
  if (length(missing) > 0) stop(sprintf("%s: missing element(s) %s", sp, paste(missing, collapse = ", ")))
  if (any(duplicated(cellhubdb[[sp]]$interaction$interaction_name))) stop(sprintf("%s: duplicated interaction names", sp))
}

summary_df <- do.call(rbind, lapply(names(cellhubdb), function(sp) summarize_species(sp, cellhubdb[[sp]])))
contrib_df <- do.call(rbind, lapply(names(cellhubdb), function(sp) source_contribution(sp, cellhubdb[[sp]])))

# Provenance stored as attributes so the list structure (db[["mouse"]]) stays unchanged for users
attr(cellhubdb, "version")    <- DB_VERSION
attr(cellhubdb, "build_date") <- format(Sys.Date())
attr(cellhubdb, "sources")    <- sort(unique(contrib_df$source))
attr(cellhubdb, "summary")    <- summary_df

saveRDS(cellhubdb, OUTPUT_RDS)
write.csv(summary_df, "CellHubDB_species_summary.csv", row.names = FALSE)
write.csv(contrib_df, "CellHubDB_source_contribution.csv", row.names = FALSE)

cat(sprintf("[Build] %s written (%s)\n", OUTPUT_RDS, format(structure(file.size(OUTPUT_RDS), class = "object_size"), units = "auto")))
cat(sprintf("[Build] Version %s | %d species | %d sources | %s interactions\n",
            DB_VERSION, length(cellhubdb), length(attr(cellhubdb, "sources")), format(sum(summary_df$interactions), big.mark = ",")))
print(summary_df[, c("species", "interactions", "n_sources", "multi_source", "pct_multi_source")], row.names = FALSE)
