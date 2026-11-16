# patch_release_2.1.1.R — CellHubDB 2.1.1: CellChat-native `annotation` column, nothing else changed.
# CellChat filters on `annotation` (subsetDB(search = "Secreted Signaling")); in 2.1 the canonical interaction type lives in
# `interaction_type` and `annotation` holds free text for most rows. This patch writes the canonical type into `annotation`
# and keeps the former text in `annotation_detail`, for CellHubDB.rds and every app_databases file.
# Usage: Rscript patch_release_2.1.1.R <release_2.1_dir> <output_2.1.1_dir>
args <- commandArgs(trailingOnly = TRUE)
in_dir <- args[1]; out_dir <- args[2]
if (dir.exists(out_dir)) stop("Output folder already exists: ", out_dir)
canonical <- c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor", "Non-protein Signaling")
log_rows <- list()

# Canonical annotation for one CellChat database; returns the patched database
patch_db <- function(db, label) {
  it <- db$interaction
  if (!"interaction_type" %in% colnames(it)) { log_rows[[length(log_rows) + 1]] <<- data.frame(file = label, n = nrow(it), changed = 0, note = "no interaction_type column"); return(db) }
  old <- if ("annotation" %in% colnames(it)) as.character(it$annotation) else rep(NA_character_, nrow(it))
  new <- ifelse(it$interaction_type %in% canonical, as.character(it$interaction_type), NA_character_)
  if (!"annotation_detail" %in% colnames(it)) it$annotation_detail <- old
  it$annotation <- new
  log_rows[[length(log_rows) + 1]] <<- data.frame(file = label, n = nrow(it), changed = sum(xor(is.na(old), is.na(new)) | (!is.na(old) & !is.na(new) & old != new)),
                                                  note = sprintf("canonical annotation: %.1f%% -> %.1f%%", 100 * mean(old %in% canonical), 100 * mean(new %in% canonical)))
  db$interaction <- it
  db
}

# Same database except annotation / annotation_detail: interactions, tiers, sources and every other column identical
check_same <- function(a, b, label) {
  ia <- a$interaction; ib <- b$interaction
  cols <- setdiff(colnames(ia), "annotation")
  ok <- nrow(ia) == nrow(ib) && isTRUE(all.equal(ia[, cols, drop = FALSE], ib[, cols, drop = FALSE], check.attributes = FALSE)) &&
    identical(a$complex, b$complex) && identical(a$cofactor, b$cofactor) && identical(a$geneInfo, b$geneInfo)
  if (!ok) stop("Content changed beyond the annotation column: ", label)
}

dir.create(out_dir, recursive = TRUE)
file.copy(list.files(in_dir, full.names = TRUE), out_dir, recursive = TRUE)
# Per-species app files
for (f in list.files(file.path(out_dir, "app_databases"), pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)) {
  db <- readRDS(f); p <- patch_db(db, sub(paste0("^", out_dir, "/?"), "", f)); check_same(db, p, f)
  attr(p, "version") <- "2.1.1"
  saveRDS(p, f)
}
# Release object (list by species), attributes kept, version updated
rel <- readRDS(file.path(out_dir, "CellHubDB.rds"))
rel_attr <- attributes(rel)
patched <- lapply(names(rel), function(sp) { p <- patch_db(rel[[sp]], paste0("CellHubDB.rds$", sp)); check_same(rel[[sp]], p, sp); p })
names(patched) <- names(rel)
for (a in setdiff(names(rel_attr), "names")) attr(patched, a) <- rel_attr[[a]]
attr(patched, "version") <- "2.1.1"
attr(patched, "build_date") <- format(Sys.Date())
saveRDS(patched, file.path(out_dir, "CellHubDB.rds"))
# App files and release object still identical, species by species
for (sp in names(patched)) {
  app <- readRDS(file.path(out_dir, "app_databases", sp, paste0("CellHubDB_", sp, "_merged.rds")))
  if (!isTRUE(all.equal(app$interaction, patched[[sp]]$interaction, check.attributes = FALSE))) stop("App merged file and CellHubDB.rds differ for ", sp)
}
log_df <- do.call(rbind, log_rows)
write.csv(log_df, file.path(out_dir, "reports", "patch_2.1.1_annotation.csv"), row.names = FALSE)
writeLines(c("# CellHubDB 2.1.1", "", "Patch of release 2.1: the CellChat `annotation` column now holds the canonical interaction type",
             "(Secreted Signaling, Cell-Cell Contact, ECM-Receptor, Non-protein Signaling; NA when unknown), copied from `interaction_type`.",
             "The former free-text annotation is kept in `annotation_detail`. Interactions, tiers, sources, complexes, cofactors and geneInfo are unchanged",
             "(checked file by file). Standalone CellChat users can therefore filter with subsetDB(db, search = \"Secreted Signaling\")."),
           file.path(out_dir, "CHANGES_2.1.1.md"))
# MD5 checksums of every file of the release
files <- list.files(out_dir, recursive = TRUE)
files <- files[files != "MD5SUMS"]
writeLines(paste(unname(tools::md5sum(file.path(out_dir, files))), files, sep = "  "), file.path(out_dir, "MD5SUMS"))
print(log_df[grepl("merged|CellHubDB.rds", log_df$file), ], row.names = FALSE)
message("CellHubDB 2.1.1 written to ", out_dir)
