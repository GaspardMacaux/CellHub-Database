# Step 02b (Reviewer 2, other 1 + other 6): pathway labels inherited from the human / mouse pair, "Unassigned" for the rest,
# then CellChat smoke test of the corrected database on all species
source("review/revision_helpers.R")
suppressPackageStartupMessages(library(CellChat))
future::plan("sequential")
db_path <- "review/02_geneinfo/CellHubDB_geneinfo.rds"; out_dir <- "review/02b_pathways"; placeholder <- "Unassigned"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
db <- readRDS(db_path)

# Case-insensitive ligand|receptor key, shared by orthologous pairs that keep the same symbol across species
pair_key <- function(it) tolower(paste(it$ligand, it$receptor, sep = "|"))

# Lookup table key -> pathway from one reference species (first non-missing label per key)
pathway_lookup <- function(it) {
  ok <- !is.na(it$pathway_name) & it$pathway_name != ""
  keys <- pair_key(it)[ok]; pw <- it$pathway_name[ok]
  return(pw[!duplicated(keys)] |> setNames(keys[!duplicated(keys)]))
}

# Fill missing pathway labels from the human then the mouse pair, record the origin, label the rest with the placeholder
propagate_pathways <- function(sp_db, lookups, placeholder) {
  it <- sp_db$interaction
  missing <- is.na(it$pathway_name) | it$pathway_name == ""
  it$pathway_origin <- ifelse(missing, NA_character_, "source")
  keys <- pair_key(it)
  for (ref in names(lookups)) {
    fill <- missing & is.na(it$pathway_origin) & keys %in% names(lookups[[ref]])
    it$pathway_name[fill] <- unname(lookups[[ref]][keys[fill]])
    it$pathway_origin[fill] <- paste0("inherited_", ref)
  }
  rest <- is.na(it$pathway_origin)
  it$pathway_name[rest] <- placeholder; it$pathway_origin[rest] <- "unassigned"
  sp_db$interaction <- it
  return(sp_db)
}

lookups <- list(human = pathway_lookup(db$human$interaction), mouse = pathway_lookup(db$mouse$interaction))
db <- lapply(db, propagate_pathways, lookups = lookups, placeholder = placeholder)

# Pathway origin per species
origin <- do.call(rbind, lapply(names(db), function(sp) {
  o <- table(factor(db[[sp]]$interaction$pathway_origin, levels = c("source", "inherited_human", "inherited_mouse", "unassigned")))
  pw <- db[[sp]]$interaction$pathway_name
  data.frame(species = sp, n_interactions = length(pw), t(as.integer(o)) |> as.data.frame() |> setNames(names(o)),
             pct_unassigned = round(100 * o[["unassigned"]] / length(pw), 1), n_pathways = length(unique(pw[pw != placeholder])), stringsAsFactors = FALSE)
}))
write.csv(origin, file.path(out_dir, "pathway_origin.csv"), row.names = FALSE)
print(origin, row.names = FALSE)

# CellChat smoke test of the corrected database (schema + geneInfo + pathways)
smoke <- do.call(rbind, lapply(names(db), function(sp) { message("Smoke test (corrected): ", sp); run_cellchat_smoke(db[[sp]], sp) }))
write.csv(smoke, file.path(out_dir, "smoke_corrected.csv"), row.names = FALSE)
print(smoke[, c("species", "computeCommunProbPathway", "n_significant", "first_error")], row.names = FALSE)
saveRDS(db, file.path(out_dir, "CellHubDB_pathways.rds"))
