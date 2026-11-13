# Step 03 (Reviewer 1 reliability + complex flag; Reviewer 2 major 3, 4, 5): confidence annotation from the evidence of each pair
# (curated vs orthology-inferred / predicted support per source, predicted_only, is_complex, role corroboration, tier), EGFR check,
# Direct/Inferred summary. Pair-level rules: SUPPORT_RULES in revision_helpers.R, read from the per-source files of src_root.
source("review/revision_helpers.R")
db_path <- "review/02b_pathways/CellHubDB_pathways.rds"; out_dir <- "review/03_confidence"; src_root <- "databases"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
db <- readRDS(db_path)
# Sources computationally inferred for every pair (decision to state in the manuscript)
predicted_sources <- c("MultiNicheNet")

# Curated support of every pair, per species and source, from the per-source files
support_of_species <- function(sp) {
  files <- list.files(file.path(src_root, sp), pattern = paste0("^CellHubDB_", sp, "_.*\\.rds$"))
  files <- files[!grepl("_merged\\.rds$", files)]
  srcs <- sub(paste0("^CellHubDB_", sp, "_(.*)\\.rds$"), "\\1", files)
  return(setNames(lapply(seq_along(files), function(k) source_support_status(readRDS(file.path(src_root, sp, files[k])), srcs[k])), srcs))
}
support <- setNames(lapply(names(db), support_of_species), names(db))

# Pair-level support per source (how many pairs of each source keep curated support)
support_tab <- do.call(rbind, lapply(names(support), function(sp) do.call(rbind, lapply(names(support[[sp]]), function(s) {
  v <- support[[sp]][[s]]
  data.frame(species = sp, source = s, pairs = length(v), curated = sum(v), inferred_or_predicted = sum(!v), pct_curated = round(100 * mean(v), 1))
}))))
write.csv(support_tab, file.path(out_dir, "pair_support_per_source.csv"), row.names = FALSE)
print(support_tab, row.names = FALSE)

db <- setNames(lapply(names(db), function(sp) annotate_confidence(db[[sp]], predicted_sources, support[[sp]])), names(db))
# Merged pairs not found in their source file (key differences): source default used
unmatched <- do.call(rbind, lapply(names(db), function(sp) { u <- attr(db[[sp]], "unmatched"); if (length(u) == 0) NULL else data.frame(species = sp, source = names(u), pairs_not_found = as.integer(u)) }))
if (!is.null(unmatched)) { write.csv(unmatched, file.path(out_dir, "pairs_not_found_in_source_file.csv"), row.names = FALSE); message("Pairs not found in their source file (source default used):"); print(unmatched, row.names = FALSE) }

# Tier distribution per species
tiers <- do.call(rbind, lapply(names(db), function(sp) {
  it <- db[[sp]]$interaction
  data.frame(species = sp, n = nrow(it), single_source = sum(it$n_sources == 1), multi_curated = sum(it$n_curated_sources >= 2), predicted_only = sum(it$predicted_only),
             complex = sum(it$is_complex), high = sum(it$confidence_tier == "high"), medium = sum(it$confidence_tier == "medium"), low = sum(it$confidence_tier == "low"),
             medium_role_corroborated = sum(it$confidence_tier == "medium" & it$role_corroborated %in% TRUE),
             medium_role_uncorroborated = sum(it$confidence_tier == "medium" & it$role_corroborated %in% FALSE))
}))
write.csv(tiers, file.path(out_dir, "confidence_tiers.csv"), row.names = FALSE)
print(tiers, row.names = FALSE)

# Reviewer 2 example: EGFR partners in human (CALM2, ARF4, GSTP1 expected with role_corroborated = FALSE)
hi <- db$human$interaction
egfr <- hi[hi$ligand == "EGFR" | hi$receptor == "EGFR" | vapply(hi$receptor, function(p) "EGFR" %in% partner_subunits(p, db$human$complex), logical(1)), ]
egfr_cols <- intersect(c("ligand", "receptor", "sources", "curated_sources", "n_sources", "n_curated_sources", "predicted_only", "ligand_role_sources", "receptor_role_sources", "role_corroborated", "confidence_tier"), colnames(egfr))
write.csv(egfr[, egfr_cols], file.path(out_dir, "EGFR_partners_human.csv"), row.names = FALSE)
print(table(tier = egfr$confidence_tier, role_corroborated = egfr$role_corroborated, useNA = "ifany"))
print(egfr[egfr$ligand %in% c("CALM2", "ARF4", "GSTP1", "HSP90AA1"), egfr_cols], row.names = FALSE)

# Direct vs Inferred evidence (Reviewer 2, major 5), from the connectomeDB2025 evidence_type column
ev <- do.call(rbind, lapply(names(db), function(sp) {
  e <- db[[sp]]$interaction$evidence_type
  if (is.null(e)) e <- rep(NA_character_, nrow(db[[sp]]$interaction))
  data.frame(species = sp, n = length(e), direct = sum(e %in% "Direct"), inferred = sum(e %in% "Inferred"), not_connectomeDB = sum(is.na(e)),
             pct_inferred_of_connectomeDB = round(100 * sum(e %in% "Inferred") / max(sum(!is.na(e)), 1), 1))
}))
write.csv(ev, file.path(out_dir, "direct_vs_inferred.csv"), row.names = FALSE)
print(ev, row.names = FALSE)

saveRDS(db, file.path(out_dir, "CellHubDB_annotated.rds"))
