# Step 13 (Reviewer 2, major 5): cross-species conservation of ligand-receptor pairs, separating curated (Direct) from orthology-inferred support.
# A pair is "curated" in a species when at least one source gives it curated support (n_curated_sources >= 1, confidence_tier rules of step 03 /
# finalize: connectomeDB2025 Direct, CellChatDB, CellPhoneDB, CellTalkDB, NeuronChatDB; connectomeDB2025 Inferred and NicheNet are not curated).
# Pairs are matched across species by upper-case "ligand|receptor" symbols: one-to-one orthologs keeping the same symbol match, renamed orthologs
# and teleost paralogs (zebrafish a/b copies) do not, so shared counts are a lower bound for distant species.
# Usage: Rscript revision/13_conservation_direct_inferred.R [release/2.1/CellHubDB.rds] [revision/results/13_conservation_v2.1]
args <- commandArgs(trailingOnly = TRUE)
db_path <- if (length(args) >= 1) args[1] else "release/2.1/CellHubDB.rds"
out_dir <- if (length(args) >= 2) args[2] else "revision/results/13_conservation_v2.1"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
db <- readRDS(db_path)

# Distinct upper-case pair keys of one species, with their curated status (a key is curated if any of its rows is)
species_keys <- function(it) {
  if (!"n_curated_sources" %in% colnames(it)) stop("n_curated_sources missing: run finalize_cellhubdb.R / step 03 first")
  keys <- toupper(paste(it$ligand, it$receptor, sep = "|"))
  cur <- tapply(it$n_curated_sources >= 1, keys, any)
  return(setNames(as.logical(cur), names(cur)))
}

# Vertebrate species covered by connectomeDB2025 (Direct / Inferred evidence present), plus human and mouse as references
has_ev <- vapply(db, function(x) "evidence_type" %in% colnames(x$interaction) && any(x$interaction$evidence_type %in% c("Direct", "Inferred")), logical(1))
vert <- names(db)[has_ev]
if (!all(c("human", "mouse") %in% vert)) stop("human and mouse must carry connectomeDB2025 evidence")
keys <- setNames(lapply(vert, function(sp) species_keys(db[[sp]]$interaction)), vert)
message(sprintf("%d vertebrate species: %s", length(vert), paste(vert, collapse = ", ")))

# 1. Per species: pairs, curated pairs, connectomeDB2025 Direct / Inferred rows
per_species <- do.call(rbind, lapply(vert, function(sp) {
  it <- db[[sp]]$interaction; k <- keys[[sp]]
  data.frame(species = sp, pairs = length(k), curated_pairs = sum(k), inferred_or_predicted_pairs = sum(!k), pct_curated = round(100 * mean(k), 1),
             connectomeDB_direct_rows = sum(it$evidence_type %in% "Direct"), connectomeDB_inferred_rows = sum(it$evidence_type %in% "Inferred"))
}))
write.csv(per_species, file.path(out_dir, "pairs_per_species.csv"), row.names = FALSE)
print(per_species, row.names = FALSE)

# 2. Conservation against each reference species (human, mouse): all pairs vs pairs curated in both species
conservation_vs <- function(sp, ref) {
  k <- keys[[sp]]; r <- keys[[ref]]
  shared <- intersect(names(k), names(r)); shared_cur <- intersect(names(k)[k], names(r)[r])
  return(data.frame(species = sp, reference = ref, pairs = length(k), shared_all = length(shared), pct_species_shared_all = round(100 * length(shared) / length(k), 1),
                    pct_reference_shared_all = round(100 * length(shared) / length(r), 1), curated_pairs = sum(k), shared_curated_both = length(shared_cur),
                    pct_species_curated_shared = if (sum(k) > 0) round(100 * length(shared_cur) / sum(k), 1) else NA_real_,
                    shared_inferred_in_species = sum(!k[shared])))
}
cons <- do.call(rbind, lapply(c("human", "mouse"), function(ref) do.call(rbind, lapply(setdiff(vert, ref), conservation_vs, ref = ref))))
write.csv(cons, file.path(out_dir, "conservation_vs_human_mouse.csv"), row.names = FALSE)
print(cons, row.names = FALSE)

# 3. Human-mouse sentence of the manuscript (v1: 7,618 shared pairs = 81.9 % of mouse, 69.6 % of human), recomputed on all and on curated pairs
hm <- cons[cons$species == "mouse" & cons$reference == "human", ]
h <- keys$human; m <- keys$mouse; shared <- intersect(names(h), names(m))
hm_split <- data.frame(category = c("curated in both", "curated in human only", "curated in mouse only", "curated in neither"),
                       shared_pairs = c(sum(h[shared] & m[shared]), sum(h[shared] & !m[shared]), sum(!h[shared] & m[shared]), sum(!h[shared] & !m[shared])))
hm_split$pct_of_shared <- round(100 * hm_split$shared_pairs / length(shared), 1)
write.csv(hm_split, file.path(out_dir, "human_mouse_shared_by_support.csv"), row.names = FALSE)
message(sprintf("Human-mouse: %d shared pairs = %.1f %% of mouse, %.1f %% of human; curated in both: %d (%.1f %% of mouse curated pairs)",
                hm$shared_all, hm$pct_species_shared_all, hm$pct_reference_shared_all, hm$shared_curated_both, hm$pct_species_curated_shared))
print(hm_split, row.names = FALSE)

# 4. Species x species Jaccard index, on all pairs and on curated pairs (supplementary heatmap)
jaccard <- function(a, b) if (length(union(a, b)) == 0) NA_real_ else length(intersect(a, b)) / length(union(a, b))
jac_matrix <- function(curated_only) {
  sets <- lapply(keys, function(k) if (curated_only) names(k)[k] else names(k))
  m <- outer(seq_along(sets), seq_along(sets), Vectorize(function(i, j) jaccard(sets[[i]], sets[[j]])))
  dimnames(m) <- list(vert, vert)
  return(round(m, 3))
}
write.csv(jac_matrix(FALSE), file.path(out_dir, "jaccard_all_pairs.csv"))
write.csv(jac_matrix(TRUE), file.path(out_dir, "jaccard_curated_pairs.csv"))

# 5. Containment matrix (asymmetric): % of the pairs of the row species found in the column species. Unlike Jaccard, it is not driven down by
# the size difference between human / mouse (~10,000 pairs) and the connectomeDB-only species (~3,000 pairs)
containment_matrix <- function(curated_only) {
  sets <- lapply(keys, function(k) if (curated_only) names(k)[k] else names(k))
  m <- outer(seq_along(sets), seq_along(sets), Vectorize(function(i, j) if (length(sets[[i]]) == 0) NA_real_ else 100 * length(intersect(sets[[i]], sets[[j]])) / length(sets[[i]])))
  dimnames(m) <- list(vert, vert)
  return(round(m, 1))
}
write.csv(containment_matrix(FALSE), file.path(out_dir, "containment_all_pairs.csv"))
write.csv(containment_matrix(TRUE), file.path(out_dir, "containment_curated_pairs.csv"))
message("Done: ", out_dir)
