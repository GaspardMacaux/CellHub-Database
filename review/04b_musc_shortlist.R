# Step 04b (Reviewer 2, major 2c + 3): MuSC calls of the merged database absent from the default CellChatDB, classified by
# novelty at the gene level (complex subunits expanded), aggregated per LR pair and direction, for manual literature checking
source("review/revision_helpers.R")
calls_path <- "review/04_benchmark_v2.1/calls_MuSC_for_manual_check.csv"; annotated_path <- "release/2.1/CellHubDB.rds"
out_dir <- "review/04_benchmark_v2.1"; focus_group <- "MuSC"; species <- "mouse"
calls <- read.csv(calls_path, stringsAsFactors = FALSE)
merged <- readRDS(annotated_path)[[species]]
ref <- CellChat::CellChatDB.mouse

# Gene sets of each partner (complex subunits expanded, case-insensitive)
genes_of <- function(p, cpx) tolower(partner_subunits(p, cpx))
ref_lig <- lapply(ref$interaction$ligand, genes_of, cpx = ref$complex); ref_rec <- lapply(ref$interaction$receptor, genes_of, cpx = ref$complex)
ref_pairs <- unique(unlist(mapply(function(l, r) as.vector(outer(l, r, paste, sep = "|")), ref_lig, ref_rec, SIMPLIFY = FALSE)))
ref_ligands <- unique(unlist(ref_lig)); ref_receptors <- unique(unlist(ref_rec))

# Novelty of one call against the default CellChatDB:
# subunit_of_default = same ligand gene and receptor gene already paired in CellChatDB (often a complex written differently);
# new_receptor_for_known_ligand / new_ligand_for_known_receptor / both_new = genuinely new pairing
classify_novelty <- function(l, r) {
  lg <- genes_of(l, merged$complex); rg <- genes_of(r, merged$complex)
  if (any(as.vector(outer(lg, rg, paste, sep = "|")) %in% ref_pairs)) return("subunit_of_default")
  kl <- any(lg %in% ref_ligands); kr <- any(rg %in% ref_receptors)
  if (kl && kr) return("new_pair_known_partners")
  if (kl) return("new_receptor_for_known_ligand")
  if (kr) return("new_ligand_for_known_receptor")
  return("both_partners_new")
}

new_calls <- calls[!calls$in_default_db, ]
new_calls$novelty <- mapply(classify_novelty, new_calls$ligand, new_calls$receptor)
new_calls$role <- ifelse(new_calls$source == focus_group, "MuSC_sender", "MuSC_receiver")
new_calls$partner <- ifelse(new_calls$source == focus_group, as.character(new_calls$target), as.character(new_calls$source))
new_calls$recommended_set <- new_calls$confidence_tier == "high" | (new_calls$confidence_tier == "medium" & new_calls$role_corroborated %in% TRUE)
print(table(novelty = new_calls$novelty, tier = new_calls$confidence_tier))

# One row per LR pair and direction, strongest probability first
key <- paste(new_calls$ligand, new_calls$receptor, new_calls$role, sep = "||")
agg <- do.call(rbind, lapply(split(seq_len(nrow(new_calls)), key), function(ix) {
  x <- new_calls[ix, ]
  data.frame(ligand = x$ligand[1], receptor = x$receptor[1], role = x$role[1], partners = paste(sort(unique(x$partner)), collapse = ", "),
             n_partners = length(unique(x$partner)), max_prob = signif(max(x$prob), 3), min_pval = min(x$pval), sources = x$sources[1],
             n_curated_sources = x$n_curated_sources[1], confidence_tier = x$confidence_tier[1], role_corroborated = x$role_corroborated[1],
             recommended_set = x$recommended_set[1], novelty = x$novelty[1], evidence = substr(x$evidence[1], 1, 200),
             literature_check = "", stringsAsFactors = FALSE)
}))
agg <- agg[order(!agg$recommended_set, agg$novelty == "subunit_of_default", agg$role, -agg$max_prob), ]
write.csv(agg, file.path(out_dir, "MuSC_new_pairs_shortlist.csv"), row.names = FALSE)
message("Distinct new LR pairs x direction: ", nrow(agg), " | in recommended set: ", sum(agg$recommended_set),
        " | subunit_of_default (complex written differently): ", sum(agg$novelty == "subunit_of_default"))
print(table(novelty = agg$novelty, recommended = agg$recommended_set))
