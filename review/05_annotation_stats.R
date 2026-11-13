# Step 05 (Reviewer 2, other 5 + 6): annotation categories (annotation vs interaction_type), pathway coverage and PMID formatting per species
source("review/revision_helpers.R")
db_path <- "review/03_confidence/CellHubDB_annotated.rds"; out_dir <- "review/05_annotation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
db <- readRDS(db_path)
canon <- c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor", "Non-protein Signaling")
placeholder <- "Unassigned"  # pathway label given in step 02b to interactions without any pathway

# Category counts of one annotation column, with Other and Unlabelled bins
count_categories <- function(x) {
  a <- ifelse(is.na(x) | x == "", "Unlabelled", x); a[!a %in% c(canon, "Unlabelled")] <- "Other"
  return(as.integer(table(factor(a, levels = c(canon, "Other", "Unlabelled")))))
}
# Evidence strings made only of PMID:<digits> references
pmid_ok <- function(ev) !is.na(ev) & grepl("^\\s*PMID:\\d+(\\s*[;,]\\s*PMID:\\d+)*\\s*$", ev)

stats <- do.call(rbind, lapply(names(db), function(sp) {
  it <- db[[sp]]$interaction
  rows <- lapply(intersect(c("annotation", "interaction_type"), colnames(it)), function(cl) {
    data.frame(species = sp, column = cl, t(setNames(count_categories(it[[cl]]), c(canon, "Other", "Unlabelled"))), check.names = FALSE)
  })
  out <- do.call(rbind, rows)
  out$pct_pathway_unassigned <- round(100 * mean(is.na(it$pathway_name) | it$pathway_name %in% c("", placeholder)), 1)
  out$pct_pathway_inherited <- if ("pathway_origin" %in% colnames(it)) round(100 * mean(grepl("^inherited", it$pathway_origin)), 1) else NA
  out$pct_evidence_na <- round(100 * mean(is.na(it$evidence) | it$evidence == ""), 1)
  out$pct_evidence_pmid <- round(100 * mean(pmid_ok(it$evidence)), 1)
  return(out)
}))
write.csv(stats, file.path(out_dir, "annotation_pathway_evidence_stats.csv"), row.names = FALSE)
print(stats)

# CellChatDB reference counts and CellChat-derived rows in the merged human database (origin of the 6 vs 424 ECM-Receptor discrepancy)
ref_h <- table(CellChat::CellChatDB.human$interaction$annotation)
hum <- db$human$interaction; from_cc <- grepl("CellChatDB", hum$sources)
print(list(CellChatDB_human = ref_h, merged_human_from_CellChatDB_annotation = table(hum$annotation[from_cc], useNA = "ifany"),
           merged_human_from_CellChatDB_interaction_type = if ("interaction_type" %in% colnames(hum)) table(hum$interaction_type[from_cc], useNA = "ifany") else NULL))

# Examples of non-PMID evidence strings per species
non_pmid <- do.call(rbind, lapply(names(db), function(sp) {
  ev <- db[[sp]]$interaction$evidence; bad <- unique(ev[!is.na(ev) & ev != "" & !pmid_ok(ev)])
  if (length(bad) == 0) return(NULL)
  data.frame(species = sp, n_distinct_non_pmid = length(bad), examples = paste(head(bad, 5), collapse = " || "))
}))
write.csv(non_pmid, file.path(out_dir, "non_pmid_evidence_examples.csv"), row.names = FALSE)
print(non_pmid)
