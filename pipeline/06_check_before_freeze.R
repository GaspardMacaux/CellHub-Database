# Diagnostics before freezing CellHubDB v2.0:
# (1) auto-generated complexes in species without a curated complex table (maize, soybean, marmoset)
# (2) reciprocal pairs (A->B and B->A) inside single-source files, collapsed by the order-insensitive merge key
root <- "data/databases"
for (sp in c("maize", "soybean", "marmoset")) {
  db <- readRDS(file.path(root, sp, sprintf("CellHubDB_%s_merged.rds", sp)))
  cx <- db$complex
  cat(sprintf("\n== %s: %d complexes ==\n", sp, nrow(cx)))
  print(head(cx[, intersect(c("subunit_1", "subunit_2", "subunit_3"), colnames(cx)), drop = FALSE], 8))
  cat("Most frequent subunit_1:\n")
  print(head(sort(table(cx$subunit_1), decreasing = TRUE), 5))
}
cat("\n== Reciprocal pairs (A->B and B->A) inside single-source files ==\n")
files <- list.files(root, pattern = "^CellHubDB_.*\\.rds$", recursive = TRUE, full.names = TRUE)
files <- files[!grepl("_merged\\.rds$", files)]
res <- do.call(rbind, lapply(files, function(f) {
  it <- readRDS(f)$interaction
  lig <- tolower(it$ligand); rec <- tolower(it$receptor)
  fwd <- paste(lig, rec, sep = "|"); rev <- paste(rec, lig, sep = "|")
  # Each reciprocal pair is counted twice (once per direction); homodimers A->A are excluded
  data.frame(file = basename(f), interactions = nrow(it), reciprocal_pairs = sum(fwd %in% rev & lig != rec) / 2)
}))
res <- res[res$reciprocal_pairs > 0, ]
print(res[order(-res$reciprocal_pairs), ], row.names = FALSE)
# Examples for the largest contributor
top <- readRDS(files[basename(files) == res$file[which.max(res$reciprocal_pairs)]])$interaction
fwd <- paste(tolower(top$ligand), tolower(top$receptor), sep = "|")
rev <- paste(tolower(top$receptor), tolower(top$ligand), sep = "|")
cat("\nExamples:\n")
print(head(top[fwd %in% rev & tolower(top$ligand) != tolower(top$receptor), c("ligand", "receptor", "interaction_type")], 10), row.names = FALSE)
