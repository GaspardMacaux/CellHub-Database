# patch_release_geneinfo.R — CellHubDB: every partner gene of an interaction is added to geneInfo, nothing else changes.
# CellChat keeps only geneInfo symbols when it extracts the signaling genes (extractGene / extractGeneSubset): a single-gene
# ligand or receptor absent from geneInfo is read as an undefined complex and silently dropped, so its interactions can never
# be scored. In 2.1 / 2.1.1, 6-16 % of the partner genes of the orthology-based species (connectomeDB2025) are missing.
# geneInfo is only a symbol whitelist: adding a symbol does not create any interaction, it only lets CellChat use the existing
# ones when the gene is present in the expression data.
# Usage: Rscript patch_release_geneinfo.R <release_dir> <output_dir> [version, default 2.1.2]
args <- commandArgs(trailingOnly = TRUE)
in_dir <- args[1]; out_dir <- args[2]; new_version <- if (length(args) >= 3) args[3] else "2.1.2"
if (dir.exists(out_dir)) stop("Output folder already exists: ", out_dir)
log_rows <- list()

# Partner genes of a CellChat-format database: single-gene partners, subunits of the complexes used, cofactor genes
partner_genes <- function(db) {
  it <- db$interaction
  partners <- unique(c(it$ligand, it$receptor))
  cx <- db$complex
  in_cx <- intersect(partners, rownames(cx))
  subunits <- if (length(in_cx) > 0) unlist(cx[in_cx, grep("^subunit", colnames(cx)), drop = FALSE], use.names = FALSE) else character(0)
  cof_ids <- unique(unlist(it[, intersect(c("agonist", "antagonist", "co_A_receptor", "co_I_receptor"), colnames(it))], use.names = FALSE))
  cf <- db$cofactor
  cof_ids <- intersect(cof_ids[!is.na(cof_ids) & cof_ids != ""], rownames(cf))
  cof_genes <- if (length(cof_ids) > 0) unlist(cf[cof_ids, grep("^cofactor", colnames(cf)), drop = FALSE], use.names = FALSE) else character(0)
  g <- unique(c(setdiff(partners, rownames(cx)), subunits, cof_genes))
  g[!is.na(g) & nzchar(trimws(g))]
}

# geneInfo completed with the missing partner genes (other columns NA, flagged in geneInfo_source); returns the patched database
patch_db <- function(db, label) {
  gi <- db$geneInfo
  genes <- partner_genes(db)
  missing <- setdiff(genes, gi$Symbol)
  if (length(missing) > 0) {
    add <- as.data.frame(matrix(NA, nrow = length(missing), ncol = ncol(gi), dimnames = list(NULL, colnames(gi))), stringsAsFactors = FALSE)
    add$Symbol <- missing
    if (!"geneInfo_source" %in% colnames(gi)) gi$geneInfo_source <- "reference"
    add$geneInfo_source <- "interaction partner"
    for (cl in colnames(gi)) if (is.factor(gi[[cl]])) gi[[cl]] <- as.character(gi[[cl]])
    gi <- rbind(gi, add[, colnames(gi), drop = FALSE])
    rownames(gi) <- NULL
    db$geneInfo <- gi
  }
  log_rows[[length(log_rows) + 1]] <<- data.frame(file = label, partner_genes = length(genes), added = length(missing),
                                                  pct_added = round(100 * length(missing) / max(1, length(genes)), 2),
                                                  examples = paste(head(missing, 10), collapse = ", "), stringsAsFactors = FALSE)
  db
}

# Same database except geneInfo rows added: interactions, complexes and cofactors identical, original geneInfo rows kept
check_same <- function(a, b, label) {
  ok <- identical(a$interaction, b$interaction) && identical(a$complex, b$complex) && identical(a$cofactor, b$cofactor) &&
    all(a$geneInfo$Symbol %in% b$geneInfo$Symbol) && nrow(b$geneInfo) >= nrow(a$geneInfo)
  if (!ok) stop("Content changed beyond geneInfo additions: ", label)
}

dir.create(out_dir, recursive = TRUE)
invisible(file.copy(list.files(in_dir, full.names = TRUE), out_dir, recursive = TRUE))
for (f in list.files(file.path(out_dir, "app_databases"), pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)) {
  db <- readRDS(f); p <- patch_db(db, sub(paste0("^", out_dir, "/?"), "", f)); check_same(db, p, f)
  for (a in setdiff(names(attributes(db)), "names")) attr(p, a) <- attr(db, a)
  attr(p, "version") <- new_version
  saveRDS(p, f)
}
rel <- readRDS(file.path(out_dir, "CellHubDB.rds"))
rel_attr <- attributes(rel)
patched <- lapply(names(rel), function(sp) { p <- patch_db(rel[[sp]], paste0("CellHubDB.rds$", sp)); check_same(rel[[sp]], p, sp); p })
names(patched) <- names(rel)
for (a in setdiff(names(rel_attr), "names")) attr(patched, a) <- rel_attr[[a]]
attr(patched, "version") <- new_version
attr(patched, "build_date") <- format(Sys.Date())
saveRDS(patched, file.path(out_dir, "CellHubDB.rds"))
# App merged files and release object still identical species by species (interactions and geneInfo)
for (sp in names(patched)) {
  app <- readRDS(file.path(out_dir, "app_databases", sp, paste0("CellHubDB_", sp, "_merged.rds")))
  if (!isTRUE(all.equal(app$interaction, patched[[sp]]$interaction, check.attributes = FALSE)) || !setequal(app$geneInfo$Symbol, patched[[sp]]$geneInfo$Symbol)) {
    stop("App merged file and CellHubDB.rds differ for ", sp)
  }
}
log_df <- do.call(rbind, log_rows)
write.csv(log_df, file.path(out_dir, "reports", paste0("patch_", new_version, "_geneinfo.csv")), row.names = FALSE)
writeLines(c(paste("# CellHubDB", new_version), "",
             "geneInfo completed with every interaction partner gene (single-gene ligands and receptors, subunits of the complexes used,",
             "cofactor genes) missing from it; added rows carry geneInfo_source = \"interaction partner\". CellChat keeps only geneInfo",
             "symbols when it extracts the signaling genes, so interactions of these genes were silently dropped before (6-16 % of the",
             "partner genes of the orthology-based species). Interactions, tiers, sources, complexes and cofactors are unchanged (checked file by file)."),
           file.path(out_dir, paste0("CHANGES_", new_version, ".md")))
files <- list.files(out_dir, recursive = TRUE)
files <- files[files != "MD5SUMS"]
writeLines(paste(unname(tools::md5sum(file.path(out_dir, files))), files, sep = "  "), file.path(out_dir, "MD5SUMS"))
merged_rows <- log_df[grepl("_merged\\.rds$", log_df$file), c("file", "partner_genes", "added", "pct_added")]
print(merged_rows[order(-merged_rows$pct_added), ], row.names = FALSE)
message("CellHubDB ", new_version, " written to ", out_dir)
