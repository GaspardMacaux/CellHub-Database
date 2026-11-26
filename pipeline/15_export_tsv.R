# Step 15: tab-separated export of a CellHubDB release, for use outside R (Python: pandas, LIANA+, pycellchat).
# Usage: Rscript pipeline/15_export_tsv.R [release CellHubDB rds] [output folder]
#   default release/2.1.2/CellHubDB_2.1.2.rds (or release/2.1.2/CellHubDB.rds) -> release/2.1.2/tsv/
# Output, one folder per species: interaction.tsv, complex.tsv, cofactor.tsv, geneInfo.tsv (the four CellChat tables, complex
# and cofactor names as a first column), and resource_liana.tsv (one row per interaction, partners written as their gene
# symbols joined by "_", the complex notation of LIANA+). Plus README.md, MD5SUMS and a zip of the whole folder.
args <- commandArgs(trailingOnly = TRUE)
release_file <- if (length(args) >= 1) args[1] else {
  cand <- c("release/2.1.2/CellHubDB_2.1.2.rds", "release/2.1.2/CellHubDB.rds")
  cand[file.exists(cand)][1]
}
if (is.na(release_file) || !file.exists(release_file)) stop("Release file not found: give its path as first argument")
db <- readRDS(release_file)
version <- if (!is.null(attr(db, "version"))) attr(db, "version") else "unversioned"
out_dir <- if (length(args) >= 2) args[2] else file.path(dirname(release_file), "tsv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Table written with row names as a named first column (complex and cofactor tables are keyed by row name in CellChat)
write_tsv <- function(df, file, rowname_col = NULL) {
  if (is.null(df)) df <- data.frame()
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df[] <- lapply(df, function(x) if (is.factor(x)) as.character(x) else x)
  if (!is.null(rowname_col)) df <- data.frame(setNames(list(rownames(df)), rowname_col), df, check.names = FALSE, stringsAsFactors = FALSE)
  utils::write.table(df, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

# Gene symbols of a partner: subunits of the complex table when the partner is a complex, otherwise the partner itself
partner_genes <- function(partners, complex) {
  has_cpx <- !is.null(complex) && nrow(complex) > 0
  vapply(as.character(partners), function(p) {
    if (has_cpx && p %in% rownames(complex)) {
      g <- unlist(complex[p, ], use.names = FALSE)
      g <- g[!is.na(g) & nzchar(g)]
      if (length(g) > 0) return(paste(g, collapse = "_"))
    }
    p
  }, character(1), USE.NAMES = FALSE)
}

# Characters that would break a tab-separated file are replaced in text columns
clean_text <- function(df) {
  df[] <- lapply(df, function(x) if (is.character(x)) gsub("[\t\r\n]+", " ", x) else x)
  df
}

summary_rows <- list()
for (sp in names(db)) {
  x <- db[[sp]]; sp_dir <- file.path(out_dir, sp); dir.create(sp_dir, showWarnings = FALSE)
  it <- clean_text(as.data.frame(x$interaction, stringsAsFactors = FALSE))
  it[] <- lapply(it, function(v) if (is.factor(v)) as.character(v) else v)
  it$ligand_genes <- partner_genes(it$ligand, x$complex)
  it$receptor_genes <- partner_genes(it$receptor, x$complex)
  rownames(it) <- NULL
  write_tsv(it, file.path(sp_dir, "interaction.tsv"))
  write_tsv(clean_text(x$complex), file.path(sp_dir, "complex.tsv"), rowname_col = "complex_name")
  write_tsv(clean_text(x$cofactor), file.path(sp_dir, "cofactor.tsv"), rowname_col = "cofactor_name")
  write_tsv(clean_text(x$geneInfo), file.path(sp_dir, "geneInfo.tsv"))
  # LIANA+ resource: columns ligand and receptor (complex subunits joined by "_"), with the CellHubDB annotation kept
  keep <- intersect(c("interaction_name", "pathway_name", "annotation", "confidence_tier", "role_corroborated", "sources", "evidence"), colnames(it))
  liana <- data.frame(ligand = it$ligand_genes, receptor = it$receptor_genes, it[, keep, drop = FALSE], stringsAsFactors = FALSE, check.names = FALSE)
  # Several CellChat interactions can share the same gene pair (complexes written differently, agonist / cofactor variants):
  # one row per pair is kept, the best-supported one (tier high > medium > low, then corroborated roles)
  if (nrow(liana) > 0) {
    tier_rank <- if ("confidence_tier" %in% colnames(liana)) match(tolower(liana$confidence_tier), c("high", "medium", "low"), nomatch = 4L) else rep(4L, nrow(liana))
    role_rank <- if ("role_corroborated" %in% colnames(liana)) ifelse(as.logical(liana$role_corroborated) %in% TRUE, 0L, 1L) else rep(1L, nrow(liana))
    liana <- liana[order(tier_rank, role_rank, seq_len(nrow(liana))), , drop = FALSE]
    liana <- liana[!duplicated(toupper(paste(liana$ligand, liana$receptor))), , drop = FALSE]
    rownames(liana) <- NULL
  }
  write_tsv(liana, file.path(sp_dir, "resource_liana.tsv"))
  summary_rows[[sp]] <- data.frame(species = sp, interactions = nrow(it), complexes = if (is.null(x$complex)) 0 else nrow(x$complex),
                                   cofactors = if (is.null(x$cofactor)) 0 else nrow(x$cofactor), genes = if (is.null(x$geneInfo)) 0 else nrow(x$geneInfo),
                                   liana_pairs = nrow(liana), stringsAsFactors = FALSE)
}
summ <- do.call(rbind, summary_rows)
write_tsv(summ, file.path(out_dir, "species_summary.tsv"))

readme <- c(
  sprintf("# CellHub Database %s: tab-separated tables", version), "",
  sprintf("Exported from %s (MD5 %s) by pipeline/15_export_tsv.R. Same content as the R object; UTF-8, tab-separated, empty cells for missing values.",
          basename(release_file), unname(tools::md5sum(release_file))), "",
  "One folder per species:", "",
  "- interaction.tsv: the CellChat interaction table (ligand, receptor, pathway_name, annotation, evidence, sources, confidence_tier, role_corroborated, ...),",
  "  plus ligand_genes and receptor_genes: the gene symbols of each partner, complex subunits joined by \"_\".",
  "- complex.tsv: complex_name and its subunits; cofactor.tsv: cofactor_name and its genes; geneInfo.tsv: gene reference table.",
  "- resource_liana.tsv: one row per ligand-receptor pair with columns ligand and receptor in the LIANA+ complex notation, and the CellHubDB annotation.", "",
  "Python example (pandas, LIANA+):", "",
  "```python",
  "import pandas as pd",
  "res = pd.read_csv(\"mouse/resource_liana.tsv\", sep=\"\\t\")",
  "# recommended set for human and mouse: high tier + medium tier with corroborated roles",
  "rec = res[(res.confidence_tier == \"high\") | ((res.confidence_tier == \"medium\") & (res.role_corroborated == True))]",
  "import liana as li",
  "li.mt.rank_aggregate(adata, groupby=\"cell_type\", resource=rec[[\"ligand\", \"receptor\"]], use_raw=False)",
  "```", "",
  "Species summary: species_summary.tsv. Checksums: MD5SUMS."
)
writeLines(readme, file.path(out_dir, "README.md"))
files <- list.files(out_dir, recursive = TRUE)
files <- files[files != "MD5SUMS"]
writeLines(paste(unname(tools::md5sum(file.path(out_dir, files))), files, sep = "  "), file.path(out_dir, "MD5SUMS"))
zip_file <- file.path(dirname(out_dir), sprintf("CellHubDB_%s_tsv.zip", version))
unlink(zip_file)
old <- setwd(dirname(out_dir)); utils::zip(basename(zip_file), basename(out_dir), flags = "-rq9X"); setwd(old)
print(summ, row.names = FALSE)
message("Written: ", normalizePath(out_dir), " and ", normalizePath(zip_file), " (MD5 ", unname(tools::md5sum(zip_file)), ")")
