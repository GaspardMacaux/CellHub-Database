# Step 11: CellHubDB numbers report — every figure computed from the released file (single source of truth), plus source versions
# and dates (manifest), Ensembl releases (GTF manifest), repairs and release gate result when those files exist
# Output: review/11_report/<version>/CellHubDB_<version>_numbers.xlsx (one sheet per table) and CellHubDB_<version>_report.md
# Usage: Rscript review/11_database_report.R [release rds] [version label]
suppressPackageStartupMessages(library(openxlsx))
args <- commandArgs(trailingOnly = TRUE)
release_file <- if (length(args) >= 1) args[1] else "release_v2.1/CellHubDB.rds"
db <- readRDS(release_file)
version <- if (length(args) >= 2) args[2] else if (!is.null(attr(db, "version"))) attr(db, "version") else "unversioned"
out_dir <- file.path("review/11_report", version); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
manifest_file <- "review/00_sources/sources_manifest.tsv"
# Pipeline outputs first (geneinfo/, reports/ from finalize_cellhubdb.R), revision-script outputs otherwise
gtf_file <- if (file.exists("geneinfo/gtf_manifest.tsv")) "geneinfo/gtf_manifest.tsv" else "review/02_geneinfo/gtf_manifest.tsv"
removed_file <- if (file.exists("reports/finalize_removed_interactions.csv")) "reports/finalize_removed_interactions.csv" else "review/07b_repair/removed_interactions.csv"; gate_file <- file.path("review/09_gate", paste0("v", sub("^v", "", version)), "gate_report.csv")
canon <- c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor", "Non-protein Signaling")
placeholder <- "Unassigned"

# Genes of a partner: complex subunits when the partner is a defined complex, otherwise the partner itself
partner_genes <- function(p, cpx) {
  if (!is.null(cpx) && p %in% rownames(cpx)) { g <- as.character(unlist(cpx[p, , drop = TRUE])); return(g[!is.na(g) & g != ""]) }
  return(p)
}

# Column of the interaction table, or NA of the right length when absent
col_or_na <- function(it, cl) if (cl %in% colnames(it)) it[[cl]] else rep(NA, nrow(it))

# Source names of each interaction (comma-separated "sources" column, falling back on source_db)
source_lists <- function(it) {
  s <- if ("sources" %in% colnames(it)) it$sources else it$source_db
  return(lapply(strsplit(ifelse(is.na(s), "", s), ",\\s*"), function(x) x[x != ""]))
}

# Percentage rounded to one decimal, 0 when the denominator is 0
pct <- function(x, n) if (n == 0) 0 else round(100 * x / n, 1)

# One row of figures per species
species_row <- function(sp, x) {
  it <- x$interaction; cpx <- x$complex; gi <- x$geneInfo$Symbol; n <- nrow(it)
  lig_genes <- lapply(it$ligand, partner_genes, cpx = cpx); rec_genes <- lapply(it$receptor, partner_genes, cpx = cpx)
  scorable <- mapply(function(l, r) all(c(l, r) %in% gi), lig_genes, rec_genes)
  srcs <- source_lists(it); pw <- col_or_na(it, "pathway_name"); type <- col_or_na(it, "interaction_type")
  tier <- col_or_na(it, "confidence_tier"); ev <- col_or_na(it, "evidence_type"); origin <- col_or_na(it, "pathway_origin")
  partners <- unique(c(it$ligand, it$receptor))
  data.frame(species = sp, interactions = n, unique_ligands = length(unique(it$ligand)), unique_receptors = length(unique(it$receptor)),
             unique_genes = length(unique(unlist(c(lig_genes, rec_genes)))), complexes_used = sum(partners %in% rownames(cpx)),
             cofactors_defined = if (is.null(x$cofactor)) 0L else nrow(x$cofactor), geneinfo_genes = length(gi),
             pct_scorable = pct(sum(scorable), n), sources = length(unique(unlist(srcs))),
             multi_source = sum(lengths(srcs) >= 2), pct_multi_source = pct(sum(lengths(srcs) >= 2), n),
             pathways = length(unique(pw[!is.na(pw) & pw != placeholder])), pct_pathway_unassigned = pct(sum(is.na(pw) | pw == placeholder), n),
             pathway_from_source = sum(origin %in% "source"), pathway_inherited = sum(grepl("^inherited", origin)),
             secreted = sum(type %in% "Secreted Signaling"), cell_contact = sum(type %in% "Cell-Cell Contact"), ecm_receptor = sum(type %in% "ECM-Receptor"),
             non_protein = sum(type %in% "Non-protein Signaling"), type_unlabelled = sum(!type %in% canon),
             tier_high = sum(tier %in% "high"), tier_medium = sum(tier %in% "medium"), tier_low = sum(tier %in% "low"),
             role_corroborated = sum(col_or_na(it, "role_corroborated") %in% TRUE), predicted_only = sum(col_or_na(it, "predicted_only") %in% TRUE),
             complex_interactions = sum(col_or_na(it, "is_complex") %in% TRUE),
             evidence_direct = sum(ev %in% "Direct"), evidence_inferred = sum(ev %in% "Inferred"), stringsAsFactors = FALSE)
}

# Interactions per species and source: supported by the source, and supported by this source only
source_rows <- function(sp, x) {
  srcs <- source_lists(x$interaction); all_src <- sort(unique(unlist(srcs)))
  do.call(rbind, lapply(all_src, function(s) {
    has <- vapply(srcs, function(v) s %in% v, logical(1))
    data.frame(species = sp, source = s, interactions = sum(has), only_this_source = sum(has & lengths(srcs) == 1),
               pct_of_species = pct(sum(has), length(srcs)), stringsAsFactors = FALSE)
  }))
}

species_tab <- do.call(rbind, lapply(names(db), function(sp) species_row(sp, db[[sp]])))
# Totals only for additive counts (sources, pathways and percentages are not additive across species)
num_cols <- setdiff(names(species_tab), c("species", "sources", "pathways", grep("^pct_", names(species_tab), value = TRUE)))
totals <- species_tab[1, ]; totals[] <- NA; totals$species <- "TOTAL (sum over species)"; totals[num_cols] <- colSums(species_tab[num_cols])
species_tab <- rbind(species_tab, totals)
source_tab <- do.call(rbind, lapply(names(db), function(sp) source_rows(sp, db[[sp]])))
source_totals <- aggregate(cbind(interactions, only_this_source) ~ source, data = source_tab, FUN = sum)
source_totals$species <- vapply(source_totals$source, function(s) paste(source_tab$species[source_tab$source == s], collapse = ", "), character(1))
names(source_totals)[names(source_totals) == "species"] <- "species_covered"

# Source versions and dates from the hand-filled manifest (one row per source)
versions_tab <- NULL
if (file.exists(manifest_file)) {
  man <- read.delim(manifest_file, stringsAsFactors = FALSE, colClasses = "character")
  man$source[is.na(man$source) | man$source == ""] <- "(not filled)"
  versions_tab <- do.call(rbind, lapply(split(man, man$source), function(m) {
    uniq <- function(v) paste(unique(v[!is.na(v) & v != ""]), collapse = "; ")
    data.frame(source = m$source[1], version = uniq(m$version), url = uniq(m$url), n_files = nrow(m),
               file_dates = paste(range(m$modified), collapse = " to "), total_size_mb = round(sum(as.numeric(m$size_bytes)) / 1e6, 1), stringsAsFactors = FALSE)
  }))
}
# R packages from which sources are extracted, as installed in this R session
pkg_tab <- do.call(rbind, lapply(c("CellChat", "NeuronChat", "nichenetr", "multinichenetr", "Seurat"), function(p) {
  if (!requireNamespace(p, quietly = TRUE)) return(data.frame(package = p, version = "not installed", commit = "", stringsAsFactors = FALSE))
  d <- utils::packageDescription(p)
  data.frame(package = p, version = d$Version, commit = if (is.null(d$RemoteSha)) "" else substr(d$RemoteSha, 1, 10), stringsAsFactors = FALSE)
}))
gtf_tab <- if (file.exists(gtf_file)) read.delim(gtf_file, stringsAsFactors = FALSE) else NULL
removed_tab <- NULL
if (file.exists(removed_file)) {
  rem <- read.csv(removed_file, stringsAsFactors = FALSE)
  # The finalize log lists every file; only the merged databases are counted
  if ("file" %in% colnames(rem)) rem <- rem[grepl("_merged\\.rds$", rem$file), ]
  if (nrow(rem) > 0) { removed_tab <- as.data.frame(table(species = rem$species, reason = rem$reason), stringsAsFactors = FALSE); removed_tab <- removed_tab[removed_tab$Freq > 0, ]; names(removed_tab)[3] <- "interactions_removed" }
}
gate_tab <- if (file.exists(gate_file)) read.csv(gate_file, stringsAsFactors = FALSE) else NULL

# Excel workbook: one sheet per table, header bold, columns sized
wb <- createWorkbook()
add_sheet <- function(name, df) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  addWorksheet(wb, name); writeData(wb, name, df, headerStyle = createStyle(textDecoration = "bold"))
  setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto"); freezePane(wb, name, firstRow = TRUE)
  return(invisible(NULL))
}
info <- data.frame(item = c("Release file", "Version", "MD5", "Report date", "Species", "Interactions (sum over species)"),
                   value = c(release_file, version, unname(tools::md5sum(release_file)), format(Sys.Date()), length(db), sum(species_tab$interactions[species_tab$species != "TOTAL (sum over species)"])), stringsAsFactors = FALSE)
add_sheet("Release", info); add_sheet("Per species", species_tab); add_sheet("Per species x source", source_tab); add_sheet("Per source", source_totals)
add_sheet("Source versions", versions_tab); add_sheet("R packages", pkg_tab); add_sheet("Ensembl GTF", gtf_tab); add_sheet("Removed (repairs)", removed_tab)
if (!is.null(gate_tab)) add_sheet("Release gate", gate_tab[gate_tab$level != "PASS", ])
xlsx_path <- file.path(out_dir, paste0("CellHubDB_", version, "_numbers.xlsx")); saveWorkbook(wb, xlsx_path, overwrite = TRUE)

# Markdown report (convert to Word with: pandoc report.md -o report.docx)
md_table <- function(df) {
  if (is.null(df) || nrow(df) == 0) return("_No data._\n")
  df[] <- lapply(df, function(v) ifelse(is.na(v), "", if (is.numeric(v)) format(v, big.mark = ",", scientific = FALSE, trim = TRUE) else as.character(v)))
  return(paste0(c(paste0("| ", paste(names(df), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
                  apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))), collapse = "\n"))
}
main_cols <- c("species", "interactions", "unique_genes", "sources", "multi_source", "pathways", "pct_pathway_unassigned", "pct_scorable", "tier_high", "tier_medium", "tier_low")
type_cols <- c("species", "secreted", "cell_contact", "ecm_receptor", "non_protein", "type_unlabelled", "evidence_direct", "evidence_inferred")
gate_line <- if (is.null(gate_tab)) "Release gate report not found." else paste0("Release gate: ", paste(names(table(gate_tab$level)), table(gate_tab$level), collapse = ", "), ".")
md <- c(paste("# CellHub Database", version, "— numbers report"), "",
        paste0("Release file: `", release_file, "` (MD5 ", info$value[3], "). Report generated ", format(Sys.Date()), ". ", gate_line), "",
        "## Overview per species", "", md_table(species_tab[, intersect(main_cols, names(species_tab))]), "",
        "## Interaction types and evidence", "", md_table(species_tab[, intersect(type_cols, names(species_tab))]), "",
        "## Interactions per source (all species)", "", md_table(source_totals), "",
        "## Source versions and dates", "", if (is.null(versions_tab)) "_Manifest not found._" else md_table(versions_tab[, c("source", "version", "file_dates", "n_files", "url")]), "",
        "## R packages used to extract sources", "", md_table(pkg_tab), "",
        "## Ensembl gene references (geneInfo)", "", if (is.null(gtf_tab)) "_GTF manifest not found._" else md_table(gtf_tab[, intersect(c("species", "release", "file"), names(gtf_tab))]), "",
        "## Interactions removed by the release repairs", "", md_table(removed_tab), "",
        "## Release gate: failures and warnings", "", if (is.null(gate_tab)) "_Not found._" else md_table(gate_tab[gate_tab$level != "PASS", intersect(c("level", "species", "file", "check", "value"), names(gate_tab))]), "",
        "Full tables (per species x source, all gate checks): ", basename(xlsx_path))
md_path <- file.path(out_dir, paste0("CellHubDB_", version, "_report.md")); writeLines(md, md_path)
print(species_tab[, intersect(main_cols, names(species_tab))], row.names = FALSE)
message("Written: ", xlsx_path, "\n         ", md_path)
