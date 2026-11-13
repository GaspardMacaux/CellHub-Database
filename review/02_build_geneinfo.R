# Step 02 (Reviewer 2, other 3): geneInfo for species without a reference table, built from versioned Ensembl GTF files
# (FTP instead of BioMart: reproducible files with MD5, no dependency on BioMart availability) and symbol coverage report
source("review/revision_helpers.R")
db_path <- "review/01_audit/CellHubDB_schema_fixed.rds"; out_dir <- "review/02_geneinfo"
ensembl_release <- 112; plants_release <- 59  # Ensembl Genomes release matching Ensembl 112
dir.create(file.path(out_dir, "gtf"), recursive = TRUE, showWarnings = FALSE)
db <- readRDS(db_path)

# Ensembl FTP species directories (frog = X. tropicalis and sheep = Rambouillet are assumptions to check against connectomeDB2025)
ens_map <- data.frame(species = c("chicken", "chimp", "cow", "dog", "frog", "horse", "macaque", "marmoset", "pig", "rat", "sheep", "zebrafish", "drosophila", "celegans",
                                  "arabidopsis", "maize", "rice", "soybean", "tomato"),
                      ftp_dir = c("gallus_gallus", "pan_troglodytes", "bos_taurus", "canis_lupus_familiaris", "xenopus_tropicalis", "equus_caballus", "macaca_mulatta",
                                  "callithrix_jacchus", "sus_scrofa", "rattus_norvegicus", "ovis_aries_rambouillet", "danio_rerio", "drosophila_melanogaster",
                                  "caenorhabditis_elegans", "arabidopsis_thaliana", "zea_mays", "oryza_sativa", "glycine_max", "solanum_lycopersicum"),
                      division = c(rep("vertebrates", 14), rep("plants", 5)), stringsAsFactors = FALSE)

# Directory URL of the GTF files for one species and division
gtf_dir_url <- function(ftp_dir, division, ensembl_release, plants_release) {
  if (division == "plants") return(sprintf("https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-%d/gtf/%s/", plants_release, ftp_dir))
  return(sprintf("https://ftp.ensembl.org/pub/release-%d/gtf/%s/", ensembl_release, ftp_dir))
}

# Download the main GTF of a species (not the chr/abinitio/patch variants) once; returns the local path
download_species_gtf <- function(dir_url, release, dest_dir) {
  options(timeout = max(1200, getOption("timeout")))
  index <- paste(readLines(dir_url, warn = FALSE), collapse = "\n")
  files <- unique(regmatches(index, gregexpr("[A-Za-z0-9_.\\-]+\\.gtf\\.gz", index))[[1]])
  files <- files[grepl(paste0("\\.", release, "\\.gtf\\.gz$"), files) & !grepl("chr|abinitio|patch|scaff", files)]
  if (length(files) != 1) stop("Expected one main GTF in ", dir_url, ", found: ", paste(files, collapse = ", "))
  dest <- file.path(dest_dir, files)
  if (!file.exists(dest)) utils::download.file(paste0(dir_url, files), dest, mode = "wb", quiet = TRUE)
  return(dest)
}

# Gene table (gene_id, gene_name, gene_biotype) from the 'gene' lines of a GTF
read_gtf_genes <- function(gtf_path) {
  lines <- data.table::fread(cmd = paste("zcat", shQuote(gtf_path), "| awk -F'\\t' '$3==\"gene\"'"), sep = "\t", header = FALSE, quote = "", select = 9)[[1]]
  attr_value <- function(key) { v <- regmatches(lines, regexpr(paste0(key, " \"[^\"]*\""), lines)); out <- rep(NA_character_, length(lines)); hit <- grepl(paste0(key, " \""), lines); out[hit] <- sub(paste0(key, " \"([^\"]*)\""), "\\1", v); out }
  genes <- data.frame(ensembl_gene_id = attr_value("gene_id"), external_gene_name = attr_value("gene_name"), gene_biotype = attr_value("gene_biotype"), stringsAsFactors = FALSE)
  genes$description <- NA_character_  # GTF files carry no gene description
  return(genes)
}

# Reference identifiers rewritten to the CellHubDB format when the pipeline changed it
# (soybean: Ensembl GLYMA_07G056900 -> Glyma.07G056900, because CellChat reads "_" as a complex separator)
normalize_reference_ids <- function(ens, species) {
  if (species == "soybean") ens$ensembl_gene_id <- sub("^GLYMA_", "Glyma.", ens$ensembl_gene_id)
  return(ens)
}

# geneInfo keyed on the identifier type (symbol or gene ID) that best matches the database tokens
build_geneinfo <- function(ens, db_genes) {
  cov_symbol <- mean(toupper(db_genes) %in% toupper(ens$external_gene_name))
  cov_id <- mean(toupper(db_genes) %in% toupper(ens$ensembl_gene_id))
  key <- if (cov_id > cov_symbol) "ensembl_gene_id" else "external_gene_name"
  sym <- ens[[key]]; keep <- !is.na(sym) & sym != "" & !duplicated(sym)
  gi <- data.frame(Symbol = sym[keep], Name = ens$description[keep], Ensembl.Gene.ID = ens$ensembl_gene_id[keep], Biotype = ens$gene_biotype[keep], stringsAsFactors = FALSE)
  matched_exact <- db_genes %in% gi$Symbol
  matched_ci <- toupper(db_genes) %in% toupper(gi$Symbol)
  return(list(geneInfo = gi, key = key, coverage_symbol = cov_symbol, coverage_id = cov_id, coverage_exact = mean(matched_exact),
              unmatched = db_genes[!matched_ci], case_only = db_genes[matched_ci & !matched_exact]))
}

report <- list(); gtf_manifest <- list()
for (i in seq_len(nrow(ens_map))) {
  sp <- ens_map$species[i]
  if (!sp %in% names(db)) next
  message("geneInfo: ", sp)
  db_genes <- extract_partner_genes(db[[sp]])
  rel <- if (ens_map$division[i] == "plants") plants_release else ensembl_release
  dir_url <- gtf_dir_url(ens_map$ftp_dir[i], ens_map$division[i], ensembl_release, plants_release)
  gtf <- tryCatch(download_species_gtf(dir_url, rel, file.path(out_dir, "gtf")), error = function(e) { message("  GTF download failed: ", conditionMessage(e)[1]); NULL })
  if (is.null(gtf)) { report[[sp]] <- data.frame(species = sp, status = "gtf_failed"); next }
  gtf_manifest[[sp]] <- data.frame(species = sp, file = basename(gtf), url = paste0(dir_url, basename(gtf)), release = rel, md5 = unname(tools::md5sum(gtf)), stringsAsFactors = FALSE)
  gi <- build_geneinfo(normalize_reference_ids(read_gtf_genes(gtf), sp), db_genes)
  writeLines(gi$unmatched, file.path(out_dir, paste0("unmatched_", sp, ".txt")))
  writeLines(gi$case_only, file.path(out_dir, paste0("case_mismatch_", sp, ".txt")))
  # Empty geneInfo tables are replaced by the Ensembl table; existing ones (zebrafish from CellChatDB) are completed with the Ensembl
  # symbols they lack, because the CellChatDB zebrafish table misses valid ZFIN genes used by connectomeDB2025 (b2m, c1qb, areg...)
  n_before <- nrow(db[[sp]]$geneInfo)
  if (n_before == 0) {
    db[[sp]]$geneInfo <- gi$geneInfo; mode <- "replaced"
  } else {
    old <- db[[sp]]$geneInfo; add <- gi$geneInfo[!gi$geneInfo$Symbol %in% old$Symbol, , drop = FALSE]
    for (cl in setdiff(colnames(old), colnames(add))) add[[cl]] <- NA
    for (cl in setdiff(colnames(add), colnames(old))) old[[cl]] <- NA
    db[[sp]]$geneInfo <- rbind(old, add[, colnames(old), drop = FALSE]); mode <- "completed"
  }
  report[[sp]] <- data.frame(species = sp, status = "ok", key = gi$key, n_db_genes = length(db_genes), pct_symbol = round(100 * gi$coverage_symbol, 1),
                             pct_id = round(100 * gi$coverage_id, 1), pct_exact_case = round(100 * gi$coverage_exact, 1), n_unmatched = length(gi$unmatched),
                             n_case_mismatch = length(gi$case_only), geneinfo_mode = mode, n_geneinfo_before = n_before, n_geneinfo_after = nrow(db[[sp]]$geneInfo))
}
report_cols <- c("species", "status", "key", "n_db_genes", "pct_symbol", "pct_id", "pct_exact_case", "n_unmatched", "n_case_mismatch", "geneinfo_mode", "n_geneinfo_before", "n_geneinfo_after")
report_df <- do.call(rbind, lapply(report, function(r) { for (cl in setdiff(report_cols, names(r))) r[[cl]] <- NA; r[, report_cols] }))
write.csv(report_df, file.path(out_dir, "geneinfo_coverage.csv"), row.names = FALSE)
if (length(gtf_manifest) > 0) write.table(do.call(rbind, gtf_manifest), file.path(out_dir, "gtf_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(report_df)
# Interactions that CellChat can score: every gene of both partners present in the final geneInfo (exact case), all species
inter_cov <- do.call(rbind, lapply(names(db), function(sp) {
  it <- db[[sp]]$interaction; cpx <- db[[sp]]$complex; ref <- db[[sp]]$geneInfo$Symbol
  scorable <- vapply(seq_len(nrow(it)), function(i) all(c(partner_subunits(it$ligand[i], cpx), partner_subunits(it$receptor[i], cpx)) %in% ref), logical(1))
  missing_tokens <- setdiff(extract_partner_genes(db[[sp]]), ref)
  writeLines(missing_tokens, file.path(out_dir, paste0("not_in_geneinfo_", sp, ".txt")))
  by_src <- table(it$source_db[!scorable])
  data.frame(species = sp, n_interactions = nrow(it), n_scorable = sum(scorable), pct_scorable = round(100 * mean(scorable), 1),
             n_tokens_missing = length(missing_tokens), unscorable_by_source = paste(names(by_src), by_src, sep = ":", collapse = " | "),
             examples_missing = paste(head(missing_tokens, 8), collapse = " | "), stringsAsFactors = FALSE)
}))
write.csv(inter_cov, file.path(out_dir, "interaction_coverage.csv"), row.names = FALSE)
print(inter_cov[, c("species", "n_interactions", "n_scorable", "pct_scorable", "n_tokens_missing")], row.names = FALSE)
saveRDS(db, file.path(out_dir, "CellHubDB_geneinfo.rds"))
