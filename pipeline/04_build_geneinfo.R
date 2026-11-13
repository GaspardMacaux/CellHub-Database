################################################################################
#        CellHub Database — genome-wide geneInfo tables from Ensembl GTF       #
#                                                                              #
# Run after the convert_*.R scripts and before curate_CellHub_DB.R.            #
# Builds one geneInfo table per species from versioned Ensembl GTF files       #
# (Ensembl 112, Ensembl Plants 59). Existing tables (CellChatDB geneInfo for   #
# human, mouse, zebrafish) are completed with the Ensembl symbols they lack,   #
# so that symbol validation in curate_CellHub_DB.R (step 6B) uses a            #
# genome-wide reference (DEFB4A, GPR15LG are absent from CellChatDB geneInfo). #
#                                                                              #
# Input:  converted/*.rds, geneinfo/geneinfo_<species>.rds (optional)          #
# Output: geneinfo/geneinfo_<species>.rds, geneinfo/gtf/*.gtf.gz,              #
#         geneinfo/gtf_manifest.tsv, geneinfo/geneinfo_coverage.csv            #
################################################################################

CONVERTED_DIR <- "data/converted"
GENEINFO_DIR <- "data/geneinfo"
GTF_DIR <- file.path(GENEINFO_DIR, "gtf")
ENSEMBL_RELEASE <- 112
PLANTS_RELEASE <- 59  # Ensembl Genomes release matching Ensembl 112
dir.create(GTF_DIR, recursive = TRUE, showWarnings = FALSE)

# Ensembl FTP species directories (frog = X. tropicalis, sheep = Rambouillet)
ENS_MAP <- data.frame(
  species = c("human", "mouse", "chicken", "chimp", "cow", "dog", "frog", "horse", "macaque", "marmoset", "pig", "rat", "sheep",
              "zebrafish", "drosophila", "celegans", "arabidopsis", "maize", "rice", "soybean", "tomato"),
  ftp_dir = c("homo_sapiens", "mus_musculus", "gallus_gallus", "pan_troglodytes", "bos_taurus", "canis_lupus_familiaris", "xenopus_tropicalis",
              "equus_caballus", "macaca_mulatta", "callithrix_jacchus", "sus_scrofa", "rattus_norvegicus", "ovis_aries_rambouillet", "danio_rerio",
              "drosophila_melanogaster", "caenorhabditis_elegans", "arabidopsis_thaliana", "zea_mays", "oryza_sativa", "glycine_max", "solanum_lycopersicum"),
  division = c(rep("vertebrates", 16), rep("plants", 5)), stringsAsFactors = FALSE)

# Print a timestamped log message
log_msg <- function(msg) {
  cat(sprintf("[%s] GENEINFO | %s\n", format(Sys.time(), "%H:%M:%S"), msg))
  return(invisible(NULL))
}

# Directory URL of the GTF files for one species and division
gtf_dir_url <- function(ftp_dir, division) {
  if (division == "plants") return(sprintf("https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-%d/gtf/%s/", PLANTS_RELEASE, ftp_dir))
  return(sprintf("https://ftp.ensembl.org/pub/release-%d/gtf/%s/", ENSEMBL_RELEASE, ftp_dir))
}

# Download the main GTF of a species (not the chr/abinitio/patch/scaffold variants) once; returns the local path
download_species_gtf <- function(dir_url, release, dest_dir) {
  # Cached GTF of this species and release used first (no network needed to rebuild; the MD5 goes to gtf_manifest.tsv)
  ftp_dir <- basename(sub("/$", "", dir_url))
  prefix <- paste0(toupper(substr(ftp_dir, 1, 1)), substring(ftp_dir, 2))
  cached <- list.files(dest_dir, pattern = paste0("^", prefix, "\\..*\\.", release, "\\.gtf\\.gz$"))
  cached <- cached[!grepl("chr|abinitio|patch|scaff", cached)]
  if (length(cached) == 1) return(file.path(dest_dir, cached))
  if (length(cached) > 1) stop("Several cached GTFs for ", ftp_dir, " release ", release, ": ", paste(cached, collapse = ", "))
  # Not cached: list the Ensembl directory and download the main GTF
  options(timeout = max(1200, getOption("timeout")))
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  index <- paste(readLines(dir_url, warn = FALSE), collapse = "\n")
  files <- unique(regmatches(index, gregexpr("[A-Za-z0-9_.\\-]+\\.gtf\\.gz", index))[[1]])
  files <- files[grepl(paste0("\\.", release, "\\.gtf\\.gz$"), files) & !grepl("chr|abinitio|patch|scaff", files)]
  if (length(files) != 1) stop("Expected one main GTF in ", dir_url, ", found: ", paste(files, collapse = ", "))
  dest <- file.path(dest_dir, files)
  utils::download.file(paste0(dir_url, files), dest, mode = "wb", quiet = TRUE)
  return(dest)
}

# Gene table (gene_id, gene_name, gene_biotype) from the 'gene' lines of a GTF
read_gtf_genes <- function(gtf_path) {
  lines <- data.table::fread(cmd = paste("zcat", shQuote(gtf_path), "| awk -F'\\t' '$3==\"gene\"'"), sep = "\t", header = FALSE, quote = "", select = 9)[[1]]
  attr_value <- function(key) {
    out <- rep(NA_character_, length(lines))
    hit <- grepl(paste0(key, " \""), lines)
    out[hit] <- sub(paste0(".*", key, " \"([^\"]*)\".*"), "\\1", lines[hit])
    return(out)
  }
  genes <- data.frame(ensembl_gene_id = attr_value("gene_id"), external_gene_name = attr_value("gene_name"), gene_biotype = attr_value("gene_biotype"), stringsAsFactors = FALSE)
  return(genes)
}

# Identifier prefixes rewritten exactly as in curate_CellHub_DB.R (step 1-ID), so that coverage is measured on the final tokens;
# the maize v4 -> v5 mapping is applied by the curation step and is not repeated here
normalize_ids <- function(x, species) {
  if (species == "maize") x <- sub("^ZEAMMB73_", "", x)
  if (species == "soybean") x <- sub("^GLYMA_", "Glyma.", x)
  return(x)
}

# Gene tokens of all converted sources for one species (partners split on the complex separator "_")
species_tokens <- function(conv, species) {
  ids <- trimws(c(conv$ligand[conv$species == species], conv$receptor[conv$species == species]))
  ids <- normalize_ids(ids[!is.na(ids) & ids != ""], species)
  tokens <- unique(unlist(strsplit(ids, "_"), use.names = FALSE))
  return(tokens[tokens != ""])
}

# Load the ligand / receptor / species columns of every converted source
conv_files <- list.files(CONVERTED_DIR, pattern = "\\.rds$", full.names = TRUE)
if (length(conv_files) == 0) stop("No .rds files found in ", CONVERTED_DIR, "/. Run the convert_*.R scripts first.")
conv <- do.call(rbind, lapply(conv_files, function(f) {
  df <- readRDS(f)
  return(data.frame(ligand = as.character(df$ligand), receptor = as.character(df$receptor), species = as.character(df$species), stringsAsFactors = FALSE))
}))
log_msg(sprintf("%d converted files, %d entries, species: %s", length(conv_files), nrow(conv), paste(sort(unique(conv$species)), collapse = ", ")))

report <- list(); manifest <- list()
for (i in seq_len(nrow(ENS_MAP))) {
  sp <- ENS_MAP$species[i]
  if (!sp %in% conv$species) { log_msg(sprintf("%s: not in converted sources, skipped", sp)); next }
  release <- if (ENS_MAP$division[i] == "plants") PLANTS_RELEASE else ENSEMBL_RELEASE
  dir_url <- gtf_dir_url(ENS_MAP$ftp_dir[i], ENS_MAP$division[i])
  gtf <- tryCatch(download_species_gtf(dir_url, release, GTF_DIR), error = function(e) { log_msg(sprintf("%s: GTF download failed: %s", sp, conditionMessage(e)[1])); return(NULL) })
  if (is.null(gtf)) stop("GTF missing for ", sp, ": geneInfo would be incomplete, stopping")
  manifest[[sp]] <- data.frame(species = sp, file = basename(gtf), url = paste0(dir_url, basename(gtf)), release = release, md5 = unname(tools::md5sum(gtf)), stringsAsFactors = FALSE)
  ens <- read_gtf_genes(gtf)
  ens$ensembl_gene_id <- normalize_ids(ens$ensembl_gene_id, sp)
  tokens <- species_tokens(conv, sp)
  gi_path <- file.path(GENEINFO_DIR, sprintf("geneinfo_%s.rds", sp))
  existing <- if (file.exists(gi_path)) readRDS(gi_path) else NULL
  # Key: gene symbol when a symbol-keyed table already exists (CellChatDB); otherwise the identifier type that best matches the database tokens
  cov_symbol <- mean(toupper(tokens) %in% toupper(ens$external_gene_name))
  cov_id <- mean(toupper(tokens) %in% toupper(ens$ensembl_gene_id))
  key <- if (is.null(existing) && cov_id > cov_symbol) "ensembl_gene_id" else "external_gene_name"
  sym <- ens[[key]]
  keep <- !is.na(sym) & sym != "" & !duplicated(sym)
  gi_ens <- data.frame(Symbol = sym[keep], Ensembl.Gene.ID = ens$ensembl_gene_id[keep], Biotype = ens$gene_biotype[keep], stringsAsFactors = FALSE)
  if (is.null(existing) || nrow(existing) == 0) {
    gi <- gi_ens; mode <- "ensembl"; n_before <- 0L
  } else {
    # Complete the existing table with the Ensembl symbols it lacks; existing rows and columns are kept unchanged
    existing[] <- lapply(existing, function(x) if (is.factor(x)) as.character(x) else x)
    add <- gi_ens[!gi_ens$Symbol %in% existing$Symbol, , drop = FALSE]
    # rep(NA, nrow()) keeps the assignment valid when add is empty (geneInfo already completed by a previous run)
    for (cl in setdiff(colnames(existing), colnames(add))) add[[cl]] <- rep(NA, nrow(add))
    for (cl in setdiff(colnames(add), colnames(existing))) existing[[cl]] <- rep(NA, nrow(existing))
    gi <- rbind(existing, add[, colnames(existing), drop = FALSE]); mode <- "completed"; n_before <- nrow(existing)
  }
  rownames(gi) <- NULL
  saveRDS(gi, gi_path)
  matched <- toupper(tokens) %in% toupper(gi$Symbol)
  writeLines(sort(tokens[!matched]), file.path(GENEINFO_DIR, sprintf("unmatched_tokens_%s.txt", sp)))
  report[[sp]] <- data.frame(species = sp, gtf = basename(gtf), key = key, mode = mode, n_geneinfo_before = n_before, n_geneinfo_after = nrow(gi),
                             n_tokens = length(tokens), pct_tokens_exact = round(100 * mean(tokens %in% gi$Symbol), 1),
                             pct_tokens_any_case = round(100 * mean(matched), 1), stringsAsFactors = FALSE)
  log_msg(sprintf("%s: %s, key %s, %d -> %d genes, %.1f %% of %d tokens matched (exact case)", sp, mode, key, n_before, nrow(gi),
                  report[[sp]]$pct_tokens_exact, length(tokens)))
}

report_df <- do.call(rbind, report)
write.csv(report_df, file.path(GENEINFO_DIR, "geneinfo_coverage.csv"), row.names = FALSE)
write.table(do.call(rbind, manifest), file.path(GENEINFO_DIR, "gtf_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(report_df, row.names = FALSE)
# Reviewer sentinel symbols must now be in the human reference
if ("human" %in% names(report)) {
  hs <- readRDS(file.path(GENEINFO_DIR, "geneinfo_human.rds"))$Symbol
  sentinel <- c("DEFB4A", "GPR15LG", "S100A8", "S100A9", "P2RY12")
  log_msg(sprintf("human sentinel symbols in reference: %s", paste(sprintf("%s=%s", sentinel, sentinel %in% hs), collapse = " ")))
}
log_msg("Done")
