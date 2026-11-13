# Step 00 (Reviewer 2, major 1): manifest of raw source files (size, date, MD5) and checksum-verified re-download
raw_dir <- "data/raw"; manifest_file <- "config/sources_manifest.tsv"
if (!dir.exists(raw_dir)) stop("Raw source directory not found: ", normalizePath(raw_dir, mustWork = FALSE), " (set raw_dir to the folder holding the downloaded source files)")
if (length(list.files(raw_dir, recursive = TRUE)) == 0) stop("Raw source directory is empty: ", raw_dir)
dir.create(dirname(manifest_file), recursive = TRUE, showWarnings = FALSE)

# Hash every raw file; source / version / url columns are filled by hand once, then committed
write_source_manifest <- function(raw_dir, out_file) {
  files <- list.files(raw_dir, recursive = TRUE, full.names = TRUE)
  info <- file.info(files)
  man <- data.frame(file = substring(normalizePath(files), nchar(normalizePath(raw_dir)) + 2), size_bytes = info$size,
                    modified = format(info$mtime, "%Y-%m-%d"), md5 = unname(tools::md5sum(files)), source = "", version = "", url = "", stringsAsFactors = FALSE)
  write.table(man, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  return(man)
}

# Download every file that has a URL and check it against the recorded MD5
fetch_sources <- function(manifest_file, raw_dir) {
  man <- read.delim(manifest_file, stringsAsFactors = FALSE, colClasses = "character")
  options(timeout = max(1200, getOption("timeout")))
  man$status <- vapply(seq_len(nrow(man)), function(i) {
    if (is.na(man$url[i]) || man$url[i] == "") return("no_url (extracted from R package)")
    dest <- file.path(raw_dir, man$file[i]); dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    ok <- tryCatch({ utils::download.file(man$url[i], dest, mode = "wb", quiet = TRUE); TRUE }, error = function(e) FALSE)
    if (!ok) return("download_failed")
    if (unname(tools::md5sum(dest)) != man$md5[i]) return("md5_mismatch")
    return("ok")
  }, character(1))
  return(man)
}

# Pinned R packages from which CellChatDB, NeuronChatDB and NicheNet are extracted (refs as stated in the manuscript;
# NeuronChat repository name and nichenetr / multinichenetr commits to be checked with packageDescription(...)$RemoteSha)
install_pinned_sources <- function(nichenetr_ref, multinichenetr_ref) {
  remotes::install_github("jinworks/CellChat", ref = "75253cd", upgrade = "never")
  remotes::install_github("Wei-BioMath/NeuronChat", ref = "36e06b2", upgrade = "never")
  remotes::install_github("saeyslab/nichenetr", ref = nichenetr_ref, upgrade = "never")
  remotes::install_github("saeyslab/multinichenetr", ref = multinichenetr_ref, upgrade = "never")
  return(invisible(TRUE))
}

if (!file.exists(manifest_file)) {
  man <- write_source_manifest(raw_dir, manifest_file)
  message(nrow(man), " raw files hashed -> fill source/version/url in ", manifest_file)
} else {
  check <- fetch_sources(manifest_file, file.path(tempdir(), "refetch"))
  print(table(check$status))
  write.table(check, sub("\\.tsv$", "_refetch_check.tsv", manifest_file), sep = "\t", quote = FALSE, row.names = FALSE)
}
