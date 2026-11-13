################################################################################
#   00_download_sources.R — re-obtain every raw source file from a pinned     #
#   public origin and check it against the recorded MD5 or content.           #
#                                                                              #
# Table: config/source_retrieval.tsv (file, md5, method, location, member, note)
# Methods:                                                                     #
#   url         download a pinned URL (GitHub commit, tag or Zenodo record)    #
#   zip_member  download a pinned zip and take one member                      #
#   package     extract datasets from an installed R package at a pinned commit#
#   rds_to_tsv  download a pinned .rds table and write it as tab-separated text#
#   derived     build the file from another raw file (function named in member)#
#   archived    no versioned public download: the file comes from the CellHubDB#
#               raw-data deposit and is only checked against its MD5          #
#   unused      kept in data/raw/ but not read by the pipeline (not checked)   #
# Check: MD5 identical, or (rds / tables) identical content when only the      #
# serialisation differs. Report: output/reports/source_retrieval_check.csv    #
# Usage: Rscript pipeline/00_download_sources.R            verify only        #
#        Rscript pipeline/00_download_sources.R --install  also write missing  #
#                                                          files to data/raw/  #
################################################################################

RAW_DIR <- "data/raw"
TABLE <- "config/source_retrieval.tsv"
REPORT <- "output/reports/source_retrieval_check.csv"
install <- "--install" %in% commandArgs(trailingOnly = TRUE)
stage <- file.path(tempdir(), "cellhubdb_sources")
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(REPORT), recursive = TRUE, showWarnings = FALSE)
options(timeout = max(1200, getOption("timeout")))
sources <- read.delim(TABLE, stringsAsFactors = FALSE, colClasses = "character", na.strings = "")

# Download one URL into the staging folder (cached per URL, so a zip is fetched once for all its members)
download_once <- function(url) {
  dest <- file.path(stage, gsub("[^A-Za-z0-9._-]", "_", sub("^https?://", "", url)))
  if (!file.exists(dest)) tryCatch(utils::download.file(url, dest, mode = "wb", quiet = TRUE), error = function(e) { unlink(dest); stop("download failed: ", url, " (", conditionMessage(e), ")") })
  return(dest)
}

# Installed package must match the pinned commit (location = "Package@commitprefix")
check_package_commit <- function(location) {
  pkg <- sub("@.*$", "", location); pin <- sub("^.*@", "", location)
  if (!requireNamespace(pkg, quietly = TRUE)) stop(pkg, " not installed (pinned commit ", pin, ")")
  sha <- packageDescription(pkg)$RemoteSha
  if (is.null(sha) || !startsWith(sha, pin)) stop(pkg, " installed at commit ", if (is.null(sha)) "unknown" else substr(sha, 1, 10), ", pinned ", pin)
  return(pkg)
}

# Load datasets of a package; member = "obj" (single object) or "name=obj;name=obj" (named list)
package_object <- function(pkg, member) {
  specs <- strsplit(strsplit(member, ";", fixed = TRUE)[[1]], "=", fixed = TRUE)
  load_one <- function(obj) { e <- new.env(); utils::data(list = obj, package = pkg, envir = e); if (!exists(obj, envir = e)) stop(obj, " not found in package ", pkg); return(e[[obj]]) }
  if (length(specs) == 1 && length(specs[[1]]) == 1) return(load_one(specs[[1]][1]))
  return(setNames(lapply(specs, function(s) load_one(s[2])), vapply(specs, `[`, character(1), 1)))
}

# One row per v4 / v5 pair of the MaizeGDB cross-reference (second column holds comma-separated v5 models)
explode_maize_xref <- function(input, output) {
  x <- read.delim(input, header = FALSE, stringsAsFactors = FALSE, colClasses = "character")
  v5 <- strsplit(x[[2]], ",", fixed = TRUE)
  long <- data.frame(v4 = rep(x[[1]], lengths(v5)), v5 = trimws(unlist(v5)), stringsAsFactors = FALSE)
  long <- long[!is.na(long$v5) & long$v5 != "", ]
  write.table(long, output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  return(invisible(output))
}

# Produce one source file in the staging folder; returns its path
retrieve <- function(s) {
  out <- file.path(stage, "files", s$file); dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  if (s$method == "url") file.copy(download_once(s$location), out, overwrite = TRUE)
  else if (s$method == "zip_member") {
    ex <- file.path(stage, "unzip", basename(s$location)); utils::unzip(download_once(s$location), files = s$member, exdir = ex, junkpaths = TRUE)
    file.copy(file.path(ex, basename(s$member)), out, overwrite = TRUE)
  } else if (s$method == "package") saveRDS(package_object(check_package_commit(s$location), s$member), out)
  else if (s$method == "rds_to_tsv") write.table(as.data.frame(readRDS(download_once(s$location))), out, sep = "\t", quote = FALSE, row.names = FALSE)
  else if (s$method == "derived") {
    input <- file.path(RAW_DIR, s$location)
    if (!file.exists(input)) input <- file.path(stage, "files", s$location)
    if (!file.exists(input)) stop("input of the derived file not available: ", s$location)
    match.fun(s$member)(input, out)
  } else stop("unknown method: ", s$method)
  return(out)
}

# Content of a file for comparison when only the serialisation differs (rds objects, delimited tables)
read_content <- function(path) {
  if (grepl("\\.rds$", path)) return(readRDS(path))
  if (grepl("\\.(txt|tsv)$", path)) return(read.delim(path, colClasses = "character", quote = "", fileEncoding = "latin1"))
  if (grepl("\\.csv$", path)) return(read.csv(path, colClasses = "character"))
  return(NULL)
}

# Compare a retrieved file with the expected MD5, then with the local copy's content
compare <- function(fetched, s) {
  md5 <- unname(tools::md5sum(fetched))
  if (identical(md5, s$md5)) return("ok: md5 identical")
  local <- file.path(RAW_DIR, s$file)
  if (file.exists(local)) {
    a <- read_content(fetched); b <- read_content(local)
    if (!is.null(a) && isTRUE(all.equal(a, b, check.attributes = FALSE))) return("ok: content identical (md5 differs, serialisation only)")
  }
  return(paste("MISMATCH: md5", md5, "expected", s$md5))
}

res <- vector("list", nrow(sources))
# Derived files after their inputs
order_rows <- order(sources$method == "derived")
for (i in order_rows) {
  s <- as.list(sources[i, ]); s[is.na(s)] <- ""
  local <- file.path(RAW_DIR, s$file)
  status <- tryCatch({
    if (s$method == "unused") "skipped: not used by the pipeline"
    else if (s$method == "archived") {
      if (!file.exists(local)) "MISSING: archived file absent from data/raw (take it from the CellHubDB raw-data deposit)"
      else if (identical(unname(tools::md5sum(local)), s$md5)) "ok: archived file, md5 identical"
      else paste("MISMATCH: archived file md5", unname(tools::md5sum(local)))
    } else {
      fetched <- retrieve(s); st <- compare(fetched, s)
      if (install && startsWith(st, "ok") && !file.exists(local)) { dir.create(dirname(local), recursive = TRUE, showWarnings = FALSE); file.copy(fetched, local); st <- paste(st, "| installed") }
      st
    }
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  res[[i]] <- data.frame(file = s$file, method = s$method, status = status, location = s$location, stringsAsFactors = FALSE)
  cat(sprintf("%-60s %-11s %s\n", s$file, s$method, status))
}
report <- do.call(rbind, res)
write.csv(report, REPORT, row.names = FALSE)
cat("\n"); print(table(sub(":.*$", "", report$status)))
bad <- !grepl("^(ok|skipped)", report$status)
if (any(bad)) stop(sum(bad), " source file(s) not reproduced, see ", REPORT)
cat("All source files reproduced or verified. Report:", REPORT, "\n")
