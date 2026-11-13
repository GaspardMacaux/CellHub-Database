################################################################################
#   collect_source_provenance.R — inventory of the raw source files used to    #
#   build CellHub Database (path, size, modification date, MD5) and of the R   #
#   packages that ship some of them. Run from the CellHub_DB root.              #
#                                                                              #
# Usage:  Rscript Script/collect_source_provenance.R [raw_dir] [output_csv]     #
# Output: source_provenance.csv + package_versions.csv                         #
################################################################################

`%||%` <- function(a, b) if (is.null(a)) b else a

# Raw file folder per source, relative to raw_dir
SOURCE_FOLDERS <- c(
  CellChatDB_v2      = "CellChat",
  CellPhoneDB        = "cellphonedb-data-5.0.0",
  CellTalkDB         = "CellTalkDB",
  connectomeDB2025   = "ConnectomeDB2025",
  NicheNet_v2        = "Multinichenet",
  NeuronChatDB       = "NeuronChat",
  FlyPhoneDB2        = "Flyphone",
  cell2cell_Celegans = "Celegans_LR",
  PlantPhoneDB       = "Plants-2022-Xu-LR-pairs",
  PlantCellChatDB    = "PlantCellChatDB"
)

# One row per raw file: source, relative path, size, last modification time, MD5
inventory_source_files <- function(raw_dir) {
  raw_root <- normalizePath(raw_dir)
  rows <- lapply(names(SOURCE_FOLDERS), function(src) {
    folder <- file.path(raw_root, SOURCE_FOLDERS[[src]])
    files <- if (dir.exists(folder)) list.files(folder, full.names = TRUE, recursive = TRUE) else character(0)
    if (length(files) == 0) {
      return(data.frame(source = src, file = NA_character_, size_bytes = NA_real_, modified = NA_character_,
                        md5 = NA_character_, status = "no file found", stringsAsFactors = FALSE))
    }
    info <- file.info(files)
    data.frame(
      source     = src,
      file       = sub(paste0("^", raw_root, "/?"), "", normalizePath(files)),
      size_bytes = info$size,
      modified   = format(info$mtime, "%Y-%m-%d %H:%M"),
      md5        = unname(tools::md5sum(files)),
      status     = "ok",
      stringsAsFactors = FALSE
    )
  })
  return(do.call(rbind, rows))
}

# Versions of the R packages that distribute or describe some of the sources
collect_package_versions <- function() {
  packages <- c("CellChat", "NeuronChat", "nichenetr", "multinichenetr", "Seurat")
  versions <- vapply(packages, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "not installed"
  }, character(1))
  remote <- vapply(packages, function(p) {
    desc <- tryCatch(utils::packageDescription(p), error = function(e) NULL)
    if (is.null(desc) || is.null(desc$RemoteSha)) "" else paste0(desc$RemoteRepo %||% "", "@", substr(desc$RemoteSha, 1, 7))
  }, character(1))
  return(data.frame(package = packages, version = versions, github = remote, stringsAsFactors = FALSE))
}

args <- commandArgs(trailingOnly = TRUE)
raw_dir <- if (length(args) >= 1) args[1] else "data/raw"
out_csv <- if (length(args) >= 2) args[2] else "output/source_provenance.csv"
if (!dir.exists(raw_dir)) stop("Raw source folder not found: ", raw_dir)

provenance <- inventory_source_files(raw_dir)
write.csv(provenance, out_csv, row.names = FALSE)
message(sprintf("%d raw files inventoried -> %s", sum(provenance$status == "ok"), out_csv))
missing <- provenance$source[provenance$status != "ok"]
if (length(missing) > 0) message("No raw file found for: ", paste(missing, collapse = ", "))

pkg <- collect_package_versions()
write.csv(pkg, "output/package_versions.csv", row.names = FALSE)
print(pkg, row.names = FALSE)