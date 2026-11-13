# Pure R helpers for the CellHubDB revision (NARGAB-2026-310); sourced by the revision scripts and the release gate.
# Schema, support rules and confidence annotation come from R/cellhubdb_functions.R (single definition shared with the pipeline).
source("R/cellhubdb_functions.R")

# All gene tokens used by the interaction partners of a species database (complexes expanded to subunits)
extract_partner_genes <- function(sp_db) {
  partners <- unique(c(sp_db$interaction$ligand, sp_db$interaction$receptor))
  partners <- partners[!is.na(partners) & partners != ""]
  genes <- unique(unlist(lapply(partners, partner_subunits, complex_df = sp_db$complex), use.names = FALSE))
  return(genes)
}

# One-row structural audit of a species database
audit_species_db <- function(sp_db, species, schema = get_cellchat_schema()) {
  inter <- sp_db$interaction; cpx <- sp_db$complex; cof <- sp_db$cofactor; gi <- sp_db$geneInfo
  nr <- function(x) if (is.null(x)) 0L else nrow(x)
  nc <- function(x) if (is.null(x)) 0L else ncol(x)
  partners <- unique(c(inter$ligand, inter$receptor)); partners <- partners[!is.na(partners) & partners != ""]
  multi <- partners[grepl("_", partners)]
  cof_ids <- unlist(inter[, intersect(schema$cofactor_refs, colnames(inter)), drop = FALSE], use.names = FALSE)
  cof_ids <- unique(cof_ids[!is.na(cof_ids) & cof_ids != ""])
  return(data.frame(species = species, n_interaction = nr(inter), n_complex = nr(cpx), n_cofactor = nr(cof), n_geneinfo = nr(gi),
                    complex_ncol = nc(cpx), cofactor_ncol = nc(cof), geneinfo_has_symbol = !is.null(gi) && "Symbol" %in% colnames(gi),
                    missing_core_cols = paste(setdiff(schema$core_interaction, colnames(inter)), collapse = ";"),
                    n_partners = length(partners), n_multisubunit = length(multi), n_multisubunit_undefined = sum(!multi %in% rownames(cpx)),
                    n_cofactor_refs = length(cof_ids), n_cofactor_refs_undefined = sum(!cof_ids %in% rownames(cof)),
                    pct_pathway_na = round(100 * mean(is.na(inter$pathway_name) | inter$pathway_name == ""), 1), stringsAsFactors = FALSE))
}

# Synthetic log-normalized dataset containing every database gene, with group-specific means so that over-expression is detected
make_synthetic_data <- function(genes, n_groups = 3, cells_per_group = 60, n_background = 300, seed = 42) {
  set.seed(seed)
  all_genes <- c(genes, paste0("bgGene", seq_len(n_background)))
  groups <- rep(paste0("G", seq_len(n_groups)), each = cells_per_group)
  lambda <- matrix(sample(c(0.05, 0.5, 3), length(all_genes) * n_groups, replace = TRUE), nrow = length(all_genes))
  counts <- do.call(cbind, lapply(seq_along(groups), function(i) rpois(length(all_genes), lambda[, match(groups[i], unique(groups))])))
  rownames(counts) <- all_genes; colnames(counts) <- paste0("cell", seq_along(groups))
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  lib <- Matrix::colSums(counts); lib[lib == 0] <- 1
  data_norm <- as(log1p(Matrix::t(Matrix::t(counts) / lib) * 1e4), "CsparseMatrix")
  meta <- data.frame(labels = factor(groups), row.names = colnames(counts))
  return(list(data = data_norm, meta = meta))
}

# Run the CellChat steps used by Cell-Hub on synthetic data; records the first failing step and its error
run_cellchat_smoke <- function(sp_db, species, seed = 42) {
  steps <- c("createCellChat", "subsetData", "identifyOverExpressedGenes", "identifyOverExpressedInteractions",
             "computeCommunProb", "filterCommunication", "computeCommunProbPathway")
  status <- setNames(rep("not_run", length(steps)), steps); first_error <- NA_character_; n_sig <- NA_integer_
  syn <- make_synthetic_data(extract_partner_genes(sp_db), seed = seed)
  cc <- NULL
  for (st in steps) {
    res <- tryCatch({
      cc <- switch(st,
                   createCellChat = { x <- CellChat::createCellChat(object = syn$data, meta = syn$meta, group.by = "labels"); x@DB <- sp_db; x },
                   subsetData = CellChat::subsetData(cc),
                   identifyOverExpressedGenes = CellChat::identifyOverExpressedGenes(cc),
                   identifyOverExpressedInteractions = CellChat::identifyOverExpressedInteractions(cc),
                   computeCommunProb = CellChat::computeCommunProb(cc, type = "triMean", nboot = 20, seed.use = seed),
                   filterCommunication = CellChat::filterCommunication(cc, min.cells = 10),
                   computeCommunProbPathway = CellChat::computeCommunProbPathway(cc))
      "ok"
    }, error = function(e) conditionMessage(e)[1])
    if (res != "ok") { status[st] <- "error"; first_error <- paste0(st, ": ", res); break }
    status[st] <- "ok"
    if (st == "filterCommunication") n_sig <- tryCatch(nrow(CellChat::subsetCommunication(cc, thresh = 0.05)), error = function(e) NA_integer_)
  }
  return(data.frame(species = species, t(status), n_significant = n_sig, first_error = first_error, stringsAsFactors = FALSE, check.names = FALSE))
}

# Keep a subset of interactions (complex and cofactor tables are kept whole, unused rows are harmless for CellChat)
subset_cellchat_db <- function(sp_db, keep) {
  keep[is.na(keep)] <- FALSE
  sp_db$interaction <- sp_db$interaction[keep, , drop = FALSE]
  return(sp_db)
}

# Run CellChat with the Cell-Hub parameters and return the significant communications
run_cellchat_db <- function(data_norm, meta, db, seed = 1, nboot = 100) {
  cc <- CellChat::createCellChat(object = data_norm, meta = meta, group.by = "labels")
  cc@DB <- db
  cc <- CellChat::subsetData(cc)
  cc <- CellChat::identifyOverExpressedGenes(cc)
  cc <- CellChat::identifyOverExpressedInteractions(cc)
  cc <- CellChat::computeCommunProb(cc, type = "triMean", nboot = nboot, seed.use = seed)
  cc <- CellChat::filterCommunication(cc, min.cells = 10)
  res <- tryCatch(CellChat::subsetCommunication(cc, thresh = 0.05), error = function(e) NULL)
  return(res)
}

# Counts of significant calls, distinct LR pairs and pathways in a CellChat result
summarise_calls <- function(res) {
  if (is.null(res) || nrow(res) == 0) return(c(n_calls = 0, n_lr = 0, n_pathways = 0))
  return(c(n_calls = nrow(res), n_lr = length(unique(tolower(paste(res$ligand, res$receptor)))), n_pathways = length(unique(res$pathway_name))))
}

# Download the released CellHubDB file once, verify its checksum, return all species or one
get_cellhub_db <- function(url, md5, species = NULL, file_name = "CellHubDB.rds", cache_dir = tools::R_user_dir("CellHubDB", which = "cache")) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(cache_dir, file_name)
  if (!file.exists(dest) || unname(tools::md5sum(dest)) != md5) {
    options(timeout = max(600, getOption("timeout")))
    utils::download.file(url, dest, mode = "wb")
    if (unname(tools::md5sum(dest)) != md5) stop("Checksum mismatch for ", dest)
  }
  db <- readRDS(dest)
  if (is.null(species)) return(db)
  if (!species %in% names(db)) stop("Unknown species: ", species, ". Available: ", paste(names(db), collapse = ", "))
  return(db[[species]])
}
