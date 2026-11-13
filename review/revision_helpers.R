# Pure R helpers for the CellHubDB revision (NARGAB-2026-310); sourced by the revision scripts

# CellChat database schema, taken from the CellChat mouse release used as reference
get_cellchat_schema <- function() {
  ref <- CellChat::CellChatDB.mouse
  core <- c("interaction_name", "pathway_name", "ligand", "receptor", "agonist", "antagonist", "co_A_receptor", "co_I_receptor", "annotation", "interaction_name_2")
  return(list(core_interaction = core, complex = colnames(ref$complex), cofactor = colnames(ref$cofactor),
              cofactor_refs = c("agonist", "antagonist", "co_A_receptor", "co_I_receptor")))
}

# Force a species database into the CellChat schema without adding or removing any interaction:
# factors -> character, NA cofactor references -> "", missing columns added empty, empty tables keep every expected column
standardize_cellchat_schema <- function(sp_db, schema = get_cellchat_schema()) {
  fix_table <- function(tab, cols) {
    if (is.null(tab) || ncol(tab) == 0) tab <- data.frame(row.names = character(0))
    tab[] <- lapply(tab, as.character)
    for (cl in setdiff(cols, colnames(tab))) tab[[cl]] <- rep("", nrow(tab))
    tab <- tab[, c(cols, setdiff(colnames(tab), cols)), drop = FALSE]
    tab[is.na(tab)] <- ""
    return(tab)
  }
  inter <- sp_db$interaction
  inter[] <- lapply(inter, function(x) if (is.factor(x)) as.character(x) else x)
  for (cl in schema$cofactor_refs) {
    if (!cl %in% colnames(inter)) inter[[cl]] <- rep("", nrow(inter))
    inter[[cl]][is.na(inter[[cl]])] <- ""
  }
  if (!"interaction_name_2" %in% colnames(inter)) inter$interaction_name_2 <- paste(inter$ligand, "-", inter$receptor)
  complex_cols <- union(schema$complex, grep("^subunit_", colnames(sp_db$complex), value = TRUE))
  cofactor_cols <- union(schema$cofactor, grep("^cofactor", colnames(sp_db$cofactor), value = TRUE))
  gene_info <- sp_db$geneInfo
  if (is.null(gene_info) || !"Symbol" %in% colnames(gene_info)) gene_info <- data.frame(Symbol = character(0), stringsAsFactors = FALSE)
  gene_info[] <- lapply(gene_info, function(x) if (is.factor(x)) as.character(x) else x)
  return(list(interaction = inter, complex = fix_table(sp_db$complex, complex_cols), cofactor = fix_table(sp_db$cofactor, cofactor_cols), geneInfo = gene_info))
}

# Genes behind one interaction partner: complex subunits when the partner is a defined complex, the partner itself otherwise
partner_subunits <- function(partner, complex_df) {
  if (!is.null(complex_df) && partner %in% rownames(complex_df)) {
    g <- unlist(complex_df[partner, grep("^subunit_", colnames(complex_df)), drop = FALSE], use.names = FALSE)
    return(g[!is.na(g) & g != ""])
  }
  return(partner)
}

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

# Number of curated sources using each gene in a given role (ligand or receptor), complexes expanded to subunits
role_source_counts <- function(partners, src_list, cpx, curated) {
  genes <- lapply(partners, partner_subunits, complex_df = cpx)
  idx <- rep(seq_along(partners), lengths(genes)); g <- unlist(genes, use.names = FALSE)
  by_gene <- tapply(idx, g, function(ix) length(intersect(unique(unlist(src_list[ix])), curated)))
  return(by_gene)
}

# Pair-level evidence rules, as published by each source: TRUE = curated (literature or experimental) support,
# FALSE = orthology-inferred or computationally predicted support. Sources not listed are curated for every pair.
# connectomeDB2025: Direct vs Inferred (orthology); FlyPhoneDB2: high = literature-curated, moderate/low = AlphaFold-Multimer / ortholog predictions;
# cell2cell (Armingol 2022): "Literature supported" vs "Putative"; PlantPhoneDB: literature and experimental PPI resources vs orthologs / STRING;
# PlantCellChatDB: "experimentally validated" vs "STRING-predicted"; NicheNet v2 (MultiNicheNet): computational aggregation.
# CellChatDB_v2 is curated for human and mouse; its zebrafish database is an orthology projection of the mouse database
# (99.9 % of pairs with a mouse equivalent carry the mouse evidence, remaining pairs cite mouse KEGG maps such as mmu04310)
CELLCHATDB_PROJECTED_SPECIES <- c("zebrafish")
SUPPORT_RULES <- list(
  connectomeDB2025 = function(it) it$evidence_type %in% "Direct",
  CellChatDB_v2 = function(it) {
    if (!"species" %in% colnames(it)) stop("CellChatDB_v2 support rule needs the species column")
    return(!it$species %in% CELLCHATDB_PROJECTED_SPECIES)
  },
  FlyPhoneDB2 = function(it) it$flyphone_confidence %in% "high",
  cell2cell_Armingol = function(it) it$pair_confidence %in% "Literature supported",
  PlantPhoneDB = function(it) it$plantphonedb_source %in% c("literature", "IntAct", "BIOGRID", "interactome2.0", "plant.MAP"),
  PlantCellChatDB = function(it) it$annotation %in% "experimentally validated",
  MultiNicheNet = function(it) rep(FALSE, nrow(it)))
# Sources that do not encode the ligand -> receptor direction: their pairs are matched in both orientations
NON_ORIENTED_SOURCES <- c("CellPhoneDB", "PlantCellChatDB")

# Curated support of every pair of one per-source database, as a named logical vector keyed by "ligand|receptor" (lower case);
# a pair listed several times is curated if any of its rows is
source_support_status <- function(src_db, source, rules = SUPPORT_RULES) {
  it <- src_db$interaction
  status <- if (source %in% names(rules)) rules[[source]](it) else rep(TRUE, nrow(it))
  keys <- tolower(paste(it$ligand, it$receptor, sep = "|"))
  if (source %in% NON_ORIENTED_SOURCES) { keys <- c(keys, tolower(paste(it$receptor, it$ligand, sep = "|"))); status <- c(status, status) }
  out <- tapply(status, keys, any)
  return(setNames(as.logical(out), names(out)))
}

# Confidence annotation, two independent dimensions:
# (1) confidence_tier from the evidence of each supporting source for this pair: high = >= 2 sources with curated support,
#     medium = 1, low = only orthology-inferred or predicted support. With support = NULL the evidence is taken at the source level
#     (all sources curated except predicted_sources). curated_sources lists the sources giving curated support.
# (2) role_corroborated flag: the ligand is used as a ligand, and the receptor as a receptor, with curated support from >= 2 sources
#     across all pairs of the species (gene-level check, NA when the species has fewer than two such sources). It removes non-canonical
#     ligands supported by a single resource but does not validate the pair itself.
# support: named list (one element per source) of named logical vectors from source_support_status(); a pair missing from its
# source file falls back to the source default (curated unless the source is in predicted_sources) and is counted in attr "unmatched"
annotate_confidence <- function(sp_db, predicted_sources, support = NULL) {
  inter <- sp_db$interaction; cpx <- sp_db$complex
  src <- strsplit(ifelse(is.na(inter$sources), "", inter$sources), ",\\s*")
  src <- lapply(src, function(s) s[s != ""])
  keys <- tolower(paste(inter$ligand, inter$receptor, sep = "|"))
  unmatched <- setNames(integer(0), character(0))
  # Sources with curated support, pair by pair
  src_cur <- lapply(seq_along(src), function(i) {
    s <- src[[i]]
    if (length(s) == 0) return(character(0))
    ok <- vapply(s, function(x) {
      if (!is.null(support) && x %in% names(support)) {
        v <- support[[x]][keys[i]]
        if (!is.na(v)) return(unname(v))
        unmatched[x] <<- if (is.na(unmatched[x])) 1L else unmatched[x] + 1L
      }
      return(!x %in% predicted_sources)
    }, logical(1))
    return(s[ok])
  })
  inter$n_sources <- lengths(src)
  inter$n_curated_sources <- lengths(src_cur)
  inter$curated_sources <- vapply(src_cur, paste, character(1), collapse = ", ")
  inter$predicted_only <- inter$n_sources > 0 & inter$n_curated_sources == 0
  inter$is_complex <- inter$ligand %in% rownames(cpx) | inter$receptor %in% rownames(cpx) | grepl("_", inter$ligand) | grepl("_", inter$receptor)
  curated_any <- unique(unlist(src_cur))
  if (length(curated_any) >= 2) {
    lig_n <- role_source_counts(inter$ligand, src_cur, cpx, curated_any); rec_n <- role_source_counts(inter$receptor, src_cur, cpx, curated_any)
    best <- function(p, counts) max(c(0, counts[intersect(partner_subunits(p, cpx), names(counts))]), na.rm = TRUE)
    inter$ligand_role_sources <- vapply(inter$ligand, best, numeric(1), counts = lig_n, USE.NAMES = FALSE)
    inter$receptor_role_sources <- vapply(inter$receptor, best, numeric(1), counts = rec_n, USE.NAMES = FALSE)
    inter$role_corroborated <- inter$ligand_role_sources >= 2 & inter$receptor_role_sources >= 2
  } else {
    inter$ligand_role_sources <- NA_real_; inter$receptor_role_sources <- NA_real_; inter$role_corroborated <- NA
  }
  inter$confidence_tier <- ifelse(inter$n_curated_sources >= 2, "high", ifelse(inter$predicted_only, "low", "medium"))
  sp_db$interaction <- inter
  attr(sp_db, "unmatched") <- unmatched
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
