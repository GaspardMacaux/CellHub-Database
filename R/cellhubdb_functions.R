# Shared functions of the CellHubDB pipeline (finalize step, release gate): CellChat schema, pair-level confidence,
# release repairs, pathway propagation and interaction types. Pure R, no side effects; sourced by finalize_cellhubdb.R.

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

cofactor_ref_cols <- c("agonist", "antagonist", "co_A_receptor", "co_I_receptor")

# One gene vector per row of a complex or cofactor table (factors converted, empty cells removed)
table_to_sets <- function(tab) {
  if (is.null(tab) || nrow(tab) == 0 || ncol(tab) == 0) return(list())
  tab[] <- lapply(tab, as.character)
  sets <- lapply(seq_len(nrow(tab)), function(i) { v <- unname(unlist(tab[i, , drop = TRUE])); v[!is.na(v) & v != ""] })
  return(setNames(sets, rownames(tab)))
}

# Table with columns <prefix>1..n (n >= min_cols) built from named gene vectors, padded with ""
sets_to_table <- function(sets, prefix, min_cols) {
  n <- max(min_cols, if (length(sets) > 0) max(lengths(sets)) else 0, 1)
  m <- matrix("", nrow = length(sets), ncol = n, dimnames = list(names(sets), paste0(prefix, seq_len(n))))
  for (i in seq_along(sets)) if (length(sets[[i]]) > 0) m[i, seq_along(sets[[i]])] <- sets[[i]]
  return(as.data.frame(m, stringsAsFactors = FALSE))
}

# Escape a gene symbol for use inside a regular expression
escape_regex <- function(x) gsub("([][.|()^{}+$*?\\\\])", "\\\\\\1", x, perl = TRUE)

# CellChat database of the same species, used as the reference for definitions missing from CellHubDB (NULL when CellChat has none)
cellchat_reference <- function(species) {
  ref_name <- switch(species, human = "CellChatDB.human", mouse = "CellChatDB.mouse", zebrafish = "CellChatDB.zebrafish", NULL)
  if (is.null(ref_name)) return(NULL)
  return(getExportedValue("CellChat", ref_name))
}

# Complexes whose definition contains a purely numeric subunit: redefined from the reference (case-insensitive name, then name without the _<n> suffix)
restore_complex_definitions <- function(db, ref_complex) {
  cpx <- db$complex; cpx[] <- lapply(cpx, as.character)
  broken <- rownames(cpx)[apply(cpx, 1, function(r) any(grepl("^[0-9]+$", r)))]
  log <- data.frame(complex = character(0), action = character(0), subunits = character(0), stringsAsFactors = FALSE)
  if (length(broken) == 0) return(list(db = db, unrestored = character(0), log = log))
  ref_sets <- table_to_sets(ref_complex)
  ref_keys <- toupper(names(ref_sets))
  unrestored <- character(0)
  for (name in broken) {
    k <- match(toupper(name), ref_keys)
    if (is.na(k)) k <- match(toupper(sub("_[0-9]+$", "", name)), ref_keys)
    if (is.na(k)) { unrestored <- c(unrestored, name); next }
    genes <- ref_sets[[k]]
    need <- paste0("subunit_", seq_along(genes))
    for (cl in setdiff(need, colnames(cpx))) cpx[[cl]] <- rep("", nrow(cpx))
    cpx[name, grepl("^subunit_", colnames(cpx))] <- ""
    cpx[name, need] <- genes
    log[nrow(log) + 1, ] <- c(name, paste("restored from reference complex", names(ref_sets)[k]), paste(genes, collapse = "+"))
  }
  db$complex <- cpx
  return(list(db = db, unrestored = unrestored, log = log))
}

# Cofactors: undefined references defined (gene identifier -> single-gene cofactor; numbered source rows -> union) or removed,
# cofactor genes absent from geneInfo replaced by the species symbol (zebrafish a/b paralogs, otherwise case-insensitive match)
repair_cofactors <- function(db, ref_cofactor, species) {
  cof <- db$cofactor
  cof_cols <- if (!is.null(cof)) grep("^cofactor[0-9]+$", colnames(cof), value = TRUE) else character(0)
  sets <- if (length(cof_cols) > 0) table_to_sets(cof[, cof_cols, drop = FALSE]) else list()
  ref_sets <- table_to_sets(ref_cofactor)
  gi <- db$geneInfo$Symbol
  log <- data.frame(item = character(0), action = character(0), stringsAsFactors = FALSE)
  cols <- intersect(cofactor_ref_cols, colnames(db$interaction))
  ids <- unique(unlist(db$interaction[, cols, drop = FALSE])); ids <- ids[!is.na(ids) & ids != ""]
  for (id in setdiff(ids, names(sets))) {
    pool <- c(ref_sets, sets)
    numbered <- pool[grepl("_[0-9]+$", names(pool)) & sub("_[0-9]+$", "", names(pool)) == id]
    numbered <- numbered[!duplicated(names(numbered))]
    if (id %in% gi) {
      sets[[id]] <- id
      log[nrow(log) + 1, ] <- c(id, "cofactor defined as the single gene of the same name")
    } else if (length(numbered) > 0) {
      sets[[id]] <- unique(unlist(numbered))
      log[nrow(log) + 1, ] <- c(id, paste("cofactor defined as the union of", paste(names(numbered), collapse = " + ")))
    } else {
      for (cl in cols) db$interaction[[cl]][db$interaction[[cl]] %in% id] <- ""
      log[nrow(log) + 1, ] <- c(id, "reference removed from interactions (no definition in any source)")
    }
  }
  # Gene symbols of another species' case: replaced by the symbols of this species
  fix_gene <- function(g) {
    if (g %in% gi) return(g)
    if (species == "zebrafish") {
      paralogs <- sort(gi[grepl(paste0("^", escape_regex(tolower(g)), "[ab]$"), gi)])
      if (length(paralogs) > 0) return(paralogs)
    }
    same <- gi[tolower(gi) == tolower(g)]
    if (length(same) > 0) return(same[1])
    return(g)
  }
  for (id in names(sets)) {
    old <- sets[[id]]
    new <- unique(unlist(lapply(old, fix_gene)))
    if (!identical(old, new)) { sets[[id]] <- new; log[nrow(log) + 1, ] <- c(id, paste("genes", paste(setdiff(old, new), collapse = ","), "->", paste(setdiff(new, old), collapse = ","))) }
  }
  db$cofactor <- sets_to_table(sets, "cofactor", max(16, length(cof_cols)))
  return(list(db = db, log = log))
}

# Interactions whose partner CellChat can never resolve, removed with the reason
remove_unresolvable_partners <- function(db, species, extra = character(0), extra_reason = "") {
  it <- db$interaction; gi <- db$geneInfo$Symbol
  p <- unique(c(it$ligand, it$receptor))
  unres <- p[!p %in% gi & !p %in% rownames(db$complex)]
  reason <- setNames(rep(NA_character_, length(unres)), unres)
  reason[grepl(" ", unres)] <- "named complex without definition"
  if (species != "mouse") reason[is.na(reason) & grepl("^h2-", unres, ignore.case = TRUE)] <- "mouse MHC (H2) symbol in a non-mouse species"
  fused <- vapply(strsplit(unres, "-", fixed = TRUE), function(x) length(x) == 2 && all(x %in% gi), logical(1))
  reason[is.na(reason) & fused] <- "two gene symbols fused with '-'"
  reason <- c(reason[!is.na(reason)], setNames(rep(extra_reason, length(extra)), extra))
  hit_l <- it$ligand %in% names(reason); hit_r <- it$receptor %in% names(reason); hit <- hit_l | hit_r
  src <- if ("sources" %in% colnames(it)) it$sources[hit] else rep(NA_character_, sum(hit))
  removed <- data.frame(species = rep(species, sum(hit)), interaction_name = it$interaction_name[hit], ligand = it$ligand[hit], receptor = it$receptor[hit], sources = src,
                        reason = unname(ifelse(hit_l[hit], reason[it$ligand[hit]], reason[it$receptor[hit]])), stringsAsFactors = FALSE)
  db$interaction <- it[!hit, , drop = FALSE]
  return(list(db = db, removed = removed))
}

# All repairs on one species database; returns the repaired database, the repair log and the removed interactions
repair_species_db <- function(db, species, ref = cellchat_reference(species)) {
  cpx <- restore_complex_definitions(db, ref$complex)
  cof <- repair_cofactors(cpx$db, ref$cofactor, species)
  rem <- remove_unresolvable_partners(cof$db, species, extra = cpx$unrestored, extra_reason = "complex with numeric subunit, no definition in the source")
  out <- rem$db
  # Complex rows still holding a numeric subunit are no longer used: dropped
  if (length(cpx$unrestored) > 0) out$complex <- out$complex[!rownames(out$complex) %in% cpx$unrestored, , drop = FALSE]
  repairs <- rbind(
    if (nrow(cpx$log) > 0) data.frame(species = species, item = cpx$log$complex, action = paste(cpx$log$action, "->", cpx$log$subunits), stringsAsFactors = FALSE),
    if (nrow(cof$log) > 0) data.frame(species = species, item = cof$log$item, action = cof$log$action, stringsAsFactors = FALSE))
  return(list(db = out, repairs = repairs, removed = rem$removed))
}

# Case-insensitive ligand|receptor key, shared by orthologous pairs that keep the same symbol across species
pair_key <- function(it) tolower(paste(it$ligand, it$receptor, sep = "|"))

# Lookup table key -> pathway from one reference species (first non-missing label per key)
pathway_lookup <- function(it) {
  ok <- !is.na(it$pathway_name) & it$pathway_name != ""
  keys <- pair_key(it)[ok]; pw <- it$pathway_name[ok]
  return(pw[!duplicated(keys)] |> setNames(keys[!duplicated(keys)]))
}

# Fill missing pathway labels from the human then the mouse pair, record the origin, label the rest with the placeholder
propagate_pathways <- function(sp_db, lookups, placeholder) {
  it <- sp_db$interaction
  missing <- is.na(it$pathway_name) | it$pathway_name == ""
  it$pathway_origin <- ifelse(missing, NA_character_, "source")
  keys <- pair_key(it)
  for (ref in names(lookups)) {
    fill <- missing & is.na(it$pathway_origin) & keys %in% names(lookups[[ref]])
    it$pathway_name[fill] <- unname(lookups[[ref]][keys[fill]])
    it$pathway_origin[fill] <- paste0("inherited_", ref)
  }
  rest <- is.na(it$pathway_origin)
  it$pathway_name[rest] <- placeholder; it$pathway_origin[rest] <- "unassigned"
  sp_db$interaction <- it
  return(sp_db)
}

# Canonical CellChat interaction categories
CANONICAL_TYPES <- c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor", "Non-protein Signaling")

# interaction_type taken from the CellChat annotation column when it holds a canonical category (the conversion lost these labels)
harmonize_interaction_type <- function(sp_db) {
  it <- sp_db$interaction
  if (!"interaction_type" %in% colnames(it)) it$interaction_type <- NA_character_
  if ("annotation" %in% colnames(it)) {
    use_ann <- !is.na(it$annotation) & it$annotation %in% CANONICAL_TYPES
    it$interaction_type[use_ann] <- it$annotation[use_ann]
  }
  sp_db$interaction <- it
  return(sp_db)
}
