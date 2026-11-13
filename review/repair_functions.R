# Repairs of CellHubDB species databases for the defects found by the release gate (steps 10 / 10b):
# complexes generated from suffixed names (numeric subunits), undefined cofactor references, cofactor genes in another species' case,
# partners that CellChat can never resolve (named complexes without definition, mouse MHC in another species, two symbols fused with "-")
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
