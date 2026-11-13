# Step 10b: origin of the gate warnings that step 10 could not attribute (zebrafish cofactors and numeric complexes, mouse unresolvable partners, arabidopsis cofactor)
# Usage: Rscript review/10b_trace_gate_warnings.R [app_databases_dir]
suppressPackageStartupMessages(library(CellChat))
args <- commandArgs(trailingOnly = TRUE)
db_dir <- if (length(args) >= 1) args[1] else "app_databases_v2.0.1"
out_dir <- "review/10_gate_warnings"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cof_cols_all <- c("agonist", "antagonist", "co_A_receptor", "co_I_receptor")

# Cofactor identifiers referenced by interactions but missing from the cofactor table
undefined_cofactors <- function(db) {
  cols <- intersect(cof_cols_all, colnames(db$interaction))
  ids <- setdiff(unique(unlist(db$interaction[, cols, drop = FALSE])), c("", NA))
  return(setdiff(ids, rownames(db$cofactor)))
}

# Interaction partners resolvable by CellChat: a gene of geneInfo or a defined complex
unresolvable_partners <- function(db) {
  p <- unique(c(db$interaction$ligand, db$interaction$receptor))
  return(p[!p %in% db$geneInfo$Symbol & !p %in% rownames(db$complex)])
}

# Source of the interactions that use given partners
source_of_partners <- function(db, genes) {
  it <- db$interaction
  hit <- it$ligand %in% genes | it$receptor %in% genes
  src <- if ("sources" %in% colnames(it)) it$sources[hit] else if ("source_db" %in% colnames(it)) it$source_db[hit] else rep(NA, sum(hit))
  return(data.frame(ligand = it$ligand[hit], receptor = it$receptor[hit], sources = src, stringsAsFactors = FALSE))
}

z <- readRDS(file.path(db_dir, "zebrafish", "CellHubDB_zebrafish_merged.rds"))
ref_z <- CellChat::CellChatDB.zebrafish

# 1. Are the 21 undefined cofactor references already undefined in CellChatDB.zebrafish itself?
undef_z <- undefined_cofactors(z); undef_ref <- undefined_cofactors(ref_z)
message("\n== 1. Zebrafish undefined cofactor references")
message("CellHubDB: ", length(undef_z), " | CellChatDB.zebrafish: ", length(undef_ref), " | shared: ", length(intersect(undef_z, undef_ref)))
message("Only in CellHubDB (build defect): ", paste(setdiff(undef_z, undef_ref), collapse = ", "))
# Cofactor identifiers of the source with a numbered suffix (e.g. 'TGFb antagonist_1') matching an undefined base name
message("Source cofactor ids sharing the base name of an undefined reference: ",
        paste(rownames(ref_z$cofactor)[sub("_[0-9]+$", "", rownames(ref_z$cofactor)) %in% undef_z], collapse = ", "))

# 2. Origin of the complexes with numeric subunits (acvr1_bmpr2_1 ...)
cpx_num <- rownames(z$complex)[apply(z$complex, 1, function(r) any(grepl("^[0-9]+$", r)))]
message("\n== 2. Zebrafish complexes with numeric subunits: ", length(cpx_num))
print(data.frame(complex = cpx_num, defined_in_source_complex = cpx_num %in% rownames(ref_z$complex),
                 used_in_source_interactions = cpx_num %in% c(ref_z$interaction$ligand, ref_z$interaction$receptor)), row.names = FALSE)
print(head(source_of_partners(z, cpx_num), 20), row.names = FALSE)
if (any(cpx_num %in% rownames(ref_z$complex))) { message("Source definitions:"); print(ref_z$complex[intersect(cpx_num, rownames(ref_z$complex)), , drop = FALSE]) }

# 3. Mouse-case cofactor genes of the zebrafish table: is there a zebrafish symbol differing only by case?
cof_genes <- unique(unlist(z$cofactor)); cof_genes <- cof_genes[!is.na(cof_genes) & cof_genes != "" & grepl("^[A-Z][a-z]", cof_genes)]
zf_lower <- tolower(z$geneInfo$Symbol)
message("\n== 3. Zebrafish cofactor genes in mouse case: ", length(cof_genes))
print(data.frame(gene = cof_genes, lowercase_in_geneInfo = tolower(cof_genes) %in% zf_lower,
                 paralogs_in_geneInfo = vapply(tolower(cof_genes), function(g) paste(z$geneInfo$Symbol[grepl(paste0("^", g, "[ab]?$"), zf_lower)], collapse = " "), character(1))),
      row.names = FALSE)

# 4. Mouse partners resolvable neither as gene nor as complex (the hyphen check of step 10 also caught valid H2- genes)
m <- readRDS(file.path(db_dir, "mouse", "CellHubDB_mouse_merged.rds"))
unres_m <- unresolvable_partners(m)
message("\n== 4. Mouse unresolvable partners: ", length(unres_m))
tab_m <- source_of_partners(m, unres_m)
write.csv(tab_m, file.path(out_dir, "mouse_unresolvable_partners.csv"), row.names = FALSE)
print(head(tab_m, 30), row.names = FALSE)

# 5. Arabidopsis undefined cofactor: interactions concerned
message("\n== 5. Arabidopsis interactions with an undefined cofactor")
print(read.csv(file.path(out_dir, "arabidopsis_undefined_cofactor_interactions.csv"), stringsAsFactors = FALSE), row.names = FALSE)
a <- readRDS(file.path(db_dir, "arabidopsis", "CellHubDB_arabidopsis_merged.rds"))
message("AT3G03950 in arabidopsis geneInfo: ", "AT3G03950" %in% a$geneInfo$Symbol, " | cofactor ids of the table: ", paste(head(rownames(a$cofactor), 10), collapse = ", "))
