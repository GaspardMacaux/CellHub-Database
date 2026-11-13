# Step 01 (Reviewer 2, major 7 + other 1): structural audit of CellHubDB.rds and CellChat smoke test on all species,
# before and after schema standardization
source("review/revision_helpers.R")
suppressPackageStartupMessages(library(CellChat))
future::plan("sequential")
db_path <- "CellHubDB.rds"; out_dir <- "review/01_audit"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
db <- readRDS(db_path)
schema <- get_cellchat_schema()

# Structural audit of the released file
audit <- do.call(rbind, lapply(names(db), function(sp) audit_species_db(db[[sp]], sp, schema)))
write.csv(audit, file.path(out_dir, "audit_released.csv"), row.names = FALSE)
print(audit[, c("species", "n_interaction", "n_complex", "n_cofactor", "n_geneinfo", "cofactor_ncol", "n_multisubunit", "n_multisubunit_undefined", "pct_pathway_na")])

# CellChat smoke test on the released file
smoke_released <- do.call(rbind, lapply(names(db), function(sp) { message("Smoke test (released): ", sp); run_cellchat_smoke(db[[sp]], sp) }))
write.csv(smoke_released, file.path(out_dir, "smoke_released.csv"), row.names = FALSE)

# Same test after schema standardization (no interaction added or removed)
db_fixed <- lapply(db, standardize_cellchat_schema, schema = schema)
saveRDS(db_fixed, file.path(out_dir, "CellHubDB_schema_fixed.rds"))
smoke_fixed <- do.call(rbind, lapply(names(db_fixed), function(sp) { message("Smoke test (fixed): ", sp); run_cellchat_smoke(db_fixed[[sp]], sp) }))
write.csv(smoke_fixed, file.path(out_dir, "smoke_schema_fixed.csv"), row.names = FALSE)

print(merge(smoke_released[, c("species", "first_error")], smoke_fixed[, c("species", "n_significant", "first_error")], by = "species", suffixes = c("_released", "_fixed")))
