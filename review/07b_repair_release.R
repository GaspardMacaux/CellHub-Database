# Step 07b: repairs of the release database for the defects found by the gate (steps 10 / 10b), version attribute added
# Always starts from the pre-repair copy, so it can be rerun safely. Run before 08 (app files) and 09 (gate).
# Usage: Rscript review/07b_repair_release.R <version label>
source("review/repair_functions.R")
suppressPackageStartupMessages(library(CellChat))
args <- commandArgs(trailingOnly = TRUE)
version <- if (length(args) >= 1) args[1] else "2.0.1"
release_path <- "release_v2.1/CellHubDB.rds"; backup_path <- "release_v2.1/CellHubDB_pre_repair.rds"; out_dir <- "review/07b_repair"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The pre-repair copy is written once, then always used as input
if (!file.exists(backup_path)) file.copy(release_path, backup_path)
db <- readRDS(backup_path)

results <- lapply(names(db), function(sp) { message("Repairing: ", sp); repair_species_db(db[[sp]], sp) })
names(results) <- names(db)
repaired <- lapply(results, `[[`, "db")
attr(repaired, "version") <- version

repairs <- do.call(rbind, lapply(results, `[[`, "repairs"))
removed <- do.call(rbind, lapply(results, `[[`, "removed"))
write.csv(repairs, file.path(out_dir, "repairs.csv"), row.names = FALSE)
write.csv(removed, file.path(out_dir, "removed_interactions.csv"), row.names = FALSE)

# Interactions before / after per species
counts <- data.frame(species = names(db), before = vapply(db, function(x) nrow(x$interaction), integer(1)),
                     after = vapply(repaired, function(x) nrow(x$interaction), integer(1)), row.names = NULL)
counts$removed <- counts$before - counts$after
write.csv(counts, file.path(out_dir, "interaction_counts.csv"), row.names = FALSE)

saveRDS(repaired, release_path)
print(counts[counts$removed > 0, ], row.names = FALSE)
if (!is.null(removed) && nrow(removed) > 0) print(table(species = removed$species, reason = removed$reason))
if (!is.null(repairs) && nrow(repairs) > 0) print(repairs, row.names = FALSE)
message("Release written: ", release_path, " (version ", version, ") | pre-repair copy: ", backup_path)
