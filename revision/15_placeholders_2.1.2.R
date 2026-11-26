# Step 15 (NARGAB-2026-310 revision): numbers and figures still in yellow in the manuscript and the letter, on release 2.1.2.
# Run from the CellHub_DB root: Rscript revision/15_placeholders_2.1.2.R
# Needs release/2.1.2/CellHubDB.rds and the 04 / 04b result files (found by name below this folder, ~/Bureau or Downloads). No CellChat run.
# Output: revision/results/15_placeholders_2.1.2/ (SUMMARY.txt, CSV tables, Fig 2 and Fig 3d-e PNG/PDF, raw counts for Fig 1)
suppressPackageStartupMessages({ library(ggplot2) })
release_path <- "release/2.1.2/CellHubDB.rds"
# Most recent copy of a 04 / 04b result file: searched in this folder, then in the Desktop and Downloads folders
find_result <- function(name) {
  hits <- unlist(lapply(c(".", "~/Bureau", "~/Téléchargements", "~/Downloads"), function(d) if (dir.exists(d)) list.files(d, pattern = paste0("^", name, "$"), recursive = TRUE, full.names = TRUE)))
  if (length(hits) == 0) return(NA_character_)
  return(hits[which.max(file.mtime(hits))])
}
out_dir <- "revision/results/15_placeholders_2.1.2"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(out_dir, "SUMMARY.txt"), split = TRUE)
db <- readRDS(release_path)
cat("Release:", attr(db, "version"), "| MD5", unname(tools::md5sum(release_path)), "| species:", length(db), "\n\n")
tier_levels <- c("high", "medium", "low")

# Split a partner name into its genes (complexes are written GENE1_GENE2 or listed in the complex table)
partner_genes <- function(p, cpx) {
  if (!is.null(cpx) && nrow(cpx) > 0 && p %in% rownames(cpx)) { g <- unlist(cpx[p, ], use.names = FALSE); return(toupper(g[!is.na(g) & g != ""])) }
  return(toupper(strsplit(p, "_", fixed = TRUE)[[1]]))
}

# Recommended set as described in the manuscript: high tier + medium tier with corroborated roles
recommended <- function(it) it$confidence_tier == "high" | (it$confidence_tier == "medium" & it$role_corroborated %in% TRUE)

## 1. EGFR example (Reviewer 2, major 3): interactions involving EGFR in the merged human database
hu <- db[["human"]]; it <- hu$interaction
has_egfr <- vapply(seq_len(nrow(it)), function(i) "EGFR" %in% c(partner_genes(it$ligand[i], hu$complex), partner_genes(it$receptor[i], hu$complex)), logical(1))
eg <- it[has_egfr, ]; eg$recommended <- recommended(eg)
cat("== EGFR (human, merged) ==\n")
cat("Interactions involving EGFR:", nrow(eg), "\n")
print(table(tier = factor(eg$confidence_tier, tier_levels), role_corroborated = eg$role_corroborated, useNA = "ifany"))
cat("In the recommended set:", sum(eg$recommended), "| excluded:", sum(!eg$recommended), "\n")
cat("Supported by >= 2 sources:", sum(lengths(strsplit(eg$sources, ",\\s*")) >= 2), "\n")
probe <- c("CALM2", "ARF4", "GSTP1", "HSP90AA1")
probe_rows <- eg[toupper(eg$ligand) %in% probe | toupper(eg$receptor) %in% probe, c("ligand", "receptor", "sources", "confidence_tier", "role_corroborated", "recommended")]
cat("\nPartners cited by the reviewer:\n"); print(probe_rows, row.names = FALSE)
write.csv(eg, file.path(out_dir, "EGFR_human_interactions.csv"), row.names = FALSE)

## 2. Is the 04 benchmark (run on the 2.1 build 338ba8a1) still valid for 2.1.2? Compare the mouse database used by CellChat
cat("\n== Mouse database: benchmark build vs 2.1.2 ==\n")
# Candidate files: multi-species bundles only (a list with a "mouse" element holding an interaction table)
old_path <- list.files(c("archive", "release"), pattern = "^CellHubDB(_[0-9.]+)?\\.rds$", recursive = TRUE, full.names = TRUE)
old_path <- old_path[grepl("338ba8a1", old_path)]
old_path <- Filter(function(f) { x <- tryCatch(readRDS(f), error = function(e) NULL); is.list(x) && is.data.frame(x[["mouse"]]$interaction) }, old_path)
mo <- db[["mouse"]]; key_cols <- c("ligand", "receptor", "confidence_tier", "role_corroborated", "sources")
# Interaction table reduced to the compared columns, sorted by ligand-receptor key so row order does not matter
sorted_keys <- function(it) { x <- it[, intersect(key_cols, colnames(it)), drop = FALSE]; x[order(toupper(x$ligand), toupper(x$receptor)), , drop = FALSE] }
if (length(old_path) == 0) {
  cat("Build 338ba8a1 not found under archive/ or release/: benchmark validity not checked\n")
} else {
  old <- readRDS(old_path[[1]])[["mouse"]]
  cat("Columns missing in the benchmark build:", paste(setdiff(key_cols, colnames(old$interaction)), collapse = ", "), "\n")
  same_pairs <- isTRUE(all.equal(sorted_keys(old$interaction), sorted_keys(mo$interaction), check.attributes = FALSE))
  cat("File:", old_path[[1]], "| mouse interactions: benchmark build", nrow(old$interaction), "| 2.1.2", nrow(mo$interaction), "\n")
  cat("Same mouse interactions, tiers, roles and sources:", same_pairs, "\n")
  cat("Same mouse complexes:", isTRUE(all.equal(old$complex, mo$complex, check.attributes = FALSE)), "\n")
  gene_col <- intersect(c("Symbol", "symbol", "gene"), colnames(mo$geneInfo))[1]
  partners <- unique(unlist(lapply(c(mo$interaction$ligand, mo$interaction$receptor), partner_genes, cpx = mo$complex)))
  missing_old <- setdiff(partners, toupper(old$geneInfo[[gene_col]])); missing_new <- setdiff(partners, toupper(mo$geneInfo[[gene_col]]))
  cat("Partner genes missing from geneInfo: benchmark build", length(missing_old), "| 2.1.2", length(missing_new), "\n")
  cat(if (same_pairs && identical(sort(missing_old), sort(missing_new))) "=> the benchmark and MuSC results apply to 2.1.2 unchanged\n" else "=> differences: rerun revision/04_benchmark_cellchat.R and 04b on release 2.1.2\n")
}

## 3. MuSC literature check (Reviewer 2, major 3): X supported / Y binding without muscle evidence / Z doubtful
cat("\n== MuSC shortlist ==\n")
sl_path <- find_result("MuSC_new_pairs_shortlist.csv")
if (is.na(sl_path)) { sink(); stop("MuSC_new_pairs_shortlist.csv not found: run revision/04_benchmark_cellchat.R then 04b, or copy the file into revision/results/") }
cat("File:", sl_path, "| modified", format(file.mtime(sl_path)), "\n")
sl <- read.csv(sl_path, stringsAsFactors = FALSE)
new_rec <- sl[sl$recommended_set %in% TRUE & sl$novelty != "subunit_of_default", ]
cat("New LR pair x direction:", nrow(sl), "| subunit_of_default:", sum(sl$novelty == "subunit_of_default"), "| recommended and truly new:", nrow(new_rec), "\n")
lc <- trimws(new_rec$literature_check); lc[is.na(lc) | lc == ""] <- "(empty)"
cat("literature_check values for the", nrow(new_rec), "recommended, truly new pairs:\n"); print(sort(table(lc), decreasing = TRUE))
write.csv(new_rec, file.path(out_dir, "SupplementaryTable_S3_MuSC_recommended_new.csv"), row.names = FALSE)

## 4. Figure 2: (a) interactions per species by confidence tier, (b, c) source overlap in human and mouse
tiers <- do.call(rbind, lapply(names(db), function(sp) data.frame(species = sp, tier = db[[sp]]$interaction$confidence_tier)))
tiers$tier <- factor(tiers$tier, tier_levels)
order_sp <- names(sort(table(tiers$species), decreasing = TRUE)); tiers$species <- factor(tiers$species, order_sp)
cat("\n== Figure 2a: interactions per species and tier ==\n"); print(table(tiers$species, tiers$tier))
p2a <- ggplot(tiers, aes(species, fill = tier)) + geom_bar(width = 0.75) +
  scale_fill_manual(values = c(high = "#1b7837", medium = "#7fbf7b", low = "#d9d9d9"), name = "Confidence tier") +
  labs(x = NULL, y = "Interactions") + theme_classic(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
for (ext in c("png", "pdf")) ggsave(file.path(out_dir, paste0("Fig2a_tiers_per_species.", ext)), p2a, width = 9, height = 4.5, dpi = 300, bg = "white")
if (requireNamespace("UpSetR", quietly = TRUE)) {
  for (sp in c("human", "mouse")) {
    src <- strsplit(db[[sp]]$interaction$sources, ",\\s*"); all_src <- sort(unique(unlist(src)))
    m <- as.data.frame(sapply(all_src, function(s) as.integer(vapply(src, function(x) s %in% x, logical(1)))))
    for (ext in c("png", "pdf")) {
      f <- file.path(out_dir, paste0("Fig2", ifelse(sp == "human", "b", "c"), "_upset_", sp, ".", ext))
      if (ext == "png") png(f, width = 2400, height = 1500, res = 300) else pdf(f, width = 8, height = 5)
      print(UpSetR::upset(m, nsets = length(all_src), nintersects = 25, order.by = "freq", mainbar.y.label = paste("Interactions,", sp), sets.x.label = "Interactions per source"))
      dev.off()
    }
  }
} else cat("UpSetR not installed: Fig 2b-c skipped\n")

## 5. Figure 3d-e: MuSC interactions of the recommended set absent from CellChatDB (bubble plot and pathway chord diagram)
calls_path <- find_result("calls_MuSC_for_manual_check.csv")
if (is.na(calls_path)) { sink(); stop("calls_MuSC_for_manual_check.csv not found: run revision/04_benchmark_cellchat.R") }
cat("\nCalls file:", calls_path, "\n")
calls <- read.csv(calls_path, stringsAsFactors = FALSE)
calls$role <- ifelse(calls$source == "MuSC", "MuSC_sender", "MuSC_receiver")
keep <- paste(calls$ligand, calls$receptor, calls$role) %in% paste(new_rec$ligand, new_rec$receptor, new_rec$role)
f3 <- calls[keep, ]
cat("\n== Figure 3d-e ==\nCalls kept:", nrow(f3), "| LR pairs:", length(unique(paste(f3$ligand, f3$receptor))), "| pathways:", length(unique(f3$pathway_name)), "\n")
f3$lr <- paste(f3$ligand, f3$receptor, sep = " - "); f3$cells <- paste(f3$source, "->", f3$target)
f3$tier <- factor(f3$confidence_tier, tier_levels); f3$logp <- pmin(-log10(pmax(f3$pval, 1e-300)), 3)
p3d <- ggplot(f3, aes(cells, lr)) + geom_point(aes(size = prob, colour = tier)) +
  scale_colour_manual(values = c(high = "#1b7837", medium = "#e08214", low = "#999999"), name = "Confidence tier", drop = FALSE) +
  scale_size_continuous(name = "Communication\nprobability", range = c(1, 5)) + labs(x = NULL, y = NULL) +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(face = "italic"))
h <- max(4, 0.14 * length(unique(f3$lr)) + 1.5); w <- max(5, 0.25 * length(unique(f3$cells)) + 3)
for (ext in c("png", "pdf")) ggsave(file.path(out_dir, paste0("Fig3d_MuSC_bubble.", ext)), p3d, width = w, height = h, dpi = 300, bg = "white", limitsize = FALSE)
if (requireNamespace("circlize", quietly = TRUE) && nrow(f3) > 0) {
  flows <- aggregate(prob ~ source + target, data = f3, FUN = sum)
  for (ext in c("png", "pdf")) {
    f <- file.path(out_dir, paste0("Fig3e_MuSC_chord.", ext))
    if (ext == "png") png(f, width = 2000, height = 2000, res = 300) else pdf(f, width = 7, height = 7)
    circlize::circos.clear(); circlize::chordDiagram(flows, directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow", annotationTrack = c("name", "grid"))
    title(paste0("MuSC interactions absent from CellChatDB (", length(unique(f3$pathway_name)), " pathways)"), cex.main = 0.8)
    dev.off(); circlize::circos.clear()
  }
  write.csv(aggregate(prob ~ pathway_name + source + target, data = f3, FUN = sum), file.path(out_dir, "Fig3e_pathway_flows.csv"), row.names = FALSE)
} else cat("circlize not installed or no call: Fig 3e skipped\n")

## 6. Figure 1: per-step counts from the pipeline logs (raw lines, to be read and placed on the workflow diagram)
cat("\n== Figure 1: lines with counts in the pipeline reports ==\n")
rep_files <- c(list.files("output/reports", full.names = TRUE, recursive = TRUE), list.files("release/2.1.2", pattern = "\\.(txt|md|log)$", full.names = TRUE, recursive = TRUE))
rep_files <- rep_files[grepl("\\.(txt|md|log|csv)$", rep_files)]
fig1 <- unlist(lapply(rep_files, function(f) {
  l <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
  hit <- l[grepl("[0-9]", l) & grepl("convert|remov|duplic|merge|entries|interactions|renamed|repair|dropped|kept|final|total", l, ignore.case = TRUE)]
  if (length(hit)) c(paste0("--- ", f), head(hit, 60)) else NULL
}))
writeLines(fig1, file.path(out_dir, "Fig1_counts_raw.txt"))
cat(length(rep_files), "report files scanned; lines written to Fig1_counts_raw.txt\n")
cat("Final interactions per species (2.1.2):\n"); print(vapply(db, function(x) nrow(x$interaction), integer(1)))
sink()
message("Done: ", normalizePath(out_dir))
