################################################################################
#         cellhubdb_stats.R — Complete statistics report for manuscript         #
#                                                                              #
# Outputs a formatted Markdown report + console log with all numbers needed    #
# for both papers (Application Note + Data Resource).                          #
################################################################################

DATABASES_DIR <- "databases"
OUTPUT_FILE <- "CellHubDB_Statistics_Report.md"

# Collect all markdown lines
md_lines <- c()
md <- function(...) {
  line <- paste0(...)
  md_lines <<- c(md_lines, line)
}

md("# CellHub Database — Statistics Report")
md(sprintf("*Generated: %s*\n", Sys.Date()))
md("---")
md("")

# ==============================================================================
# Load all merged databases
# ==============================================================================

species_dirs <- sort(list.dirs(DATABASES_DIR, recursive = FALSE))
species_list <- basename(species_dirs)

all_merged <- list()
grand_total <- 0

for (sp in species_list) {
  merged_file <- list.files(file.path(DATABASES_DIR, sp), pattern = "_merged\\.rds$", full.names = TRUE)
  if (length(merged_file) == 0) next
  db <- readRDS(merged_file[1])
  all_merged[[sp]] <- db$interaction
  grand_total <- grand_total + nrow(db$interaction)
}

# ==============================================================================
# 1. GRAND TOTALS
# ==============================================================================

md("## 1. Grand Totals")
md("")
md(sprintf("- **Total species:** %d", length(species_list)))
md(sprintf("- **Total interactions (merged):** %s", format(grand_total, big.mark = ",")))
md("")
md("### Per-species summary")
md("")
md("| Species | Interactions | Sources | Multi-source | Multi% | Pathways | Path% | Evid% |")
md("|---------|-------------|---------|-------------|--------|----------|-------|-------|")

species_stats <- data.frame(species = character(), n = integer(), stringsAsFactors = FALSE)

for (sp in species_list) {
  if (!sp %in% names(all_merged)) next
  inter <- all_merged[[sp]]
  n <- nrow(inter)
  
  # Count individual source files
  ind_files <- list.files(file.path(DATABASES_DIR, sp), pattern = "\\.rds$")
  ind_files <- ind_files[!grepl("_merged\\.rds$", ind_files)]
  n_src <- length(ind_files)
  
  # Multi-source
  n_multi <- if ("sources" %in% colnames(inter)) sum(grepl(",", inter$sources), na.rm = TRUE) else 0
  pct_multi <- round(100 * n_multi / n, 1)
  
  # Pathways
  n_pw <- sum(!is.na(inter$pathway_name))
  n_unique_pw <- length(unique(inter$pathway_name[!is.na(inter$pathway_name)]))
  pct_pw <- round(100 * n_pw / n, 1)
  
  # Evidence
  n_ev <- if ("evidence" %in% colnames(inter)) sum(!is.na(inter$evidence)) else 0
  pct_ev <- round(100 * n_ev / n, 1)
  
  md(sprintf("| %s | %s | %d | %s | %s%% | %d | %s%% | %s%% |",
    sp, format(n, big.mark = ","), n_src, format(n_multi, big.mark = ","),
    pct_multi, n_unique_pw, pct_pw, pct_ev))
  
  species_stats <- rbind(species_stats, data.frame(species = sp, n = n, stringsAsFactors = FALSE))
}

md(sprintf("| **TOTAL** | **%s** | | | | | | |", format(grand_total, big.mark = ",")))
md("")

# ==============================================================================
# 2. SOURCE × SPECIES MATRIX
# ==============================================================================

md("## 2. Source × Species Matrix")
md("")

# Build matrix
source_data <- list()
all_sources <- c()

for (sp in species_list) {
  sp_dir <- file.path(DATABASES_DIR, sp)
  ind_files <- list.files(sp_dir, pattern = "\\.rds$", full.names = TRUE)
  ind_files <- ind_files[!grepl("_merged\\.rds$", ind_files)]
  
  for (f in ind_files) {
    db <- readRDS(f)
    src <- unique(db$interaction$source_db)[1]
    if (is.na(src)) next
    all_sources <- unique(c(all_sources, src))
    if (is.null(source_data[[src]])) source_data[[src]] <- list()
    source_data[[src]][[sp]] <- nrow(db$interaction)
  }
}
all_sources <- sort(all_sources)

# Header
header <- paste(c("| Source", sapply(species_list, function(s) substr(s, 1, 8)), "Total |"), collapse = " | ")
md(header)
sep <- paste(c("|---", rep("---", length(species_list)), "---|"), collapse = "|")
md(sep)

source_totals <- c()
for (src in all_sources) {
  vals <- c()
  src_total <- 0
  for (sp in species_list) {
    n <- source_data[[src]][[sp]]
    if (is.null(n)) {
      vals <- c(vals, "-")
    } else {
      vals <- c(vals, format(n, big.mark = ","))
      src_total <- src_total + n
    }
  }
  md(sprintf("| %s | %s | %s |", src, paste(vals, collapse = " | "), format(src_total, big.mark = ",")))
  source_totals[src] <- src_total
}
md("")

md("### Per-source summary")
md("")
md("| Source | Total interactions | Species covered |")
md("|--------|-------------------|-----------------|")
for (src in names(sort(source_totals, decreasing = TRUE))) {
  n_sp <- length(source_data[[src]])
  md(sprintf("| %s | %s | %d |", src, format(source_totals[src], big.mark = ","), n_sp))
}
md(sprintf("| **Total sources** | | **%d** |", length(all_sources)))
md("")

# ==============================================================================
# 3. CROSS-SOURCE OVERLAP
# ==============================================================================

md("## 3. Cross-source Overlap")
md("")

for (sp in c("human", "mouse")) {
  if (!sp %in% names(all_merged)) next
  inter <- all_merged[[sp]]
  if (!"sources" %in% colnames(inter)) next
  
  md(sprintf("### %s (%s interactions)", toupper(sp), format(nrow(inter), big.mark = ",")))
  md("")
  
  # By number of sources
  n_sources_per <- sapply(inter$sources, function(s) length(strsplit(s, ", ")[[1]]))
  
  md("**By number of contributing sources:**")
  md("")
  md("| Sources | Interactions | % |")
  md("|---------|-------------|---|")
  for (n in sort(unique(n_sources_per))) {
    count <- sum(n_sources_per == n)
    md(sprintf("| %d source(s) | %s | %.1f%% |", n, format(count, big.mark = ","), 100 * count / nrow(inter)))
  }
  md("")
  
  # Top combinations
  multi <- inter[n_sources_per > 1, ]
  if (nrow(multi) > 0) {
    combos <- sort(table(multi$sources), decreasing = TRUE)
    md(sprintf("**Top 10 source combinations (%s multi-source entries):**", format(nrow(multi), big.mark = ",")))
    md("")
    md("| Combination | Count |")
    md("|-------------|-------|")
    for (combo in head(names(combos), 10)) {
      md(sprintf("| %s | %d |", combo, combos[combo]))
    }
    md("")
  }
  
  # Per-source contribution
  md("**Per-source contribution:**")
  md("")
  md("| Source | Winning | Present in | Exclusive |")
  md("|--------|---------|-----------|-----------|")
  all_src_sp <- sort(unique(unlist(strsplit(inter$sources, ", "))))
  for (src in all_src_sp) {
    n_winner <- sum(inter$source_db == src, na.rm = TRUE)
    n_present <- sum(grepl(src, inter$sources, fixed = TRUE))
    n_exclusive <- sum(inter$sources == src, na.rm = TRUE)
    md(sprintf("| %s | %s | %s | %s |", src,
      format(n_winner, big.mark = ","), format(n_present, big.mark = ","),
      format(n_exclusive, big.mark = ",")))
  }
  md("")
}

# ==============================================================================
# 4. INTERACTION TYPE DISTRIBUTION
# ==============================================================================

md("## 4. Interaction Type Distribution")
md("")

for (sp in c("human", "mouse", "ALL")) {
  if (sp == "ALL") {
    inter_list <- do.call(rbind, all_merged)
    label <- "All species"
    total <- nrow(inter_list)
  } else {
    if (!sp %in% names(all_merged)) next
    inter_list <- all_merged[[sp]]
    label <- toupper(sp)
    total <- nrow(inter_list)
  }
  
  md(sprintf("### %s", label))
  md("")
  md("| Type | Count | % |")
  md("|------|-------|---|")
  type_table <- table(inter_list$interaction_type, useNA = "always")
  type_names <- names(type_table)
  for (i in seq_along(type_table)) {
    # The NA bin has an NA name: index by position, not by name
    tname <- if (is.na(type_names[i])) "Not annotated" else type_names[i]
    n_type <- as.integer(type_table[i])
    md(sprintf("| %s | %s | %.1f%% |", tname, format(n_type, big.mark = ","), 100 * n_type / total))
  }
  md("")
}

# ==============================================================================
# 5. PATHWAY COVERAGE
# ==============================================================================

md("## 5. Pathway Coverage")
md("")

for (sp in c("human", "mouse")) {
  if (!sp %in% names(all_merged)) next
  inter <- all_merged[[sp]]
  
  n_pw <- sum(!is.na(inter$pathway_name))
  n_unique <- length(unique(inter$pathway_name[!is.na(inter$pathway_name)]))
  
  md(sprintf("### %s", toupper(sp)))
  md("")
  md(sprintf("- Interactions with pathway: %s / %s (%.1f%%)",
    format(n_pw, big.mark = ","), format(nrow(inter), big.mark = ","), 100 * n_pw / nrow(inter)))
  md(sprintf("- Unique pathways: %d", n_unique))
  md("")
  
  pw_table <- sort(table(inter$pathway_name), decreasing = TRUE)
  md("**Top 15 pathways:**")
  md("")
  md("| Pathway | Count |")
  md("|---------|-------|")
  for (pw in head(names(pw_table), 15)) {
    md(sprintf("| %s | %d |", pw, pw_table[pw]))
  }
  md("")
}

# ==============================================================================
# 6. HUMAN-MOUSE OVERLAP
# ==============================================================================

md("## 6. Human-Mouse Ortholog Overlap")
md("")

if (all(c("human", "mouse") %in% names(all_merged))) {
  human <- all_merged[["human"]]
  mouse <- all_merged[["mouse"]]
  
  human_pairs <- unique(paste(toupper(human$ligand), toupper(human$receptor), sep = "|||"))
  mouse_pairs <- unique(paste(toupper(mouse$ligand), toupper(mouse$receptor), sep = "|||"))
  
  overlap <- sum(mouse_pairs %in% human_pairs)
  
  md(sprintf("- Human unique L-R pairs: %s", format(length(human_pairs), big.mark = ",")))
  md(sprintf("- Mouse unique L-R pairs: %s", format(length(mouse_pairs), big.mark = ",")))
  md(sprintf("- Overlap: %s pairs (%.1f%% of mouse, %.1f%% of human)",
    format(overlap, big.mark = ","),
    100 * overlap / length(mouse_pairs),
    100 * overlap / length(human_pairs)))
  md("")
}

# ==============================================================================
# 7. GENE STATISTICS
# ==============================================================================

md("## 7. Gene Statistics")
md("")

md("| Species | Unique ligands | Unique receptors | Unique genes | Both L&R |")
md("|---------|---------------|-----------------|-------------|----------|")

for (sp in c("human", "mouse")) {
  if (!sp %in% names(all_merged)) next
  inter <- all_merged[[sp]]
  
  all_lig <- unique(inter$ligand)
  all_rec <- unique(inter$receptor)
  all_genes <- unique(c(unlist(strsplit(inter$ligand, "_")), unlist(strsplit(inter$receptor, "_"))))
  both <- length(intersect(all_lig, all_rec))
  
  md(sprintf("| %s | %s | %s | %s | %d |",
    sp, format(length(all_lig), big.mark = ","), format(length(all_rec), big.mark = ","),
    format(length(all_genes), big.mark = ","), both))
}
md("")

# ==============================================================================
# 8. ORGANISM CATEGORIES
# ==============================================================================

md("## 8. Organism Categories")
md("")

mammals <- c("human", "mouse", "rat", "chimp", "macaque", "marmoset",
             "cow", "dog", "horse", "pig", "sheep")
vertebrates_other <- c("chicken", "frog", "zebrafish")
invertebrates <- c("drosophila", "celegans")
plants <- c("arabidopsis", "maize", "rice", "soybean", "tomato")

present_mammals <- species_list[species_list %in% mammals]
present_vert <- species_list[species_list %in% vertebrates_other]
present_inv <- species_list[species_list %in% invertebrates]
present_plants <- species_list[species_list %in% plants]

md("| Category | Count | Species |")
md("|----------|-------|---------|")
md(sprintf("| Mammals | %d | %s |", length(present_mammals), paste(present_mammals, collapse = ", ")))
md(sprintf("| Other vertebrates | %d | %s |", length(present_vert), paste(present_vert, collapse = ", ")))
md(sprintf("| Invertebrates | %d | %s |", length(present_inv), paste(present_inv, collapse = ", ")))
md(sprintf("| Plants | %d | %s |", length(present_plants), paste(present_plants, collapse = ", ")))
md(sprintf("| **Total** | **%d** | |", length(species_list)))
md("")

# ==============================================================================
# 9. CURATION IMPACT
# ==============================================================================

md("## 9. Curation Impact")
md("")

removed_files <- list.files(DATABASES_DIR, pattern = "removed_.*\\.csv$", full.names = TRUE)
for (f in removed_files) {
  removed <- read.csv(f, stringsAsFactors = FALSE)
  md(sprintf("### %s (%s entries removed)", basename(f), format(nrow(removed), big.mark = ",")))
  md("")
  md("| Source | Removed |")
  md("|--------|---------|")
  src_table <- table(removed$source_db)
  for (src in names(src_table)) {
    md(sprintf("| %s | %d |", src, src_table[src]))
  }
  md("")
}

# ==============================================================================
# 10. KEY NUMBERS FOR MANUSCRIPT
# ==============================================================================

md("## 10. Key Numbers for Manuscript")
md("")
md("```")
md(sprintf("Total curated interactions:  %s", format(grand_total, big.mark = ",")))
md(sprintf("Total species:              %d", length(species_list)))
md(sprintf("Total source databases:     %d", length(all_sources)))

if ("human" %in% names(all_merged)) {
  h <- all_merged[["human"]]
  n_multi_h <- sum(grepl(",", h$sources))
  md(sprintf("Human interactions:         %s (%s multi-source, %.0f%%)",
    format(nrow(h), big.mark = ","), format(n_multi_h, big.mark = ","), 100 * n_multi_h / nrow(h)))
}
if ("mouse" %in% names(all_merged)) {
  m <- all_merged[["mouse"]]
  n_multi_m <- sum(grepl(",", m$sources))
  md(sprintf("Mouse interactions:         %s (%s multi-source, %.0f%%)",
    format(nrow(m), big.mark = ","), format(n_multi_m, big.mark = ","), 100 * n_multi_m / nrow(m)))
}

md(sprintf("Mammals:                    %d species", length(present_mammals)))
md(sprintf("Other vertebrates:          %d species", length(present_vert)))
md(sprintf("Invertebrates:              %d species", length(present_inv)))
md(sprintf("Plants:                     %d species", length(present_plants)))

if ("drosophila" %in% names(all_merged))
  md(sprintf("Drosophila (first):         %s interactions", format(nrow(all_merged[["drosophila"]]), big.mark = ",")))
if ("celegans" %in% names(all_merged))
  md(sprintf("C. elegans (first):         %s interactions", format(nrow(all_merged[["celegans"]]), big.mark = ",")))

plant_total <- sum(vapply(present_plants, function(sp) {
  if (sp %in% names(all_merged)) nrow(all_merged[[sp]]) else 0L
}, integer(1)))
md(sprintf("Plants total (first):       %s interactions across %d species",
  format(plant_total, big.mark = ","), length(present_plants)))
md("```")
md("")
# ==============================================================================
# UNIQUE BIOLOGICAL PAIRS — global and per species
# ==============================================================================

md("## 11. Unique Biological Pairs")
md("")

# Global: normalize all pairs to uppercase for cross-species comparison
all_pairs_global <- c()
per_species_unique <- list()

for (sp in names(all_merged)) {
  inter <- all_merged[[sp]]
  pairs <- unique(paste(toupper(inter$ligand), toupper(inter$receptor), sep = "|||"))
  per_species_unique[[sp]] <- length(pairs)
  all_pairs_global <- c(all_pairs_global, pairs)
}

n_unique_global <- length(unique(all_pairs_global))

md(sprintf("- **Total species-specific interactions:** %s (sum of merged per species)",
           format(grand_total, big.mark = ",")))
md(sprintf("- **Unique biological L-R pairs (cross-species):** %s",
           format(n_unique_global, big.mark = ",")))
md(sprintf("- **Redundancy ratio:** %.1fx (same pair counted across multiple species)",
           grand_total / n_unique_global))
md("")

# Per species: unique pairs within that species
md("### Unique pairs per species")
md("")
md("| Species | Merged interactions | Unique L-R pairs | Deduplicated |")
md("|---------|-------------------|-----------------|-------------|")
for (sp in names(all_merged)) {
  n_merged <- nrow(all_merged[[sp]])
  n_unique <- per_species_unique[[sp]]
  n_dedup <- n_merged - n_unique
  md(sprintf("| %s | %s | %s | %d |",
             sp, format(n_merged, big.mark = ","), format(n_unique, big.mark = ","), n_dedup))
}
md("")

# ==============================================================================
# PER-SOURCE UNIQUE CONTRIBUTIONS
# ==============================================================================

md("## 12. Per-Source Unique Biological Pairs")
md("")

# Collect pairs per source (across all species, uppercase for comparison)
source_pairs <- list()
for (sp in names(all_merged)) {
  inter <- all_merged[[sp]]
  if (!"sources" %in% colnames(inter)) next
  all_src <- unique(unlist(strsplit(inter$sources, ", ")))
  for (src in all_src) {
    src_rows <- inter[grepl(src, inter$sources, fixed = TRUE), ]
    pairs <- unique(paste(toupper(src_rows$ligand), toupper(src_rows$receptor), sep = "|||"))
    source_pairs[[src]] <- unique(c(source_pairs[[src]], pairs))
  }
}

md("| Source | Unique pairs | Exclusive pairs | Exclusive% |")
md("|--------|-------------|----------------|-----------|")

for (src in sort(names(source_pairs))) {
  my_pairs <- source_pairs[[src]]
  other_pairs <- unique(unlist(source_pairs[names(source_pairs) != src]))
  exclusive <- my_pairs[!my_pairs %in% other_pairs]
  md(sprintf("| %s | %s | %s | %.1f%% |",
             src, format(length(my_pairs), big.mark = ","),
             format(length(exclusive), big.mark = ","),
             100 * length(exclusive) / length(my_pairs)))
}
md("")

# ==============================================================================
# PAIRWISE OVERLAP MATRIX
# ==============================================================================

md("## 13. Pairwise Overlap Between Sources")
md("")
md("Number of shared unique L-R pairs between each pair of sources (diagonal = total unique pairs in source).")
md("")

src_names <- sort(names(source_pairs))

# Header
header_line <- paste(c("| Source", sapply(src_names, function(s) substr(s, 1, 10)), "|"), collapse = " | ")
md(header_line)
sep_line <- paste(c("|---", rep("---", length(src_names)), "|"), collapse = "|")
md(sep_line)

for (s1 in src_names) {
  vals <- c()
  for (s2 in src_names) {
    if (s1 == s2) {
      vals <- c(vals, format(length(source_pairs[[s1]]), big.mark = ","))
    } else {
      overlap <- sum(source_pairs[[s1]] %in% source_pairs[[s2]])
      vals <- c(vals, format(overlap, big.mark = ","))
    }
  }
  md(sprintf("| %s | %s |", s1, paste(vals, collapse = " | ")))
}
md("")
# ==============================================================================
# SAVE MARKDOWN FILE
# ==============================================================================

writeLines(md_lines, OUTPUT_FILE)
cat(sprintf("\nReport saved to: %s\n", OUTPUT_FILE))
cat("Done!\n")
