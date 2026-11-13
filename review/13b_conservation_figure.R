# Step 13b (Reviewer 2, major 5): supplementary figure of cross-species conservation, from the tables of step 13.
# Panel A: curated vs orthology-inferred / predicted pairs per vertebrate species.
# Panel B: containment (% of the pairs of the row species found in the column species), on all pairs and on pairs curated in both species.
# Containment is used instead of Jaccard, which is driven down by the size gap between human / mouse (~10,000 pairs) and the other species (~3,000).
# Usage: Rscript review/13b_conservation_figure.R [review/13_conservation_v2.1]
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) >= 1) args[1] else "review/13_conservation_v2.1"
# Species order: primates, rodents, other mammals, birds, amphibians, fish (species absent from the tables are dropped)
species_order <- c("human", "chimp", "macaque", "marmoset", "mouse", "rat", "cow", "sheep", "pig", "horse", "dog", "chicken", "frog", "zebrafish")
# Species with fewer curated pairs are greyed in the curated heatmap (a percentage of a handful of pairs is not interpretable)
MIN_CURATED <- 20

# Long table of one species x species matrix, with species as ordered factors (rows from top to bottom, columns from left to right)
read_matrix <- function(file, label, order) {
  m <- as.matrix(read.csv(file.path(in_dir, file), row.names = 1, check.names = FALSE))
  sp <- c(intersect(order, rownames(m)), setdiff(rownames(m), order))
  m <- m[sp, sp]
  return(data.frame(row_species = factor(rep(sp, times = length(sp)), levels = rev(sp)), col_species = factor(rep(sp, each = length(sp)), levels = sp),
                    value = as.vector(m), set = label))
}

# Panel A: pairs per species by support
ps <- read.csv(file.path(in_dir, "pairs_per_species.csv"), stringsAsFactors = FALSE)
sp_levels <- c(intersect(species_order, ps$species), setdiff(ps$species, species_order))
pa <- rbind(data.frame(species = ps$species, support = "Curated", pairs = ps$curated_pairs),
            data.frame(species = ps$species, support = "Orthology-inferred or predicted", pairs = ps$inferred_or_predicted_pairs))
pa$species <- factor(pa$species, levels = rev(sp_levels))
pa$support <- factor(pa$support, levels = c("Orthology-inferred or predicted", "Curated"))
lab <- data.frame(species = factor(ps$species, levels = rev(sp_levels)), total = ps$pairs, label = sprintf("%.1f %%", ps$pct_curated))
p_a <- ggplot(pa, aes(x = pairs, y = species, fill = support)) + geom_col(width = 0.75) +
  geom_text(data = lab, aes(x = total, y = species, label = label), inherit.aes = FALSE, hjust = -0.1, size = 3) +
  scale_fill_manual(values = c("Curated" = "#2166ac", "Orthology-inferred or predicted" = "#d1d1d1"), breaks = c("Curated", "Orthology-inferred or predicted")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Ligand-receptor pairs", y = NULL, fill = NULL, title = "A. Curated support per species (label: % curated)") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", panel.grid.major.y = element_blank())

# Panel B: containment heatmaps, all pairs vs curated pairs (species with too few curated pairs greyed in the curated panel)
cur_label <- sprintf("Pairs curated in both species (species with >= %d curated pairs)", MIN_CURATED)
mc <- read_matrix("containment_curated_pairs.csv", cur_label, species_order)
few <- ps$species[ps$curated_pairs < MIN_CURATED]
mc$value[mc$row_species %in% few | mc$col_species %in% few] <- NA
mb <- rbind(read_matrix("containment_all_pairs.csv", "All pairs", species_order), mc)
mb$set <- factor(mb$set, levels = c("All pairs", cur_label))
mb$label <- ifelse(is.na(mb$value), "", sprintf("%.0f", mb$value))
p_b <- ggplot(mb, aes(x = col_species, y = row_species, fill = value)) + geom_tile(colour = "white") +
  geom_text(aes(label = label, colour = value > 60), size = 2.2, show.legend = FALSE) + scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white"), na.value = "black") +
  facet_wrap(~set) + scale_fill_gradient(low = "#f7fbff", high = "#08306b", limits = c(0, 100), na.value = "grey85", name = "% of row\nspecies pairs") +
  labs(x = "Found in species", y = "Pairs of species", title = "B. Pairs of each species found in the other species") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

# Save each panel (PDF for the manuscript, PNG for checking)
save_panel <- function(p, name, w, h) {
  for (ext in c("pdf", "png")) ggsave(file.path(in_dir, paste0(name, ".", ext)), p, width = w, height = h, dpi = 300)
  return(invisible(file.path(in_dir, paste0(name, ".pdf"))))
}
save_panel(p_a, "FigSX_A_curated_per_species", 7, 4.5)
save_panel(p_b, "FigSX_B_containment_species", 12, 6)
message("Figures written to ", in_dir)
