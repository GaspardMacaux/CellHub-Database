################################################################################
#         upset_plot.R — Source overlap visualization for CellHub Database      #
#                                                                              #
# Generates UpSet plots showing how L-R interactions overlap between source    #
# databases for human and mouse. Output: PNG files for manuscript figures.     #
################################################################################

# install.packages("UpSetR")  # Uncomment if not installed
library(UpSetR)

DATABASES_DIR <- "databases"
OUTPUT_DIR <- "figures"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

for (species in c("human", "mouse")) {
  cat(sprintf("=== Generating UpSet plot for %s ===\n", species))
  
  merged_file <- list.files(file.path(DATABASES_DIR, species),
    pattern = "_merged\\.rds$", full.names = TRUE)
  if (length(merged_file) == 0) {
    cat(sprintf("  No merged file found for %s, skipping\n", species))
    next
  }
  
  db <- readRDS(merged_file[1])
  inter <- db$interaction
  
  if (!"sources" %in% colnames(inter)) {
    cat("  No 'sources' column, skipping\n")
    next
  }
  
  # Build binary membership matrix: each row = interaction, each col = source
  all_sources <- sort(unique(unlist(strsplit(inter$sources, ", "))))
  cat(sprintf("  Sources: %s\n", paste(all_sources, collapse = ", ")))
  cat(sprintf("  Interactions: %d\n", nrow(inter)))
  
  membership <- data.frame(matrix(0, nrow = nrow(inter), ncol = length(all_sources)))
  colnames(membership) <- all_sources
  
  for (i in seq_len(nrow(inter))) {
    srcs <- strsplit(inter$sources[i], ", ")[[1]]
    for (src in srcs) {
      if (src %in% colnames(membership)) {
        membership[i, src] <- 1
      }
    }
  }
  
  # Generate plot
  out_file <- file.path(OUTPUT_DIR, sprintf("upset_%s.png", species))
  
  png(out_file, width = 1200, height = 800, res = 150)
  
  print(upset(membership,
    sets = all_sources,
    order.by = "freq",
    nsets = length(all_sources),
    nintersects = 20,
    mb.ratio = c(0.6, 0.4),
    point.size = 3,
    line.size = 1.2,
    mainbar.y.label = "Number of interactions",
    sets.x.label = "Interactions per source",
    text.scale = c(1.3, 1.1, 1.1, 1.0, 1.3, 1.0),
    main.bar.color = "#2F5496",
    sets.bar.color = "#C55A11",
    matrix.color = "#7030A0",
    shade.color = "#E2D1F0",
    show.numbers = "yes"
  ))
  
  dev.off()
  cat(sprintf("  Saved: %s\n", out_file))
}

cat("\nDone!\n")
library(ggplot2)

DATABASES_DIR <- "databases"

species_data <- data.frame(species = character(), n = integer(), 
                           category = character(), stringsAsFactors = FALSE)

mammals <- c("human","mouse","rat","chimp","macaque","marmoset",
             "cow","dog","horse","pig","sheep")
vertebrates <- c("chicken","frog","zebrafish")
invertebrates <- c("drosophila","celegans")
plants <- c("arabidopsis","maize","rice","soybean","tomato")

species_dirs <- list.dirs(DATABASES_DIR, recursive = FALSE)
for (sp_dir in species_dirs) {
  sp <- basename(sp_dir)
  mf <- list.files(sp_dir, pattern = "_merged\\.rds$", full.names = TRUE)
  if (length(mf) == 0) next
  db <- readRDS(mf[1])
  
  cat_label <- if (sp %in% mammals) "Mammals" else
    if (sp %in% vertebrates) "Other vertebrates" else
      if (sp %in% invertebrates) "Invertebrates" else "Plants"
  
  species_data <- rbind(species_data, data.frame(
    species = sp, n = nrow(db$interaction), category = cat_label,
    stringsAsFactors = FALSE))
}

cat(sprintf("Loaded %d species\n", nrow(species_data)))
print(species_data)

species_data$category <- factor(species_data$category, 
                                levels = c("Mammals","Other vertebrates","Invertebrates","Plants"))
species_data <- species_data[order(species_data$category, -species_data$n), ]
species_data$species <- factor(species_data$species, levels = species_data$species)

p <- ggplot(species_data, aes(x = species, y = n, fill = category)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c(
    "Mammals" = "#2F5496",
    "Other vertebrates" = "#548235",
    "Invertebrates" = "#C55A11",
    "Plants" = "#7030A0")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Number of interactions", fill = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "top",
    panel.grid.major.y = element_line(color = "gray90"))

ggsave("figures/barplot_species.png", p, width = 10, height = 5, dpi = 300, bg = "white")
cat("Saved!\n")