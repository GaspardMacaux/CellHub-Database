################################################################################
#          convert_neuronchatdb.R — NeuronChatDB → CellHub Database                    #
#                                                                              #
# Input:  raw_databases/NeuronChat/NeuronChatDB_all.rds                        #
# Output: converted/neuronchatdb_mouse.rds, converted/neuronchatdb_human.rds   #
#                                                                              #
# NeuronChatDB contains neuro-specific interactions in 5 categories:           #
#   - Neuropeptides (Vip_Vipr1) — already proper gene symbols                 #
#   - Synaptic adhesion (Nrxn1_Nlgn1) — already proper gene symbols           #
#   - Gap junctions (Gja1_Gja1) — already proper gene symbols                 #
#   - Neurotransmitters (Glu_Grin1) — ligand is a metabolite name             #
#   - Gas (CO_Gucy1b1) — ligand is a gas molecule name                        #
#                                                                              #
# For neurotransmitter/gas interactions, the ligand is not a gene but a        #
# metabolite (Glu, GABA, 5-HT, etc.). NeuronChat provides a lig_contributor   #
# column listing the biosynthetic enzymes/transporters as gene proxies         #
# (e.g. Gls-Slc17a6-Slc17a7 for glutamate). A cell "produces" the metabolite  #
# if it expresses ANY contributor (OR logic).                                  #
#                                                                              #
# This converter expands each metabolite interaction into N interactions,      #
# one per contributor gene, so CellChat can match them against the expression  #
# matrix. The original metabolite name is preserved in the annotation column.  #
#                                                                              #
# Reference: Zhao et al. (2023) Nat Methods. PMID:36959322                    #
################################################################################

INPUT_DIR   <- "raw_databases/NeuronChat"
OUTPUT_DIR  <- "converted"
SOURCE_NAME <- "NeuronChatDB"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

cat("[NeuronChatDB] Loading...\n")
db <- readRDS(file.path(INPUT_DIR, "NeuronChatDB_all.rds"))

# Map NeuronChat ligand_type to CellChat interaction_type
map_type <- function(ligand_type) {
  switch(ligand_type,
         "Gap junction protein"        = "Cell-Cell Contact",
         "Synaptic adhesion molecule"  = "Cell-Cell Contact",
         "Gas"                         = "Non-protein Signaling",
         "Neurotransmitter"            = "Non-protein Signaling",
         "Neuropeptide"                = "Secreted Signaling",
         NA_character_
  )
}

for (species in c("mouse", "human")) {
  df_name <- sprintf("df_interactionDB_%s", species)
  if (!df_name %in% names(db)) {
    cat(sprintf("[NeuronChatDB] SKIPPED %s — %s not found\n", species, df_name))
    next
  }
  
  df <- db[[df_name]]
  cat(sprintf("[NeuronChatDB] Processing %s (%d raw interactions)...\n", species, nrow(df)))
  
  results <- list()
  
  for (i in seq_len(nrow(df))) {
    ligand_type <- df$ligand_type[i]
    int_type <- map_type(ligand_type)
    receptor <- strsplit(df$interaction_name[i], "_")[[1]]
    # Receptor is everything after the first underscore
    receptor_name <- paste(receptor[-1], collapse = "_")
    
    if (ligand_type %in% c("Neurotransmitter", "Gas")) {
      # Metabolite/gas interaction — expand using lig_contributor
      # lig_contributor format: "Gls-Slc17a6-Slc17a7-Slc17a8"
      metabolite_name <- receptor[1]  # e.g. "Glu", "GABA", "CO", "NO"
      contributors <- strsplit(df$lig_contributor[i], "-")[[1]]
      contributors <- trimws(contributors)
      contributors <- contributors[contributors != ""]
      
      for (contrib in contributors) {
        int_name <- paste(contrib, receptor_name, sep = "_")
        int_name_2 <- paste(contrib, "-", receptor_name)
        
        results[[length(results) + 1]] <- data.frame(
          interaction_name       = int_name,
          interaction_name_2     = int_name_2,
          ligand                 = contrib,
          receptor               = receptor_name,
          species                = species,
          source_db              = SOURCE_NAME,
          pathway_name           = NA_character_,
          interaction_type       = int_type,
          annotation             = paste0(ligand_type, " (", metabolite_name, ")"),
          evidence               = NA_character_,
          neurochat_ligand_type  = ligand_type,
          neurochat_metabolite   = metabolite_name,
          neurochat_original_name = df$interaction_name[i],
          stringsAsFactors       = FALSE
        )
      }
      
    } else {
      # Neuropeptide, synaptic adhesion, gap junction — already proper gene symbols
      ligand_name <- receptor[1]  # First part of interaction_name
      
      int_name <- df$interaction_name[i]
      int_name_2 <- paste(ligand_name, "-", receptor_name)
      
      results[[length(results) + 1]] <- data.frame(
        interaction_name       = int_name,
        interaction_name_2     = int_name_2,
        ligand                 = ligand_name,
        receptor               = receptor_name,
        species                = species,
        source_db              = SOURCE_NAME,
        pathway_name           = NA_character_,
        interaction_type       = int_type,
        annotation             = ligand_type,
        evidence               = NA_character_,
        neurochat_ligand_type  = ligand_type,
        neurochat_metabolite   = NA_character_,
        neurochat_original_name = df$interaction_name[i],
        stringsAsFactors       = FALSE
      )
    }
  }
  
  converted <- do.call(rbind, results)
  
  # Remove duplicates (same contributor-receptor pair can appear if multiple
  # metabolites share a contributor gene)
  n_before <- nrow(converted)
  converted <- converted[!duplicated(converted$interaction_name), ]
  n_after <- nrow(converted)
  
  out_path <- file.path(OUTPUT_DIR, sprintf("neuronchatdb_%s.rds", species))
  saveRDS(converted, out_path)
  
  # Stats
  n_direct <- sum(converted$neurochat_ligand_type %in% c("Neuropeptide", "Synaptic adhesion molecule", "Gap junction protein"))
  n_expanded <- sum(converted$neurochat_ligand_type %in% c("Neurotransmitter", "Gas"))
  
  cat(sprintf("  %s: %d -> %d interactions (-%d dups)\n", species, n_before, n_after, n_before - n_after))
  cat(sprintf("    Direct gene symbols: %d (Neuropeptide=%d, Adhesion=%d, Gap=%d)\n",
              n_direct,
              sum(converted$neurochat_ligand_type == "Neuropeptide"),
              sum(converted$neurochat_ligand_type == "Synaptic adhesion molecule"),
              sum(converted$neurochat_ligand_type == "Gap junction protein")))
  cat(sprintf("    Expanded from metabolites: %d (Neurotransmitter=%d, Gas=%d)\n",
              n_expanded,
              sum(converted$neurochat_ligand_type == "Neurotransmitter"),
              sum(converted$neurochat_ligand_type == "Gas")))
  cat(sprintf("    -> %s\n", out_path))
}

cat("[NeuronChatDB] Done!\n")