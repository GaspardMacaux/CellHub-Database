################################################################################
#              GaspouDB — Extraction de toutes les sources brutes              #
#                                                                              #
# Ce script extrait les données LR de chaque package R et les sauvegarde       #
# en RDS bruts dans data/raw/. Les sources web (CSV/Excel) doivent être   #
# téléchargées manuellement et placées dans data/raw/ avant de lancer     #
# les scripts de conversion.                                                   #
#                                                                              #
# Prérequis packages R (installer AVANT de lancer ce script):                  #
#   pak::pak("jinworks/CellChat")              # CellChatDB v1 + v2           #
#   pak::pak("saeyslab/multinichenetr")         # MultiNicheNet                #
#   pak::pak("Wei-BioMath/NeuronChat")          # NeuronChatDB                #
#   BiocManager::install("SingleCellSignalR")   # LRdb                        #
#   pak::pak("soumelis-lab/ICELLNET")           # ICELLNET                     #
#   pak::pak("FullStackGoogler/FlyPhoneDB2")    # FlyPhoneDB2 (pkg = FlyPhone)#
#                                                                              #
# Prérequis téléchargements manuels:                                           #
#   - connectomeDB2025 : https://connectomedb.org -> all_species.csv           #
#   - CellPhoneDB v5   : https://github.com/ventolab/cellphonedb-data         #
#   - CellTalkDB       : https://github.com/ZJUFanLab/CellTalkDB              #
#   - MRCLinkdb        : http://ccc.mrclink.com -> download                   #
#   - cell2cell C.elegans : wget du Excel (voir section 7)                     #
#   - PlantPhoneDB     : https://github.com/plant-phonedb/PlantPhoneDB        #
################################################################################

OUTPUT_DIR <- "data/raw"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

cat("============================================================\n")
cat("  GaspouDB — Extraction des sources brutes\n")
cat("============================================================\n\n")

# ==============================================================================
# SOURCE 1: CellChatDB (v1 + v2) — package CellChat
# ==============================================================================
# Package: CellChat
# Install: pak::pak("jinworks/CellChat")
# Contains: mouse, human, zebrafish (v2 only for zebrafish)
# Format: CellChat-native list (interaction, complex, cofactor, geneInfo)

cat("=== [1/6] CellChatDB ===\n")
if (requireNamespace("CellChat", quietly = TRUE)) {

  # Mouse
  db <- CellChat::CellChatDB.mouse
  saveRDS(db, file.path(OUTPUT_DIR, "CellChatDB_v2_mouse.rds"))
  cat(sprintf("  ✓ CellChatDB v2 mouse: %d interactions\n", nrow(db$interaction)))

  # Human
  db <- CellChat::CellChatDB.human
  saveRDS(db, file.path(OUTPUT_DIR, "CellChatDB_v2_human.rds"))
  cat(sprintf("  ✓ CellChatDB v2 human: %d interactions\n", nrow(db$interaction)))

  # Zebrafish (v2 only)
  if (exists("CellChatDB.zebrafish", where = asNamespace("CellChat"))) {
    db <- CellChat::CellChatDB.zebrafish
    saveRDS(db, file.path(OUTPUT_DIR, "CellChatDB_v2_zebrafish.rds"))
    cat(sprintf("  ✓ CellChatDB v2 zebrafish: %d interactions\n", nrow(db$interaction)))
  } else {
    cat("  ⚠ CellChatDB zebrafish not found in this CellChat version\n")
  }

} else {
  cat("  ✗ SKIPPED — CellChat not installed\n")
  cat("    Install: pak::pak('jinworks/CellChat')\n")
}

# ==============================================================================
# SOURCE 2: MultiNicheNet — package multinichenetr
# ==============================================================================
# Package: multinichenetr
# Install: pak::pak("saeyslab/multinichenetr")
# Contains: human + mouse LR network (from NicheNet + extensions)

cat("\n=== [2/6] MultiNicheNet ===\n")
if (requireNamespace("multinichenetr", quietly = TRUE)) {

  lr_network <- multinichenetr::lr_network
  saveRDS(lr_network, file.path(OUTPUT_DIR, "MultiNicheNet_lr_network.rds"))
  cat(sprintf("  ✓ MultiNicheNet LR network: %d rows\n", nrow(lr_network)))
  cat(sprintf("    Columns: %s\n", paste(colnames(lr_network), collapse = ", ")))

} else {
  cat("  ✗ SKIPPED — multinichenetr not installed\n")
  cat("    Install: pak::pak('saeyslab/multinichenetr')\n")
}

# ==============================================================================
# SOURCE 3: NeuronChatDB — package NeuronChat
# ==============================================================================
# Package: NeuronChat
# Install: pak::pak("Wei-BioMath/NeuronChat")
# Contains: ~220 neuro-specific interactions (neurotransmitters, neuropeptides)

cat("\n=== [3/6] NeuronChatDB ===\n")
if (requireNamespace("NeuronChat", quietly = TRUE)) {

  # NeuronChat stores its DB as an internal data object
  if (exists("NeuronChat_DB", where = asNamespace("NeuronChat"))) {
    db <- NeuronChat::NeuronChat_DB
    saveRDS(db, file.path(OUTPUT_DIR, "NeuronChatDB.rds"))
    n <- if (is.data.frame(db)) nrow(db) else if (is.list(db)) length(db) else NA
    cat(sprintf("  ✓ NeuronChatDB: %s (class: %s)\n", n, paste(class(db), collapse = ", ")))
  } else {
    # Try loading via data()
    tryCatch({
      e <- new.env()
      data("NeuronChat_DB", package = "NeuronChat", envir = e)
      db <- e$NeuronChat_DB
      saveRDS(db, file.path(OUTPUT_DIR, "NeuronChatDB.rds"))
      cat(sprintf("  ✓ NeuronChatDB loaded via data()\n"))
    }, error = function(err) {
      cat(sprintf("  ⚠ Could not find DB object. Available objects:\n"))
      cat(sprintf("    %s\n", paste(ls(asNamespace("NeuronChat")), collapse = ", ")))
      cat("    Check manually and extract.\n")
    })
  }

} else {
  cat("  ✗ SKIPPED — NeuronChat not installed\n")
  cat("    Install: pak::pak('Wei-BioMath/NeuronChat')\n")
}

# ==============================================================================
# SOURCE 4: SingleCellSignalR / LRdb
# ==============================================================================
# Package: SingleCellSignalR
# Install: BiocManager::install("SingleCellSignalR")
# Contains: 3191 human LR pairs (compiled from FANTOM5/HPRD/IUPHAR)

cat("\n=== [4/6] SingleCellSignalR / LRdb ===\n")
if (requireNamespace("SingleCellSignalR", quietly = TRUE)) {

  db <- SingleCellSignalR::LRdb
  saveRDS(db, file.path(OUTPUT_DIR, "SingleCellSignalR_LRdb.rds"))
  cat(sprintf("  ✓ LRdb: %d interactions\n", nrow(db)))
  cat(sprintf("    Columns: %s\n", paste(colnames(db), collapse = ", ")))

} else {
  cat("  ✗ SKIPPED — SingleCellSignalR not installed\n")
  cat("    Install: BiocManager::install('SingleCellSignalR')\n")
}

# ==============================================================================
# SOURCE 5: ICELLNET
# ==============================================================================
# Package: icellnet
# Install: pak::pak("soumelis-lab/ICELLNET")
# Contains: ~380 curated human LR pairs, enriched in immune interactions

cat("\n=== [5/6] ICELLNET ===\n")
if (requireNamespace("icellnet", quietly = TRUE)) {

  # ICELLNET stores its DB as internal data — name varies across versions
  ns_objects <- ls(asNamespace("icellnet"))
  db_candidates <- intersect(ns_objects, c("db", "LR.DB", "LR_DB", "icellnet_db", "DB"))

  if (length(db_candidates) > 0) {
    db <- get(db_candidates[1], envir = asNamespace("icellnet"))
    saveRDS(db, file.path(OUTPUT_DIR, "ICELLNET_db.rds"))
    n <- if (is.data.frame(db)) nrow(db) else length(db)
    cat(sprintf("  ✓ ICELLNET DB (%s): %d entries\n", db_candidates[1], n))
  } else {
    # Try loading via data()
    tryCatch({
      e <- new.env()
      # Try common names
      for (nm in c("icellnet_db", "db", "LR.DB")) {
        tryCatch({
          data(list = nm, package = "icellnet", envir = e)
          db <- e[[nm]]
          saveRDS(db, file.path(OUTPUT_DIR, "ICELLNET_db.rds"))
          cat(sprintf("  ✓ ICELLNET DB loaded via data('%s')\n", nm))
          break
        }, error = function(err) NULL)
      }
      if (!file.exists(file.path(OUTPUT_DIR, "ICELLNET_db.rds"))) {
        cat("  ⚠ Could not find DB object. Available:\n")
        cat(sprintf("    Namespace: %s\n", paste(head(ns_objects, 20), collapse = ", ")))
        # Also check data()
        avail_data <- data(package = "icellnet")$results[, "Item"]
        cat(sprintf("    data(): %s\n", paste(avail_data, collapse = ", ")))
      }
    }, error = function(err) {
      cat(sprintf("  ⚠ Error: %s\n", conditionMessage(err)))
    })
  }

} else {
  cat("  ✗ SKIPPED — icellnet not installed\n")
  cat("    Install: pak::pak('soumelis-lab/ICELLNET')\n")
}

# ==============================================================================
# SOURCE 6: FlyPhoneDB2 — package FlyPhone
# ==============================================================================
# Package name: FlyPhone (NOT FlyPhoneDB2)
# Install: pak::pak("FullStackGoogler/FlyPhoneDB2")
# Contains: Drosophila melanogaster LR pairs in 3 confidence tiers:
#   V2A  = 1804 pairs (All: high + moderate + low confidence)
#   V2HM = 967 pairs  (High + Moderate confidence only)
#   V2H  = 316 pairs  (High confidence only)
#   V1   = legacy version (2021)
# Gene symbols: FlyBase convention (Burs, sNPF, dpp, Toll)
# Columns: Pair_ID, Gene_secreted, Gene_receptor, Pathway_receptor,
#   rank (confidence level), FBgn IDs, mammalian orthologs, annotations

cat("\n=== [6/6] FlyPhoneDB2 (package: FlyPhone) ===\n")
if (requireNamespace("FlyPhone", quietly = TRUE)) {

  # Load all 4 dataset versions
  e <- new.env()
  data("V1", package = "FlyPhone", envir = e)
  data("V2A", package = "FlyPhone", envir = e)
  data("V2H", package = "FlyPhone", envir = e)
  data("V2HM", package = "FlyPhone", envir = e)

  # Also load pathway components (useful for annotation)
  tryCatch({
    data("pathway_components_v2", package = "FlyPhone", envir = e)
  }, error = function(err) NULL)

  # Save all versions in a single list
  fly_db <- list(
    V2A  = e$V2A,   # All confidence levels (1804 pairs)
    V2HM = e$V2HM,  # High + Moderate (967 pairs)
    V2H  = e$V2H,   # High only (316 pairs)
    V1   = e$V1     # Legacy (2021)
  )
  if (exists("pathway_components_v2", envir = e)) {
    fly_db$pathway_components <- e$pathway_components_v2
  }

  saveRDS(fly_db, file.path(OUTPUT_DIR, "FlyPhoneDB2_all.rds"))

  cat(sprintf("  ✓ V2A (all):            %d interactions\n", nrow(e$V2A)))
  cat(sprintf("  ✓ V2HM (high+moderate): %d interactions\n", nrow(e$V2HM)))
  cat(sprintf("  ✓ V2H (high only):      %d interactions\n", nrow(e$V2H)))
  cat(sprintf("  ✓ V1 (legacy 2021):     %d interactions\n", nrow(e$V1)))
  cat(sprintf("  Columns: %s\n", paste(colnames(e$V2A), collapse = ", ")))
  cat("  Saved to: data/raw/FlyPhoneDB2_all.rds\n")

} else {
  cat("  ✗ SKIPPED — FlyPhone not installed\n")
  cat("    Install: pak::pak('FullStackGoogler/FlyPhoneDB2')\n")
  cat("    Note: the R package name is 'FlyPhone', not 'FlyPhoneDB2'\n")
}

# ==============================================================================
# SOURCES NÉCESSITANT UN TÉLÉCHARGEMENT MANUEL
# ==============================================================================

cat("\n============================================================\n")
cat("  SOURCES À TÉLÉCHARGER MANUELLEMENT\n")
cat("============================================================\n")
cat("  Placer les fichiers dans:", OUTPUT_DIR, "/\n\n")

cat("  1. connectomeDB2025 (PRIORITÉ HAUTE — 14 espèces)\n")
cat("     URL: https://connectomedb.org\n")
cat("     Télécharger: all_species.csv\n")
cat("     -> data/raw/all_species.csv\n\n")

cat("  2. CellPhoneDB v5 (2911 interactions humain)\n")
cat("     URL: https://github.com/ventolab/cellphonedb-data\n")
cat("     Télécharger: data_interaction_input.csv, data_gene_input.csv,\n")
cat("                  data_complex_input.csv, data_protein_input.csv\n")
cat("     -> data/raw/cellphonedb_*.csv\n\n")

cat("  3. CellTalkDB (humain + souris)\n")
cat("     wget https://raw.githubusercontent.com/ZJUFanLab/CellTalkDB/master/database/human_lr_pair.txt \\\n")
cat("       -O data/raw/CellTalkDB_human.txt\n")
cat("     wget https://raw.githubusercontent.com/ZJUFanLab/CellTalkDB/master/database/mouse_lr_pair.txt \\\n")
cat("       -O data/raw/CellTalkDB_mouse.txt\n\n")

cat("  4. MRCLinkdb (métabolite-récepteur, humain + souris)\n")
cat("     URL: http://ccc.mrclink.com -> section download\n")
cat("     -> data/raw/MRCLinkdb_human.csv\n")
cat("     -> data/raw/MRCLinkdb_mouse.csv\n\n")

cat("  5. C. elegans (cell2cell / Armingol et al.)\n")
cat("     wget https://github.com/LewisLabUCSD/Celegans-cell2cell/raw/master/data/PPI-Networks/Celegans-Curated-LR-pairs.xlsx \\\n")
cat("       -O data/raw/Celegans_Curated_LR_pairs.xlsx\n\n")

cat("  6. PlantPhoneDB (plantes / Arabidopsis)\n")
cat("     URL: https://github.com/plant-phonedb/PlantPhoneDB\n")
cat("     -> data/raw/PlantPhoneDB_LR.csv\n\n")



# ==============================================================================
# SOURCE 7: PlantCellChat — package PlantCellChat
# ==============================================================================
# Package: PlantCellChat
# Install: pak::pak("mrliuw/PlantCellChat")
# Contains: LR interactions for 5 plant species:
#   Arabidopsis thaliana, Oryza sativa (rice), Solanum lycopersicum (tomato),
#   Zea mays (maize), Glycine max (soybean)
# PlantCellChatDB: 12,712 LR pairs across 5 species
# Includes hormone-receptor interactions and receptor-like kinases
# Reference: Liu et al. (2026) Plant J. PMID:42119054

cat("\n=== [7/7] PlantCellChat (PlantCellChatDB) ===\n")
if (requireNamespace("PlantCellChat", quietly = TRUE)) {
  
  # List available data objects
  avail_data <- data(package = "PlantCellChat")$results[, "Item"]
  cat(sprintf("  Available data objects: %s\n", paste(avail_data, collapse = ", ")))
  
  # Try to load all DB-related objects
  e <- new.env()
  loaded <- character(0)
  for (nm in avail_data) {
    tryCatch({
      data(list = nm, package = "PlantCellChat", envir = e)
      loaded <- c(loaded, nm)
    }, error = function(err) NULL)
  }
  cat(sprintf("  Successfully loaded: %s\n", paste(loaded, collapse = ", ")))
  
  # Inspect each loaded object
  for (nm in loaded) {
    obj <- e[[nm]]
    if (is.data.frame(obj)) {
      cat(sprintf("  %s: %d rows x %d cols\n", nm, nrow(obj), ncol(obj)))
      cat(sprintf("    Columns: %s\n", paste(colnames(obj), collapse = ", ")))
    } else if (is.list(obj)) {
      cat(sprintf("  %s: list with %d elements (%s)\n", nm, length(obj),
                  paste(head(names(obj), 10), collapse = ", ")))
      # If it's a CellChat-style DB list, inspect interaction table
      if ("interaction" %in% names(obj)) {
        cat(sprintf("    interaction: %d rows\n", nrow(obj$interaction)))
        cat(sprintf("    Columns: %s\n", paste(colnames(obj$interaction), collapse = ", ")))
      }
    } else {
      cat(sprintf("  %s: class=%s\n", nm, paste(class(obj), collapse = ", ")))
    }
  }
  
  # Also check namespace for unexported objects
  ns_objects <- grep("DB\\.?|db\\.?|LR|lr|interact|database",
                     ls(asNamespace("PlantCellChat")),
                     value = TRUE, ignore.case = TRUE)
  if (length(ns_objects) > 0) {
    cat(sprintf("  Namespace DB-related objects: %s\n", paste(ns_objects, collapse = ", ")))
  }
  
  # Save everything in one list
  plantcellchat_db <- list()
  for (nm in loaded) {
    plantcellchat_db[[nm]] <- e[[nm]]
  }
  saveRDS(plantcellchat_db, file.path(OUTPUT_DIR, "PlantCellChatDB_all.rds"))
  cat(sprintf("  ✓ Saved to: %s/PlantCellChatDB_all.rds\n", OUTPUT_DIR))
  
} else {
  cat("  ✗ SKIPPED — PlantCellChat not installed\n")
  cat("    Install: pak::pak('mrliuw/PlantCellChat')\n")
  cat("    GitHub: https://github.com/mrliuw/PlantCellChat\n")
  cat("    Web DB: https://bioinformatics.fafu.edu.cn/PlantCellChat/\n")
}
# ==============================================================================
# RÉSUMÉ
# ==============================================================================

cat("============================================================\n")
cat("  RÉSUMÉ\n")
cat("============================================================\n")

# Check extracted R package files
extracted <- list.files(OUTPUT_DIR, pattern = "\\.rds$")
cat(sprintf("\nFichiers extraits depuis packages R: %d\n", length(extracted)))
for (f in extracted) {
  fsize <- file.info(file.path(OUTPUT_DIR, f))$size
  cat(sprintf("  ✓ %-45s (%5.0f KB)\n", f, fsize / 1024))
}

# Check manual downloads
manual_files <- c(
  "all_species.csv",
  "cellphonedb_interaction.csv" = "data_interaction_input.csv",
  "CellTalkDB_human.txt",
  "CellTalkDB_mouse.txt",
  "MRCLinkdb_human.csv",
  "MRCLinkdb_mouse.csv",
  "Celegans_Curated_LR_pairs.xlsx",
  "PlantPhoneDB_LR.csv"
)

present <- manual_files[file.exists(file.path(OUTPUT_DIR, manual_files))]
absent <- manual_files[!file.exists(file.path(OUTPUT_DIR, manual_files))]

if (length(present) > 0) {
  cat(sprintf("\nTéléchargements manuels présents: %d\n", length(present)))
  for (f in present) cat(sprintf("  ✓ %s\n", f))
}
if (length(absent) > 0) {
  cat(sprintf("\nTéléchargements manuels manquants: %d\n", length(absent)))
  for (f in absent) cat(sprintf("  ✗ %s\n", f))
}

cat("\n")
