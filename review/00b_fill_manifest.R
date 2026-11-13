# Step 00b: pre-fill source and version columns of the raw source manifest (URLs are added by hand)
manifest_file <- "review/00_sources/sources_manifest.tsv"
man <- read.delim(manifest_file, stringsAsFactors = FALSE, colClasses = "character")
man[is.na(man)] <- ""
# Folder prefix -> source name and version; "package:" marks files extracted from an R package by extract_all_databases.R
lookup <- data.frame(prefix = c("Celegans_LR/", "CellChat/", "cellphonedb-data-5.0.0/", "CellTalkDB/", "ConnectomeDB2025/", "Flyphone/",
                                "Multinichenet/", "NeuronChat/", "PlantCellChatDB/", "Plants-2022-Xu-LR-pairs/", "extract_all_databases.R"),
                     source = c("cell2cell_Armingol", "CellChatDB_v2", "CellPhoneDB", "CellTalkDB", "connectomeDB2025", "FlyPhoneDB2",
                                "MultiNicheNet", "NeuronChatDB", "PlantCellChatDB", "PlantPhoneDB", "extraction_script"),
                     version = c("Armingol et al. 2022", "package:CellChat 2.2.0.9001 (75253cd)", "v5.0.0 data release", "", "2025", "",
                                 "package:nichenetr 2.2.1.1 / multinichenetr 2.1.0", "package:NeuronChat 1.0.0 (36e06b2)", "", "Xu et al. 2022", "pipeline"),
                     stringsAsFactors = FALSE)
for (i in seq_len(nrow(lookup))) {
  hit <- startsWith(man$file, lookup$prefix[i])
  man$source[hit & man$source == ""] <- lookup$source[i]
  man$version[hit & man$version == ""] <- lookup$version[i]
}
write.table(man, manifest_file, sep = "\t", quote = FALSE, row.names = FALSE)
print(man[, c("file", "source", "version")])
message("Still empty version: ", sum(man$version == ""), " | still empty url: ", sum(man$url == ""))