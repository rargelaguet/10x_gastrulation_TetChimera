io <- list()

io$min.counts.per.gene <- 50
io$k <- 30
io$min.mean <- 1e-3
io$npcs <- 50

io$atlas.RDS <- NULL # Choose between NULL and path to atlas SCE RDS
io$atlas.metadata <- NULL # Choose between NULL and path to atlas metadata txt
io$corrected.atlas.RDS <- NULL # Choose between NULL and path to atlas SCE RDS
io$corrected.atlas.metadata <- NULL  # Choose between NULL and path to atlas metadata txt
io$atlas_stages <- NULL # Choose between NULL and atlas stages
io$testing <- NULL # Choose between NULL and TRUE

io$subset_celltypes <- list(
    "Haematoendothelial" = c("Haematoendothelial progenitors", "Blood progenitors 1", "Blood progenitors 2", "Erythroid1", "Erythroid2", "Erythroid3", "Endothelium", "Cardiomyocytes"),
    "Endoderm" = c("ExE endoderm", "Visceral endoderm", "Gut", "Def. endoderm", "Notochord"),
    "EpiblastPSNeuro" = c("Epiblast", "Primitive Streak", "Rostral neurectoderm", "Spinal cord", "Surface ectoderm", "Nascent mesoderm", "NMP", "Neural crest", "Caudal neurectoderm", "Caudal epiblast", "Anterior Primitive Streak", "Forebrain/Midbrain/Hindbrain", "PGC"),
    "Mesoderm" = c("Nascent mesoderm", "Mesenchyme", "Mixed mesoderm", "ExE mesoderm, Intermediate mesoderm", "Pharyngeal mesoderm", "Paraxial mesoderm", "Somitic mesoderm", "Caudal Mesoderm", "Allantois")
    # Leaving out Parietal Endoderm and ExE ectoderm because they should already be well identified.
)

io$order <- list()
io$order[["WT"]] <- c("ATLAS", "SIGAG3_E8.5_hashing_Host-WT_L007", "SIGAC3_E8.5_pool2_Host-WT_L003", "SIGAA3_E8.5_pool1_Host-WT_L001", "SIGAE3_E7.5_pool1_Host-WT_L005")
io$order[["TKO"]] <- c("ATLAS", "SIGAH3_E8.5_hasting_TET-TKO_L008", "SIGAB3_E8.5_pool1_TET-TKO_L002", "SIGAD3_E8.5_pool2_TET-TKO_L004", "SIGAF3_E7.5_pool1_TET-TKO_L006")