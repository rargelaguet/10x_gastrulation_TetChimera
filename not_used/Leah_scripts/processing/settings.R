io <- list()

io$samples <- list()

io$samples$second_batch <- c(
  "E75_TET_TKO_L002",
  "E75_WT_Host_L001",
  "E85_Rep1_TET_TKO_L004",
  "E85_Rep2_TET_TKO_L006",
  "E85_Rep1_WT_Host_L003",
  "E85_Rep2_WT_Host_L005"
)

io$samples$third_batch <- c(
  "SIGAE4_E105_3_TET123_Chimera_Host_L005",
  "SIGAF4_E105_3_TET123_Chimera_TKO_L006",
  "SIGAG4_E105_5_TET123_Chimera_Host_L007",
  "SIGAH4_E105_5_TET123_Chimera_TKO_L008"
)

io$subset.proteincoding <- NULL # Choose between NULL and path to proteincoding genes
io$qc <- TRUE # Choose between NULL and TRUE

io$min_nFeature_RNA <- list()
io$min_nCount_RNA <- list()
io$max_percent.mt <- list()
io$max_doublet_score <- list()
io$min_nFeature_RNA$second_batch <- c(
  "E75_TET_TKO_L002" = 3300,
  "E75_WT_Host_L001" = 2000,
  "E85_Rep1_TET_TKO_L004" = 2500,
  "E85_Rep2_TET_TKO_L006" = 2000,
  "E85_Rep1_WT_Host_L003" = 2000,
  "E85_Rep2_WT_Host_L005" = 2000
)
io$min_nCount_RNA$second_batch <- c(
  "E75_TET_TKO_L002" = 14000,
  "E75_WT_Host_L001" = 5000,
  "E85_Rep1_TET_TKO_L004" = 6000,
  "E85_Rep2_TET_TKO_L006" = 5000,
  "E85_Rep1_WT_Host_L003" = 5000,
  "E85_Rep2_WT_Host_L005" = 5000
)
io$max_percent.mt$second_batch <- c(
  "E75_TET_TKO_L002" = 10,
  "E75_WT_Host_L001" = 10,
  "E85_Rep1_TET_TKO_L004" = 10,
  "E85_Rep2_TET_TKO_L006" = 10,
  "E85_Rep1_WT_Host_L003" = 10,
  "E85_Rep2_WT_Host_L005" = 10
)
io$max_doublet_score$second_batch <- c(
  "E75_TET_TKO_L002" = 1000,
  "E75_WT_Host_L001" = 1000,
  "E85_Rep1_TET_TKO_L004" = 1000,
  "E85_Rep2_TET_TKO_L006" = 1000,
  "E85_Rep1_WT_Host_L003" = 1000,
  "E85_Rep2_WT_Host_L005" = 1000
)

io$min_nFeature_RNA$third_batch <- c(
  "SIGAE4_E105_3_TET123_Chimera_Host_L005" = 3500,
  "SIGAF4_E105_3_TET123_Chimera_TKO_L006" = 2000,
  "SIGAG4_E105_5_TET123_Chimera_Host_L007" = 3500,
  "SIGAH4_E105_5_TET123_Chimera_TKO_L008" = 2000
)
io$min_nCount_RNA$third_batch <- c(
  "SIGAE4_E105_3_TET123_Chimera_Host_L005" = 10000,
  "SIGAF4_E105_3_TET123_Chimera_TKO_L006" = 8000,
  "SIGAG4_E105_5_TET123_Chimera_Host_L007" = 8000,
  "SIGAH4_E105_5_TET123_Chimera_TKO_L008" = 8000
)
io$max_percent.mt$third_batch <- c(
  "SIGAE4_E105_3_TET123_Chimera_Host_L005" = 10,
  "SIGAF4_E105_3_TET123_Chimera_TKO_L006" = 10,
  "SIGAG4_E105_5_TET123_Chimera_Host_L007" = 10,
  "SIGAH4_E105_5_TET123_Chimera_TKO_L008" = 10
)
io$max_doublet_score$third_batch <- c(
  "SIGAE4_E105_3_TET123_Chimera_Host_L005" = 1000,
  "SIGAF4_E105_3_TET123_Chimera_TKO_L006" = 1000,
  "SIGAG4_E105_5_TET123_Chimera_Host_L007" = 1000,
  "SIGAH4_E105_5_TET123_Chimera_TKO_L008" = 1000
)