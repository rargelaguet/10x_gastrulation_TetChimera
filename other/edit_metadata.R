
source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")

io$output.metadata <- "/Users/ricard/data/10x_gastrulation_TetChimera/sample_metadata.txt.gz"
sample_metadata <- fread(io$metadata)

# Assign stage
# sample_metadata <- fread(io$metadata) %>%
#   .[,stage:=sapply(stringr::str_split(batch,"_"),"[[",1)] %>%
#   .[stage=="E75",stage:="E7.5"] %>% .[stage%in%c("E8","E85"),stage:="E8.5"]

# Assign embryos
opts$batch2embryo <- c(
  "E75_TET_TKO_L002" = "E7.5_embryo1",
  "E75_WT_Host_L001" = "E7.5_embryo1",
  "E85_Rep1_TET_TKO_L004" = "E8.5_embryo1",
  "E85_Rep1_WT_Host_L003" = "E8.5_embryo1",
  "E85_Rep2_TET_TKO_L006" = "E8.5_embryo2",
  "E85_Rep2_WT_Host_L005" = "E8.5_embryo2",
  "E8_5_TET_WT_rep1_SIGAG8" = "E8.5_embryo3",
  "E8_5_TET_WT_rep2_SIGAH8" = "E8.5_embryo4"
  # "SIGAE4_E105_3_TET123_Chimera_Host_L005" = "E10.5_embryo1",
  # "SIGAF4_E105_3_TET123_Chimera_TKO_L006" = "E10.5_embryo1",
  # "SIGAG4_E105_5_TET123_Chimera_Host_L007" = "E10.5_embryo2",
  # "SIGAH4_E105_5_TET123_Chimera_TKO_L008" = "E10.5_embryo2"
)

opts$batch2class <- c(
  "E75_TET_TKO_L002" = "E7.5_TET_TKO",
  "E75_WT_Host_L001" = "E7.5_Host",
  "E85_Rep1_TET_TKO_L004" = "E8.5_TET_TKO",
  "E85_Rep1_WT_Host_L003" = "E8.5_Host",
  "E85_Rep2_TET_TKO_L006" = "E8.5_TET_TKO",
  "E85_Rep2_WT_Host_L005" = "E8.5_Host",
  "E8_5_TET_WT_rep1_SIGAG8" = "E8.5_WT",
  "E8_5_TET_WT_rep2_SIGAH8" = "E8.5_WT"
  # "SIGAE4_E105_3_TET123_Chimera_Host_L005" = "E10.5_embryo1",
  # "SIGAF4_E105_3_TET123_Chimera_TKO_L006" = "E10.5_embryo1",
  # "SIGAG4_E105_5_TET123_Chimera_Host_L007" = "E10.5_embryo2",
  # "SIGAH4_E105_5_TET123_Chimera_TKO_L008" = "E10.5_embryo2"
)

sample_metadata[,embryo:=stringr::str_replace_all(batch,opts$batch2embryo)]
sample_metadata[,class:=stringr::str_replace_all(batch,opts$batch2class)]

sample_metadata[,.N,c("embryo","class","batch")]

fwrite(sample_metadata, io$output.metadata, sep="\t", na="NA", quote=F)

