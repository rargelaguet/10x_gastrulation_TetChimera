#########
## I/O ##
#########

source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")

io$mapping.dir <- paste0(io$basedir,"/results/mapping")
# io$metadata <- paste0(io$basedir,"/processed/sample_metadata.txt.gz")

opts$batches <- c(
  # Second batch
  "E75_TET_TKO_L002",
  "E75_WT_Host_L001",
  "E85_Rep1_TET_TKO_L004",
  "E85_Rep1_WT_Host_L003",
  "E85_Rep2_TET_TKO_L006",
  "E85_Rep2_WT_Host_L005",
  
  # Fifth batch
  "E8_5_TET_WT_rep1_SIGAG8",
  "E8_5_TET_WT_rep2_SIGAH8"
)

###############
## Load data ##
###############

# Load metadata
metadata <- fread(io$metadata) %>% 
  .[,c("celltype.mapped","celltype.score","stage.mapped","closest.cell"):=NULL]

# Load mapping results
mapping.dt <- opts$batches %>% map(function(x) {
  file <- sprintf("%s/mapping_mnn_%s.rds",io$mapping.dir,x)
  if (file.exists(file)) {
    readRDS(file)$mapping %>% as.data.table %>%
      .[,c("cell","celltype.mapped","celltype.score","stage.mapped","closest.cell")] %>%
      .[,batch:=x]
  }
}
) %>% rbindlist# %>% .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,"_"," ")]

table(mapping.dt$celltype.mapped)
table(mapping.dt$stage.mapped)
unique(mapping.dt$batch)

#########################################
## Quick comparisons to other mappings ##
#########################################

# mapping.old <- fread("/Users/ricard/data/10x_gastrulation_TetChimera/backups/sample_metadata_6No.txt.gz") %>%
#   .[,c("cell","batch","celltype.mapped","celltype.score")]
# foo <- merge(mapping.dt,mapping.old, by=c("cell","batch"))
# foo[celltype.mapped.x==celltype.mapped.y]

###########
## Merge ##
###########

metadata <- metadata %>% 
  merge(mapping.dt,by=c("cell","batch"),all.x=TRUE)

head(metadata)

#################
## Save output ##
#################

fwrite(metadata, io$metadata, sep="\t", na="NA", quote=F)


