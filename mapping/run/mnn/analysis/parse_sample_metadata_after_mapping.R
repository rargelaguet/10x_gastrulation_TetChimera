#########
## I/O ##
#########

source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")

io$mapping.dir <- paste0(io$basedir,"/results/fourth_batch/mapping")
io$metadata <- paste0(io$basedir,"/processed/fourth_batch/sample_metadata.txt.gz")

opts$batches <- c(
  "SIGAC2_TET_TKO_E9_5_Head1_L002",
  "SIGAD2_TET_TKO_E9_5_Trunk1_L002",
  "SIGAE2_TET_TKO_E9_5_Tail1_L002",
  "SIGAE6_TET_TKO_E9_5_Head2_L003",
  "SIGAF2_TET_TKO_E9_5_YS1_L002",
  "SIGAF6_TET_TKO_E9_5_Trunk2_L003",
  "SIGAG6_TET_TKO_E9_5_Tail2_L003",
  "SIGAH6_TET_TKO_E9_5_YS2_L003"
)

###############
## Load data ##
###############

# Load metadata
metadata <- fread(io$metadata) %>% .[,c("celltype.mapped","celltype.score"):=NULL]

# Load mapping results
mapping.dt <- opts$batches %>% map(function(x) 
  readRDS(sprintf("%s/mapping_mnn_%s.rds",io$mapping.dir,x))$mapping %>% .[,c("cell","celltype.mapped","celltype.score","closest.cell")] %>% as.data.table
) %>% rbindlist# %>% .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,"_"," ")]

unique(mapping.dt$celltype.mapped)

###########
## Merge ##
###########

metadata <- metadata %>% merge(mapping.dt,by="cell",all.x=TRUE)

#################
## Save output ##
#################

fwrite(metadata, io$metadata, sep="\t", na="NA", quote=F)


