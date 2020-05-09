library(Seurat)
library(batchelor)
library(SingleCellExperiment)
library(data.table)
library(purrr)

################
## Define I/O ##
################

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_TetChimera/mapping/run/mapping_functions.R")
  io$path2atlas <- "/Users/ricard/data/gastrulation10x"
  io$path2query <- "/Users/ricard/data/10x_gastrulation_TetChimera"
  io$outdir <- "/Users/ricard/data/10x_gastrulation_TetChimera/results/mapping"
} else {
  source("/homes/ricard/10x_gastrulation_TetChimera/mapping/run/mapping_functions.R")
  io$path2atlas <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x"
  io$path2query <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_TetChimera"
  io$outdir <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_TetChimera/results/mapping"
}

####################
## Define options ##
####################

opts <- list()
opts$atlas.stages <- c(
  # "E6.5",
  # "E6.75",
  # "E7.0",
  # "E7.25",
  # "E7.5",
  # "E7.75",
  "E8.0",
  "E8.25",
  "E8.5"
  # "mixed_gastrulation"
)

opts$query.batches <- c(
  # "E75_TET_TKO_L002", 
  # "E75_WT_Host_L001", 
  "E85_Rep1_TET_TKO_L004", 
  "E85_Rep2_TET_TKO_L006", 
  "E85_Rep1_WT_Host_L003", 
  "E85_Rep2_WT_Host_L005"
  # "E125_DNMT3A_HET_A_L001", 
  # "E125_DNMT3A_HET_A_L003",
  # "E125_DNMT3A_KO_B_L002", 
  # "E125_DNMT3A_KO_E_L004"
)

################
## Load atlas ##
################

sce_atlas  <- readRDS(paste0(io$path2atlas,"/processed/SingleCellExperiment.rds"))
meta_atlas <- fread(paste0(io$path2atlas,"/sample_metadata.txt.gz")) %>%
  .[stripped==F & doublet==F]

# Filter
meta_atlas <- meta_atlas[!meta_atlas$stage%in%opts$atlas.stages,]
sce_atlas <- sce_atlas[,meta_atlas$cell] 

################
## Load query ##
################

# Load Seurat object
sce_query <- readRDS(paste0(io$path2query,"/processed/second_batch/SingleCellExperiment.rds"))
meta_query <- fread(paste0(io$path2query,"/processed/second_batch/sample_metadata.txt.gz")) %>%
  .[pass_QC==T & batch%in%opts$query.batches & cell%in%colnames(sce_query)]

# table(meta_query$cell%in%colnames(sce_query))
# foo <- meta_query[!cell%in%colnames(sce_query)]

# Filter
sce_query <- sce_query[,meta_query$cell]

#############
## Prepare ## 
#############

genes.intersect <- intersect(rownames(sce_query), rownames(sce_atlas))
sce_query  <- sce_query[genes.intersect,]
sce_atlas <- sce_atlas[genes.intersect,]

# meta_query_list <- list()
# meta_query_list$cell <- meta_query$cell[match(colnames(sce_query), meta_query$cell)]
# meta_query_list$stage <- meta_query$stage[match(colnames(sce_query), meta_query$cell)]

#########
## Map ##
#########

marker_genes.dt <- fread("/Users/ricard/data/gastrulation10x/results/marker_genes/marker_genes.txt.gz")
marker_genes.dt <- marker_genes.dt[,head(.SD,n=50),by="celltype"]
marker_genes <- marker_genes.dt$ens_id
marker_genes <- marker_genes[marker_genes%in%genes.intersect]

mapping  <- mapWrap(
  atlas_sce = sce_atlas, 
  atlas_meta = meta_atlas,
  map_sce = sce_query, 
  map_meta = meta_query, 
  genes = marker_genes,
  order = NULL,
  npcs = 50,
  k = 25,
  return.list = FALSE
)

##########
## Save ##
##########

# save mapping results as an .rds file
saveRDS(mapping, paste0(io$outdir,"/mapping10x_mnn.rds"))
