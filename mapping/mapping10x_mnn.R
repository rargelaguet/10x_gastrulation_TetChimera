library(Seurat)
library(SingleCellExperiment)
library(data.table)
library(purrr)


io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_TetChimera/mapping/mapping_settings.R")
  io$path2atlas <- "/Users/ricard/data/gastrulation10x/processed"
  io$path2query <- "/Users/ricard/data/10x_gastrulation_TetChimera/processed"
  io$outdir <- "/Users/ricard/data/10x_gastrulation_TetChimera/mapping"
} else {
  source("/homes/ricard/10x_gastrulation_TetChimera/mapping/mapping_settings.R")
  io$path2atlas <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/processed"
  io$path2query <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_TetChimera/processed"
  io$outdir <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_TetChimera/mapping"
}

################
## Load atlas ##
################

sce_atlas  <- readRDS(paste0(io$path2atlas,"/SingleCellExperiment.rds"))
meta_atlas <- readRDS(paste0(io$path2atlas,"/sample_metadata.rds"))

# Filter
meta_atlas <- meta_atlas[!meta_atlas$stage%in%c("E6.5","E6.75","E7.0","E7.25"),]
meta_atlas <- meta_atlas[!meta_atlas$celltype%in%c("PGC"),]
sce_atlas <- sce_atlas[,meta_atlas$cell] 

################
## Load query ##
################

# Load Seurat object
seurat_query  <- readRDS(paste0(io$path2query, "/seurat.rds"))
meta_query <- seurat_query@meta.data

# Filter
meta_query <- meta_query[meta_query$stage=="E8.5",]
seurat_query <- seurat_query[,meta_query$cells]

# Convert from Seurat to SCE
sce_query <- Seurat::as.SingleCellExperiment(seurat_query)

#############
## Prepare ## 
#############

genes <- intersect(rownames(sce_query), rownames(sce_atlas))
sce_query  <- sce_query[genes,]
sce_atlas <- sce_atlas[genes,]

meta_query_list <- list()
meta_query_list$cell <- meta_query$cell[match(colnames(sce_query), meta_query$cell)]
meta_query_list$stage <- meta_query$stage[match(colnames(sce_query), meta_query$cell)]

#########
## Map ##
#########

# TO-DO: BATCH EFFECT CORRECTION IN THE QUERY?

mapping  <- mapWrap(
  atlas_sce = sce_atlas, atlas_meta = meta_atlas,
  map_sce = sce_query, map_meta = meta_query_list, 
  k = 25
)

##########
## Save ##
##########

# save mapping results as an .rds file
saveRDS(mapping, paste0(io$outdir,"/mapping10x_mnn_v2.rds"))
