# library(SingleCellExperiment)
library(data.table)
library(purrr)
library(ggplot2)
library(Seurat)

io <- list()
io$path2atlas <- "/Users/ricard/data/gastrulation10x/processed"
io$path2query <- "/Users/ricard/data/10x_gastrulation_TetChimera/processed"
io$outdir <- "/Users/ricard/data/10x_gastrulation_TetChimera/mapping"

################
## Load atlas ##
################

# Load SingleCellExperiment object
sce_atlas  <- readRDS(paste0(io$path2atlas, "/SingleCellExperiment.rds"))
meta_atlas <- readRDS(paste0(io$path2atlas, "/sample_metadata.rds"))

# Filter
meta_atlas <- meta_atlas[meta_atlas$stage%in%c("E6.75"),]
sce_atlas <- sce_atlas[,meta_atlas$cell] 

# Convert to Seurat object
seurat_atlas   <- CreateSeuratObject(as.matrix(counts(sce_atlas)), project = "ATLAS", min.cells = 25)
seurat_atlas@meta.data$map <- "ATLAS"
seurat_atlas   <- NormalizeData(seurat_atlas)
rownames(meta_atlas) <- meta_atlas$cell
seurat_atlas   <- AddMetaData(object = seurat_atlas, metadata = meta_atlas)
seurat_atlas   <- ScaleData(seurat_atlas)

################
## Load query ##
################

# Load Seurat object
seurat_query  <- readRDS(paste0(io$path2query, "/seurat.rds"))

# Filter
seurat_query <- seurat_query[,seurat_query@meta.data$genotype=="WT"]
meta_query <- seurat_query@meta.data

#############
## Prepare ## 
#############

genes <- intersect(rownames(seurat_query), rownames(seurat_atlas))
seurat_query  <- seurat_query[genes,]
seurat_atlas <- seurat_atlas[genes,]

#########
## Map ##
#########

# Gene selection for input to CCA
seurat_atlas <- FindVariableFeatures(seurat_atlas,selection.method = "vst")
seurat_query <- FindVariableFeatures(seurat_query,selection.method = "vst")

# Mapping
anchors <- FindTransferAnchors(reference = seurat_atlas, query = seurat_query, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = seurat_atlas$celltype, dims = 1:30)
seurat_query <- AddMetaData(object = seurat_query, metadata = predictions)

# table(seurat_query$predicted.id)

##########
## Save ##
##########

# save mapping results as an .rds file
saveRDS(mapping, paste0(io$outdir,"/mapping10x_seurat.rds"))
