library(Seurat)

samples <- c(
  "SIGAA3_E8.5_pool1_Host-WT_L001", 
  # "SIGAB3_E8.5_pool1_TET-TKO_L002", 
  "SIGAC3_E8.5_pool2_Host-WT_L003", 
  # "SIGAD3_E8.5_pool2_TET-TKO_L004", 
  # "SIGAE3_E7.5_pool1_Host-WT_L005", 
  # "SIGAF3_E7.5_pool1_TET-TKO_L006", 
  "SIGAG3_E8.5_hashing_Host-WT_L007"
  # "SIGAH3_E8.5_hasting_TET-TKO_L008"
)

########################
## Load Seurat object ##
########################

srat <- readRDS("/Users/ricard/data/10x_gastrulation_TetChimera/processed/seurat.rds")

srat <- srat[,srat@meta.data$genotype=="WT" & srat@meta.data$batch%in%samples]
# meta.data <- srat@meta.data

############################
## Batch effect correcton ##
############################

# Split seurat object by batch
srat.split <- SplitObject(srat, split.by = "batch")
for (i in 1:length(srat.split)) {
  srat.split[[i]] <- NormalizeData(srat.split[[i]])
  srat.split[[i]] <- FindVariableFeatures(srat.split[[i]], nfeatures=1000)
  srat.split[[i]] <- ScaleData(srat.split[[i]])
}

# Find common HVGs
genes.use <- c()
for (i in 1:length(srat.split)) {
  genes.use <- c(genes.use, srat.split[[i]]@assays$RNA@var.features)
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(srat.split)) {
  srat.split[[i]] <- srat.split[[i]][genes.use,]
}

# Find anchors
# The returned object will contain a new Assay, which holds an integrated expression matrix for all cells
anchors <- FindIntegrationAnchors(srat.split, dims = 1:30)

# Do batch correction
# After running IntegrateData, the Seurat object will contain a new Assay with the integrated expression matrix. 
# Note that the original uncorrected values are still stored in the object in the “RNA” assay
srat.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

# Re-scale
srat.integrated <- ScaleData(srat.integrated, verbose = FALSE)

srat <- srat.integrated

# Change default assay
DefaultAssay(srat.integrated) <- "RNA"

#############
## Run PCA ##
#############

srat <- RunPCA(srat, pcs.compute = 15)

###################
## Run UMAP/TSNE ##
###################

# srat <- RunTSNE(srat)
srat <- RunUMAP(srat, dims = 1:15)

###################################
## Plot dimensionality reduction ##
###################################

DimPlot(srat, reduction="pca", group.by = "batch", dims = c(5,6))
DimPlot(srat, reduction="umap", group.by = "batch")

# Nascent mesoderm
FeaturePlot(srat, features = c("ENSMUSG00000006574"))

# ExE endoderm
FeaturePlot(srat, features = c("ENSMUSG00000031883","ENSMUSG00000028307"))

# Forebrain/Midbrain/Hindbrain (none?????)
FeaturePlot(srat, features = c("ENSMUSG00000058665","ENSMUSG00000039419"))

# ExE endoderm
FeaturePlot(srat, features = c("",""))

# Erithroyd3
FeaturePlot(srat, features = c("ENSMUSG00000006574","ENSMUSG00000028332"))

# Notochord
FeaturePlot(srat, features = c("ENSMUSG00000010136","ENSMUSG00000039676"))

# ExE ectoderm
FeaturePlot(srat, features = c("ENSMUSG00000070473","ENSMUSG00000061186"))

FeaturePlot(srat, features = c("ENSMUSG00000027562","ENSMUSG00000032698"))
FeaturePlot(srat, features = c("",""))
FeaturePlot(srat, features = c("",""))
FeaturePlot(srat, features = c("",""))

###############################################################
## Add dimensionality reduction information to the meta data ##
###############################################################

##########
## Save ##
##########

# saveRDS(srat, )


for (i in head(srat@assays$integrated@var.features,n=1)) {
  FeaturePlot(srat, features = i) %>% print
}
