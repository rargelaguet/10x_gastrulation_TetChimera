library(Seurat)

srat <- readRDS("/Users/ricard/data/10x_gastrulation_TetChimera/processed/seurat.rds")
meta.data <- srat@meta.data

# Scale data and regress out technical effects
srat <- ScaleData(srat, 
  # vars.to.regress=c("sample"),
  # model.use = "linear",
  do.scale = FALSE, do.center = TRUE
)

srat <- FindVariableFeatures(srat)

# Run
srat <- RunPCA(srat, pcs.compute = 15)
# srat <- RunTSNE(srat)
srat <- RunUMAP(srat, dims = 1:15)

# Plot
# DimPlot(srat, reduction.use="pca", group.by = "celltype", dim.1 = 2, dim.2 = 3)
# DimPlot(srat, reduction.use="tsne", group.by = "celltype")
DimPlot(srat, reduction.use="umap", group.by = "genotype")


