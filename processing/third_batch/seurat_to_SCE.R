library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scran)

# Load default settings
source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
io$outfile <- paste0(io$basedir,"/processed/third_batch/SingleCellExperiment.rds")

# Define options
opts$batches <- c(
  "SIGAA6_E85_2_Dnmt3aKO_Dnmt3b_WT_L001",
  "SIGAB6_E85_3_Dnmt3aWT_Dnmt3b_WT_L002",
  "SIGAC6_E85_5_Dnmt3aKO_Dnmt3b_Het_L003",
  "SIGAD6_E85_8_Dnmt3aHet_Dnmt3b_KO_L004",
  "SIGAE4_E105_3_TET123_Chimera_Host_L005",
  "SIGAF4_E105_3_TET123_Chimera_TKO_L006",
  "SIGAG4_E105_5_TET123_Chimera_Host_L007",
  "SIGAH4_E105_5_TET123_Chimera_TKO_L008"
)

# Update sample metadata
sample_metadata <- sample_metadata %>%
  .[batch%in%opts$batches & pass_QC==T]

# Load Seurat object
seurat <- readRDS(io$seurat)[,sample_metadata$cell]
dim(seurat)

# Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(seurat)

# Add metadata
colData(sce) <- sample_metadata %>% as.data.frame %>% tibble::column_to_rownames("cell") %>%
  .[colnames(sce),] %>% DataFrame()

# Normalise using scran
clusts = as.numeric(quickCluster(sce, method = "igraph", min.size = 100, BPPARAM = mcparam))
min.clust = min(table(clusts))/2
new_sizes = c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce = computeSumFactors(sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

to.plot <- data.frame(X = Matrix::colSums(counts(sce)), Y = sizeFactors(sce))
ggplot(to.plot, mapping = aes(x = X, y = Y)) +
  geom_point() +
  labs(x = "Number of UMIs", y = "Size Factor") +
  theme_classic()

saveRDS(sce, io$outfile)