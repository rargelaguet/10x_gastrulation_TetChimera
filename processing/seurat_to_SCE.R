library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scran)

# Load default settings
# source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
source("/homes/ricard/10x_gastrulation_TetChimera/settings.R")

# Define I/O
# io$outfile <- paste0(io$basedir,"/processed/SingleCellExperiment.rds")
# io$outfile <- "/Users/ricard/data/10x_gastrulation_TetChimera/processed/first_batch/SingleCellExperiment.rds"
io$outfile <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_TetChimera/processed/first_batch/SingleCellExperiment.rds"

# Load Seurat object
# io$seurat <- "/Users/ricard/data/10x_gastrulation_TetChimera/processed/first_batch/seurat.rds"
io$seurat <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_TetChimera/processed/first_batch/seurat.rds"
seurat <- readRDS(io$seurat)

# Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(seurat)

# Add metadata
stopifnot(sample_metadata$cell%in%colnames(sce))
stopifnot(colnames(sce)%in%sample_metadata$cell)
sample_metadata <- sample_metadata %>% setkey(cell) %>% .[colnames(sce)]
stopifnot(sample_metadata$cell == colnames(sce))
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

sce <- logNormCounts(sce)

saveRDS(sce, io$outfile)
