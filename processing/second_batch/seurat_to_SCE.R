library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scran)

# Load default settings
source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
io$outfile <- paste0(io$basedir,"/processed/second_batch/SingleCellExperiment.rds")

# Define options
opts$batches <- c(
	"E75_TET_TKO_L002",
	"E75_WT_Host_L001",
	"E85_Rep1_TET_TKO_L004",
	"E85_Rep2_TET_TKO_L006",
	"E85_Rep1_WT_Host_L003",
	"E85_Rep2_WT_Host_L005"
	# "E125_DNMT3A_HET_A_L001",
	# "E125_DNMT3A_HET_A_L003",
	# "E125_DNMT3A_KO_B_L002",
	# "E125_DNMT3A_KO_E_L004"
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