library(MOFA2)
library(SingleCellExperiment)
library(Seurat)
library(RColorBrewer)
# library(scater)
# library(scran)

#####################
## Define settings ##
#####################

source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
io$outfile <- paste0(io$basedir,"/mofa/model.hdf5")

opts$batches <- c(
  "E75_TET_TKO_L002",
  "E75_WT_Host_L001"
  # "E85_Rep1_TET_TKO_L004",
  # "E85_Rep2_TET_TKO_L006",
  # "E85_Rep1_WT_Host_L003",
  # "E85_Rep2_WT_Host_L005"
  # "E125_DNMT3A_HET_A_L001",
  # "E125_DNMT3A_HET_A_L003",
  # "E125_DNMT3A_KO_B_L002",
  # "E125_DNMT3A_KO_E_L004"
)


# Update sample metadata
sample_metadata <- sample_metadata %>%
  .[batch%in%opts$batches & pass_QC==T]
table(sample_metadata$batch)

###############
## Load data ##
###############

# Load RNA expression data as Seurat object
seurat <- readRDS(io$seurat)[,sample_metadata$cell]
dim(seurat)

# Normalise
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize")

# Scale
# seurat <- ScaleData(seurat)

################
## Parse data ##
################

# Change gene names from ENSEMBL ID to gene symbols
gene_metadata <- fread(io$gene_metadata) %>%  .[,c("symbol","ens_id")]
seurat <- subset(seurat, features = rownames(seurat)[rownames(seurat)%in%gene_metadata$ens_id])
gene_metadata <- gene_metadata %>% .[ens_id%in%rownames(seurat)] %>% setkey(ens_id) %>% .[rownames(seurat)]
rownames(seurat@assays$RNA@counts) <- rownames(seurat@assays$RNA@data) <- gene_metadata$symbol
seurat[["RNA"]]@meta.features <- data.frame(row.names=gene_metadata$symbol)


# Select HVG (per batch)
hvg_list <- list()
for (i in unique(sample_metadata$batch)) {
  seurat.subset <- seurat[,sample_metadata[batch==i,cell]]
  seurat.subset <- FindVariableFeatures(seurat.subset, selection.method = 'vst', nfeatures = 1500)
  hvg_list[[i]] <- seurat.subset@assays$RNA@var.features
}
hvgs <- unique(unlist(hvg_list))

# Plot venn diagram with the overlap of HVG
foo <- VennDiagram::venn.diagram(
  x = hvgs,
  filename = NULL,
  col = "transparent", 
  fill = brewer.pal(n=length(hvg_list), "Spectral")[1:length(hvg_list)],
  alpha = 0.60, 
  cex = 1.5, 
  fontfamily = "serif", fontface = "bold"
)
grid.draw(foo)

# pdf(file=sprintf("%s/venn_all.pdf",io$outdir))
# dev.off()

#######################
# Create MOFA object ##
#######################

MOFAobject <- create_mofa(seurat, groups="batch", features=hvgs)

# Visualise data structure
# plot_data_overview(MOFAobject)

####################
## Define options ##
####################

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 15

# Training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "fast"
train_opts$seed <- 42
train_opts$maxiter <- 5

#########################
## Prepare MOFA object ##
#########################

MOFAobject <- prepare_mofa(MOFAobject,
  model_options = model_opts,
  training_options = train_opts
)

#####################
## Train the model ##
#####################

model <- run_mofa(MOFAobject, io$outfile)
