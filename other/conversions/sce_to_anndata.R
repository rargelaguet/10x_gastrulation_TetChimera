library(SingleCellExperiment)
library(scran)
library(zellkonverter)


#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
  source("/Users/ricard/10x_gastrulation_TetChimera/utils.R")
} else {
  source("/homes/ricard/10x_gastrulation_TetChimera/settings.R")
  source("/homes/ricard/10x_gastrulation_TetChimera/utils.R")
}

setwd("/Users/ricard/10x_gastrulation_TetChimera/scanpy")

# Define I/O
# io$metadata <- paste0(io$basedir,"/results/rna/doublets/sample_metadata_after_doublets.txt.gz")
io$outfile <- io$anndata

# Define options
opts$samples <- c(
  "E75_TET_TKO_L002",
  "E75_WT_Host_L001",
  "E85_Rep1_TET_TKO_L004",
  "E85_Rep1_WT_Host_L003",
  "E85_Rep2_TET_TKO_L006",
  "E85_Rep2_WT_Host_L005",
  "E8_5_TET_WT_rep1_SIGAG8",
  "E8_5_TET_WT_rep2_SIGAH8"
)

###############
## Load data ##
###############

# Load cell metadata
sample_metadata <- fread(io$metadata) %>%
  .[pass_QC==TRUE & sample%in%opts$samples]

# Load SingleCellExperiment
sce <- load_SingleCellExperiment(io$sce,  cells=sample_metadata$cell, remove_non_expressed_genes = T)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

################
## Conversion ##
################

adata <- SCE2AnnData(sce, X_name = "counts", skip_assays = FALSE)

##########
## Save ##
##########

# writeH5AD(anndata.object, file = io$outfile)
adata$write_h5ad(io$outfile)

##########
## Test ##
##########

# library(DelayedArray)
# library(reticulate)

# R.utils::sourceDirectory("/Users/ricard/git/zellkonverter/R", verbose=T, modifiedOnly=FALSE)

# sce.backup <- sce
# sce <- sce[1:1000,1:1000]
# X_name = "counts"
# skip_assays = FALSE
