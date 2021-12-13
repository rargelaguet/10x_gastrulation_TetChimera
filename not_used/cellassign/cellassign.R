# Set up reticulate connection
suppressPackageStartupMessages(library(reticulate))
if (grepl("ricard",Sys.info()['nodename'])) {
  use_python("/Users/ricard/anaconda3/envs/gpflow_v2/bin/python", required = TRUE)
} else if(grepl("ebi",Sys.info()['nodename'])){
  use_python("/nfs/research1/stegle/users/ricard/conda-envs/gpflow2/bin/python", required = TRUE)
}  
py_config()

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(tensorflow))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(cellassign))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--batch',    type="character",   nargs='+',   help='Batch')
p$add_argument('--outdir',    type="character",               help='Output directory')
p$add_argument('--test',      action = "store_true",          help='Testing mode')
args <- p$parse_args(commandArgs(TRUE))

## START TEST 
args$batch <- c("E75_TET_TKO_L002")
args$outdir <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_TetChimera/results/second_batch/cellassign"
args$test <- TRUE
## END TEST

#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/10x_gastrulation_TetChimera/settings.R")
}  

# I/O
dir.create(args$outdir, showWarnings = F); dir.create(paste0(args$outdir,"/pdf"), showWarnings = F)

# Cell types to use
opts$celltypes <- opts$celltypes.1

# Maximum of M marker genes per cell type (sorted according to marker score)
opts$max.genes <- 50

# Update metadata
sample_metadata <- sample_metadata %>% 
  .[batch%in%args$batch]


# Acitvate test mode
if (isTRUE(args$test)) {
  print("Testing mode, subsetting cells...")
  sample_metadata <- sample_metadata %>% head(n=100)
}

###############
## Load data ##
###############

# Load RNA expression data
sce <- readRDS(io$sce)[,sample_metadata$cell]

#######################
## Load marker genes ##
#######################

marker_genes.dt <- fread(io$atlas.marker_genes) %>%
  .[celltype%in%opts$celltypes] %>%
  setorder(celltype,-score)

# Sanity check
stopifnot(all(unique(marker_genes.dt$group)%in%opts$celltypes))

# Maximum number of marker genes per cell type
marker_genes.dt <- marker_genes.dt[,head(.SD,n=opts$max.genes),by="celltype"]

######################
## Print statistics ##
######################

# print("Number of marker genes per cell type")
# print(marker_genes.dt[,.N,by="celltype"])

print("Total number of marker genes")
print(length(unique(marker_genes.dt$ens_id)))

print("Total number of cells")
print(ncol(sce))

################
## cellassign ##
################

# Create binary membership matrix
marker_gene_list <- split(marker_genes.dt, by="celltype") %>% map(~ .$ens_id)
bmat <- marker_list_to_mat(marker_gene_list)

# Subset SingleCellExperiment
sce <- sce[rownames(bmat),]

# Extract size factors
s <- sizeFactors(sce)

# Run
fit <- cellassign(sce, 
  marker_gene_info = bmat, 
  min_delta = 1,
  s = s, 
  # learning_rate = 1e-2, 
  shrinkage = TRUE,
  verbose = FALSE
)
print(fit)

##################
## Query output ##
##################

# Plot heatmap of cell type probabilities
outfile <- sprintf("%s/pdf/heatmap_probabilities_%s.pdf",args$outdir,paste(args$batch,collapse="-"))
pdf(outfile, width = 8, height = 8)
pheatmap::pheatmap(cellprobs(fit))
dev.off()

# Compare to ground truth
# foo <- table(sce$group, celltypes(fit))
# pdf(sprintf("%s/pdf/heatmap_assignments.pdf",args$outdir), width = 8, height = 8)
# pheatmap::pheatmap(foo, cluster_rows = F, cluster_cols = F)
# dev.off()

##########
## Save ##
##########

outfile <- sprintf("%s/cellassign_fit_%s.rds",args$outdir,paste(args$batch,collapse="-"))
print(outfile)
saveRDS(fit, outfile)