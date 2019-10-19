library(data.table)
library(purrr)
library(Seurat)

################
## Define I/O ##
################

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$directory <- "/Users/ricard/data/10x_gastrulation_TetChimera/original/test"
  io$genes <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  io$outfile <- "/Users/ricard/data/10x_gastrulation_TetChimera/processed/seurat.rds"
} else {
  io$directory <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_TetChimera/original/test"
  io$genes <- "/hps/nobackup/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  io$outfile <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_TetChimera/processed/seurat.rds"
}

##############################
## Load and merge data sets ##
##############################

samples <- c("SIGAA3_E8.5_pool1_Host-WT_L001", "SIGAB3_E8.5_pool1_TET-TKO_L002", "SIGAC3_E8.5_pool2_Host-WT_L003", "SIGAD3_E8.5_pool2_TET-TKO_L004", "SIGAE3_E7.5_pool1_Host-WT_L005", "SIGAF3_E7.5_pool1_TET-TKO_L006", "SIGAG3_E8.5_hashing_Host-WT_L007", "SIGAH3_E8.5_hasting_TET-TKO_L008")
genotype <- c("WT", "TKO", "WT", "TKO", "WT", "TKO", "WT", "TKO"); names(genotype) <- samples

mtx <- list()
cell.info <- list()
gene.info <- list()
for (i in samples) {
  
  # Load cell metadata
  barcode.loc <- sprintf("%s/%s_barcodes.tsv",io$directory,i)
  cell.info[[i]] <- read.delim(barcode.loc, header = FALSE, colClasses = "character", stringsAsFactors = FALSE, sep="\t")
  colnames(cell.info[[i]]) <- c("barcode")
  cell.info[[i]]$genotype <- genotype[[i]]
  cell.info[[i]]$batch <- i
  
  # Load gene metadata (note we could just load this once)
  gene.loc <- sprintf("%s/%s_features.tsv",io$directory,i)
  gene.info[[i]] <- read.delim(gene.loc, header = FALSE, colClasses = "character", stringsAsFactors = FALSE, sep="\t")[,c(1,2)]
  colnames(gene.info[[i]]) <- c("ens_id","gene_id")
  
  # Load matrix  
  matrix.loc <- sprintf("%s/%s_matrix.mtx",io$directory,i)
  mtx[[i]] <- Matrix::readMM(matrix.loc)
  rownames(mtx[[i]]) <- gene.info[[i]]$ens_id
  colnames(mtx[[i]]) <- cell.info[[i]]$V1
}
gene.info <- gene.info[[1]]

# Concatenate cell  metadata
cell.info <- do.call("rbind",cell.info)
cell.info$cell <- paste("cell",1:nrow(cell.info),sep="_")
rownames(cell.info) <- cell.info$cell

# Concatenate matrices
mtx <- do.call("cbind",mtx)
colnames(mtx) <- cell.info$cell

################
## Processing ##
################

# Remove cells with few UMIs mapping to the mouse genome
mouse.genes <- grep("mm10",rownames(mtx))
foo <- Matrix::colSums( mtx[mouse.genes,] )
cells.to.keep <- colnames(mtx)[which(foo>1000)]
mtx <- mtx[mouse.genes,cells.to.keep]

# Rename genes
mouse.genes <- rownames(mtx)
mouse.genes <- substr(mouse.genes,8,nchar(mouse.genes))
rownames(mtx) <- mouse.genes

# Remove duplicated genes (BECAUSE OF GENE SYMBOLS MATCHING TO MULTIPLE ESNEMBL IDS, TO-DO...)
# rownames(mtx)[duplicated(rownames(mtx))]

# Subset protein-coding genes
genes <- fread(io$genes)[,ens_id]
genes <- genes[genes %in% mouse.genes]
mouse.genes <- mouse.genes[mouse.genes %in% genes]
mtx <- mtx[mouse.genes,]

# Subset cell metadata
cell.info <- cell.info[colnames(mtx),]

# Subset gene metadata
gene.info <- gene.info[grep("mm10",gene.info$gene_id),]
gene.info$gene_id <- substr(gene.info$gene_id,8,nchar(gene.info$gene_id))
gene.info <- gene.info[rownames(mtx),]

############
## Seurat ##
############

# Create seurat object
srat <- CreateSeuratObject(mtx,
  min.features = 1, min.cells = 100,
  meta.data = cell.info
)

# Normalize data
srat <- NormalizeData(srat)

# Scale data
# srat <- ScaleData(srat,
#   vars.to.regress="nUMI",
#   model.use = "linear",
#   do.scale = FALSE, do.center = TRUE,
# )

##########################
## SingleCellExperiment ##
##########################

# sce <- as.SingleCellExperiment(srat)

##########
## Save ##
##########

saveRDS(srat, io$outfile)

