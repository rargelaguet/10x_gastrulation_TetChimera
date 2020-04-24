library(data.table)
library(purrr)
library(Seurat)

################
## Define I/O ##
################

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$directory <- "/Users/ricard/data/10x_gastrulation_TetChimera/original/second_batch"
  io$genes <- "¡/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  io$outfile <- "/Users/ricard/data/10x_gastrulation_TetChimera/processed/second_batch/seurat.rds"
} else {
  # io$directory <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_TetChimera/original/test"
  # io$genes <- "/hps/nobackup2/research/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  # io$outfile <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_TetChimera/processed/seurat.rds"
}

##############################
## Load and merge data sets ##
##############################

samples <- c(
  "SIGAA3_E8.5_pool1_Host-WT_L001", 
  "SIGAB3_E8.5_pool1_TET-TKO_L002", 
  "SIGAC3_E8.5_pool2_Host-WT_L003", 
  "SIGAD3_E8.5_pool2_TET-TKO_L004", 
  "SIGAE3_E7.5_pool1_Host-WT_L005", 
  "SIGAF3_E7.5_pool1_TET-TKO_L006", 
  "SIGAG3_E8.5_hashing_Host-WT_L007", 
  "SIGAH3_E8.5_hasting_TET-TKO_L008"
)
genotype <- c("WT", "TKO", "WT", "TKO", "WT", "TKO", "WT", "TKO"); names(genotype) <- samples

mtx <- list()
cell.info <- list()
gene.info <- list()

for (i in samples) {
# for (i in head(samples,n=2)) {
  
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
rownames(gene.info) = gene.info$ens_id

# Concatenate cell  metadata
cell.info <- do.call("rbind",cell.info)
cell.info$cell <- paste("cell",1:nrow(cell.info),sep="_")
rownames(cell.info) <- cell.info$cell

# Add stage information to the sample metadata
cell.info$stage <- sapply(strsplit(cell.info$batch,"_"), "[[", 2)

# Concatenate matrices
mtx <- do.call("cbind",mtx)
colnames(mtx) <- cell.info$cell

################
## Processing ##
################

# Remove human genes
mouse.genes <- grep("mm10",rownames(mtx))
mtx <- mtx[mouse.genes,]

# Remove prefix from genes
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
gene.info$ens_id <- substr(gene.info$ens_id,8,nchar(gene.info$ens_id))
rownames(gene.info) <- substr(rownames(gene.info),8,nchar(rownames(gene.info)))
gene.info <- gene.info[rownames(mtx),]

##############
## QC cells ##
##############

# Remove cells with too few or too many UMIs
foo <- Matrix::colSums(mtx)
hist(foo)
cells.to.keep <- colnames(mtx)[which(foo>1000 & foo<5000)]
mtx <- mtx[,cells.to.keep]
hist(Matrix::colSums(mtx))

# Filter cells by mithocondrial gene expression
# (Note) Atlas: Cells with mitochondrial gene-expression fractions greater than 2.37% were excluded.
mt.genes <- gene.info$ens_id[ grep("mt-",gene.info$gene_id) ]
foo <- Matrix::colSums(mtx[mt.genes,]) / Matrix::colSums(mtx)
cells.to.keep <- colnames(mtx)[which(foo<0.15)]
mtx <- mtx[,cells.to.keep]
hist(Matrix::colSums(mtx[mt.genes,]))

# Subset sample metadata
cell.info <- cell.info[cell.info$cell%in%colnames(mtx),]
cell.info <- cell.info[colnames(mtx),]

##############
## QC genes ##
##############

# Remove genes with high dropout rate
foo <- Matrix::rowSums(mtx>0)
mtx <- mtx[foo>10,]

############
## Seurat ##
############

stopifnot(colnames(mtx) == cell.info$cell)

# Create seurat object
srat <- CreateSeuratObject(mtx, meta.data = cell.info)

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

max(srat@assays$RNA@data)
srat@assays$RNA@counts[1:5,1:5]


