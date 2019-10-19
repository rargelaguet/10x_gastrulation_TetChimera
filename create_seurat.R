library(data.table)
library(purrr)
# library(scater)
library(Seurat)

directory <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/original"

# Load cell metadata
barcode.loc <- file.path(directory, "meta.tab")
cell.info <- read.delim(barcode.loc, header = TRUE, colClasses = "character", stringsAsFactors = FALSE, sep="\t")
rownames(cell.info) <- cell.info$cell
  
# Load gene metadata
gene.loc <- file.path(directory, "genes.tsv")
gene.info <- read.delim(gene.loc, header = FALSE, colClasses = "character", stringsAsFactors = FALSE, sep="\t")
colnames(gene.info) <- c("ens_id","gene_id")

# Load matrix  
matrix.loc <- file.path(directory, "raw_counts.mtx")
data_mat <- Matrix::readMM(matrix.loc)
rownames(data_mat) <- gene.info$gene_id
colnames(data_mat) <- cell.info$cell

# Subset  
data_mat <- data_mat[,cell.info$doublet==F & cell.info$stripped==F]
cell.info <- cell.info[cell.info$doublet==F & cell.info$stripped==F,]

# Create seurat object
srat <- CreateSeuratObject(data_mat,
  min.features = 1, min.cells = 100,
  meta.data = cell.info
)

# Normalize data
# srat <- NormalizeData(srat)

# Scale data
# srat <- ScaleData(srat, 
#   vars.to.regress="nUMI",
#   model.use = "linear",
#   do.scale = FALSE, do.center = TRUE, 
# )


# Save
saveRDS(srat, "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/processed/seurat/seurat.rds")