suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))

#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/10x_gastrulation_TetChimera"
  io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  io$basedir <- "/hps/nobackup2/stegle/users/ricard/10x_gastrulation_TetChimera"
} else {
  stop("Computer not recognised")
}

io$metadata <- paste0(io$basedir,"/processed/second_batch/sample_metadata.txt.gz")
# io$gene_metadata <- paste0(io$basedir,"/features/genes/Mmusculus_genes_BioMart.87.txt")

io$rna.seurat <- paste0(io$basedir,"/processed/second_batch/seurat.rds")
# io$rna.sce <- paste0(io$basedir,"/rna/SingleCellExperiment.rds")
# io$rna.counts <- paste0(io$basedir,"/rna/counts_datatable.tsv.gz")

#############
## Options ##
#############

opts <- list()


##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) 
  
  
