library(data.table)
library(purrr)
library(Seurat)

seurat <- readRDS("/Users/ricard/data/10x_gastrulation_TetChimera/processed/second_batch/seurat.rds")
metadata <- fread("/Users/ricard/data/10x_gastrulation_TetChimera/processed/second_batch/sample_metadata.txt.gz")
metadata.filt <- metadata %>% .[pass_QC==T]

dim(seurat)
dim(metadata)
dim(metadata.filt)

nrow(metadata.filt[!cell %in% colnames(seurat)]) # 2131
length(colnames(seurat)[!colnames(seurat)%in%metadata.filt$cell]) # 2196
