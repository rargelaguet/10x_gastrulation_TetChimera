library(data.table)
library(purrr)
library(Seurat)

seurat <- readRDS("/Users/ricard/data/10x_gastrulation_TetChimera/processed/second_batch/seurat.rds")
metadata <- fread("/Users/ricard/data/10x_gastrulation_TetChimera/processed/second_batch/sample_metadata.txt.gz")