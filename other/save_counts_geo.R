library(Seurat)
outdir <- "/Users/argelagr/data/10x_gastrulation_TetChimera/processed/geo"
seurat <- readRDS("/Users/argelagr/data/10x_gastrulation_TetChimera/processed_all/seurat.rds")

sample2alias <- c(
  "E75_TET_TKO_L002"           = "E7.5_Tet_TKO#",
  "E75_WT_Host_L001"           = "E7.5_WT_tdTomato-_1#",
  "E85_Rep1_TET_TKO_L004"      = "E8.5_Tet_TKO_1#",
  "E85_Rep1_WT_Host_L003"      = "E8.5_WT_tdTomato-_1#",
  "E85_Rep2_TET_TKO_L006"      = "E8.5_Tet_TKO_2#",
  "E85_Rep2_WT_Host_L005"      = "E8.5_WT_tdTomato-_2#",
  "E8_5_TET_WT_rep1_SIGAG8"    = "E8.5_Tet_TKO_3#",
  "E8_5_TET_WT_rep2_SIGAH8"    = "E8.5_Tet_TKO_4#"
)

seurat <- seurat[,seurat$sample%in%names(sample2alias)]

# rename samples
cellnames <- colnames(seurat@assays$RNA@counts)
cellnames <- stringr::str_replace_all(cellnames, sample2alias)
cellnames <- stringr::str_replace_all(cellnames, "#_","#")

# mtx <- 
tmp <- stringr::str_split(cellnames,"#") %>% map_chr(1)
table(tmp)

# thsi doesnt work
# length(cellnames)==ncol(seurat)
# colnames(seurat) <- cellnames

Matrix::writeMM(seurat@assays$RNA@counts, file.path(outdir,"counts.mtx"))
write.table(cellnames, file.path(outdir,"barcodes.txt"), quote = F, row.names = F, col.names = F)
write.table(rownames(seurat), file.path(outdir,"features.txt"), quote = F, row.names = F, col.names = F)

to.save <- data.table(table(tmp)) %>% setnames(c("sample","ncells"))
fwrite(to.save, file.path(outdir,"stats.txt"), sep="\t")
