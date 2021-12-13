library(scater)
library(data.table)
library(purrr)

io <- list()
io$vel_file <- "D:/GASTRULATION_DATA/rna/velocyto/velocyto.rds"
io$sample_metadata <- "D:/GASTRULATION_DATA/sample_metadata_scNMT.txt"
io$out_folder <- "D:/GASTRULATION_DATA/rna/velocyto/sce"

vel <- readRDS(io$vel_file)
meta <- fread(io$sample_metadata ) %>% 
  .[id_rna %in% colnames(vel[[1]])] %>% 
  as.data.frame() %>% 
  tibble::column_to_rownames("id_rna") %>% 
  .[colnames(vel[[1]]),]

sce <- map(vel, ~SingleCellExperiment(., colData = meta))

dir.create(io$out_folder)

filenames <- paste0(io$out_folder, "/", c("exons", "introns", "exon_spanning"), "_sce.rds")

walk2(sce, filenames, saveRDS)
