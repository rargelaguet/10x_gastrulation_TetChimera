library(data.table)
library(purrr)
library(rhdf5)
library(Matrix)


indir <- "/Users/clarks/data/tet_tko/processed/velocyto/"


looms <- dir(indir, pattern = ".loom$", full = TRUE)
loom=looms[[1]]

walk(looms, function(loom){
  h5ls(loom)
  cells <- h5read(loom, "col_attrs")
  genes <- h5read(loom, "row_attrs")
  mats <- h5read(loom, "layers")
  
  # matrices do not have dim names and need transposing to genes x cells
  mats <- map(mats, function(mat){
    mat <- t(mat)
    colnames(mat) <- cells$CellID
    rownames(mat) <- genes$Gene
    Matrix(mat, sparse = TRUE)
  })
  
  map(mats, ~.[1:5, 1:5])
  
  cell_meta <- as.data.table(cells) %>%
    .[, barcode := gsub(".*:", "", CellID) %>% gsub("x", "-1", .)]
  gene_meta <- as.data.table(genes)
  
  output <- list(cell_meta = cell_meta, 
                 gene_meta = gene_meta,
                 mats = mats)
  outname <- gsub(".loom", ".rds", loom)
  saveRDS(output, outname)
})


