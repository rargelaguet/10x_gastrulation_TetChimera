  sce_all <- SingleCellExperiment::SingleCellExperiment(
    list(counts=Matrix::Matrix(cbind(counts(atlas_sce),counts(map_sce)),sparse=TRUE)))
  big_sce <- scater::normalize(sce_all)
  message("Done\n")