io <- list()

io$samples <- list()
io$samples[["TettKOchimera"]] <- c(
  "WT",
  "TKO"
)
#io$samples <- NULL

io$subset.proteincoding <- NULL # Choose between NULL and path to proteincoding genes

io$minUMIs <- 500
io$min_nFeature_RNA <- list()
io$min_nCount_RNA <- list()
io$min_nFeature_RNA[["TettKOchimera"]] <- c(
    'WT' = 675,
    'TKO' = 675
)
io$min_nCount_RNA[["TettKOchimera"]] <- c(
    'WT' = 1250,
    'TKO' = 1250
)