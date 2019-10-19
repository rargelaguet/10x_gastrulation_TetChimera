library(SingleCellExperiment)
library(scater)

sce <- readRDS("/Users/ricard/data/10x_gastrulation_TetChimera/processed/sce.rds")

sce2 <- scater::normalize(sce)

logcounts(sce)[1:10,1:10]
logcounts(sce2)[1:10,1:10]

foo <- as.numeric(logcounts(sce))
bar <- as.numeric(logcounts(sce2))

cor(foo,bar)
