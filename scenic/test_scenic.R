library(Seurat)
library(SCENIC)
library(purrr)

seurat <- readRDS("/Users/ricard/data/10x_rna_atac/seurat.rds")

seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", assay = "RNA")
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 5000, assay = "RNA", verbose = FALSE)

org="hgnc"
dbDir="/Users/ricard/data/hg38_regulation/cistarget"
myDatasetTitle="SCENIC example"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
dbs[["10kb"]] <- "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
dbs[["500bp"]] <- "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=2) 
# scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
# scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
# scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr"="mm9-tss-centered-5kb-10species.mc8nr.feather")
# scenicOptions@settings$db_mcVersion <- "v8"

# Gene selection
exprMat <- seurat@assays$RNA@data %>% as.matrix
# genesKept <- geneFiltering(exprMat, scenicOptions, minCountsPerGene=3*.01*ncol(exprMat), minSamples=ncol(exprMat)*.01)
# exprMat_filtered <- exprMat[genesKept, ]

# Calculate correlations between genes
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 

# Genie
runGenie3(exprMat_filtered_log, scenicOptions)

# SCENIC
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) #** Only for toy run!!
runSCENIC_3_scoreCells(scenicOptions, logMat)