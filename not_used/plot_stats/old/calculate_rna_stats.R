library(SingleCellExperiment)
library(scater)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation10x/settings.R")
}

io$outfile <- paste0(io$basedir,"/results/general_stats/general_stats.txt.gz")

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- readRDS(io$rna.sce)[,as.character(sample_metadata$cell)]

############################
## Extract RNA statistics ##
############################

qc.metrics <- perCellQCMetrics(sce) %>% as.data.table %>%
  .[,c("detected","total")] %>%
  .[,cell:=colnames(sce)] %>%
  setnames(c("nFeature_RNA","nCount_RNA"))

##########
## Save ##
##########

fwrite(qc.metrics, io$outfile, sep="\t")

################################
## Merge with sample metadata ##
################################

# metadata <- fread(io$metadata)
# metadata <- metadata %>% merge(qc.metrics,by="cell",all.x=T)
# metadata[,c("celltype2","celltype3"):=NULL]
# 
# fwrite(metadata, io$metadata, sep="\t", quote=F, na="NA")
