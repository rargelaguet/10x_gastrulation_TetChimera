matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--embryo',   type="character",       help = 'Embryo name')
p$add_argument('--test',     action = "store_true",  help = 'Testing mode')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##  
# args$embryo <- "embryo1"
# args$test <- TRUE
## END TEST ##  

#####################
## Define settings ##
#####################

# Load default settings
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_spatial/settings.R")
  source("/Users/ricard/gastrulation_spatial/iterative_mapping_10x/run/utils.R")
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation_spatial/settings.R")
  source("/homes/ricard/gastrulation_spatial/iterative_mapping_10x/run/utils.R")
}  
io$outdir <- paste0(io$basedir,"/iterative_mapping"); dir.create(io$outdir, showWarnings = F)

###################
## Load metadata ##
###################

# Load spatial metadata
sample_metadata_query <- fread(io$metadata) %>%
  .[embryo%in%args$embryo] 

# Load atlas metadata
sample_metadata_atlas <- fread(io$atlas.metadata) %>%
  .[stripped==F & doublet==F] %>%
  .[stage%in%c("E8.5")] %>%
  .[,celltype:=stringr::str_replace_all(celltype,"_", " ")] %>%
  .[!celltype%in%c("Erythroid1","Erythroid3","Blood progenitors 1")]
  
# Merge some lineages in the atlas
# sample_metadata_atlas[celltype%in%c("Erythroid1","Erythroid2","Erythroid3"),celltype:="Erythroid"]
# sample_metadata_atlas[celltype%in%c("Blood progenitors 1","Blood progenitors 2"),celltype:="Blood progenitors"]

# Subset cells
if (isTRUE(args$test)) {
  sample_metadata_query <- sample_metadata_query %>% head(n=1000)
  sample_metadata_atlas <- sample_metadata_atlas %>% split(.$celltype) %>% map(~ head(.,n=100)) %>% rbindlist
}

# Filter cell types with low numbers
opts$celltypes <- which(table(sample_metadata_atlas$celltype)>25) %>% names
sample_metadata_atlas <- sample_metadata_atlas[celltype%in%opts$celltypes]

# print statistics
table(sample_metadata_atlas$celltype)

#########################
## Load RNA expression ##
#########################

# Load RNA expression data for the query
sce.query <- readRDS(io$sce)[,sample_metadata_query$cell]
sce.query <- sce.query[rownames(sce.query)!="Xist",]

# sce.query$embryo_pos_z = factor(paste0(sce.query$embryo,"_", sce.query$pos, "_", sce.query$z))

# Load RNA expression data for the atlas
sce.atlas <- readRDS(io$atlas.sce)[,sample_metadata_atlas$cell]
colData(sce.atlas) <- sample_metadata_atlas %>% tibble::column_to_rownames("cell") %>% DataFrame

################
## Parse data ##
################

# Change gene names from ENSEMBL to symbols
gene_metadata <- fread(io$gene_metadata) %>% .[,c("ens_id","symbol")] %>%
  .[symbol!="" & ens_id%in%rownames(sce.atlas)]
sce.atlas <- sce.atlas[rownames(sce.atlas)%in%gene_metadata$ens_id,]
foo <- gene_metadata$symbol; names(foo) <- gene_metadata$ens_id
new.names <- foo[rownames(sce.atlas)]
rownames(sce.atlas) <- new.names
rownames(sce.atlas)[rownames(sce.atlas) == "Prkcdbp"] <- "Cavin3"

# remove non-expressed genes
sce.query <- sce.query[rowMeans(counts(sce.query))>0,]
sce.atlas <- sce.atlas[rowMeans(counts(sce.atlas))>0,]

# intersect genes
genes <- intersect(rownames(sce.query), rownames(sce.atlas))
sce.query  <- sce.query[genes,]
sce.atlas <- sce.atlas[genes,]
dim(sce.query); dim(sce.atlas)

#############
## Run MNN ##
#############

mapping.dt <- mnn.fn(sce.query, sce.atlas, npcs = 30, k = 25, cosineNorm = TRUE)
  
##########
## Save ##
##########

fwrite(mapping.dt, sprintf("%s/%s_standard_mnn.txt.gz",io$outdir,args$embryo), sep="\t")

