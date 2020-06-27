suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--atlas_stages',    type="character",   nargs='+',  help='Atlas stage(s)')
p$add_argument('--query_batches',   type="character",   nargs='+',  help='Query batch(es)')
p$add_argument('--test',            action = "store_true",          help='Testing mode')
p$add_argument('--outdir',          type="character",               help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
args$atlas_stages <- c(
  "E6.5",
  "E6.75",
  "E7.0",
  "E7.25",
  "E7.5",
  "E7.75",
  "E8.0",
  "E8.25",
  "E8.5",
  "mixed_gastrulation"
)

args$query_batches <- c(
  # "SIGAA3_E8.5_pool1_Host-WT_L001"
  # "SIGAB3_E8.5_pool1_TET-TKO_L002",
  "SIGAC3_E8.5_pool2_Host-WT_L003",
  # "SIGAD3_E8.5_pool2_TET-TKO_L004", 
  # "SIGAE3_E7.5_pool1_Host-WT_L005", 
  # "SIGAF3_E7.5_pool1_TET-TKO_L006", 
  "SIGAG3_E8.5_hashing_Host-WT_L007",
  "SIGAH3_E8.5_hasting_TET-TKO_L008"
)

args$test <- FALSE
## END TEST ##

################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
  source("/Users/ricard/10x_gastrulation_TetChimera/mapping/run/mapping_functions.R")
} else {
  source("/homes/ricard/10x_gastrulation_TetChimera/settings.R")
  source("/homes/ricard/10x_gastrulation_TetChimera/mapping/run/mapping_functions.R")
}
io$path2atlas <- io$atlas.basedir
io$path2query <- io$basedir

if (isTRUE(args$test)) print("Test mode activated...")

################
## Load atlas ##
################

# Load SingleCellExperiment
sce_atlas  <- readRDS(io$atlas.sce)

# Load cell metadata
meta_atlas <- fread(io$atlas.metadata) %>%
  .[stripped==F & doublet==F & stage%in%args$atlas_stages]

# Filter
if (isTRUE(args$test)) meta_atlas <- head(meta_atlas,n=1000)
sce_atlas <- sce_atlas[,meta_atlas$cell] 

################
## Load query ##
################

# Load SingleCellExperiment
io$sce <- "/Users/ricard/data/10x_gastrulation_TetChimera/processed/first_batch/SingleCellExperiment.rds"
sce_query <- readRDS(io$sce)

# Load cell metadata
io$metadata <- "/Users/ricard/data/10x_gastrulation_TetChimera/processed/first_batch/sample_metadata.txt.gz"
# meta_query <- fread(io$metadata) %>% .[pass_QC==T & batch%in%args$query_batches]
# meta_query <- fread(io$metadata) %>% .[nFeature_RNA>1000 & batch%in%args$query_batches]
meta_query <- fread(io$metadata) %>% .[nCount_RNA>1250 & batch%in%args$query_batches]
table(meta_query$batch)

# Filter
if (isTRUE(args$test)) meta_query <- head(meta_query,n=1000)
sce_query <- sce_query[,meta_query$cell]

#############
## Prepare ## 
#############

# Filter out non-expressed genes
sce_query <- sce_query[rowMeans(counts(sce_query))>1e-3,]
sce_atlas <- sce_atlas[rowMeans(counts(sce_atlas))>1e-3,]

# Intersect genes
genes.intersect <- intersect(rownames(sce_query), rownames(sce_atlas))
sce_query  <- sce_query[genes.intersect,]
sce_atlas <- sce_atlas[genes.intersect,]

# Load gene markers to be used as HVGs
# marker_genes.dt <- fread(io$atlas.marker_genes)
# marker_genes.dt <- marker_genes.dt[,head(.SD,n=50),by="celltype"]
# marker_genes <- unique(marker_genes.dt$ens_id)
# marker_genes <- marker_genes[marker_genes%in%genes.intersect]
# print(marker_genes.dt[,.N,by="celltype"])
# stopifnot(all(marker_genes%in%rownames(sce_atlas)))

#########
## Map ##
#########

mapping  <- mapWrap(
  atlas_sce = sce_atlas, atlas_meta = meta_atlas,
  map_sce = sce_query, map_meta = meta_query, 
  genes = NULL, npcs = 50, k = 25
)

##########
## Save ##
##########

outfile <- sprintf("%s/mapping_mnn_%s.rds",args$outdir,paste(args$query_batches,collapse="-"))
saveRDS(mapping, outfile)
