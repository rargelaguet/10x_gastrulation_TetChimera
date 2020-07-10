# Load libraries
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(batchelor))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(edgeR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--atlas_stages',    type="character",   nargs='+',  help='Atlas stage(s)')
p$add_argument('--query_batches',   type="character",   nargs='+',  help='Query batch(es)')
p$add_argument('--test',            action = "store_true",  help = 'Testing mode')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$atlas_stages <- c(
#   # "E6.5"
#   # "E6.75",
#   # "E7.0",
#   # "E7.25",
#   # "E7.5",
#   # "E7.75",
#   # "E8.0",
#   "E8.25",
#   "E8.5"
#   # "mixed_gastrulation"
# )
# 
# args$query_batches <- c(
#   # "E75_TET_TKO_L002",
#   # "E75_WT_Host_L001",
#   # "E85_Rep1_TET_TKO_L004"
#   # "E85_Rep2_TET_TKO_L006",
#   # "E85_Rep1_WT_Host_L003"
#   "E85_Rep2_WT_Host_L005"
#   # "E125_DNMT3A_HET_A_L001",
#   # "E125_DNMT3A_HET_A_L003",
#   # "E125_DNMT3A_KO_B_L002",
#   # "E125_DNMT3A_KO_E_L004"
# )
# 
# args$test <- FALSE
## END TEST ##

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
  source("/Users/ricard/10x_gastrulation_TetChimera/iterative_mapping/run/utils.R")
  io$atlas.marker_genes <- "/Users/ricard/data/gastrulation10x/results/marker_genes/E8.5/marker_genes.txt.gz"
  io$script_load_data <- "/Users/ricard/10x_gastrulation_TetChimera/iterative_mapping/run/load_data.R"
} else {
  source("/homes/ricard/10x_gastrulation_TetChimera/settings.R")
  source("/homes/ricard/10x_gastrulation_TetChimera/iterative_mapping/run/utils.R")
  io$atlas.marker_genes <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/results/marker_genes/E8.5/marker_genes.txt.gz"
  io$script_load_data <- "/homes/ricard/10x_gastrulation_TetChimera/iterative_mapping/run/load_data.R"
}
io$path2atlas <- io$atlas.basedir
io$path2query <- io$basedir
io$outdir <- paste0(io$basedir,"/results/iterative_mapping")

if (isTRUE(args$test)) print("Test mode activated...")

###############
## Load data ##
###############

source(io$script_load_data)

########################################################
## Define distance matrix for hierarchical clustering ##
########################################################

opts$celltypes <- which(table(meta_atlas$celltype)>25) %>% names %>% stringr::str_replace_all("_", " ")
# opts$celltypes <- unique(sample_metadata_atlas$celltype) 

dist <- fread(paste0(io$atlas.basedir,"/results/phylogenetic_tree/PAGA_distances.csv.gz")) %>%
  as.data.frame %>% tibble::column_to_rownames("V1") %>% as.matrix %>%
  .[opts$celltypes,opts$celltypes] %>%
  as.dist

#######################
## Recursive mapping ##
#######################

sce_query$celltype_mapped <- paste(opts$celltypes,collapse="%")

while (any(grepl("%",sce_query$celltype_mapped))) {
  print(table(sce_query$celltype_mapped))
  mapping_dt <- recursive.fn(sce_query, sce_atlas, dist)
  ids <- match(mapping_dt$cell,colnames(sce_query))
  sce_query$celltype_mapped[ids] <- mapping_dt$celltype_mapped
  sce_query$mapping_score[ids] <- mapping_dt$mapping_score
}

table(sce_query$celltype_mapped)

##########
## Save ##
##########

mapping_dt <- data.table(
  cell = colnames(sce_query), 
  celltype_mapped = sce_query$celltype_mapped,
  mapping_score = sce_query$mapping_score
)
# foo <- mapping_dt %>% merge(meta_query[,c("cell","celltype.mapped","celltype.score")] %>% setnames(c("cell","celltype_old","score_old")))

fwrite(mapping_dt, sprintf("%s/%s_iterative_mnn.txt.gz",io$outdir,paste(args$query_batches,collapse="_")), sep="\t")

