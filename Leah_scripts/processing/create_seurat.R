#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library(optparse))

###################
## Parse options ##
###################
 
option_list = list(
    make_option(c("-i", "--inputseurat"), type="character", default=NULL, 
                help="input seurat file", metavar="character"),
    make_option(c("-S", "--settings"), type="character", default=NULL, 
                help="path to settings file", metavar="character"),
    make_option(c("-O", "--output.seurat"), type="character", default=NULL, 
                help="path to output seurat file", metavar="character"),
    make_option(c("-o", "--output.metadata"), type="character", default=NULL, 
                help="path to output sample metadata file", metavar="character"),
    make_option(c("-e", "--experiment"), type="character", default=NULL, 
                help="specify experiment (currently either 10x_old_EBsvGastruloids or 10x_new_EBsvMPs)", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);

if (is.null(opts$inputseurat)){
    print_help(opt_parser)
    stop("An input seurat file must be supplied.n", call.=FALSE)
} else if (is.null(opts$settings)) {
    print_help(opt_parser)
    stop("A settings file must be supplied.n", call.=FALSE)
} else if (is.null(opts$output.seurat)) {
    print_help(opt_parser)
    stop("A directory to a Seurat output file must be supplied.n", call.=FALSE)
} else if (is.null(opts$output.metadata)) {
    print_help(opt_parser)
    stop("A directory to a sample metadata output file must be supplied.n", call.=FALSE)
} else if (is.null(opts$experiment)) {
    print_help(opt_parser)
    stop("The experiment must be specified.n", call.=FALSE)
}

##########################
## Source settings file ##
##########################

source(opts$settings)


####################
## Load data sets ##
####################

message("Loading dataset")

srat.processed <- readRDS(paste0(opts$inputseurat))

srat <- CreateSeuratObject(srat.processed@assays$RNA@counts, meta.data = srat.processed@meta.data); rm(srat.processed)


########
## QC ##
########

mtx_sub_list <- list()
cell.info_sub_list <- list()
celltypes <- unique(srat@meta.data$genotype)
for (c in celltypes) {
    foo <- io$min_nFeature_RNA[[opts$experiment]][[c]]
    bar <- io$min_nCount_RNA[[opts$experiment]][[c]]
    tmp <- subset(srat, subset = nFeature_RNA > foo & nCount_RNA > bar, cells = rownames(srat@meta.data)[which(srat@meta.data$genotype == c)])
    cell.info_sub_list[[c]] <- tmp@meta.data
    cell.info_sub_list[[c]]$origrn <- rownames(cell.info_sub_list[[c]])
    mtx_sub_list[[c]] <- tmp@assays$RNA@counts
}
mtx <- do.call("cbind", mtx_sub_list)
cell.info <- do.call("rbind", cell.info_sub_list)
rownames(cell.info) <- cell.info$origrn

############
## Seurat ##
############

# Create seurat object

message("Creating seurat...")

srat <- CreateSeuratObject(mtx, meta.data = cell.info)


##########################
## SingleCellExperiment ##
##########################

# sce <- as.SingleCellExperiment(srat)

##########
## Save ##
##########

saveRDS(srat, opts$output.seurat)
fwrite(cell.info, opts$output.metadata, quote=F, na="NA", sep="\t")
