here::i_am("processing/1_create_seurat_rna.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(Seurat))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--inputdir',        type="character",                    help='Input directory')
p$add_argument('--outdir',       type="character",                    help='Output directory')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
p$add_argument('--test',            action="store_true",                 help='Testing mode')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args <- list()
# args$inputdir <- paste0(io$basedir,"/original/td_tomato/all_batches_symbolic")
# args$outdir <- paste0(io$basedir,"/processed_new")
# args$samples <- opts$samples[1:2]
# args$test <- FALSE
## END TEST ##

# These settings are important for consistency with ArchR, which provides little flexibility to edit cell names
opts$trim.barcode <- FALSE
opts$sample_cell_separator <- "_"

##############################
## Load and merge data sets ##
##############################

stopifnot(args$samples%in%opts$samples)
if (args$test) args$samples <- head(args$samples,n=2)

count_mtx <- list()
cell.info <- list()

for (i in args$samples) {
  print(i)
    
  # Load cell metadata
  barcode.loc <- sprintf("%s/%s/barcodes.tsv.gz",args$inputdir,i)
  cell.info[[i]] <- fread(barcode.loc, header=F) %>%
    setnames("barcode") %>%
    .[,barcode:=ifelse(rep(opts$trim.barcode,.N),gsub("-1","",barcode),barcode)] %>%
    .[,c("sample","cell"):=list(i,sprintf("%s%s%s",i,opts$sample_cell_separator,barcode))]
  dim(cell.info[[i]])
  
  # Load matrix  
  count_mtx[[i]] <- Read10X(sprintf("%s/%s",args$inputdir,i))
}

print(lapply(count_mtx,dim))

#######################
## Keep common genes ##
#######################

genes <- Reduce("intersect",lapply(count_mtx,rownames))
for (i in 1:length(count_mtx)) {
  count_mtx[[i]] <- count_mtx[[i]][genes,]
}

stopifnot(length(unique(lapply(count_mtx,nrow)))==1)
stopifnot(length(unique(lapply(count_mtx,rownames)))==1)

#################
## Concatenate ##
#################

# Concatenate cell metadata
cell.info <- rbindlist(cell.info)
rownames(cell.info) <- cell.info$cell

# Concatenate matrices
count_mtx <- do.call("cbind",count_mtx)
colnames(count_mtx) <- cell.info$cell

##################
## Filter genes ##
##################

# Remove duplicated genes
count_mtx <- count_mtx[!duplicated(rownames(count_mtx)),]

# Sanity checks
stopifnot(sum(duplicated(rownames(count_mtx)))==0)
stopifnot(sum(duplicated(colnames(count_mtx)))==0)
stopifnot("tomato-td"%in%rownames(count_mtx))

##########################
## Create Seurat object ##
##########################

cell.info.to.seurat <- cell.info[cell%in%colnames(count_mtx)] %>% setkey(cell) %>% .[colnames(count_mtx)] %>% as.data.frame
rownames(cell.info.to.seurat) <- cell.info.to.seurat$cell
stopifnot(rownames(cell.info.to.seurat)==colnames(count_mtx))
stopifnot(sum(is.na(rownames(cell.info.to.seurat$cell)))==0)

seurat <- CreateSeuratObject(count_mtx, meta.data = cell.info.to.seurat)

# Add mit percenatge
seurat[["mit_percent_RNA"]] <- PercentageFeatureSet(seurat, pattern = "^mt-") %>% round(2)

# Add rib RNA content
ribo.genes <- grep(pattern = "^Rp[l|s]", x = rownames(seurat), value = TRUE)
seurat[["rib_percent_RNA"]] <- PercentageFeatureSet(seurat, features = ribo.genes) %>% round(2)

#####################
## Create metadata ##
#####################

metadata <- seurat@meta.data %>% as.data.table %>% .[,orig.ident:=NULL] %>%
  .[,c("cell","barcode","sample","nFeature_RNA","nCount_RNA","mit_percent_RNA","rib_percent_RNA")]

# Add stage information
metadata[,stage:=as.character(NA)] %>%
  .[grepl("E75",sample),stage:="E7.5"] %>%
  .[grepl("E8",sample),stage:="E8.5"] %>%
  .[grepl("E9_5",sample),stage:="E9.5"]

print(table(metadata$stage))
stopifnot(!is.na(metadata$stage))

##########
## Save ##
##########

fwrite(metadata, file.path(args$outdir,"metadata.txt.gz"), quote=F, na="NA", sep="\t")
# fwrite(cell.info, paste0(args$outdir,"/cell_info.txt.gz"), quote=F, na="NA", sep="\t")
# fwrite(gene.info, paste0(args$outdir,"/gene_info.txt.gz"), quote=F, na="NA", sep="\t")
saveRDS(seurat, file.path(args$outdir,"seurat.rds"))
