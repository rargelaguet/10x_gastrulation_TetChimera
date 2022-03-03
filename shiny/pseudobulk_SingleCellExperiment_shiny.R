here::i_am("pseudobulk/pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))
source(here::here("pseudobulk/utils.R"))

# I/O
io$metadata <- "/Users/argelagr/shiny_tet/data/cell_metadata.txt.gz"
io$sce <- file.path(io$basedir,"processed/SingleCellExperiment.rds")
io$outdir <- "/Users/argelagr/shiny_tet/data"

# Options
opts$group_by <- "genotype_sample_celltype"
opts$normalisation_method <- "cpm"
opts$rename_celltypes <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Anterior_Primitive_Streak" = "Primitive_Streak"
)

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  .[,celltype:=stringr::str_replace_all(celltype,opts$rename_celltypes)] %>%
  .[,genotype_sample_celltype:=sprintf("%s-%s-%s",genotype,sample,celltype)] %>%
  .[!is.na(eval(as.name(opts$group_by)))]

# unique(sample_metadata[,c("genotype","sample")])

######################################################
## Calculate pseudovulk stats and do some filtering ##
######################################################

pseudobulk_stats.dt <- sample_metadata[,.N,by=c("genotype","sample","celltype","genotype_sample_celltype")]

# Filter each instance by minimum number of cells
pseudobulk_stats.dt <- pseudobulk_stats.dt[N>=30]

# Select celltypes that are measured in at least 3 WT samples
celltypes.to.use <- pseudobulk_stats.dt[genotype=="WT",.N,by=c("celltype")] %>% .[N>=3,celltype] 
pseudobulk_stats.dt <- pseudobulk_stats.dt[celltype%in%celltypes.to.use]

# For each genotype and celltype combination, require at least 3 samples
tmp <- pseudobulk_stats.dt[,.N,by=c("celltype","genotype")] %>% .[N>=3] %>% .[,N:=NULL]
pseudobulk_stats.dt <- pseudobulk_stats.dt %>% merge(tmp,by=c("celltype","genotype"))

# Update metadata
sample_metadata <- sample_metadata[genotype_sample_celltype%in%pseudobulk_stats.dt$genotype_sample_celltype]

##############################
## Load RNA expression data ##
##############################

sce <- load_SingleCellExperiment(io$sce, cells=sample_metadata$cell)
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

################
## Pseudobulk ##
################

sce_pseudobulk <- pseudobulk_sce_fn(
  x = sce,
  assay = "counts",
  by = opts$group_by,
  fun = "sum",
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
)

assayNames(sce_pseudobulk) <- "counts"

###############
## Parse SCE ##
###############

# Add metadata
# sce_pseudobulk$class <- stringr::str_split(colnames(sce_pseudobulk), pattern = "-") %>% map_chr(1)
sce_pseudobulk$genotype <- stringr::str_split(colnames(sce_pseudobulk), pattern = "-") %>% map_chr(1)
sce_pseudobulk$sample <- stringr::str_split(colnames(sce_pseudobulk), pattern = "-") %>% map_chr(2)
sce_pseudobulk$celltype <- stringr::str_split(colnames(sce_pseudobulk), pattern = "-") %>% map_chr(3)

# Filter genes
genes <- fread("/Users/argelagr/data/shiny_dnmt_tet/genes.txt", header = F)[[1]]
sce_pseudobulk <- sce_pseudobulk[genes,]

# Sanity checks
# stopifnot(unique(sce_pseudobulk$class)%in%opts$classes)
stopifnot(unique(sce_pseudobulk$sample)%in%opts$samples)
stopifnot(unique(sce_pseudobulk$celltype)%in%c(opts$celltypes,c("Erythroid", "Blood_progenitors")))

###################
## Normalisation ##
###################

if (opts$normalisation_method=="deseq2") {
  
  suppressPackageStartupMessages(library(DESeq2))
  dds <- DESeqDataSet(sce_pseudobulk, design=~1)
  dds <- varianceStabilizingTransformation(dds)
  logcounts(sce_pseudobulk) <- assay(dds)
  
} else if (opts$normalisation_method=="cpm") {
  
  logcounts(sce_pseudobulk) <- log2(1e6*(sweep(counts(sce_pseudobulk),2,colSums(counts(sce_pseudobulk)),"/"))+1)
  # logcounts(sce_pseudobulk) <- edgeR::cpm(counts(sce_pseudobulk), log=TRUE, prior.count = 1)
  
} else {
  stop("Normalisation method not recognised")
}

# Remove counts assay
assays(sce_pseudobulk)["counts"] <- NULL

# Save
saveRDS(sce_pseudobulk, file.path(io$outdir,"SingleCellExperiment_shiny.rds"))

#####################
## Save statistics ##
#####################

# stats.dt <- data.table(table(sce[[opts$group_by]])) %>% setnames(c("sample","N"))
fwrite(pseudobulk_stats.dt, file.path(io$outdir,"pseudobulk_stats_shiny.txt.gz"), sep="\t")
