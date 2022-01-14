source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
# io$sce.pseudobulk <- file.path(io$basedir,"results_new/pseudobulk/SingleCellExperiment_pseudobulk_class_celltype.rds")
io$sce.pseudobulk <- file.path(io$basedir,"results_new/pseudobulk/SingleCellExperiment_pseudobulk_class_celltype.rds")
io$outdir <- file.path(io$basedir,"results_new/differential/pseudobulk"); dir.create(io$outdir, showWarnings = F)

# Options
opts$min.cells <- 25

opts$ko.classes <- c(
  "E7.5_TET_TKO", 
  "E8.5_TET_TKO",
  "E9.5_TET_TKO"
)

# opts$wt.classes <- c("E7.5_WT","E8.5_WT")
opts$wt.class <- "WT"

##############################
## Load pseudobulk RNA data ##
##############################

sce <- readRDS(io$sce.pseudobulk)

# Add metadata
sce$class <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(1)
sce$celltype <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(2)

# Filter by minimum number of cells
sce <- sce[,names(which(metadata(sce)$n_cells>=opts$min.cells))]

######################
## Subset WT and KO ##
######################
stop("FIX THIS CODE")

# Subset WT samples from the SingleCellExperiment
sce.wt <- sce[,sce$class=="WT"]
metadata(sce.wt)$n_cells <- metadata(sce.wt)$n_cells[grep(opts$wt.class, names(metadata(sce.wt)$n_cells), value=T)]
names(metadata(sce.wt)$n_cells) <- gsub(paste0(opts$wt.class,"-"),"",names(metadata(sce.wt)$n_cells))
colnames(sce.wt) <- sce.wt$celltype

# Select celltypes with sufficient number of cells
celltypes.wt <- names(which(metadata(sce.wt)$n_cells >= opts$min.cells))


# Subset KO samples from the SingleCellExperiment
sce.ko <- sce[,sce$class=="TET_TKO"]
metadata(sce.ko)$n_cells <- metadata(sce.ko)$n_cells[grep("TET_TKO", names(metadata(sce.ko)$n_cells), value=T)]
names(metadata(sce.ko)$n_cells) <- gsub(paste0("TET_TKO","-"),"",names(metadata(sce.ko)$n_cells))
colnames(sce.ko) <- sce.ko$celltype

# Select celltypes with sufficient number of cells
celltypes.ko <- names(which(metadata(sce.ko)$n_cells >= opts$min.cells))

####################################################
## Differential expression per class and celltype ##
####################################################

celltypes.to.use <- intersect(celltypes.wt,celltypes.ko)

for (j in celltypes.to.use) {
    
  tmp <- data.table(
    gene = rownames(sce.ko),
    expr_ko = logcounts(sce.ko[,j])[,1] %>% round(2),
    expr_wt = logcounts(sce.wt[,j])[,1] %>% round(2)
  ) %>% .[,diff:=round(expr_ko-expr_wt,2)] %>% sort.abs("diff") 

  # save      
  fwrite(tmp, sprintf("%s/%s_WT_vs_TET_TKO.txt.gz",io$outdir,j), sep="\t")
}
