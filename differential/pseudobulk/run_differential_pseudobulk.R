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

####################################################
## Differential expression per class and celltype ##
####################################################

# Subset WT samples from the SingleCellExperiment
sce.wt <- sce[,sce$class==opts$wt.class]
metadata(sce.wt)$n_cells <- metadata(sce.wt)$n_cells[grep(opts$wt.class, names(metadata(sce.wt)$n_cells), value=T)]
names(metadata(sce.wt)$n_cells) <- gsub(paste0(opts$wt.class,"-"),"",names(metadata(sce.wt)$n_cells))
colnames(sce.wt) <- sce.wt$celltype

# Select celltypes with sufficient number of cells
celltypes.wt <- names(which(metadata(sce.wt)$n_cells >= opts$min.cells))

# i <- "E8.5_Dnmt3aKO_Dnmt3bKO"
for (i in opts$ko.classes) {
  
  # Subset KO samples from the SingleCellExperiment
  sce.ko <- sce[,sce$class==i]
  metadata(sce.ko)$n_cells <- metadata(sce.ko)$n_cells[grep(i, names(metadata(sce.ko)$n_cells), value=T)]
  names(metadata(sce.ko)$n_cells) <- gsub(paste0(i,"-"),"",names(metadata(sce.ko)$n_cells))
  colnames(sce.ko) <- sce.ko$celltype
  
  # Select celltypes with sufficient number of cells
  celltypes.ko <- names(which(metadata(sce.ko)$n_cells >= opts$min.cells))
  celltypes.to.use <- intersect(celltypes.wt,celltypes.ko)
  
  for (j in celltypes.to.use) {
      
    tmp <- data.table(
      gene = rownames(sce.ko),
      expr_ko = logcounts(sce.ko[,j])[,1] %>% round(2),
      expr_wt = logcounts(sce.ko[,j])[,1] %>% round(2)
    ) %>% .[,diff:=round(expr_ko-expr_wt,2)] %>% sort.abs("diff") 

    # save      
    outfile <- sprintf("%s/%s/%s_%s_vs_%s.txt.gz", io$outdir,i,j,opts$wt.class,i); dir.create(dirname(outfile), showWarnings = F)
    fwrite(tmp, outfile, sep="\t")
  }
}

#######################################
## Differential expression per class ##
#######################################

# Subset WT samples from the SingleCellExperiment
sce.wt <- sce[,sce$class==opts$wt.class]

for (i in opts$ko.classes) {
  
  # Subset KO samples from the SingleCellExperiment
  sce.ko <- sce[,sce$class==i]
  
  tmp <- data.table(
    gene = rownames(sce.ko),
    expr_ko = logcounts(sce.ko)[,1] %>% round(2),
    expr_wt = logcounts(sce.wt)[,1] %>% round(2)
  ) %>% .[,diff:=round(expr_ko-expr_wt,2)] %>% sort.abs("diff") 
    
  # save      
  outfile <- sprintf("%s/%s/%s_vs_%s.txt.gz", io$outdir,i,opts$wt.class,i); dir.create(dirname(outfile), showWarnings = F)
  fwrite(tmp, outfile, sep="\t")
}