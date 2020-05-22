library(ggpubr)

##############
## Settings ##
##############

source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
io$outdir <- paste0(io$basedir,"/results/celltype_affinity"); dir.create(io$outdir, showWarnings = F)

opts$scale <- FALSE

opts$classes <- c(
  # "E7.5_Host", 
  # "E7.5_TET_TKO", 
  "E8.5_Host", 
  "E8.5_TET_TKO"
  # "E10.5_Host", 
  # "E10.5_TET_TKO"
)

###############
## Load data ##
###############

# Load metadata
sample_metadata <- sample_metadata %>% 
  .[class%in%opts$classes] %>%
  .[!is.na(celltype.mapped)] %>%
  .[,c("cell","class","celltype.mapped")]

# Load gene markers
marker_genes.dt <- fread(io$atlas.marker_genes)

# Load SingleCellExperiment
sce <- readRDS(io$sce)[,sample_metadata$cell]

# Add relevant information to the colData slot
sce$celltype.mapped <- sample_metadata$celltype.mapped
sce$class <- sample_metadata$class

################
## Parse data ##
################

# Subset genes
sce <- sce[unique(marker_genes.dt$ens_id),]

# if (isTRUE(opts$scale)) {
#   stop("TO-DO")
# } else {
#   expr.matrix <- as.matrix(logcounts(sce))
# }


##################
## Computations ##
##################

# Calculate atlas-based cell type "score" for each predicted cell type in the query
dt <- unique(sce$celltype.mapped) %>% map(function(i) {
# dt <- unique(sce$celltype.mapped) %>% head(n=5) %>% map(function(i) {
  opts$classes %>% map(function(j) {
    logcounts(sce[,sce$celltype.mapped==i & sce$class==j]) %>% 
      rowMeans %>% as.data.table(keep.rownames = T) %>% 
      setnames(c("ens_id","expr")) %>%
      merge(marker_genes.dt[,c("celltype","ens_id")], by="ens_id", allow.cartesian=TRUE) %>%
      .[,.(score=mean(expr)),by="celltype"] %>%
      .[,c("celltype.mapped","class"):=list(i,j)]
  }) %>% rbindlist 
}) %>% rbindlist

# TO-DO Calculate cell type "score" for each single cell

##########
## Plot ##
##########

colors <- opts$celltype.colors
names(colors) <- names(colors) %>% stringr::str_replace_all("_"," ")

# Plot number of marker genes per cell types
for (i in unique(dt$celltype.mapped)) {
  
  to.plot <- dt[celltype.mapped==i] %>%
    .[,celltype:=stringr::str_replace_all(celltype,"_"," ")] %>%
    .[,celltype:=factor(celltype,levels=names(colors))]
  
  p <- ggbarplot(to.plot, x="celltype", y="score", fill="celltype") +
    facet_wrap(~class, ncol=1) +
    scale_fill_manual(values=colors) +
    labs(x="", y="Cell type affinity") +
    theme(
      axis.text.y = element_text(size=rel(0.75)),
      axis.text.x = element_text(colour="black",size=rel(0.8), angle=90, hjust=1, vjust=0.5),
      axis.ticks.x = element_blank(),
      legend.position = "none"
  )
  
  pdf(sprintf("%s/%s_affinity.pdf",io$outdir,i), width = 9, height = 5)
  print(p)
  dev.off()
}

##########
## Save ##
##########

# fwrite(dt.filt, paste0(io$outdir,"/marker_genes.txt.gz"))
