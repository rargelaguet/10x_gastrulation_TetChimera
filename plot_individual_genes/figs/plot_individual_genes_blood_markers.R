here::i_am("plot_individual_genes/plot_individual_genes.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O ##
io$outdir <- file.path(io$basedir,"results_new/individual_genes/figs"); dir.create(io$outdir, showWarnings = F)

## Define options ##

# Define cell types to plot
opts$celltypes = c(
	"Epiblast",
	"Primitive_Streak",
	"Caudal_epiblast",
	# "PGC",
	"Anterior_Primitive_Streak",
	"Notochord",
	"Def._endoderm",
	"Gut",
	"Nascent_mesoderm",
	"Mixed_mesoderm",
	"Intermediate_mesoderm",
	"Caudal_Mesoderm",
	"Paraxial_mesoderm",
	"Somitic_mesoderm",
	"Pharyngeal_mesoderm",
	"Cardiomyocytes",
	"Allantois",
	"ExE_mesoderm",
	"Mesenchyme",
	"Haematoendothelial_progenitors",
	"Endothelium",
	# "Blood_progenitors_1",
	# "Blood_progenitors_2",
	"Blood_progenitors",
	# "Erythroid1",
	# "Erythroid2",
	# "Erythroid3",
	"Erythroid",
	"NMP",
	"Rostral_neurectoderm",
	"Caudal_neurectoderm",
	"Neural_crest",
	"Forebrain_Midbrain_Hindbrain",
	"Spinal_cord",
	"Surface_ectoderm"
	# "Visceral_endoderm",
	# "ExE_endoderm",
	# "ExE_ectoderm",
	# "Parietal_endoderm"
)

# Define classes to plot
opts$stage_classes <- c(
  "E7.5_WT", 
  "E7.5_TET_TKO", 
  "E8.5_WT",
  "E8.5_TET_TKO"
)

# Define colors
opts$stage_class_colors <- c(
  # "E7.5_WT" = "#CCCCCC", 
  "E7.5_WT" = "#5CACEE", 
  "E7.5_TET_TKO" = "#FF7F50", 
  "E8.5_WT" = "#436EEE", 
  # "E8.5_WT" = "#B0B0B0", 
  "E8.5_TET_TKO" = "#EE4000"
)

opts$min.cells <- 25

opts$rename_celltypes <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors"
  # "Intermediate_mesoderm" = "Mixed_mesoderm",
  # "Paraxial_mesoderm" = "Mixed_mesoderm",
  # "Nascent_mesoderm" = "Mixed_mesoderm",
  # "Pharyngeal_mesoderm" = "Mixed_mesoderm"
  # "Visceral_endoderm" = "ExE_endoderm"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE] %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$rename_celltypes)] %>% .[celltype.mapped%in%opts$celltypes] %>% .[,celltype.mapped:=factor(celltype.mapped, levels=opts$celltypes)] %>%
  .[,stage_class:=sprintf("%s_%s",stage,class)] %>% .[stage_class%in%opts$stage_classes] %>% .[,stage_class:=factor(stage_class,levels=opts$stage_classes)]
  
table(sample_metadata$stage_class)
table(sample_metadata$sample)
table(sample_metadata$celltype.mapped)

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

##############################
## Define plotting function ##
##############################

plot_fn <- function(to.plot, to.plot.subset, gene.to.plot) {
  p <- ggplot(to.plot, aes(x=stage_class, y=expr, fill=stage_class)) +
    geom_violin(scale = "width", alpha=0.75) +
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.75) +
    geom_jitter(width=0.05, alpha=0.25, size=0.5, shape=21, data=to.plot.subset) +
    stat_compare_means(comparisons = list(c("E8.5_WT", "E8.5_TET_TKO")), aes(label = paste0("p = ", ..p.format..)), size=3, method="t.test") +
    # stat_summary(fun.data = give.n, geom = "text", size=3) +
    scale_fill_manual(values=opts$stage_class_colors) +
    facet_wrap(~celltype, scales="fixed", nrow=1) +
    theme_classic() +
    labs(x="",y=sprintf("%s expression",gene.to.plot)) +
    guides(x = guide_axis(angle = 90)) +
    theme(
      strip.text = element_text(size=rel(0.85)),
      axis.text.x = element_text(colour="black",size=rel(0.75)),
      axis.text.y = element_text(colour="black",size=rel(0.9)),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(colour="black",size=rel(1.0)),
      legend.position = "none",
      legend.title = element_blank(),
      legend.text = element_text(size=rel(0.85))
    )
  
  pdf(sprintf("%s/%s_boxplots_single_cells.pdf",io$outdir,gene.to.plot), width=8, height=4)
  print(p)
  dev.off()
}

##########
## Plot ##
##########

gene.to.plot <- "Hba-x"
celltypes.to.plot <- c("Haematoendothelial_progenitors", "Endothelium", "Blood_progenitors", "Erythroid")

to.plot <- data.table(
  cell = colnames(sce),
  expr = logcounts(sce)[gene.to.plot,],
  class = sce$class,
  stage_class = sce$stage_class,
  celltype = sce$celltype.mapped
) %>% .[celltype%in%celltypes.to.plot] %>% .[,N:=.N,by=c("class","celltype")] %>% .[N>=opts$min.cells]
    
# to.plot[expr==0,expr:=0.1]
to.plot.subset <- to.plot %>% .[sample.int(n=nrow(.), size=nrow(.)/5)] %>% .[expr>0]

plot_fn(to.plot, to.plot.subset, gene.to.plot)

##########
## Plot ##
##########

gene.to.plot <- "Klf1"
celltypes.to.plot <- c("Haematoendothelial_progenitors", "Endothelium", "Blood_progenitors", "Erythroid")

to.plot <- data.table(
  cell = colnames(sce),
  expr = logcounts(sce)[gene.to.plot,],
  class = sce$class,
  stage_class = sce$stage_class,
  celltype = sce$celltype.mapped
) %>% .[celltype%in%celltypes.to.plot] %>% .[,N:=.N,by=c("class","celltype")] %>% .[N>=opts$min.cells]

# to.plot[expr==0,expr:=0.1]
to.plot[expr>4,expr:=4]
to.plot.subset <- to.plot %>% .[expr>0] %>% .[sample.int(n=nrow(.), size=nrow(.)/2)]

plot_fn(to.plot, to.plot.subset, gene.to.plot)

##########
## Plot ##
##########

gene.to.plot <- "Lefty2"
celltypes.to.plot <- c("Primitive_Streak","Nascent_mesoderm", "Mixed_mesoderm", "Somitic_mesoderm")

to.plot <- data.table(
  cell = colnames(sce),
  expr = logcounts(sce)[gene.to.plot,],
  class = sce$class,
  stage_class = sce$stage_class,
  celltype = sce$celltype.mapped
) %>% .[celltype%in%celltypes.to.plot] %>% .[,N:=.N,by=c("class","celltype")] %>% .[N>=opts$min.cells]

# to.plot[expr==0,expr:=0.1]
to.plot[expr>5,expr:=5]
to.plot.subset <- to.plot %>% .[expr>0] %>% .[sample.int(n=nrow(.), size=nrow(.)/2)]

plot_fn(to.plot, to.plot.subset, gene.to.plot)


##########
## Plot ##
##########

gene.to.plot <- "Tspan33"
celltypes.to.plot <- c("Haematoendothelial_progenitors", "Endothelium", "Blood_progenitors", "Erythroid")

to.plot <- data.table(
  cell = colnames(sce),
  expr = logcounts(sce)[gene.to.plot,],
  class = sce$class,
  stage_class = sce$stage_class,
  celltype = sce$celltype.mapped
) %>% .[celltype%in%celltypes.to.plot] %>% .[,N:=.N,by=c("class","celltype")] %>% .[N>=opts$min.cells]

# to.plot[expr==0,expr:=0.1]
to.plot[expr>4,expr:=4]
to.plot.subset <- to.plot %>% .[expr>0] %>% .[sample.int(n=nrow(.), size=nrow(.)/2)]

plot_fn(to.plot, to.plot.subset, gene.to.plot)


##########
## Plot ##
##########

gene.to.plot <- "Fgf3"
celltypes.to.plot <- c("ExE_mesoderm","Haematoendothelial_progenitors", "Blood_progenitors", "Erythroid")

to.plot <- data.table(
  cell = colnames(sce),
  expr = logcounts(sce)[gene.to.plot,],
  class = sce$class,
  stage_class = sce$stage_class,
  celltype = sce$celltype.mapped
) %>% .[celltype%in%celltypes.to.plot] %>% .[,N:=.N,by=c("class","celltype")] %>% .[N>=opts$min.cells]

# to.plot[expr==0,expr:=0.1]
to.plot[expr>5,expr:=5]
to.plot.subset <- to.plot %>% .[expr>0] %>% .[sample.int(n=nrow(.), size=nrow(.)/2)]

plot_fn(to.plot, to.plot.subset, gene.to.plot)



##########
## Plot ##
##########

gene.to.plot <- "Fgf8"
celltypes.to.plot <- c("Intermediate_mesoderm","ExE_mesoderm", "Mixed_mesoderm", "Somitic_mesoderm")

to.plot <- data.table(
  cell = colnames(sce),
  expr = logcounts(sce)[gene.to.plot,],
  class = sce$class,
  stage_class = sce$stage_class,
  celltype = sce$celltype.mapped
) %>% .[celltype%in%celltypes.to.plot] %>% .[,N:=.N,by=c("class","celltype")] %>% .[N>=opts$min.cells]

# to.plot[expr==0,expr:=0.1]
to.plot[expr>4.75,expr:=4.75]
to.plot.subset <- to.plot %>% .[expr>0] %>% .[sample.int(n=nrow(.), size=nrow(.)/2)]

plot_fn(to.plot, to.plot.subset, gene.to.plot)


