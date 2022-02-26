here::i_am("plot_individual_genes/plot_individual_genes.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O ##
io$outdir <- file.path(io$basedir,"results_all/individual_genes/figs"); dir.create(io$outdir, showWarnings = F)

# Define classes to plot
opts$stage_classes <- c(
  "E7.5_WT", 
  "E7.5_TET_TKO", 
  "E8.5_WT",
  "E8.5_TET_TKO"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE] %>%
  .[,stage_class:=sprintf("%s_%s",stage,class2)] %>% .[stage_class%in%opts$stage_classes] %>% .[,stage_class:=factor(stage_class,levels=opts$stage_classes)]
  
table(sample_metadata$stage_class)

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=sample_metadata$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

##########
## Plot ##
##########

gene <- "tomato-td"

rna.dt <- data.table(
  cell = colnames(sce),
  expr = counts(sce)[gene,],
  class = factor(sce$class, levels=opts$classes),
  stage = sce$stage,
  stage_class = factor(sce$stage_class, levels=opts$stage_classes)
)

to.plot <- rna.dt[,mean(expr>0),by=c("stage","class")]

p <- ggplot(to.plot, aes(x=class, y=V1, fill=class)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values=opts$class_colors) +
  facet_wrap(~stage) +
  theme_classic() +
  labs(x="",y=sprintf("%s expression",gene)) +
  # guides(x = guide_axis(angle = 90)) +
  theme(
    strip.text = element_text(size=rel(0.85)),
    strip.background = element_blank(),
    axis.text.x = element_text(colour="black",size=rel(0.75)),
    # axis.text.x = element_blank(),
    axis.text.y = element_text(colour="black",size=rel(0.9)),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(colour="black",size=rel(1.0)),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size=rel(0.85))
  )

pdf(sprintf("%s/tdTomato_boxplots_single_cells.pdf",io$outdir), width=8.5, height=3)
print(p)
dev.off()
