here::i_am("plot_individual_genes/pseudobulk/plot_individual_genes_pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results_all/individual_genes/pseudobulk/tdTomato"); dir.create(io$outdir, showWarnings = F, recursive = T)
io$sce.pseudobulk <- file.path(io$basedir,"results_all/pseudobulk/SingleCellExperiment_pseudobulk_stage_class2_celltype.rds")

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- readRDS(io$sce.pseudobulk)

# Parse metadata
sce$stage <- stringr::str_split(colnames(sce), pattern = "/") %>% map_chr(1)
sce$class <- stringr::str_split(colnames(sce), pattern = "/") %>% map_chr(2)
sce$celltype <- stringr::str_split(colnames(sce), pattern = "/") %>% map_chr(3)
sce$stage_class <- paste(sce$stage,sce$class,sep="_")

######################
## Plot Hemoglobins ##
######################

gene <- "Hba-x"

to.plot <- data.table(
  sample = colnames(sce),
  expr = logcounts(sce)[gene,],
  stage = sce$stage,
  stage_class = sce$stage_class,
  class = sce$class,
  celltype = sce$celltype
) %>% .[,class:=factor(class, levels=opts$classes)]

p <- ggplot(to.plot, aes(x=class, y=expr, fill=class)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values=opts$class_colors) +
  facet_wrap(~stage, scales="fixed") +
  theme_classic() +
  labs(x="",y=sprintf("%s expression",gene)) +
  # guides(x = guide_axis(angle = 90)) +
  theme(
    strip.text = element_text(size=rel(0.85)),
    axis.text.x = element_text(colour="black",size=rel(0.9)),
    axis.text.y = element_text(colour="black",size=rel(0.9)),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(colour="black",size=rel(1.0)),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size=rel(0.85))
  )

pdf(sprintf("%s/%s_barplot_pseudobulk.pdf",io$outdir,gene), width=8, height=5)
print(p)
dev.off()
    
