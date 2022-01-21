here::i_am("plot_individual_genes/plot_individual_genes.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O ##
io$metadata <- file.path(io$basedir,"results_all/mapping/sample_metadata_after_mapping.txt.gz")
io$sce <- file.path(io$basedir,"processed_all/SingleCellExperiment.rds")
io$outdir <- file.path(io$basedir,"results_all/individual_genes"); dir.create(io$outdir, showWarnings = F)

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
	"Blood_progenitors_1",
	"Blood_progenitors_2",
	"Erythroid1",
	"Erythroid2",
	"Erythroid3",
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
# opts$stage_classes <- c(
#   "E7.5_WT", 
#   "E7.5_TET_TKO", 
#   "E8.5_WT",
#   "E8.5_TET_TKO"
#   # "E9.5_TET_TKO"
# )

opts$stage_classes <- c(
  "E7.5_WT_tdTomato+", 
  "E7.5_WT_tdTomato-", 
  "E7.5_TET_TKO", 
  "E8.5_WT_tdTomato+", 
  "E8.5_WT_tdTomato-", 
  "E8.5_TET_TKO"
  # "E9.5_TET_TKO"
)



# Define colors
# opts$stage_class_colors <- c(
#   "E7.5_WT" = "#CCCCCC", 
#   "E7.5_TET_TKO" = "#FF7F50", 
#   "E8.5_WT" = "#B0B0B0", 
#   "E8.5_TET_TKO" = "#EE4000",
#   "E8.5_WT" = "#B0B0B0"
#   # "E9.5_TET_TKO" = "#B22222"
# )

opts$stage_class_colors <- c(
  "E7.5_WT_tdTomato++" = "#CCCCCC", 
  "E7.5_WT_tdTomato-" = "#CCCCCC", 
  "E7.5_TET_TKO"  = "#FF7F50", 
  "E8.5_WT_tdTomato++" = "#B0B0B0", 
  "E8.5_WT_tdTomato-" = "#B0B0B0", 
  "E8.5_TET_TKO" = "#EE4000"
  # "E9.5_TET_TKO"
)

opts$min.cells <- 25

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_rnaQC==TRUE & celltype.mapped%in%opts$celltypes] %>%
  .[,stage_class:=sprintf("%s_%s",stage,class)] %>% .[stage_class%in%opts$stage_classes] %>% .[,stage_class:=factor(stage_class,levels=opts$stage_classes)] %>%
  .[,celltype.mapped:=factor(celltype.mapped, levels=opts$celltypes)]

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

##########
## Plot ##
##########

genes.to.plot <- c("Hba-x","Dppa4","Fgf3","Fgf8","Klf1","Klf2","Kdr","Cd34","Lefty2","Hemgn","Slc4a1","Mrap","Tspan33")
# genes.to.plot <- c("Lefty1","Cd34","Tmsb4x","Fgf3","Spata7","Cer1","Spink1","Dppa4","Dppa5a","Prc1","Lefty2","Ube2c","Hba-x","Hbb-y","Hba-a1","Hbb-bh1")
# genes.to.plot <- c("Vegfa","Vegfb","Vegfc","Vegfd","Kdr","Flt1","Tal1","Runx1","Etv2)
# genes.to.plot <- c("Tet1","Tet2","Tet3","Dnmt1","Dnmt3a","Dnmt3b","Dnmt3l")
# genes.to.plot <- rownames(sce)[grep("tomato",rownames(sce))]
# genes.to.plot <- fread(io$atlas.marker_genes)$gene %>% unique %>% .[!grepl("Rik$",.)]
# genes.to.plot <- fread(io$atlas.marker_genes) %>% .[grep("Erythroid",celltype),gene] %>% unique
# genes.to.plot <- fread("/Users/ricard/data/gastrulation10x/results/differential/celltypes/E8.5/Neural_crest_vs_Forebrain_Midbrain_Hindbrain.txt.gz") %>% .[sig==T & logFC<(-2),gene]

for (i in 1:length(genes.to.plot)) {
  
  gene <- genes.to.plot[i]
  
  if (gene %in% rownames(sce)) {
    print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
    outfile <- sprintf("%s/%s_boxplots_single_cells.pdf",io$outdir,gene)
    
    to.plot <- data.table(
      cell = colnames(sce),
      expr = logcounts(sce)[gene,],
      class = sce$class,
      stage_class = sce$stage_class,
      celltype = sce$celltype.mapped
    ) %>% .[,N:=.N,by=c("class","celltype")] %>% .[N>=opts$min.cells]
    
  p <- ggplot(to.plot, aes(x=stage_class, y=expr, fill=stage_class)) +
      geom_violin(scale = "width", alpha=0.8) +
      geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
      stat_summary(fun.data = give.n, geom = "text", size=3) +
      # geom_jitter(size=2, shape=21, stroke=0.2, alpha=0.5) +
      # scale_fill_manual(values=opts$colors) +
      # scale_fill_brewer(palette="Dark2") +
      scale_fill_manual(values=opts$stage_class_colors) +
      facet_wrap(~celltype, scales="fixed") +
      theme_classic() +
      labs(x="",y=sprintf("%s expression",gene)) +
      guides(x = guide_axis(angle = 90)) +
      theme(
        strip.text = element_text(size=rel(0.85)),
        axis.text.x = element_text(colour="black",size=rel(0.9)),
        axis.text.y = element_text(colour="black",size=rel(0.9)),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(colour="black",size=rel(1.0)),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=rel(0.85))
      )
    
    pdf(outfile, width=10, height=9)
    print(p)
    dev.off()

  } else {
    print(sprintf("%s not found",gene))
  }
}


##########
## TEST ##
##########

# # genes.to.plot <- c("Lefty1","Cd34","Tmsb4x","Fgf3","Spata7","Cer1","Spink1","Dppa4","Dppa5a","Prc1","Lefty2","Ube2c","Hba-x","Hbb-y","Hba-a1","Hbb-bh1")
# # genes.to.plot <- c("Vegfa","Vegfb","Vegfc","Vegfd","Kdr","Flt1","Tal1","Runx1","Etv2)
# # genes.to.plot <- c("Tet1","Tet2","Tet3","Dnmt1","Dnmt3a","Dnmt3b","Dnmt3l")
# # genes.to.plot <- rownames(sce)[grep("tomato",rownames(sce))]
# # genes.to.plot <- fread(io$atlas.marker_genes)$gene %>% unique %>% .[!grepl("Rik$",.)]
# # genes.to.plot <- fread(io$atlas.marker_genes) %>% .[grep("Erythroid",celltype),gene] %>% unique
# # genes.to.plot <- fread("/Users/ricard/data/gastrulation10x/results/differential/celltypes/E8.5/Neural_crest_vs_Forebrain_Midbrain_Hindbrain.txt.gz") %>% .[sig==T & logFC<(-2),gene]
# 
# # gene <- "Cd34"
# for (gene in genes.to.plot) {
#   print(gene)
#   outfile <- sprintf("%s/%s_boxplots_single_cells.pdf",io$outdir,gene)
#   
#   to.plot <- data.table(
#     cell = colnames(sce),
#     expr = logcounts(sce)[gene,],
#     class = sce$class,
#     sample = sce$sample,
#     stage_class = sce$stage_class,
#     celltype = sce$celltype.mapped
#   ) %>% .[,N:=.N,by=c("class","celltype")] %>% .[N>=opts$min.cells]
#   
#   p <- ggplot(to.plot, aes(x=sample, y=expr, fill=stage_class)) +
#     geom_violin(scale = "width", alpha=0.8) +
#     geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
#     # stat_summary(fun.data = give.n, geom = "text", size=3) +
#     # geom_jitter(size=2, shape=21, stroke=0.2, alpha=0.5) +
#     # scale_fill_manual(values=opts$colors) +
#     # scale_fill_brewer(palette="Dark2") +
#     scale_fill_manual(values=opts$stage_class_colors) +
#     facet_wrap(~celltype, scales="fixed") +
#     theme_classic() +
#     labs(x="",y=sprintf("%s expression",gene)) +
#     guides(x = guide_axis(angle = 90)) +
#     theme(
#       strip.text = element_text(size=rel(0.85)),
#       axis.text.x = element_text(colour="black",size=rel(0.7)),
#       axis.text.y = element_text(colour="black",size=rel(0.9)),
#       axis.ticks.x = element_blank(),
#       axis.title.y = element_text(colour="black",size=rel(1.0)),
#       legend.position = "none",
#       legend.title = element_blank(),
#       legend.text = element_text(size=rel(0.85))
#     )
#   
#   pdf(outfile, width=10, height=9)
#   print(p)
#   dev.off()
# }

