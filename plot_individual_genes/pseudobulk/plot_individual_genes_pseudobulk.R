here::i_am("plot_individual_genes/pseudobulk/plot_individual_genes_pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outdir <- file.path(io$basedir,"results/individual_genes/pseudobulk")
io$sce.pseudobulk <- file.path(io$basedir,"results/pseudobulk/SingleCellExperiment_pseudobulk_stage_class2_celltype.rds")

# Define options
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
  # "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm"
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)
opts$min.cells <- 30

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- readRDS(io$sce.pseudobulk)

# Add metadata
sce$stage <- stringr::str_split(colnames(sce), pattern = "/") %>% map_chr(1)
sce$class <- stringr::str_split(colnames(sce), pattern = "/") %>% map_chr(2)
sce$celltype <- stringr::str_split(colnames(sce), pattern = "/") %>% map_chr(3)
sce$stage_class <- paste(sce$stage,sce$class,sep="_")

# Filter celltypes
# sce <- sce[,sce$celltype%in%opts$celltypes]

# Filter by minimum number of cells
sce <- sce[,names(which(metadata(sce)$n_cells>=opts$min.cells))]

# Remove NAs
sce <- sce[,!grepl("NA",colnames(sce))]

sce <- sce[,sce$celltype%in%names(which(table(sce$celltype)==2))]

#############
## Barplots ##
#############

# genes.to.plot <- c("Tex19.1","Morc1","Dppa3","Rex1","Dppa5a","Dppa4","Dppa2","Zfp981")
# genes.to.plot <- c("Pou5f1","Epcam","Fst","Lefty2","Cdkn1c","Acta1")
# genes.to.plot <- grep("^Hox",rownames(sce),value=T)
genes.to.plot <- fread(io$atlas.marker_genes)[,gene] %>% unique %>% head(n=3)
genes.to.plot <- genes.to.plot[genes.to.plot%in%rownames(sce)]

for (i in 1:length(genes.to.plot)) {
  
  gene <- genes.to.plot[i]
  
  print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
  outfile <- sprintf("%s/%s_barplot_pseudobulk.pdf",io$outdir,gene)

  to.plot <- data.table(
    sample = colnames(sce),
    expr = logcounts(sce)[gene,],
    class = sce$class,
    stage = sce$stage,
    stage_class = sce$stage_class,
    celltype = sce$celltype
  )
  
  p <- ggplot(to.plot, aes(x=stage_class, y=expr, fill=stage_class)) +
    geom_bar(stat="identity", color="black") +
    scale_fill_manual(values=opts$stage_class2_colors) +
    facet_wrap(~celltype, scales="fixed") +
    theme_classic() +
    labs(x="",y=sprintf("%s expression",gene)) +
    # guides(x = guide_axis(angle = 90)) +
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
}

