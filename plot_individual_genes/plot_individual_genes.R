suppressPackageStartupMessages(library(SingleCellExperiment))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/10x_gastrulation_TetChimera/settings.R")
}

## I/O ##

io$outdir <- paste0(io$basedir,"/results/individual_genes")

## Define options ##

# Define cell types to plot
opts$celltypes = c(
	"Epiblast",
	"Primitive_Streak",
	"Caudal_epiblast",
	# "PGC",
	# "Anterior_Primitive_Streak",
	# "Notochord",
	"Def._endoderm",
	"Gut",
	"Nascent_mesoderm",
	"Mixed_mesoderm",
	"Intermediate_mesoderm",
	"Caudal_Mesoderm",
	"Paraxial_mesoderm",
	"Somitic_mesoderm",
	"Pharyngeal_mesoderm",
	# "Cardiomyocytes",
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
# opts$celltypes <- c("Blood_progenitors_1")

# Define classes to plot
opts$classes <- c(
  "E7.5_Host",
  "E7.5_TET_TKO",
  "E8.5_Host",
  "E8.5_TET_TKO"
  # "E10.5_Host", 
  # "E10.5_TET_TKO"
)

# Define colors
opts$colors <- c(
  "E7.5_Host" = "#CCCCCC", 
  "E7.5_TET_TKO" = "#FF7F50", 
  "E8.5_Host" = "#B0B0B0", 
  "E8.5_TET_TKO" = "#EE4000"
)

# Update sample metadata
sample_metadata <- fread(io$metadata) %>% 
  # .[class%in%opts$classes & celltype.mapped%in%opts$celltypes] %>%
  .[celltype.mapped%in%opts$celltypes] %>%
  .[,celltype.mapped:=factor(celltype.mapped, levels=opts$celltypes)]

table(sample_metadata$class)
table(sample_metadata$celltype.mapped)

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- readRDS(io$sce)[,sample_metadata$cell]

# Remove genes that are not expressed
sce <- sce[rowMeans(counts(sce))>0,]

# Load gene metadata
gene_metadata <- fread(io$gene_metadata) %>%
  .[symbol!="" & symbol%in%rownames(sce)]

################
## Parse data ##
################

# Rename genes
# new.names <- gene_metadata$symbol
# names(new.names) <- gene_metadata$ens_id
# sce <- sce[rownames(sce) %in% names(new.names),]
# rownames(sce) <- new.names[rownames(sce)]
# stopifnot(sum(is.na(new.names))==0)
# stopifnot(sum(duplicated(new.names))==0)

##########
## Plot ##
##########

# genes.to.plot <- c("Lefty1","Cd34","Tmsb4x","Fgf3","Spata7","Cer1","Spink1","Dppa4","Dppa5a","Prc1","Lefty2","Ube2c","Hba-x","Hbb-y","Hba-a1","Hbb-bh1")
# genes.to.plot <- c("Vegfa","Vegfb","Vegfc","Vegfd","Kdr","Flt1","Tal1","Runx1","Etv2)
# genes.to.plot <- c("Tet1","Tet2","Tet3","Dnmt1","Dnmt3a","Dnmt3b","Dnmt3l")
genes.to.plot <- rownames(sce)[grep("Tet",rownames(sce))]
genes.to.plot <- rownames(sce)[grep("tomato",rownames(sce))]

for (i in 1:length(genes.to.plot)) {
  gene <- genes.to.plot[i]
  print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
  
  # Create data.table to plot
  to.plot <- data.table(
    cell = colnames(sce),
    expr = logcounts(sce[gene,])[1,]
  ) %>% merge(sample_metadata, by="cell")

  # Plot
  # p <- ggplot(to.plot, aes(x=class, y=expr, fill=class)) +
  p <- ggplot(to.plot, aes(x=stage, y=expr, fill=batch)) +
    # geom_violin(scale = "width", alpha=0.8) +
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
    # geom_jitter(size=2, shape=21, stroke=0.2, alpha=0.5) +
    # scale_fill_manual(values=opts$colors) +
    facet_wrap(~celltype.mapped, scales="fixed") +
    theme_classic() +
    labs(title=gene, x="",y=sprintf("%s expression",gene)) +
    theme(
      strip.text = element_text(size=rel(1.0)),
      # plot.title = element_text(hjust = 0.5, size=rel(1.1), color="black"),
      plot.title = element_blank(),
      # axis.text.x = element_text(colour="black",size=rel(1.2), angle=50, hjust=1),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(colour="black",size=rel(1.0)),
      axis.title.y = element_text(colour="black",size=rel(1.2)),
      legend.position="right",
      legend.title = element_blank(),
      legend.text = element_text(size=rel(1.1))
    )
    
  # pdf(sprintf("%s/%s.pdf",io$outdir,i), width=5, height=3.5, useDingbats = F)
  # ggsave("ggtest.png", width = 3.25, height = 3.25, dpi = 1200)
  jpeg(sprintf("%s/%s.jpeg",io$outdir,gene), width = 1500, height = 800)
  # jpeg(sprintf("%s/test_%s.jpeg",io$outdir,gene), width = 900, height = 900)
  print(p)
  dev.off()
}
