suppressPackageStartupMessages(library(SingleCellExperiment))

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
  source("/Users/ricard/10x_gastrulation_TetChimera/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/10x_gastrulation_TetChimera/settings.R")
  source("/homes/ricard/10x_gastrulation_TetChimera/utils.R")
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
# opts$celltypes <- c("Blood_progenitors_1")

# Define classes to plot
opts$classes <- c(
  "E7.5_Host",
  "E7.5_TET_TKO",
  "E8.5_Host",
  "E8.5_TET_TKO",
  "E8.5_WT"
)

# Define colors
opts$colors <- c(
  "E7.5_Host" = "#CCCCCC", 
  "E7.5_TET_TKO" = "#FF7F50", 
  "E8.5_Host" = "#B0B0B0", 
  "E8.5_TET_TKO" = "#EE4000",
  "E8.5_WT" = "#B0B0B0"
)

# Update sample metadata
sample_metadata <- fread(io$metadata) %>% 
  .[pass_QC==TRUE & class%in%opts$classes & celltype.mapped%in%opts$celltypes] %>%
  .[,celltype.mapped:=factor(celltype.mapped, levels=opts$celltypes)] %>%
  .[,class:=factor(class,levels=opts$classes)]

table(sample_metadata$class)
table(sample_metadata$celltype.mapped)

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- load_SingleCellExperiment(io$sce, cells=sample_metadata$cell, normalise = FALSE)

# Add sample metadata as colData
colData(sce) <- sample_metadata %>% tibble::column_to_rownames("cell") %>% DataFrame

# Normalise
sce <- logNormCounts(sce)

##########
## Plot ##
##########

genes.to.plot <- c("Eomes","Dppa4")
# genes.to.plot <- c("Lefty1","Cd34","Tmsb4x","Fgf3","Spata7","Cer1","Spink1","Dppa4","Dppa5a","Prc1","Lefty2","Ube2c","Hba-x","Hbb-y","Hba-a1","Hbb-bh1")
# genes.to.plot <- c("Vegfa","Vegfb","Vegfc","Vegfd","Kdr","Flt1","Tal1","Runx1","Etv2)
# genes.to.plot <- c("Tet1","Tet2","Tet3","Dnmt1","Dnmt3a","Dnmt3b","Dnmt3l")
# genes.to.plot <- rownames(sce)[grep("tomato",rownames(sce))]
# genes.to.plot <- fread(io$atlas.marker_genes)$gene %>% unique %>% .[!grepl("Rik$",.)]
genes.to.plot <- fread("/Users/ricard/data/gastrulation10x/results/differential/celltypes/E8.5/Neural_crest_vs_Forebrain_Midbrain_Hindbrain.txt.gz") %>% .[sig==T & logFC<(-2),gene]

for (i in 1:length(genes.to.plot)) {
  
  gene <- genes.to.plot[i]
  
  if (gene %in% rownames(sce)) {
    print(sprintf("%s/%s: %s",i,length(genes.to.plot),gene))
    outfile <- sprintf("%s/Tim/%s.jpeg",io$outdir,gene)
    
    if (!file.exists(outfile)) {
      
      to.plot <- data.table(
        cell = colnames(sce),
        expr = logcounts(sce)[gene,]
      ) %>% merge(sample_metadata[,c("cell","sample","class","celltype.mapped")], by="cell") %>%
        .[,N:=.N,by=c("sample","celltype.mapped")] %>% .[N>=10]
      
      p <- ggplot(to.plot, aes(x=class, y=expr, fill=class)) +
        geom_violin(scale = "width", alpha=0.8) +
        geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.8) +
        # stat_summary(fun.data = give.n, geom = "text", size=2.5) + 
        # geom_jitter(size=2, shape=21, stroke=0.2, alpha=0.5) +
        scale_fill_manual(values=opts$colors) +
        facet_wrap(~celltype.mapped, scales="fixed") +
        theme_classic() +
        labs(title=gene, x="",y=sprintf("%s expression",gene)) +
        guides(x = guide_axis(angle = 90)) +
        theme(
          strip.text = element_text(size=rel(0.85)),
          # plot.title = element_text(hjust = 0.5, size=rel(1.1), color="black"),
          plot.title = element_blank(),
          axis.text.x = element_text(colour="black",size=rel(0.95)),
          # axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(colour="black",size=rel(1.0)),
          axis.title.y = element_text(colour="black",size=rel(1.0)),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size=rel(0.85))
        )
      
        # pdf(outfile, width=11, height=10)
        # png(outfile, width = 1100, height = 1000)
        jpeg(outfile, width = 700, height = 600)
        print(p)
        dev.off()
        
    } else {
      print(sprintf("%s already exists...",outfile))
    }

  } else {
    print(sprintf("%s not found",gene))
  }
}


#########
## Old ##
#########

# Load gene metadata
# gene_metadata <- fread(io$gene_metadata) %>%
#   .[symbol!="" & symbol%in%rownames(sce)]

# Rename genes
# new.names <- gene_metadata$symbol
# names(new.names) <- gene_metadata$ens_id
# sce <- sce[rownames(sce) %in% names(new.names),]
# rownames(sce) <- new.names[rownames(sce)]
# stopifnot(sum(is.na(new.names))==0)
# stopifnot(sum(duplicated(new.names))==0)