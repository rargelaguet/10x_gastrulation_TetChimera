library(ggpubr)

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/10x_gastrulation_TetChimera/settings.R")
}

opts$batches <- c(
  "E75_TET_TKO_L002", 
  "E75_WT_Host_L001", 
  "E85_Rep1_TET_TKO_L004", 
  "E85_Rep1_WT_Host_L003", 
  "E85_Rep2_TET_TKO_L006", 
  "E85_Rep2_WT_Host_L005"
  # "SIGAE4_E105_3_TET123_Chimera_Host_L005", 
  # "SIGAF4_E105_3_TET123_Chimera_TKO_L006", 
  # "SIGAG4_E105_5_TET123_Chimera_Host_L007", 
  # "SIGAH4_E105_5_TET123_Chimera_TKO_L008"
)

io$outdir <- paste0(io$basedir,"/results/general_stats/pdf")

sample_metadata <- sample_metadata %>% 
  .[batch%in%opts$batches]

sample_metadata <- fread("/Users/ricard/data/10x_gastrulation_TetChimera/processed/first_batch/sample_metadata.txt.gz")

###############################################
## Boxplots of general statistics per batch ##
###############################################

to.plot <- sample_metadata %>% 
  melt(id.vars=c("cell","batch","stage"), measure.vars=c("nCount_RNA","nFeature_RNA"))

p <- ggboxplot(to.plot, x = "batch", y = "value", outlier.shape=NA) +
  yscale("log10", .format = TRUE) +
  labs(x="", y="") +
  facet_wrap(~variable, scales="free_y", nrow=2) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    # axis.text.x = element_blank(),
    axis.text.x = element_text(colour="black",size=rel(0.9), angle=50, hjust=1),
    # axis.ticks.x = element_blank()
  )

pdf(paste0(io$outdir,"/general_stats_per_batch.pdf"), width=8, height=6, useDingbats = F)
print(p)
dev.off()

########################################
## Barplots number of cells per batch ##
########################################

to.plot <- sample_metadata[,.N,by="batch"]

p <- ggbarplot(to.plot, x = "batch", y = "N", fill="gray70") +
  labs(x="", y="Number of cells (after QC)") +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    # axis.text.x = element_blank(),
    axis.text.x = element_text(colour="black",size=rel(0.8), angle=40, hjust=1),
    # axis.ticks.x = element_blank()
  )

pdf(paste0(io$outdir,"/N_per_batch.pdf"), width=8, height=6, useDingbats = F)
print(p)
dev.off()


##################################################
## Boxplots of general statistics per cell type ##
##################################################

to.plot <- sample_metadata %>% 
  melt(id.vars=c("cell","celltype.mapped"), measure.vars=c("nCount_RNA","nFeature_RNA"))

p <- ggboxplot(to.plot, x = "celltype.mapped", y = "value", fill="celltype.mapped", outlier.shape=NA) +
  yscale("log10", .format = TRUE) +
  labs(x="", y="") +
  scale_fill_manual(values=opts$celltype.colors) +
  facet_wrap(~variable, scales="free_y") +
  guides(fill = guide_legend(override.aes = list(size=0.25), ncol=1)) +
  scale_size(guide = 'none') +
  theme(
    legend.position = "right",
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

pdf(paste0(io$outdir,"/general_stats_per_celltype.pdf"), width=16, height=10, useDingbats = F)
print(p)
dev.off()

