library(ggpubr)
library(SingleCellExperiment)
library(scran)

#####################
## Define settings ##
#####################

source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
io$seurat <- paste0(io$basedir,"/processed/fourth_batch/seurat.rds")
io$metadata <- paste0(io$basedir,"/processed/fourth_batch/metadata.txt.gz")
io$outdir <- paste0(io$basedir,"/processed/fourth_batch/qc/")

opts$nFeature_RNA <- 2500
opts$nCount_RNA <- 5000
opts$percent.mt <- 10


###############
## Load data ##
###############

srat <- readRDS(io$seurat)

metadata <- fread(io$metadata)
# md$doublet_score <- NA

# mtx_sub_list <- list()
# cell.info_sub_list <- list()
# srat_sub_list <- list()

###########
## Plots ##
###########

to.plot <- metadata %>%
    melt(id.vars=c("batch","cell"), measure.vars=c("nCount_RNA","nFeature_RNA","percent.mt"))

tmp <- data.table(
    variable = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
    value = c(opts$nCount_RNA, opts$nFeature_RNA,opts$percent.mt)
)

p <- gghistogram(to.plot, x="value", fill="batch", bins=50) +
    geom_vline(aes(xintercept=value), linetype="dashed", data=tmp) +
    facet_wrap(~variable, scales="free") +
    theme(
        axis.text =  element_text(size=rel(0.8)),
        axis.title.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank()
    )
    
pdf(sprintf("%s/qc_metrics.pdf",io$outdir), width=14, height=8, useDingbats = F)
print(p)
dev.off()

#############################
## Calculate doublet score ##
#############################

opts$max_doublet_score <- 10000

sce <- as.SingleCellExperiment(srat)
metadata[,doublet_score:=doubletCells(sce)]

# Plot
to.plot <- metadata %>% 
    melt(id.vars=c("batch","cell"), measure.vars=c("doublet_score"))

p <- gghistogram(to.plot, x="value", fill="batch", bins=50) +
    geom_vline(xintercept=opts$max_doublet_score, linetype="dashed") +
    theme(
        axis.text =  element_text(size=rel(0.8)),
        legend.position = "right"
    )

table(metadata$doublet_score<opts$max_doublet_score)

pdf(sprintf("%s/doublet_score.pdf",io$outdir), width=9, height=5, useDingbats = F)
print(p)
dev.off()

###########
## Filter ##
###########

cells <- metadata %>%
    # .[ nFeature_RNA>opts$nFeature_RNA & nCount_RNA>opts$nCount_RNA & percent.mt<opts$percent.mt & doublet_score<opts$max_doublet_score,cell]
    .[ nFeature_RNA>opts$nFeature_RNA & nCount_RNA>opts$nCount_RNA & percent.mt<opts$percent.mt,cell]
length(cells) / nrow(srat@meta.data)

# Subset
metadata[,pass_QC:=cell%in%cells]
srat <- srat[,cells]

##########
## Save ##
##########

saveRDS(srat, io$seurat)
fwrite(metadata, io$metadata, quote=F, na="NA", sep="\t")

