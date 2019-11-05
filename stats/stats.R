library(Seurat)
library(data.table)
library(purrr)
library(ggpubr)

########################
## Load Seurat object ##
########################

# srat <- readRDS("/Users/ricard/data/10x_gastrulation_TetChimera/processed/seurat.rds")
srat <- readRDS("/Users/ricard/data/10x_gastrulation_TetChimera/processed/seurat.rds")

# srat <- srat[,srat@meta.data$genotype=="WT"]
# meta.data <- srat@meta.data

# Plot QC stats
colnames(srat@meta.data)

dt <- srat@meta.data %>% as.data.table
# VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "genotype")

dt2 = dt %>% melt(id.vars=c("batch"), measure.vars=c("nFeature_RNA", "nCount_RNA"))

ggboxplot(dt2, x="batch", y="value") +
  facet_wrap(~variable, scales="free") +
  theme(
    axis.text.x = element_text(color="black", angle=30, hjust=1, vjust=1)
  )



dt <- fread("/Users/ricard/data/gastrulation10x/original/meta.tab")
head(dt)


#######

asd <- 