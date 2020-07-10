library(SingleCellExperiment)
library(scater)
library(scran)
library(ggpubr)

source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")

mapping <- fread("/Users/ricard/data/10x_gastrulation_TetChimera/results/second_batch/mapping/disagreement_leah_vs_ricard.txt") %>%
  .[,celltype.mapped.leah:=stringr::str_replace_all(celltype.mapped.leah," ","_")] %>%
  .[,celltype.mapped.leah:=stringr::str_replace_all(celltype.mapped.leah,"/","_")] %>%
  .[,celltype.mapped.ricard:=stringr::str_replace_all(celltype.mapped.ricard," ","_")]


sample_metadata <- sample_metadata %>%
  .[cell%in%mapping$cell]


sce <- readRDS(io$sce)[,mapping$cell]

sce <- logNormCounts(sce)
# foo <- fread(io$atlas.average_expression_per_celltype)
gene_markers <- fread(io$atlas.marker_genes)

foo <- mapping
for (i in head(mapping$cell,n=1000)) {
  print(i)
  # i <- "cell_35907"
  celltype.leah <- mapping[cell==i,celltype.mapped.leah] 
  celltype.ricard <- mapping[cell==i,celltype.mapped.ricard]
  
  diff <- fread(sprintf("%s/%s_vs_%s.txt.gz",io$atlas.differential,celltype.leah,celltype.ricard)) %>%
    .[sig==T & abs(logFC)>1.5]
  
  to.plot <- as.matrix(logcounts(sce[diff$ens_id,i])) %>% t %>% as.data.table(keep.rownames = "cell") %>%
  melt(id.vars = "cell", value.name = "expr", variable.name = "ens_id") %>%
    merge(diff[,c("ens_id","gene","logFC")], by="ens_id") %>%
    .[,sign:=sign(logFC)] %>%
    # merge(gene_markers[celltype%in%c(celltype.leah,celltype.ricard),c("celltype","ens_id")], by=c("ens_id")) %>%
    .[,.(expr=mean(expr)),by="sign"]
  
  
  sign <- to.plot[which.max(to.plot$expr),sign]
  if (sign==(-1)) {
    foo[cell==i,predicted:=celltype.leah]
  } else if (sign==(1)) {
    foo[cell==i,predicted:=celltype.ricard]
  }
}

foo[,correct_ricard:=celltype.mapped.ricard==predicted]
table(foo$correct_ricard)
