source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")

files <- c(
  "/Users/ricard/data/10x_gastrulation_TetChimera/results/second_batch/mapping/leah/mapping_mnn_WT.txt.gz",
  "/Users/ricard/data/10x_gastrulation_TetChimera/results/second_batch/mapping/leah/mapping_mnn_TET_TKO.txt.gz",
  "/Users/ricard/data/10x_gastrulation_TetChimera/results/third_batch/mapping/leah/mapping_mnn_TET_Host.txt.gz",
  "/Users/ricard/data/10x_gastrulation_TetChimera/results/third_batch/mapping/leah/mapping_mnn_TET_TKO.txt.gz"
)



dt <- files %>% map(~ fread(.)[!is.na(celltype.mapped),c("barcode","batch","celltype.mapped","stage.mapped","closest.cell","celltype.score","cellstage.score")]) %>% rbindlist %>%
  .[,cell2:=paste(batch,barcode,sep="_")] %>%
  .[,c("barcode","batch"):=NULL]
  # .[cell2%in%sample_metadata$cell2]

sample_metadata <- fread("/Users/ricard/data/10x_gastrulation_TetChimera/processed/merged/sample_metadata.txt.gz") %>%
  .[,c("cell2", "cell", "barcode", "batch", "stage", "class", "target", "nCount_RNA", "nFeature_RNA", "percent.mt", "doublet_score", "pass_QC")] %>%
  merge(dt,by="cell2",all.x=T)

sample_metadata %>%
  .[,cell:=NULL] %>%
  setnames("cell2","cell")

fwrite(sample_metadata, io$metadata, sep="\t", col.names=T, row.names=F, na="NA", quote=F)
