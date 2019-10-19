library(data.table)
library(purrr)

######################
## Define  settings ##
######################

# I/O
io <- list()
io$sample_metadata <- "/Users/ricard/data/scnmt_gastrulation_TetKO/sample_metadata.txt"
io$mapping <- "/Users/ricard/data/scnmt_gastrulation_TetKO/rna/mapping_10x/mapping10x_mnn.rds"
io$outdir <- "/Users/ricard/data/scnmt_gastrulation_TetKO/rna/mapping_10x"

###############
## Load data ##
###############

mapping <- readRDS(io$mapping)
sample_metadata <- fread(io$sample_metadata)

################################################
## Add lineage annotations to sample metadata ##
################################################

mapping.dt <- mapping$mapping %>% as.data.table %>% 
  setnames("cell","id_rna") %>%
  setnames("celltype.mapped","lineage10x") %>%
  .[,closest.cell:=NULL]

sample_metadata <- merge(sample_metadata,mapping.dt, by="id_rna", all.x=T)

sample_metadata %>%
  .[,lineage10x_2:=lineage10x] %>%
  
  # Mesoderm
  .[lineage10x%in%c("Nascent mesoderm","Pharyngeal mesoderm","Paraxial mesoderm","ExE mesoderm","Mesenchyme","Mixed mesoderm","Intermediate mesoderm","Somitic mesoderm","Caudal mesoderm"), lineage10x_2:="Mesoderm"] %>%
  # Blood
  .[lineage10x%in%c("Blood progenitors 1","Blood progenitors 2","Haematoendothelial progenitors"), lineage10x_2:="Blood"] %>%
  # Embryonic endoderm
  .[lineage10x%in%c("Gut","Def. endoderm","Notochord"), lineage10x_2:="Endoderm"] %>%
  # Extra-embryonic endoderm
  .[lineage10x%in%c("Parietal endoderm","ExE endoderm","Visceral endoderm"), lineage10x_2:="ExE endoderm"] %>%
  # Primitive streak
  .[lineage10x%in%c("Primitive Streak", "Caudal epiblast","Anterior Primitive Streak"), lineage10x_2:="Primitive streak"] %>%
  # Ectoderm
  .[lineage10x%in%c("Rostral neurectoderm","Surface ectoderm"), lineage10x_2:="Ectoderm"] %>%
  # PGC are not expected, they are most likely epiblast
  .[lineage10x%in%c("PGC"), c("lineage10x","lineage10x_2"):="Epiblast"]
  # At E6.5 we do not expect any mature ectoderm
  # .[stage=="E6.5" & lineage10x%in%c("Rostral neurectoderm","Surface ectoderm"), lineage10x_2:="Epiblast"]

unique(sample_metadata$lineage10x_2)

# Save
fwrite(sample_metadata, file=paste0(io$outdir,"/sample_metadata_mapping_mnn.txt"), sep="\t", row.names=F, col.names=T, na="NA", quote=F)

