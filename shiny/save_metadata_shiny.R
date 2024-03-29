source(here::here("settings.R"))

#####################
## Define settings ##
#####################

# io$metadata <- file.path(io$basedir,"results/mapping/sample_metadata_after_mapping.txt.gz")
io$outfile <- "/Users/argelagr/shiny_tet/data/cell_metadata.txt.gz"

########################
## Load cell metadata ##
########################

sample_metadata <- fread(io$metadata) %>%
  setnames("celltype.mapped","celltype") %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE] %>%
  .[,genotype:=ifelse(grepl("WT",class),"WT","Tet_TKO")] %>%
  .[,c("cell","genotype","class","stage","alias","celltype","nFeature_RNA","mit_percent_RNA","rib_percent_RNA")] %>%
  setnames("alias","sample")
  # .[,celltype:=stringr::str_replace_all(celltype,opts$rename_celltypes)] %>%

table(sample_metadata$class)
table(sample_metadata$genotype)

##########
## Save ##
##########

fwrite(sample_metadata, io$outfile, sep="\t", na="NA", quote=F)
