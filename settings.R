suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(argparse))
# suppressPackageStartupMessages(library(Seurat))

#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/10x_gastrulation_TetChimera"
  io$Rscript <- "/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript"
  io$atlas.basedir <- "/Users/ricard/data/gastrulation10x"
  # io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
  io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  io$basedir <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_TetChimera"
  io$Rscript <- "/nfs/research1/stegle/users/ricard/R-4.0.3/bin/Rscript"
  io$atlas.basedir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x"
  # io$gene_metadata <- "/hps/nobackup2/research/stegle/users/ricard/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
  io$gene_metadata <- "/hps/nobackup2/research/stegle/users/ricard/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
} else if (Sys.info()[['nodename']]=="BI2404M") {
  io$basedir <- "/Users/argelagr/data/10x_gastrulation_TetChimera"
  io$atlas.basedir <- "/Users/argelagr/data/gastrulation10x"
  io$gene_metadata <- "/Users/argelagr/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
} else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
  if (grepl("Clark", Sys.info()['effective_user'])) {
    io$basedir <- "/bi/scratch/Stephen_Clark/tet_chimera_10x/"
    io$gene_metadata <- "/bi/scratch/Stephen_Clark/annotations/Mmusculus_genes_BioMart.87.txt"
  } else if (grepl("argelag", Sys.info()['effective_user'])) {
    io$basedir <- "/bi/group/reik/ricard/data/10x_gastrulation_TetChimera"
    io$atlas.basedir <- "/bi/group/reik/ricard/data/pijuansala2019_gastrulation10x"
    io$gene_metadata <- "/bi/group/reik/ricard/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
  }
} else {
  stop("Computer not recognised")
}

io$metadata <- paste0(io$basedir,"/sample_metadata.txt.gz")
io$seurat <- paste0(io$basedir,"/processed_all/seurat.rds")
io$sce <- paste0(io$basedir,"/processed_all/SingleCellExperiment.rds")
# io$anndata <- paste0(io$basedir,"/processed_all/anndata.h5ad")

# Atlas information
io$atlas.metadata <- paste0(io$atlas.basedir,"/sample_metadata.txt.gz")
io$atlas.marker_genes <- paste0(io$atlas.basedir,"/results/marker_genes/old/all_stages/marker_genes.txt.gz")
io$atlas.differential <- paste0(io$atlas.basedir,"/results/differential")
io$atlas.average_expression_per_celltype <- paste0(io$atlas.basedir,"/results/marker_genes/all_stages/avg_expr_per_celltype_and_gene.txt.gz")
io$atlas.sce <- paste0(io$atlas.basedir,"/processed/SingleCellExperiment.rds")


#############
## Options ##
#############

opts <- list()

opts$celltypes = c(
	"Epiblast",
	"Primitive_Streak",
	"Caudal_epiblast",
	"PGC",
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
	"Surface_ectoderm",
	"Visceral_endoderm",
	"ExE_endoderm",
	"ExE_ectoderm",
	"Parietal_endoderm"
)

opts$celltype.colors = c(
	"Epiblast" = "#635547",
	"Primitive_Streak" = "#DABE99",
	"Caudal_epiblast" = "#9e6762",
	"PGC" = "#FACB12",
	"Anterior_Primitive_Streak" = "#c19f70",
	"Notochord" = "#0F4A9C",
	"Def._endoderm" = "#F397C0",
	"Gut" = "#EF5A9D",
	"Nascent_mesoderm" = "#C594BF",
	"Mixed_mesoderm" = "#DFCDE4",
	"Intermediate_mesoderm" = "#139992",
	"Caudal_Mesoderm" = "#3F84AA",
	"Paraxial_mesoderm" = "#8DB5CE",
	"Somitic_mesoderm" = "#005579",
	"Pharyngeal_mesoderm" = "#C9EBFB",
	"Cardiomyocytes" = "#B51D8D",
	"Allantois" = "#532C8A",
	"ExE_mesoderm" = "#8870ad",
	"Mesenchyme" = "#cc7818",
	"Haematoendothelial_progenitors" = "#FBBE92",
	"Endothelium" = "#ff891c",
	"Blood_progenitors" = "#c9a997",
	"Blood_progenitors_1" = "#f9decf",
	"Blood_progenitors_2" = "#c9a997",
	"Erythroid" = "#EF4E22",
	"Erythroid1" = "#C72228",
	"Erythroid2" = "#f79083",
	"Erythroid3" = "#EF4E22",
	"NMP" = "#8EC792",
	"Neurectoderm" = "#65A83E",
	"Rostral_neurectoderm" = "#65A83E",
	"Caudal_neurectoderm" = "#354E23",
	"Neural_crest" = "#C3C388",
	"Forebrain_Midbrain_Hindbrain" = "#647a4f",
	"Spinal_cord" = "#CDE088",
	"Surface_ectoderm" = "#f7f79e",
	"Visceral_endoderm" = "#F6BFCB",
	"ExE_endoderm" = "#7F6874",
	"ExE_ectoderm" = "#989898",
	"Parietal_endoderm" = "#1A1A1A"
)


opts$samples <- c(
  
  # First samples (all failed QC)
  # SIGAA3_E8.5_pool1_Host-WT_L001
  # SIGAB3_E8.5_pool1_TET-TKO_L002
  # SIGAC3_E8.5_pool2_Host-WT_L003
  # SIGAD3_E8.5_pool2_TET-TKO_L004
  # SIGAE3_E7.5_pool1_Host-WT_L005
  # SIGAF3_E7.5_pool1_TET-TKO_L006
  # SIGAG3_E8.5_hashing_Host-WT_L007
  # SIGAH3_E8.5_hasting_TET-TKO_L00
  
  # Second batch
  "E75_TET_TKO_L002",
  "E75_WT_Host_L001",
  "E85_Rep1_TET_TKO_L004",
  "E85_Rep1_WT_Host_L003",
  "E85_Rep2_TET_TKO_L006",
  "E85_Rep2_WT_Host_L005",
  
  # Third batch (all failed QC)
  # "SIGAE4_E105_3_TET123_Chimera_Host_L005", 
  # "SIGAF4_E105_3_TET123_Chimera_TKO_L006", 
  # "SIGAG4_E105_5_TET123_Chimera_Host_L007", 
  # "SIGAH4_E105_5_TET123_Chimera_TKO_L008"
  
  # Fourth batch
  # "SIGAC2_TET_TKO_E9_5_Head1",
  # "SIGAD2_TET_TKO_E9_5_Trunk1",
  # "SIGAE2_TET_TKO_E9_5_Tail1",
  # "SIGAE6_TET_TKO_E9_5_Head2",
  # "SIGAF2_TET_TKO_E9_5_YS1",
  # "SIGAF6_TET_TKO_E9_5_Trunk2",
  # "SIGAG6_TET_TKO_E9_5_Tail2",
  # "SIGAH6_TET_TKO_E9_5_YS2",
  
  # Fifth batch
  "E8_5_TET_WT_rep1_SIGAG8",
  "E8_5_TET_WT_rep2_SIGAH8",

  # WT controls from Carolina
  "E7.5_batch_1_tdTomato+",  # sample_1
  "E7.5_batch_1_tdTomato-",  # sample_2
  "E7.5_batch_2_tdTomato+",  # sample_3
  "E7.5_batch_2_tdTomato-",  # sample_4
  "E8.5_batch_3_tdTomato+",  # sample_5
  "E8.5_batch_3_tdTomato-",  # sample_6
  "E8.5_batch_4_tdTomato+",  # sample_7
  "E8.5_batch_4_tdTomato-",  # sample_8
  "E8.5_batch_5_tdTomato+",  # sample_9
  "E8.5_batch_5_tdTomato-"  # sample_10
)

opts$stage_classes <- c(
  "E7.5_WT_tdTomato+", 
  "E7.5_WT_tdTomato-", 
  "E7.5_TET_TKO", 
  "E8.5_WT_tdTomato+",
  "E8.5_WT_tdTomato-",
  "E8.5_TET_TKO"
  # "E9.5_TET_TKO"
)

opts$classes <- c("WT_tdTomato+", "WT_tdTomato-", "TET_TKO")


opts$sample2class <- c(
  "E75_TET_TKO_L002" = "TET_TKO",
  "E75_WT_Host_L001" = "WT_tdTomato-",
  "E85_Rep1_TET_TKO_L004" = "TET_TKO",
  "E85_Rep1_WT_Host_L003" = "WT_tdTomato-",
  "E85_Rep2_TET_TKO_L006" = "TET_TKO",
  "E85_Rep2_WT_Host_L005" = "WT_tdTomato-",
  # "SIGAC2_TET_TKO_E9_5_Head1" = "TET_TKO",
  # "SIGAD2_TET_TKO_E9_5_Trunk1" = "TET_TKO",
  # "SIGAE2_TET_TKO_E9_5_Tail1" = "TET_TKO",
  # "SIGAE6_TET_TKO_E9_5_Head2" = "TET_TKO",
  # "SIGAF2_TET_TKO_E9_5_YS1" = "TET_TKO",
  # "SIGAF6_TET_TKO_E9_5_Trunk2" = "TET_TKO",
  # "SIGAG6_TET_TKO_E9_5_Tail2" = "TET_TKO",
  # "SIGAH6_TET_TKO_E9_5_YS2" = "TET_TKO",
  "E8_5_TET_WT_rep1_SIGAG8" = "TET_TKO",
  "E8_5_TET_WT_rep2_SIGAH8" = "TET_TKO",

  "E7.5_batch_1_tdTomato-" = "WT_tdTomato-",
  "E7.5_batch_1_tdTomato+" = "WT_tdTomato+",
  "E7.5_batch_2_tdTomato-" = "WT_tdTomato-",
  "E7.5_batch_2_tdTomato+" = "WT_tdTomato+",
  "E8.5_batch_3_tdTomato-" = "WT_tdTomato-",
  "E8.5_batch_3_tdTomato+" = "WT_tdTomato+",
  "E8.5_batch_4_tdTomato-" = "WT_tdTomato-",
  "E8.5_batch_4_tdTomato+" = "WT_tdTomato+",
  "E8.5_batch_5_tdTomato-" = "WT_tdTomato-",
  "E8.5_batch_5_tdTomato+" = "WT_tdTomato+"
)

opts$stages <- c(
  "E7.5", 
  "E8.5"
  # "E9.5"
)

opts$atlas.stages <- c(
  "E6.5",
  "E6.75",
  "E7.0",
  "E7.25",
  "E7.5",
  "E7.75",
  "E8.0",
  "E8.25",
  "E8.5"
  # "mixed_gastrulation"
)


opts$stage.colors = c(
  "E8.5" = "#440154FF",
  "E8.25" = "#472D7BFF",
  "E8.0" = "#3B528BFF",
  "E7.75" = "#2C728EFF",
  "E7.5" = "#21908CFF",
  "E7.25" = "#27AD81FF",
  "E7.0" = "#5DC863FF",
  "E6.75" = "#AADC32FF",
  "E6.5" = "#FDE725FF"
)
# opts$stage.colors <- viridis::viridis(n=length(opts$stages))
# names(opts$stage.colors) <- rev(opts$stages)


opts$sample2alias <- c(
  "E75_TET_TKO_L002"           = "E7.5_Tet_TKO",
  "E75_WT_Host_L001"           = "E7.5_WT_tdTomato-_1",
  "E85_Rep1_TET_TKO_L004"      = "E8.5_Tet_TKO_1",
  "E85_Rep1_WT_Host_L003"      = "E8.5_WT_tdTomato-_1",
  "E85_Rep2_TET_TKO_L006"      = "E8.5_Tet_TKO_2",
  "E85_Rep2_WT_Host_L005"      = "E8.5_WT_tdTomato-_2",
  "E8_5_TET_WT_rep1_SIGAG8"    = "E8.5_Tet_TKO_3",
  "E8_5_TET_WT_rep2_SIGAH8"    = "E8.5_Tet_TKO_4",
  # "SIGAC2_TET_TKO_E9_5_Head1"  = "E9.5_Tet_TKO_Head_1",
  # "SIGAD2_TET_TKO_E9_5_Trunk1" = "E9.5_Tet_TKO_Trunk_1",
  # "SIGAE2_TET_TKO_E9_5_Tail1"  = "E9.5_Tet_TKO_Tail_1",
  # "SIGAE6_TET_TKO_E9_5_Head2"  = "E9.5_Tet_TKO_Head_2",
  # "SIGAF2_TET_TKO_E9_5_YS1"    = "E9.5_Tet_TKO_Yolk_Sac_1",
  # "SIGAF6_TET_TKO_E9_5_Trunk2" = "E9.5_Tet_TKO_Trunk_2",
  # "SIGAG6_TET_TKO_E9_5_Tail2"  = "E9.5_Tet_TKO_Tail_2",
  # "SIGAH6_TET_TKO_E9_5_YS2"    = "E9.5_Tet_TKO_Yolk_Sac_2" 
  "E7.5_batch_1_tdTomato-" = "E7.5_WT_tdTomato-_2",
  "E7.5_batch_1_tdTomato+" = "E7.5_WT_tdTomato+_1",
  "E7.5_batch_2_tdTomato-" = "E7.5_WT_tdTomato-_3",
  "E7.5_batch_2_tdTomato+" = "E7.5_WT_tdTomato+_2",
  "E8.5_batch_3_tdTomato-" = "E8.5_WT_tdTomato-_3",
  "E8.5_batch_3_tdTomato+" = "E8.5_WT_tdTomato+_1",
  "E8.5_batch_4_tdTomato-" = "E8.5_WT_tdTomato-_4",
  "E8.5_batch_4_tdTomato+" = "E8.5_WT_tdTomato+_2",
  "E8.5_batch_5_tdTomato-" = "E8.5_WT_tdTomato-_5",
  "E8.5_batch_5_tdTomato+" = "E8.5_WT_tdTomato+_3"
)

opts$stage_class_colors <- c(
  "E7.5_WT_tdTomato+" = "#5CACEE",# "#CCCCCC", 
  "E7.5_WT_tdTomato-" = "#5CACEE",#"#CCCCCC", 
  "E7.5_TET_TKO"  = "#FF7F50", 
  "E8.5_WT_tdTomato+" = "#1C86EE",#"#B0B0B0", 
  "E8.5_WT_tdTomato-" = "#1C86EE",#"#B0B0B0", 
  "E8.5_TET_TKO" = "#EE4000"
  # "E9.5_TET_TKO"
)

opts$stage_class2_colors <- c(
  "E7.5_WT" = "#5CACEE",# "#CCCCCC", 
  "E7.5_TET_TKO"  = "#FF7F50", 
  "E8.5_WT" = "#1C86EE", # "#B0B0B0", 
  "E8.5_TET_TKO" = "#EE4000"
)

opts$class_colors <- c(
  "WT_tdTomato-" = "#5CACEE", # "#CCCCCC", 
  "WT_tdTomato+" = "#5CACEE",# "#B0B0B0", 
  "TET_TKO" = "#EE4000"
)

opts$class2_colors <- c(
  "WT" = "#4F94CD",
  "TET_TKO" = "#EE4000",
  "Tet-TKO" = "#EE4000"
)

##########################
## Load sample metadata ##
##########################

# sample_metadata <- fread(io$metadata) %>% 
#   # .[pass_QC==T] %>% 
#   # .[batch%in%opts$samples] %>%
#   .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped," ","_")] %>%
#   .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,"/","_")] %>%
#   .[,celltype.mapped:=factor(celltype.mapped, levels=names(opts$celltype.colors))]

# metadata <- fread("/Users/argelagr/data/10x_gastrulation_TetChimera/results_all/mapping/sample_metadata_after_mapping.txt.gz")
# stopifnot(metadata$sample%in%names(opts$sample2alias))
# metadata[,alias:=opts$sample2alias[sample]]
# print(table(metadata$alias))
# stopifnot(!is.na(metadata$alias))
# fwrite(metadata, "/Users/argelagr/data/10x_gastrulation_TetChimera/results_all/mapping/sample_metadata_after_mapping.txt.gz", sep="\t", na="NA", quote=F)

##########
## E9.5 ##
##########

# sample_metadata[,region:=sample_metadata$batch %>% strsplit("_")  %>% map_chr(6) %>% substr(.,1,nchar(.)-1)]

