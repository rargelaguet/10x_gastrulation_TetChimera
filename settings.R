suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))

#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/10x_gastrulation_TetChimera"
  io$atlas.basedir <- "/Users/ricard/data/gastrulation10x"
  io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  io$basedir <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_TetChimera"
  io$atlas.basedir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x"
  io$gene_metadata <- "/hps/nobackup2/research/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
} else {
  stop("Computer not recognised")
}

io$metadata <- paste0(io$basedir,"/processed/second_batch/sample_metadata.txt.gz")
io$seurat <- paste0(io$basedir,"/processed/second_batch/seurat.rds")
io$sce <- paste0(io$basedir,"/processed/second_batch/SingleCellExperiment.rds")
# io$counts <- paste0(io$basedir,"/rna/counts_datatable.tsv.gz")

# Atlas information
io$atlas.metadata <- paste0(io$atlas.basedir,"/sample_metadata.txt.gz")
io$atlas.marker_genes <- paste0(io$atlas.basedir,"/results/marker_genes/marker_genes.txt.gz")
io$atlas.average_expression_per_celltype <- paste0(io$atlas.basedir,"/results/marker_genes/avg_expr_per_celltype_and_gene.txt.gz")
io$atlas.sce <- paste0(io$atlas.basedir,"/processed/SingleCellExperiment.rds")


#############
## Options ##
#############

opts <- list()

opts$celltypes.1 = c(
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

opts$colors1 = c(
	"Epiblast" = "#635547",
	"Primitive Streak" = "#DABE99",
	"Caudal epiblast" = "#9e6762",
	"PGC" = "#FACB12",
	"Anterior Primitive Streak" = "#c19f70",
	"Notochord" = "#0F4A9C",
	"Def. endoderm" = "#F397C0",
	"Gut" = "#EF5A9D",
	"Nascent mesoderm" = "#C594BF",
	"Mixed mesoderm" = "#DFCDE4",
	"Intermediate mesoderm" = "#139992",
	"Caudal Mesoderm" = "#3F84AA",
	"Paraxial mesoderm" = "#8DB5CE",
	"Somitic mesoderm" = "#005579",
	"Pharyngeal mesoderm" = "#C9EBFB",
	"Cardiomyocytes" = "#B51D8D",
	"Allantois" = "#532C8A",
	"ExE mesoderm" = "#8870ad",
	"Mesenchyme" = "#cc7818",
	"Haematoendothelial progenitors" = "#FBBE92",
	"Endothelium" = "#ff891c",
	"Blood progenitors 1" = "#f9decf",
	"Blood progenitors 2" = "#c9a997",
	"Erythroid1" = "#C72228",
	"Erythroid2" = "#f79083",
	"Erythroid3" = "#EF4E22",
	"NMP" = "#8EC792",
	"Rostral neurectoderm" = "#65A83E",
	"Caudal neurectoderm" = "#354E23",
	"Neural crest" = "#C3C388",
	"Forebrain/Midbrain/Hindbrain" = "#647a4f",
	"Spinal cord" = "#CDE088",
	"Surface ectoderm" = "#f7f79e",
	"Visceral endoderm" = "#F6BFCB",
	"ExE endoderm" = "#7F6874",
	"ExE ectoderm" = "#989898",
	"Parietal endoderm" = "#1A1A1A"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% .[pass_QC==T]
  
  
