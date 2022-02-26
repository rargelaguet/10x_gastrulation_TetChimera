here::i_am("shiny/save_expr_matrix.R")


source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$outfile <- file.path(io$basedir,"shiny/SingleCellExperiment_pseudobulk_class_celltype.rds") 
io$sce.pseudobulk <- file.path(io$basedir,"results_new/pseudobulk/SingleCellExperiment_pseudobulk_class_celltype.rds")

# Define options
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
  # "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm"
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

opts$classes <- c("WT", "TET_TKO")

# opts$rename_celltypes <- c(
#   "Erythroid3" = "Erythroid",
#   "Erythroid2" = "Erythroid",
#   "Erythroid1" = "Erythroid",
#   "Blood_progenitors_1" = "Blood_progenitors",
#   "Blood_progenitors_2" = "Blood_progenitors",
#   "Anterior_Primitive_Streak" = "Primitive_Streak"
# )

opts$min.cells <- 30

###############
## Load data ##
###############

# Load SingleCellExperiment object
sce <- readRDS(io$sce.pseudobulk)

# Add metadata
sce$class <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(1)
sce$celltype <- stringr::str_split(colnames(sce), pattern = "-") %>% map_chr(2)

# Filter by minimum number of cells
sce <- sce[,names(which(metadata(sce)$n_cells>=opts$min.cells))]

# Filter
sce <- sce[,sce$celltype%in%opts$celltypes]
sce <- sce[,sce$class%in%opts$classes]

genes <- fread("/Users/argelagr/data/10x_gastrulation_TetChimera/shiny/genes.txt", header = F)[[1]]
sce <- sce[genes,]

##########
## Save ##
##########

assays(sce)["counts"] <- NULL
saveRDS(sce, io$outfile)
