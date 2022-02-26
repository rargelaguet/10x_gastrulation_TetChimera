# TO-DO: PLOT DIFFERENTIAL PER STAGE

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$indir <- file.path(io$basedir,"results_all/differential/pseudobulk")
io$outdir <- file.path(io$basedir,"results_all/differential/pseudobulk/pdf"); dir.create(io$outdir, showWarnings = F)

# Options

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
  # "Mixed_mesoderm",
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
  "Blood_progenitors",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid",
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
  # "ExE_ectoderm"
  # "Parietal_endoderm"
)

opts$rename_celltypes <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors"
  # "Intermediate_mesoderm" = "Mixed_mesoderm",
  # "Paraxial_mesoderm" = "Mixed_mesoderm",
  # "Nascent_mesoderm" = "Mixed_mesoderm",
  # "Pharyngeal_mesoderm" = "Mixed_mesoderm"
  # "Visceral_endoderm" = "ExE_endoderm"
)

##############################
## Load precomputed results ##
##############################

# i <- opts$ko.classes[1]; j <- opts$celltypes[10]
diff.dt <- opts$celltypes %>% map(function(i) {
  file <- sprintf("%s/%s_WT_vs_TET_TKO.txt.gz", io$indir,i)
  if (file.exists(file)) {
    fread(file) %>% .[,celltype:=i]
  }
}) %>% rbindlist

# Merge celltypes 
diff.dt <- diff.dt %>% 
  .[,celltype:=stringr::str_replace_all(celltype,opts$rename_celltypes)] %>%
  .[,.(expr_ko=mean(expr_ko), expr_wt=mean(expr_wt), diff=mean(diff)), by=c("celltype","gene")]

# save
fwrite(diff.dt, file.path(io$outdir,"diff_pseudobulk.txt.gz"))

