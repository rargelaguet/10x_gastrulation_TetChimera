#####################
## Define settings ##
#####################

io <- list(); opts <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$standard.mnn.script <- "/Users/ricard/10x_gastrulation_TetChimera/iterative_mapping/run/standard_mnn.R"
  io$iterative.mnn.script <- "/Users/ricard/10x_gastrulation_TetChimera/iterative_mapping/run/iterative_mnn.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  io$standard.mnn.script <- "/homes/ricard/10x_gastrulation_TetChimera/iterative_mapping/run/standard_mnn.R"
  io$iterative.mnn.script <- "/homes/ricard/10x_gastrulation_TetChimera/iterative_mapping/run/iterative_mnn.R"
  io$tmpdir <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_TetChimera/results/iterative_mapping/tmp"
} 

# Atlas stages
opts$atlas_stages <- c(
  # "E6.5"
  # "E6.75",
  # "E7.0",
  # "E7.25",
  # "E7.5",
  # "E7.75",
  "E8.0",
  "E8.25",
  "E8.5"
  # "mixed_gastrulation"
)

# Query batches
opts$query_batches <- c(
  
  # E7.5
  # "E75_TET_TKO_L002",
  # "E75_WT_Host_L001",

  # E8.5
  "E85_Rep1_TET_TKO_L004",
  "E85_Rep2_TET_TKO_L006",
  "E85_Rep1_WT_Host_L003",
  "E85_Rep2_WT_Host_L005"

  # E12.5
  # "E125_DNMT3A_HET_A_L001",
  # "E125_DNMT3A_HET_A_L003",
  # "E125_DNMT3A_KO_B_L002",
  # "E125_DNMT3A_KO_E_L004"
)

# Test mode (subset cells)?
opts$test <- FALSE


#########
## Run ##
#########

for (i in opts$query_batches) {
  # LSF
  if (grepl("ricard",Sys.info()['nodename'])) {
    lsf <- ""
  } else {
    lsf <- sprintf("bsub -M 30000 -n 1 -q research-rh74 -o %s/%s.txt", io$tmpdir, i)
  }

  # Run standard MNN
  cmd <- sprintf("%s Rscript %s --query_batches %s --atlas_stages %s", lsf, io$standard.mnn.script, i, paste(opts$atlas_stages, collapse=" "))
  if (isTRUE(opts$test)) cmd <- paste0(cmd, " --test")
  system(cmd)

  # Run tree-guided MNN
  # cmd <- sprintf("%s Rscript %s --query_batches %s --atlas_stages %s", lsf, io$iterative.mnn.script, i, paste(opts$atlas_stages, collapse=" "))
  # if (isTRUE(opts$test)) cmd <- paste0(cmd, " --test")
  # system(cmd)
}
