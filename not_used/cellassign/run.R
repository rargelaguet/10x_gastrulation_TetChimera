
################
## Define I/O ##
################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
  io$script <- "/Users/ricard/10x_gastrulation_TetChimera/cellassign/cellassign.R"
} else {
  source("/homes/ricard/10x_gastrulation_TetChimera/settings.R")
  io$script <- "/homes/ricard/10x_gastrulation_TetChimera/cellassign/cellassign.R"
  io$tmpdir <- paste0(io$basedir,"/results/second_batch/cellassign/tmp"); dir.create(io$tmpdir)
}
io$outdir <- paste0(io$basedir,"/results/second_batch/cellassign")

####################
## Define options ##
####################

opts$batch <- c(
  "E75_TET_TKO_L002",
  "E75_WT_Host_L001"
  # "E85_Rep1_TET_TKO_L004",
  # "E85_Rep2_TET_TKO_L006",
  # "E85_Rep1_WT_Host_L003",
  # "E85_Rep2_WT_Host_L005"
  # "E125_DNMT3A_HET_A_L001",
  # "E125_DNMT3A_HET_A_L003",
  # "E125_DNMT3A_KO_B_L002",
  # "E125_DNMT3A_KO_E_L004"
)

# Test mode (subsetting cells)?
opts$test_mode <- TRUE

#########
## Run ##
#########

for (i in opts$batch) {
  
  # Define LSF command
  if (grepl("ricard",Sys.info()['nodename'])) {
    lsf <- ""
  } else if (grepl("ebi",Sys.info()['nodename'])) {
    lsf <- sprintf("bsub -M 75000 -n 1 -o %s/%s.txt", io$tmpdir,paste(i,collapse=" "))
  }
  cmd <- sprintf("%s Rscript %s --batch %s --outdir %s", lsf, io$script, paste(i,collapse=" "),io$outdir)
  if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test")
  
  # Run
  print(cmd)
  system(cmd)
}
