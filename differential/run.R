#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
  io$script <- "/Users/ricard/10x_gastrulation_TetChimera/differential/differential.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/10x_gastrulation_TetChimera/settings.R")
  io$script <- "/homes/ricard/10x_gastrulation_TetChimera/differential/differential.R"
  io$tmpdir <- "/hps/nobackup2/research/stegle/users/ricard/10x_gastrulation_TetChimera/results/second_batch/differential/tmp"; dir.create(io$tmpdir, showWarnings=F)
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/results/second_batch/differential"); dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

# Statistical test
opts$statistical.test <- "edgeR"

# Testing mode
opts$test_mode <- FALSE

# Define groups
opts$groupA <- "WT"
opts$groupB <- "TET_TKO"

#########
## Run ##
#########

for (i in head(opts$celltypes.1,n=3)) {
# for (i in opts$celltypes.1) {
    outfile <- sprintf("%s/%s_WT_vs_TKO.txt.gz", io$outdir,i)
    
    # Define LSF command
    if (grepl("ricard",Sys.info()['nodename'])) {
      lsf <- ""
    } else if (grepl("ebi",Sys.info()['nodename'])) {
      lsf <- sprintf("bsub -M 15000 -n 1 -q research-rh74 -o %s/%s_vs_%s.txt", io$tmpdir,groupA,groupB)
    }
    cmd <- sprintf("%s Rscript %s --groupA %s --groupB %s --celltype %s --test %s --outfile %s", lsf, io$script, opts$groupA, opts$groupB, i, opts$statistical.test, outfile)
    if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")
    
    # Run
    print(cmd)
    system(cmd)
}

