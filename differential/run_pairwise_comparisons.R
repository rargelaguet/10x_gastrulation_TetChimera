#########
## I/O ##
#########

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation10x/settings.R")
  io$script <- "/Users/ricard/gastrulation10x/differential/differential.R"
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation10x/settings.R")
  io$script <- "/homes/ricard/gastrulation10x/differential/differential.R"
  io$tmpdir <- "/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/results/differential/tmp"; dir.create(io$tmpdir, showWarnings=F)
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/results/differential"); dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

# Define cell types
opts$groups <- opts$celltypes.1

# Statistical test
# opts$test <- c("edgeR","t-test")
opts$statistical.test <- c("edgeR")

# Testing mode
opts$test_mode <- FALSE

#########
## Run ##
#########

# for (test in opts$statistical.test) {
#   for (i in 1:(length(opts$groups)-1)) {
#     groupA <- opts$groups[[i]]
#     for (j in 1:length(opts$groups)) {
#       if (i!=j) {
#         groupB <- opts$groups[[j]]
#         outfile <- sprintf("%s/%s_vs_%s.txt.gz", io$outdir,groupA,groupB)# %>% 
#           # stringr::str_replace_all(.," ","-")
        
#         # Define LSF command
#         if (grepl("ricard",Sys.info()['nodename'])) {
#           lsf <- ""
#         } else if (grepl("ebi",Sys.info()['nodename'])) {
#           lsf <- sprintf("bsub -M 15000 -n 1 -q research-rh74 -o %s/%s_vs_%s.txt", io$tmpdir,groupA,groupB)
#         }
#         cmd <- sprintf("%s Rscript %s --groupA %s --groupB %s --test %s --outfile %s", lsf, io$script, groupA, groupB, test, outfile)
#         if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")
        
#         # Run
#         print(cmd)
#         system(cmd)
#       }
#     }
#   }
# }


##############################
## Run selected comparisons ##
##############################

# opts$comparisons <- list(
#   c("groupA"="Erythroid1",          "groupB"="Mixed_mesoderm"),
#   c("groupA"="Blood_progenitors_2", "groupB"="Caudal_Mesoderm"),
#   c("groupA"="Allantois",           "groupB"="Haematoendothelial_progenitors"),
#   c("groupA"="Blood_progenitors_2", "groupB"="NMP"),
#   c("groupA"="Erythroid1",          "groupB"="Anterior_Primitive_Streak"),
#   c("groupA"="Erythroid1",          "groupB"="Haematoendothelial_progenitors")
# )

opts$comparisons <- list(
  c("groupA"="Allantois",           "groupB"="Haematoendothelial_progenitors"),
  c("groupA"="Blood_progenitors_2", "groupB"="Haematoendothelial_progenitors"),
  c("groupA"="Blood_progenitors_2", "groupB"="Mesenchyme"),
  c("groupA"="Cardiomyocytes",      "groupB"="ExE_ectoderm"),
  c("groupA"="Cardiomyocytes",      "groupB"="ExE_endoderm"),
  c("groupA"="Cardiomyocytes",      "groupB"="Parietal_endoderm"),
  c("groupA"="Erythroid1",          "groupB"="Parietal_endoderm"),
  c("groupA"="Caudal_neurectoderm", "groupB"="Neural_crest"),
  c("groupA"="Endothelium",          "groupB"="Epiblast"),
  c("groupA"="Endothelium",          "groupB"="Visceral_endoderm")
)



for (comparison in opts$comparisons) {
  groupA <- comparison[["groupA"]]; groupB <- comparison[["groupB"]]
  for (test in opts$statistical.test) {
    outfile <- sprintf("%s/%s_vs_%s.txt.gz", io$outdir,groupA,groupB)

    # Define LSF command
    if (grepl("ricard",Sys.info()['nodename'])) {
      lsf <- ""
    } else if (grepl("ebi",Sys.info()['nodename'])) {
      lsf <- sprintf("bsub -M 15000 -n 1 -q research-rh74 -o %s/%s_vs_%s.txt", io$tmpdir,groupA,groupB)
    }
    cmd <- sprintf("%s Rscript %s --groupA %s --groupB %s --test %s --outfile %s", lsf, io$script, groupA, groupB, test, outfile)
    if (isTRUE(opts$test_mode)) cmd <- paste0(cmd, " --test_mode")

    # Run
    print(cmd)
    system(cmd)
  }
}

