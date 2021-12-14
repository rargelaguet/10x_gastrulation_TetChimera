library(SoupX)
library(purrr)
library(data.table)
source(here::here("settings.R"))

# script to run the SoupX method to quantify and remove contamination from ambient RNA
# takes as input 10X CellRanger output files

#io <- list()
#io$basedir <- "/bi/scratch/Stephen_Clark/tet_chimera_10x/"

# output from 10x CellRanger can be found in the following folders:

io$indirs <- c("/bi/scratch/Felix/For_Stephen_Clark/200624_10X_tomato/Samples5019_26/",
               "/bi/scratch/Felix/For_Stephen_Clark/200624_10X_tomato/Samples5176_5181/",
               "/bi/scratch/Felix/For_Stephen_Clark/200624_10X_tomato/Samples5191_94/")

io$outdir <- paste0(io$basedir, "/raw_data")

# locate the 'outs' directory which contains the files we want. 
# there will be one 'outs' dir per 10X sample.

data_dirs <-  list.dirs(io$indirs, full = TRUE, recursive = TRUE) %>%
  .[grep("/outs/", .)] %>%
  gsub("(outs/).*", "\\1", .) %>%
  unique()

# iterate over the different 10X samples 

walk(data_dirs, ~{
  
  # run the SoupX commands
  sc <- load10X(.)
  sc <- autoEstCont(sc, forceAccept = TRUE)
  adj <- adjustCounts(sc)
  
  # now save the raw data and the adjusted data to a fresh output folder
  
  # find the raw data (filtered = no empty drops)
  filt <- dir(paste0(., "/filtered_feature_bc_matrix/"), full = TRUE)
  
  outdir <- paste0(io$outdir, "/", basename(dirname(.)))
  soupdir <- paste0(outdir, "/soup")
  dir.create(soupdir, recursive = TRUE)
  
  # save raw data
  file.copy(filt, paste0(outdir, "/", basename(filt)))
  
  # save adjusted data
  Matrix::writeMM(adj, paste0(soupdir, "/soupX_adjusted_matrix.mtx.gz"))
  soup <- setDT(sc$soupProfile, keep.rownames = "gene") %>%
    .[order(-rank(counts))]
  fwrite(soup, paste0(soupdir, "/soup.tsv.gz"), sep = "\t", na = "NA")
  
})
