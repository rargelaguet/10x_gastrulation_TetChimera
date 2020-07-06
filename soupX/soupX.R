library(SoupX)
library(purrr)
library(data.table)
source(here::here("settings.R"))

# script to run the SoupX method to quantify and remove contamination from ambient RNA

#io <- list()
#io$basedir <- "/bi/scratch/Stephen_Clark/tet_chimera_10x/"

io$indirs <- c("/bi/scratch/Felix/For_Stephen_Clark/200624_10X_tomato/Samples5019_26/",
               "/bi/scratch/Felix/For_Stephen_Clark/200624_10X_tomato/Samples5176_5181/",
               "/bi/scratch/Felix/For_Stephen_Clark/200624_10X_tomato/Samples5191_94/")

io$outdir <- paste0(io$basedir, "/raw_data")

data_dirs <-  list.dirs(io$indirs, full = TRUE, recursive = TRUE) %>%
  .[grep("/outs/", .)] %>%
  gsub("(outs/).*", "\\1", .) %>%
  unique()

.= data_dirs[1]

walk(data_dirs, ~{
  sc <- load10X(.)
  sc <- autoEstCont(sc, forceAccept = TRUE)
  adj <- adjustCounts(sc)
  
  filt <- dir(paste0(., "/filtered_feature_bc_matrix/"), full = TRUE)
  
  outdir <- paste0(io$outdir, "/", basename(dirname(.)))
  soupdir <- paste0(outdir, "/soup")
  dir.create(soupdir, recursive = TRUE)
  
  file.copy(filt, paste0(outdir, "/", basename(filt)))
  Matrix::writeMM(adj, paste0(soupdir, "/soupX_adjusted_matrix.mtx.gz"))
  soup <- setDT(sc$soupProfile, keep.rownames = "gene") %>%
    .[order(-rank(counts))]
  fwrite(soup, paste0(soupdir, "/soup.tsv.gz"), sep = "\t", na = "NA")
  
})
