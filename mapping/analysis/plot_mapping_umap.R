here::i_am("mapping/analysis/plot_mapping_umap.R")

source(here::here("settings.R"))
source(here::here("mapping/analysis/plot_utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--query_metadata',        type="character",                               help='Cell metadata (after mapping)')
p$add_argument('--stages',         type="character",       nargs="+",   help='Stages')
p$add_argument('--atlas_metadata',        type="character",                               help='Cell metadata (after mapping)')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
args$query_metadata <- file.path(io$basedir,"results_all/mapping/sample_metadata_after_mapping.txt.gz")
args$atlas_metadata <- file.path(io$atlas.basedir,"sample_metadata.txt.gz")
args$outdir <- file.path(io$basedir,"results_all/mapping/pdf")
args$stages <- c("E7.5","E8.5")
## END TEST ##

dir.create(args$outdir, showWarnings = F)
dir.create(file.path(args$outdir,"per_sample"), showWarnings = F)
dir.create(file.path(args$outdir,"per_class"), showWarnings = F)

#####################
## Define settings ##
#####################

# Options
opts$remove_ExE_cells <- TRUE
opts$subset_atlas <- TRUE

# Dot size
opts$size.mapped <- 0.18
opts$size.nomapped <- 0.1

# Transparency
opts$alpha.mapped <- 0.65
opts$alpha.nomapped <- 0.35

#########################
## Load query metadata ##
#########################

sample_metadata <- fread(args$query_metadata) %>%
  # .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(closest.cell)]
  .[pass_rnaQC==TRUE & !is.na(closest.cell) & stage%in%args$stages]

stopifnot("closest.cell"%in%colnames(sample_metadata))

if (opts$remove_ExE_cells) {
  print("Removing ExE cells...")
  sample_metadata <- sample_metadata %>% .[!celltype.mapped%in%c("ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

################
## Load atlas ##
################

# Load atlas cell metadata
meta_atlas <- fread(args$atlas_metadata) %>%
  # .[celltype%in%opts$celltypes] %>%
  .[stripped==F & doublet==F]

if (opts$remove_ExE_cells) {
  print("Removing ExE cells...")
  meta_atlas <- meta_atlas %>% .[!celltype%in%c("ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

if (opts$subset_atlas) {
  meta_atlas <- meta_atlas[sample.int(50000)]
}

# Extract precomputed dimensionality reduction coordinates
umap.dt <- meta_atlas %>%
  .[,c("cell","umapX","umapY","celltype")] %>%
  setnames(c("umapX","umapY"),c("V1","V2"))

####################
## Plot all cells ##
####################

# to.plot <- umap.dt %>% copy %>%
#   .[,index:=match(cell, sample_metadata[,closest.cell] )] %>% 
#   .[,mapped:=as.factor(!is.na(index))] %>% 
#   .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas cells","Query cells"))] %>%
#   setorder(mapped) 
# 
# p <- plot.dimred(to.plot, query.label = "Query cells", atlas.label = "Atlas cells")
# 
# pdf(sprintf("%s/umap_mapped_allcells.pdf",args$outdir), width=8, height=6.5)
# print(p)
# dev.off()

###############################
## Plot one sample at a time ##
###############################

samples.to.plot <- unique(sample_metadata$sample)

for (i in samples.to.plot) {
  
  to.plot <- umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata[sample==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
    setorder(mapped) 
  
  p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas")
  
  pdf(sprintf("%s/per_sample/umap_mapped_%s.pdf",args$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
}

###############################
## Plot one class at a time ##
###############################

classes.to.plot <- unique(sample_metadata$class)

for (i in classes.to.plot) {
  
  to.plot <- umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata[class==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
    setorder(mapped) 
  
  p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas") + theme(legend.position = "none")
  
  pdf(sprintf("%s/per_class/umap_mapped_%s.pdf",args$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
}


#############################
## Plot WT and KO together ##
#############################

# Subsample query cells to have the same N per class
# sample_metadata_subset <- sample_metadata %>% .[,.SD[sample.int(n=.N, size=4500)], by=c("stage","class2")]

# i <- "E7.5"
for (i in args$stages) {
  
  to.plot <- umap.dt %>% copy %>%
    .[,index.wt:=match(cell, sample_metadata[class2=="WT" & stage==i,closest.cell] )] %>%
    .[,index.ko:=match(cell, sample_metadata[class2=="TET_TKO" & stage==i,closest.cell] )] %>%
    .[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]] %>%
    .[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]] %>%
    .[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>%
    .[,mapped:=plyr::mapvalues(mapped, from = c("0","-10","10"), to = c("Atlas","WT","TET TKO"))] %>% setorder(mapped)
  
  p <- plot.dimred.wtko(to.plot, wt.label = "WT", ko.label = "TET TKO", nomapped.label = "Atlas") +
    theme(legend.position = "top", axis.line = element_blank())
  
  pdf(sprintf("%s/per_class/umap_mapped_%s_WT_and_KO.pdf",args$outdir,i), width=5.5, height=6.5)
  print(p)
  dev.off()
}
# Completion token
file.create(file.path(args$outdir,"completed.txt"))
