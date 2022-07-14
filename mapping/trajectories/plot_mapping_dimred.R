# here::i_am("mapping/analysis/plot_mapping_umap.R")

source(here::here("settings.R"))
source(here::here("mapping/trajectories/plot_utils.R"))


######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--query_metadata',        type="character",                               help='Cell metadata (after mapping)')
# p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
p$add_argument('--atlas_metadata',        type="character",                               help='Cell metadata (after mapping)')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
args$query_metadata <- file.path(io$basedir,"results/mapping/trajectories/blood/sample_metadata_after_mapping.txt.gz")
args$atlas_metadata <- file.path(io$atlas.basedir,"results/trajectories/blood_precomputed/blood_trajectory.txt.gz")
args$outdir <- file.path(io$basedir,"results/mapping/trajectories/blood/pdf/revision")
## END TEST ##

dir.create(args$outdir, showWarnings = F)
dir.create(file.path(args$outdir,"per_sample"), showWarnings = F)
dir.create(file.path(args$outdir,"per_class"), showWarnings = F)

#####################
## Define settings ##
#####################

# io$trajectory <- file.path(io$atlas.basedir,"results/trajectories/blood_scanpy/blood_trajectory.txt.gz")
io$trajectory <- file.path(io$atlas.basedir,"results/trajectories/blood_precomputed/blood_trajectory.txt.gz")

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
  .[,class:=ifelse(grepl("WT",class),"WT","Tet-TKO")] %>%
  .[!is.na(closest.cell)]

stopifnot("closest.cell"%in%colnames(sample_metadata))

################
## Load atlas ##
################

# Load atlas trajectory
# atlas_trajectory.dt <- fread(io$trajectory) %>% setnames(c("FA1","FA2"),c("V1","V2"))
atlas_trajectory.dt <- fread(io$trajectory)[,c("cell","V1","V2")]
meta_atlas <- fread(args$atlas_metadata)[,c("cell","stage","celltype")] %>% merge(atlas_trajectory.dt, by="cell")

# Load atlas cell metadata
# meta_atlas <- fread(args$atlas_metadata) %>% setnames(c("DC1","DC2"),c("V1","V2"))

stopifnot(sample_metadata$closest.cell%in%meta_atlas$cell)

###############################
## Plot one sample at a time ##
###############################

samples.to.plot <- unique(sample_metadata$sample)

for (i in samples.to.plot) {
  
  to.plot <- meta_atlas %>% copy %>%
    .[,index:=match(cell, sample_metadata[sample==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
    setorder(mapped) 
  
  p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas")
  
  pdf(sprintf("%s/per_sample/umap_mapped_%s.pdf",args$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
}

##############################
## Plot one class at a time ##
##############################

classes.to.plot <- unique(sample_metadata$class)

for (i in classes.to.plot) {
  
  to.plot <- meta_atlas %>% copy %>%
    .[,index:=match(cell, sample_metadata[class==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
    setorder(mapped) 
  
  p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas") + theme(legend.position = "none")
  
  pdf(sprintf("%s/per_class/umap_mapped_%s.pdf",args$outdir,i), width=6, height=6)
  print(p)
  dev.off()
}

#############################
## Plot WT and KO together ##
#############################

opts$size.mapped <- 1
opts$size.nomapped <- 0.08

# Transparency
opts$alpha.mapped <- 0.75
opts$alpha.nomapped <- 0.10

# Subsample query cells to have the same N per class
sample_metadata_subset <- sample_metadata %>% .[,.SD[sample.int(n=.N, size=1500)], by="class"]

to.plot <- meta_atlas %>% copy %>%
  .[,index.wt:=match(cell, sample_metadata_subset[class=="WT",closest.cell] )] %>%
  .[,index.ko:=match(cell, sample_metadata_subset[class=="Tet-TKO",closest.cell] )] %>%
  .[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]] %>%
  .[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]] %>%
  .[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>%
  # .[,mapped:=plyr::mapvalues(mapped, from = c("0","-10","10"), to = c("Atlas","WT","Tet-TKO"))] %>% setorder(mapped)
  .[,mapped:=plyr::mapvalues(mapped, from = c("10","-10","0"), to = c("Tet-TKO","WT","Atlas"))] %>% setorder(mapped)

p <- plot.dimred.wtko(to.plot, wt.label = "WT", ko.label = "Tet-TKO", nomapped.label = "Atlas") +
  scale_colour_manual(values=opts$class2_colors[c("WT","Tet-TKO")]) +
  guides(size="none", alpha="none") +
  theme(legend.position = "none", axis.line = element_blank())

pdf(sprintf("%s/per_class/umap_mapped_WT_and_KO.pdf",args$outdir), width=6, height=6)
print(p)
dev.off()

# Completion token
file.create(file.path(args$outdir,"completed.txt"))
