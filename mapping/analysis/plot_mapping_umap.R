###################################################################
## Plot dimensionality reduction of EB cells mapped to the atlas ##
###################################################################

source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
source("/Users/ricard/10x_gastrulation_TetChimera/mapping/analysis/plot_utils.R")

io$outdir <- paste0(io$basedir,"/results/mapping/pdf")

#############
## Options ##
#############

opts$classes <- c(
  "E7.5_Host", 
  "E7.5_TET_TKO", 
  "E8.5_Host", 
  "E8.5_TET_TKO", 
  "E10.5_Host", 
  "E10.5_TET_TKO"
)

# Dot size
opts$size.mapped <- 0.3
opts$size.nomapped <- 0.1

# Transparency
opts$alpha.mapped <- 0.9
opts$alpha.nomapped <- 0.35

####################
## Load 10x atlas ##
####################

# Load atlas cell metadata
meta_atlas <- fread(io$atlas.metadata) %>%
  .[stripped==F & doublet==F]

# Extract precomputed dimensionality reduction coordinates
umap.dt <- meta_atlas[,c("cell","umapX","umapY","celltype")] %>%
  setnames(c("umapX","umapY"),c("V1","V2"))

########################################################
## Plot dimensionality reduction: one batch at a time ##
########################################################

for (i in opts$classes) {
  to.plot <- umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata[class==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
    setorder(mapped) 
  
  p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas")
  
  pdf(sprintf("%s/umap_mapped_%s.pdf",io$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
}

######################################################
## Plot dimensionality reduction: WT vs KO together ##
######################################################

for (i in opts$classes) {
  to.plot <- umap.dt %>% copy %>%
    .[,index.wt:=match(cell, sample_metadata[class=="E8.5_Dnmt3aWT_Dnmt3bWT",closest.cell] )] %>%
    .[,index.ko:=match(cell, sample_metadata[class==i,closest.cell] )] %>%
    .[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]] %>%
    .[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]] %>%
    .[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>%
    .[,mapped:=plyr::mapvalues(mapped, from = c("0","-10","10"), to = c("Atlas","WT",i))] %>% setorder(mapped)
  
  p <- plot.dimred.wtko(to.plot, wt.label = "WT", ko.label = i, nomapped.label = "Atlas")
  
  pdf(sprintf("%s/umap_mapped_%s.pdf",io$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
}