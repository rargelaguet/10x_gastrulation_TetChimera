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
  "E8.5_TET_TKO"
  # "E10.5_Host", 
  # "E10.5_TET_TKO"
)


opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
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
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm"
  # "Visceral_endoderm",
  # "ExE_endoderm",
  # "ExE_ectoderm",
  # "Parietal_endoderm"
)

# Dot size
opts$size.mapped <- 0.18
opts$size.nomapped <- 0.1

# Transparency
opts$alpha.mapped <- 0.65
opts$alpha.nomapped <- 0.35

####################
## Load 10x atlas ##
####################

# Load atlas cell metadata
meta_atlas <- fread(io$atlas.metadata) %>%
  .[celltype%in%opts$celltypes] %>%
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

opts$stages <- c("E7.5","E8.5")

for (i in opts$stages) {
  
  to.plot <- umap.dt %>% copy %>%
    .[,index.wt:=match(cell, sample_metadata[class==paste0(i,"_Host"),closest.cell] )] %>%
    .[,index.ko:=match(cell, sample_metadata[class==paste0(i,"_TET_TKO"),closest.cell] )] %>%
    .[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]] %>%
    .[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]] %>%
    .[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>%
    .[,mapped:=plyr::mapvalues(mapped, from = c("0","-10","10"), to = c("Atlas","WT",i))] %>% setorder(mapped)
  
  p <- plot.dimred.wtko(to.plot, wt.label = "WT", ko.label = i, nomapped.label = "Atlas") +
    theme(legend.position = "none", axis.line = element_blank())
  
  pdf(sprintf("%s/umap_mapped_%s.pdf",io$outdir,i), width=5, height=6)
  print(p)
  dev.off()
}
