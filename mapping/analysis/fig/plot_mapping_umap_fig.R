here::i_am("mapping/analysis/plot_mapping_umap.R")

source(here::here("settings.R"))
source(here::here("mapping/analysis/plot_utils.R"))

io$query_metadata <- file.path(io$basedir,"results/mapping/sample_metadata_after_mapping.txt.gz")
io$atlas_metadata <- file.path(io$atlas.basedir,"sample_metadata.txt.gz")
io$outdir <- file.path(io$basedir,"results/mapping/pdf/fig"); dir.create(io$outdir, showWarnings = F)

#####################
## Define settings ##
#####################

# Options
opts$remove_ExE_cells <- FALSE
opts$subset_atlas <- TRUE

# Dot size
opts$size.mapped <- 0.18
opts$size.nomapped <- 0.1

# Transparency
opts$alpha.mapped <- 0.65
opts$alpha.nomapped <- 0.35

opts$samples <- c(
  "E7.5_Tet_TKO",
  "E7.5_WT_tdTomato-_1","E7.5_WT_tdTomato-_2",
  "E7.5_WT_tdTomato+_1","E7.5_WT_tdTomato+_2",
  "E8.5_WT_tdTomato+_1","E8.5_WT_tdTomato+_2","E8.5_WT_tdTomato+_3",
  "E8.5_WT_tdTomato-_1","E8.5_WT_tdTomato-_2","E8.5_WT_tdTomato-_3",
  "E8.5_Tet_TKO_2","E8.5_Tet_TKO_3","E8.5_Tet_TKO_4"
)


#########################
## Load query metadata ##
#########################

sample_metadata <- fread(io$query_metadata) %>%
  .[,alias:=factor(alias,levels=opts$samples)] %>%
  .[pass_rnaQC==TRUE & alias%in%opts$samples & !is.na(closest.cell)]

if (opts$remove_ExE_cells) {
  print("Removing ExE cells...")
  sample_metadata <- sample_metadata %>% .[!celltype.mapped%in%c("ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

################
## Load atlas ##
################

# Load atlas cell metadata
meta_atlas <- fread(io$atlas_metadata) %>%
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

##########################
## Plot E7.5 WT samples ##
##########################

opts$samples <- c(
  "E7.5_WT_tdTomato-_1", "E7.5_WT_tdTomato-_2",
  "E7.5_WT_tdTomato+_1", "E7.5_WT_tdTomato+_2"
)

# ncells.to.subset <- 2500

to.plot <- opts$samples %>% map(function(i) {
  # sample_metadata.subset <- sample_metadata[alias==i] %>% .[sample.int(min(ncells.to.subset,nrow(.)))]
  sample_metadata.subset <- sample_metadata[alias==i]
  umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata.subset[alias==i,closest.cell] )] %>%
    .[,mapped:=as.factor(!is.na(index))] %>%
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas","Query"))] %>%
    .[,sample:=factor(i,levels=opts$samples)] %>%
    setorder(mapped)
}) %>% rbindlist

p <- plot.dimred(to.plot, query.label = "Query", atlas.label = "Atlas") +
  facet_wrap(~sample, nrow=2) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(color="black", size=rel(0.8))
  )

pdf(file.path(io$outdir,"umap_mapped_E7.5_WT_samples.pdf"), width=8, height=8)
print(p)
dev.off()

#######################
## Plot E8.5 samples ##
#######################

opts$samples <- c(
  "E8.5_WT_tdTomato+_2","E8.5_WT_tdTomato+_3", # "E8.5_WT_tdTomato+_1",
  "E8.5_WT_tdTomato-_1","E8.5_WT_tdTomato-_2",#"E8.5_WT_tdTomato-_3",
  "E8.5_Tet_TKO_2","E8.5_Tet_TKO_3"#"E8.5_Tet_TKO_4"
)

# ncells.to.subset <- 2500

to.plot <- opts$samples %>% map(function(i) {
  # sample_metadata.subset <- sample_metadata[alias==i] %>% .[sample.int(min(ncells.to.subset,nrow(.)))]
  sample_metadata.subset <- sample_metadata[alias==i]
  umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata.subset[alias==i,closest.cell] )] %>%
    .[,mapped:=as.factor(!is.na(index))] %>%
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas","Query"))] %>%
    .[,sample:=factor(i,levels=opts$samples)] %>%
    setorder(mapped)
}) %>% rbindlist

p <- plot.dimred(to.plot, query.label = "Query", atlas.label = "Atlas") +
  facet_wrap(~sample, nrow=3) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(color="black", size=rel(0.8))
  )

pdf(file.path(io$outdir,"umap_mapped_E8.5_samples.pdf"), width=8, height=12)
print(p)
dev.off()

###############################
## Plot one sample at a time ##
###############################

# Subset classes to similar number of cells
# ncells.to.subset <- min(sample_metadata[,.N,by="alias"][["N"]])
# ncells.to.subset <- 5000

# for (i in head(opts$samples,n=2)) {
#   
#   # sample_metadata.subset <- sample_metadata[alias==i] %>% .[sample.int(min(ncells.to.subset,nrow(.)))]
#   sample_metadata.subset <- sample_metadata[alias==i]
#   
#   to.plot <- umap.dt %>% copy %>%
#     .[,index:=match(cell, sample_metadata.subset[,closest.cell] )] %>% 
#     .[,mapped:=as.factor(!is.na(index))] %>% 
#     .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
#     setorder(mapped) 
#   
#   p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas") + 
#     labs(title=i) +
#     theme(
#       plot.title = element_text(hjust=0.5),
#       legend.position = "none"
#     )
#   
#   pdf(sprintf("%s/umap_mapped_%s.pdf",io$outdir,i), width=6, height=5)
#   print(p)
#   dev.off()
# }

#############################
## Plot WT and KO together ##
#############################

# Subsample query cells to have the same N per class
sample_metadata_subset <- sample_metadata %>% .[,.SD[sample.int(n=.N, size=4500)], by=c("stage","class2")]

# i <- "E8.5"
for (i in opts$stages) {
  
  to.plot <- umap.dt %>% copy %>%
    .[,index.wt:=match(cell, sample_metadata[class2=="WT" & stage==i,closest.cell] )] %>%
    .[,index.ko:=match(cell, sample_metadata[class2=="TET_TKO" & stage==i,closest.cell] )] %>%
    .[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]] %>%
    .[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]] %>%
    .[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>%
    .[,mapped:=plyr::mapvalues(mapped, from = c("0","-10","10"), to = c("Atlas","TET TKO","WT"))] %>% setorder(mapped)
  
  p <- plot.dimred.wtko(to.plot, wt.label = "WT", ko.label = "TET TKO", nomapped.label = "Atlas") +
    theme(legend.position = "top", axis.line = element_blank())
  
  pdf(sprintf("%s/umap_mapped_%s_WT_and_KO.pdf",io$outdir,i), width=6.5, height=6.5)
  print(p)
  dev.off()
}
