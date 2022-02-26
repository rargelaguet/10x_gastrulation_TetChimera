here::i_am("celltype_proportions/plot_celltype_proportions.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

## START TEST ##
io$outdir <- file.path(io$basedir,"results/celltype_proportions/fig")
## END TEST ##

# I/O
dir.create(io$outdir, showWarnings = F)

opts$celltypes = c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  # "PGC",
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
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm"
  # "Parietal_endoderm"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & celltype.mapped%in%opts$celltypes] %>%
  setnames("celltype.mapped","celltype")

###################
## Plot E7.5 TKO ##
###################

opts$samples <- c(
  "E7.5_Tet_TKO"
)
stopifnot(opts$samples%in%unique(sample_metadata$alias))

to.plot <- sample_metadata %>%
  .[alias%in%opts$samples] %>%
  .[,alias:=factor(alias,levels=opts$samples)] %>%
  .[,N:=.N,by="class"] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("class","alias","celltype")] 

# Define colours and cell type order
celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$celltype)]
to.plot[,celltype:=factor(celltype, levels=rev(names(celltype.colors)))]

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill=celltype), stat="identity", color="black",) +
  scale_fill_manual(values=celltype.colors, drop=F) +
  facet_wrap(~alias, nrow=2, scales="free_x") +
  coord_flip() +
  labs(y="Number of cells") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(color="black", size=rel(0.8)),
    axis.title.x = element_text(color="black", size=rel(0.9)),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size=rel(0.75), color="black"),
    axis.text.x = element_text(size=rel(1), color="black")
  )

pdf(file.path(io$outdir,"E7.5_TET_TKO_celltype_proportions_horizontal_barplots.pdf"), width=5, height=3.75)
print(p)
dev.off()


##################
## Plot E7.5 WT ##
##################

opts$samples <- c(
  # "E7.5_Tet_TKO"
  "E7.5_WT_tdTomato-_1","E7.5_WT_tdTomato-_2",#,"E7.5_WT_tdTomato-_3",
  "E7.5_WT_tdTomato+_1","E7.5_WT_tdTomato+_2"#,"E7.5_WT_tdTomato+_3"
)
stopifnot(opts$samples%in%unique(sample_metadata$alias))

to.plot <- sample_metadata %>%
  .[alias%in%opts$samples] %>%
  .[,alias:=factor(alias,levels=opts$samples)] %>%
  .[,N:=.N,by="class"] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("class","alias","celltype")] 

# Define colours and cell type order
celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$celltype)]
to.plot[,celltype:=factor(celltype, levels=rev(names(celltype.colors)))]

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill=celltype), stat="identity", color="black",) +
  scale_fill_manual(values=celltype.colors, drop=F) +
  facet_wrap(~alias, nrow=2, scales="free_x") +
  coord_flip() +
  labs(y="Number of cells") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(color="black", size=rel(0.8)),
    axis.title.x = element_text(color="black", size=rel(0.9)),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size=rel(0.75), color="black"),
    axis.text.x = element_text(size=rel(1), color="black")
  )

# pdf(file.path(io$outdir,"E7.5_TET_TKO_celltype_proportions_horizontal_barplots.pdf"), width=6, height=5)
pdf(file.path(io$outdir,"E7.5_WT_celltype_proportions_horizontal_barplots.pdf"), width=7, height=6)
print(p)
dev.off()


###############
## Plot E8.5 ##
###############

opts$samples <- c(
  "E8.5_WT_tdTomato+_2","E8.5_WT_tdTomato+_3", # "E8.5_WT_tdTomato+_1",
  "E8.5_WT_tdTomato-_1","E8.5_WT_tdTomato-_2",#"E8.5_WT_tdTomato-_3",
  "E8.5_Tet_TKO_2","E8.5_Tet_TKO_3"#"E8.5_Tet_TKO_4"
)

to.plot <- sample_metadata %>%
  .[alias%in%opts$samples] %>%
  .[,alias:=factor(alias,levels=opts$samples)] %>%
  .[,N:=.N,by="class"] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("class","alias","celltype")] 

# Define colours and cell type order
celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$celltype)]
to.plot[,celltype:=factor(celltype, levels=rev(names(celltype.colors)))]

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill=celltype), stat="identity", color="black",) +
  scale_fill_manual(values=celltype.colors, drop=F) +
  facet_wrap(~alias, nrow=3, scales="free_x") +
  coord_flip() +
  labs(y="Number of cells") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(color="black", size=rel(0.8)),
    axis.title.x = element_text(color="black", size=rel(0.9)),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size=rel(0.75), color="black"),
    axis.text.x = element_text(size=rel(1), color="black")
  )
  
pdf(file.path(io$outdir,"E8.5_celltype_proportions_horizontal_barplots.pdf"), width=7, height=10)
print(p)
dev.off()
