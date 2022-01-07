#####################
## Define settings ##
#####################

source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")

################
## Define I/O ##
################

io$outdir <- paste0(io$basedir,"/results/celltype_proportions/barplots")
dir.create(paste0(io$outdir,"/per_sample"), showWarnings = F)
dir.create(paste0(io$outdir,"/per_class"), showWarnings = F)

####################
## Define options ##
####################

opts$classes <- c(
  "E7.5_Host", 
  "E7.5_TET_TKO", 
  "E8.5_Host", 
  "E8.5_TET_TKO",
  "E8.5_WT"
)

opts$to.merge <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Intermediate_mesoderm" = "Mixed_mesoderm",
  "Paraxial_mesoderm" = "Mixed_mesoderm",
  "Nascent_mesoderm" = "Mixed_mesoderm",
  "Pharyngeal_mesoderm" = "Mixed_mesoderm"
  # "Visceral_endoderm" = "ExE_endoderm"
)


###############
## Load data ##
###############

# mapping <- fread(io$mapping)
# sample_metadata <- sample_metadata %>% merge(mapping)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[pass_QC==TRUE & class%in%opts$classes & !is.na(celltype.mapped)] %>%
  .[class%in%opts$classes]
  
# Define cell type order
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% sample_metadata$celltype.mapped]
stopifnot(sort(unique(as.character(sample_metadata$celltype.mapped))) == sort(names(opts$celltype.colors)))
sample_metadata <- sample_metadata %>% 
  .[,celltype.mapped:=factor(celltype.mapped,levels=sort(names(opts$celltype.colors), decreasing = F))]

#####################################
## Calculate cell type proportions ##
#####################################

to.plot.sample <- sample_metadata %>% 
  .[,.N, by=c("celltype.mapped","sample","class")]# %>%
  # .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped," ", "_")] %>

to.plot.class <- sample_metadata %>%
  .[,.N, by=c("celltype.mapped","class")]


#####################
## Plot per class ##
#####################

for (i in unique(to.plot.class$class)) {
  
  p <- ggplot(to.plot.class[class==i], aes(x=celltype.mapped, y=N)) +
    geom_bar(aes(fill=celltype.mapped), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors, drop=FALSE) +
    scale_x_discrete(drop=FALSE) +
    coord_flip() +
    labs(y="Number of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(1.2)),
      axis.title.x = element_text(color="black", size=rel(1.1)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.1), color="black"),
      axis.text.x = element_text(size=rel(1.1), color="black")
    )
  
  pdf(sprintf("%s/per_class/barplots_%s.pdf",io$outdir,i), width=7, height=7)
  print(p)
  dev.off()
  
  p <- ggplot(to.plot.sample[class==i], aes(x=celltype.mapped, y=N)) +
    geom_bar(aes(fill=celltype.mapped), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors, drop=FALSE) +
    scale_x_discrete(drop=FALSE) +
    facet_wrap(~sample, nrow=1) +
    coord_flip() +
    labs(y="Number of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(0.9)),
      axis.title.x = element_text(color="black", size=rel(1.1)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.1), color="black"),
      axis.text.x = element_text(size=rel(1.1), color="black")
    )
  
  pdf(sprintf("%s/per_class/barplots_facet_sample_%s.pdf",io$outdir,i))
  print(p)
  dev.off()
}

#####################
## Plot per sample ##
#####################

for (i in unique(to.plot.sample$sample)) {
  
  p <- ggplot(to.plot.sample[sample==i], aes(x=celltype.mapped, y=N)) +
    geom_bar(aes(fill=celltype.mapped), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors, drop=FALSE) +
    scale_x_discrete(drop=FALSE) +
    coord_flip() +
    labs(y="Number of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(1.2)),
      axis.title.x = element_text(color="black", size=rel(1.1)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.1), color="black"),
      axis.text.x = element_text(size=rel(1.1), color="black")
    )
  
  pdf(sprintf("%s/per_sample/barplots_%s.pdf",io$outdir,i), width=7, height=7)
  print(p)
  dev.off()
}
