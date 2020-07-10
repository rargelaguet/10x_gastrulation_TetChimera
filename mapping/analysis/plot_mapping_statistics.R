#####################
## Define settings ##
#####################

source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
# source("/Users/ricard/10x_gastrulation_TetChimera/mapping/plot/plot_utils.R")

# io$mapping <- "/Users/ricard/data/10x_gastrulation_TetChimera/mapping/sample_metadata_mapping_mnn.txt"
io$outdir <- paste0(io$basedir,"/results/mapping/pdf")

opts$aggregated.celltypes <- c(
  "Erythroid1" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid3" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors",
  "Rostral_neurectoderm" = "Neurectoderm",
  "Caudal_neurectoderm" = "Neurectoderm",
  "Anterior_Primitive_Streak" = "Primitive_Streak"
)

###############
## Load data ##
###############

# mapping <- fread(io$mapping)
# sample_metadata <- sample_metadata %>% merge(mapping)

################
## Parse data ##
################

to.plot <- sample_metadata %>%
  .[class%in%opts$classes] %>%
  # .[!is.na(celltype.mapped),.N, by=c("stage","celltype.mapped","batch","class")] %>% 
  .[!celltype.mapped%in%c("ExE_ectoderm","ExE_endoderm","Parietal_endoderm","Visceral_endoderm")] %>%
  .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$aggregated.celltypes)] %>%
  .[!is.na(celltype.mapped),.N, by=c("stage","celltype.mapped","class")] %>% 
  droplevels() %>% .[complete.cases(.)]


##########
## Plot ##
##########

to.plot[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,"_"," ")]
names(opts$celltype.colors) <- names(opts$celltype.colors) %>% stringr::str_replace_all("_"," ")
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$celltype.mapped)]
to.plot[,celltype.mapped:=factor(celltype.mapped, levels=names(opts$celltype.colors))]

for (i in unique(to.plot$stage)) {
  p <- ggplot(to.plot[stage==i], aes(x=celltype.mapped, y=N)) +
    geom_bar(aes(fill=celltype.mapped), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~class, nrow=1, scales="fixed") +
    coord_flip() +
    labs(y="Number of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(1.2)),
      axis.title.x = element_text(color="black", size=rel(1.1)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.3), color="black"),
      axis.text.x = element_text(size=rel(1.1), color="black")
    )
  
  pdf(sprintf("%s/mapping_stats_aggregated_%s.pdf",io$outdir,i), width=7, height=6)
  print(p)
  dev.off()
}
