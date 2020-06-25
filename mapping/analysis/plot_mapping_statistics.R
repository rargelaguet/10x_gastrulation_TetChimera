#####################
## Define settings ##
#####################

source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
# source("/Users/ricard/10x_gastrulation_TetChimera/mapping/plot/plot_utils.R")

# io$mapping <- "/Users/ricard/data/10x_gastrulation_TetChimera/mapping/sample_metadata_mapping_mnn.txt"
io$outdir <- paste0(io$basedir,"/results/mapping/pdf")

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
  .[!is.na(celltype.mapped),.N, by=c("stage","celltype.mapped","batch","class")] %>% 
  .[!celltype.mapped%in%c("ExE_ectoderm","ExE_endoderm","Parietal_endoderm","Visceral_endoderm")] %>%
  # .[, celltype.mapped:=stringr::str_replace_all( celltype.mapped,"_"," ")] %>%
  # .[, celltype.mapped:=factor( celltype.mapped,levels=names(colors))] %>%
  .[complete.cases(.)]

##########
## Plot ##
##########

for (i in unique(to.plot$stage)) {
  p <- ggplot(to.plot[stage==i], aes(x=celltype.mapped, y=N)) +
    geom_bar(aes(fill=celltype.mapped), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~batch, nrow=1, scales="fixed") +
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
  
  pdf(sprintf("%s/mapping_stats_%s.pdf",io$outdir,i), width=14, height=7)
  print(p)
  dev.off()
}
