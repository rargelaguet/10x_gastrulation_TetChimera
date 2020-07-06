#####################
## Define settings ##
#####################

source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
# source("/Users/ricard/10x_gastrulation_TetChimera/mapping/plot/plot_utils.R")

# io$mapping <- "/Users/ricard/data/10x_gastrulation_TetChimera/mapping/sample_metadata_mapping_mnn.txt"
io$outdir <- paste0(io$basedir,"/results/first_batch/mapping/pdf")

###############
## Load data ##
###############

io$mapping.results <- "/Users/ricard/data/10x_gastrulation_TetChimera/results/first_batch"

opts$query_batches <- c(
  "SIGAA3_E8.5_pool1_Host-WT_L001",
  "SIGAB3_E8.5_pool1_TET-TKO_L002",
  "SIGAC3_E8.5_pool2_Host-WT_L003",
  "SIGAD3_E8.5_pool2_TET-TKO_L004", 
  "SIGAE3_E7.5_pool1_Host-WT_L005", 
  "SIGAF3_E7.5_pool1_TET-TKO_L006", 
  "SIGAG3_E8.5_hashing_Host-WT_L007",
  "SIGAH3_E8.5_hasting_TET-TKO_L008"
)

sample_metadata <- fread("/Users/ricard/data/10x_gastrulation_TetChimera/processed/first_batch/sample_metadata.txt.gz")
mapping.dt <- opts$query_batches %>% 
  map(~ fread(sprintf("%s/mapping_mnn_%s.txt.gz",io$mapping.results,.)) ) %>%
  rbindlist %>% merge(sample_metadata,by="cell")

# mapping <- fread(io$mapping)
# sample_metadata <- sample_metadata %>% merge(mapping)

################
## Parse data ##
################

to.plot <- mapping.dt %>%
  # .[class%in%opts$classes] %>%
  .[!is.na(celltype.mapped),.N, by=c("stage","celltype.mapped","batch")] %>%
  # .[!celltype.mapped%in%c("ExE_ectoderm","ExE_endoderm","Parietal_endoderm","Visceral_endoderm")] %>%
  # .[, celltype.mapped:=stringr::str_replace_all( celltype.mapped,"_"," ")] %>%
  # .[, celltype.mapped:=factor( celltype.mapped,levels=names(opts$celltype.colors))] 
  .[, celltype.mapped:=factor( celltype.mapped,levels=sort(names(opts$celltype.colors), decreasing = F))]

##########
## Plot ##
##########

for (i in unique(to.plot$stage)) {
  p <- ggplot(to.plot[stage==i], aes(x=celltype.mapped, y=N)) +
    geom_bar(aes(fill=celltype.mapped), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors) +
    facet_wrap(~batch, scales="free_x") +
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
  
  pdf(sprintf("%s/mapping_stats_%s.pdf",io$outdir,i), width=14, height=14)
  print(p)
  dev.off()
}
