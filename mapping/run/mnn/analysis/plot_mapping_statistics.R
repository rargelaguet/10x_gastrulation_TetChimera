#####################
## Define settings ##
#####################

source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
# source("/Users/ricard/10x_gastrulation_TetChimera/mapping/plot/plot_utils.R")

io$metadata <- "/Users/ricard/data/10x_gastrulation_TetChimera/processed/fourth_batch/sample_metadata.txt.gz"
io$mapping.results <- "/Users/ricard/data/10x_gastrulation_TetChimera/results/fourth_batch/mapping"
io$outdir <- paste0(io$basedir,"/results/fourth_batch/mapping/pdf")


opts$query_batches <- c(
  "SIGAC2_TET_TKO_E9_5_Head1_L002",
  "SIGAE6_TET_TKO_E9_5_Head2_L003",
  "SIGAF2_TET_TKO_E9_5_YS1_L002",
  "SIGAH6_TET_TKO_E9_5_YS2_L003",
  "SIGAE2_TET_TKO_E9_5_Tail1_L002",
  "SIGAG6_TET_TKO_E9_5_Tail2_L003",
  "SIGAD2_TET_TKO_E9_5_Trunk1_L002",
  "SIGAF6_TET_TKO_E9_5_Trunk2_L003"
)

###############
## Load data ##
###############

# Load sample metadata
sample_metadata <- fread(io$metadata) %>% .[pass_QC==TRUE]

# Load mapping results
# mapping.dt <- opts$query_batches %>% 
#   map(~ readRDS(sprintf("%s/mapping_mnn_%s.rds",io$mapping.results,.))[["mapping"]] ) %>%
#   rbindlist %>% merge(sample_metadata,by="cell")

# mapping <- fread(io$mapping)
# sample_metadata <- sample_metadata %>% merge(mapping)

################
## Parse data ##
################

to.plot <- sample_metadata %>%
  # .[class%in%opts$classes] %>%
  .[!is.na(celltype.mapped),.N, by=c("stage","celltype.mapped","batch","region")] %>%
  # .[!celltype.mapped%in%c("ExE_ectoderm","ExE_endoderm","Parietal_endoderm","Visceral_endoderm")] %>%
  # .[, celltype.mapped:=stringr::str_replace_all( celltype.mapped,"_"," ")] %>%
  # .[, celltype.mapped:=factor( celltype.mapped,levels=names(opts$celltype.colors))] 
  .[, celltype.mapped:=factor( celltype.mapped,levels=sort(names(opts$celltype.colors), decreasing = F))]

##########
## Plot ##
##########

for (i in unique(to.plot$region)) {
  p <- ggplot(to.plot[region==i], aes(x=celltype.mapped, y=N)) +
    geom_bar(aes(fill=celltype.mapped), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors, drop=F) + 
    scale_x_discrete(drop=FALSE) +
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
