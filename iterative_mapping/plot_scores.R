
#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_spatial/settings.R")
} else if(grepl("ebi",Sys.info()['nodename'])){
  source("/homes/ricard/gastrulation_spatial/settings.R")
}  
io$mapping.dir <- paste0(io$basedir,"/iterative_mapping")
io$outdir <- paste0(io$basedir,"/iterative_mapping/pdf")

opts$embryos <- c("embryo1","embryo2","embryo3")

###################
## Load metadata ##
###################

sample_metadata <- fread(io$metadata) %>%
  .[embryo%in%opts$embryos] 

################
## Load data  ##
################

mapping_dt <- opts$embryos %>% map(function(i) {
  rbind(
    fread(sprintf("%s/%s_standard_mnn.txt.gz",io$mapping.dir,i)) %>% .[,embryo:=i] %>% .[,method:="Standard MNN"],
    fread(sprintf("%s/%s_iterative_mnn.txt.gz",io$mapping.dir,i)) %>% .[,embryo:=i] %>% .[,method:="Tree-guided MNN"]
  )
}) %>% rbindlist

##########
## Plot ##
##########

for (i in opts$embryos) {
  to.plot <- mapping_dt[embryo==i] %>% 
    merge(sample_metadata[,c("cell","embryo","z")], by=c("cell","embryo")) %>% 
    .[celltype_mapped=="Forebrain Midbrain Hindbrain",celltype_mapped:="Forebrain/Midbrain/Hindbrain"]
  
  p1 <- ggplot(to.plot, aes(x=method, y=mapping_score)) +
    geom_boxplot(aes(fill=method), color="black", outlier.shape=NA) +
    facet_wrap(~celltype_mapped, scales="fixed") +
    labs(x="", y="Mapping score") +
    theme_classic() +
    theme(
      legend.title = element_blank(),
      legend.position = "top",
      axis.text.y = element_text(color="black", size=rel(0.8)),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text = element_text(color="black", size=rel(0.6))
    )
  
  p2 <- ggplot(to.plot, aes(x=method, y=mapping_score)) +
    geom_boxplot(aes(fill=method), color="black", outlier.shape=NA) +
    labs(x="", y="Mapping score") +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      # axis.text.x = element_text(color="black", size=rel(1.0))
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color="black", size=rel(1.0)),
    )
  
  p <- cowplot::plot_grid(plotlist=list(p2,p1), rel_widths=c(1/4,3/4))
  
  pdf(sprintf("%s/scores_%s.pdf",io$outdir,i), width=14, height=9, useDingbats = F)
  print(p)
  dev.off()
}
