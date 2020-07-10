##############
## Settings ##
##############

source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")
source("/Users/ricard/10x_gastrulation_TetChimera/differential/analysis/utils.R")
io$diff.dir <- paste0(io$basedir,"/results/differential")
io$outdir <- paste0(io$basedir,"/results/differential/pdf")

# Define groups
opts$comparisons <- c(
  "E7.5_Host" = "E7.5_TET_TKO", 
  "E8.5_Host" = "E8.5_TET_TKO"
)

opts$min.cells.per.group <- 50

#############################################
## Load results from differential analysis ##
#############################################

dt <- names(opts$comparisons) %>% map(function(i) { opts$celltypes %>% map(function(j) {
  file <- sprintf("%s/%s_%s_vs_%s.txt.gz", io$diff.dir,j,i,opts$comparisons[[i]])
  if (file.exists(file)) {
    fread(file) %>% 
      .[,c(1,2,4,6,7,10,11)] %>% 
      setnames(c(sprintf("N_%s",i),sprintf("N_%s",opts$comparisons[[i]])),c("N_A","N_B")) %>%
      .[N_A>opts$min.cells.per.group & N_B>opts$min.cells.per.group] %>%
      .[,c("celltype","groupA","groupB"):=list(j,opts$comparisons[[i]],i)]
  }
}) %>% rbindlist }) %>% rbindlist

unique(dt$groupA)
unique(dt$groupB)
unique(dt$celltype)

# Remove some hits
dt <- dt[!gene%in%c("Xist", "Hbb-y", "Hbb-bh1", "Hbb-bs", "Hba-x", "Hba-a2", "Hba-a1", "Hbq1b", "Hbb-bt", "Hbb-bh2")]

# Remove some cell types
dt <- dt[!celltype%in%c("ExE_ectoderm","ExE_endoderm","Parietal_endoderm","Visceral_endoderm")]

##########
## Plot ##
##########

to.plot <- dt %>% copy %>%
  .[,.(N=sum(sig,na.rm=T)) ,by=c("celltype","groupA","groupB")]

ggplot(to.plot, aes(x=factor(celltype), y=N)) +
  geom_bar(aes(fill = celltype, color=celltype), stat = 'identity') + 
  # geom_polygon(color="black", fill=NA, alpha=0.5, linetype="dashed", data=foo) +
  facet_wrap(~groupA, nrow=2) +
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  scale_color_manual(values=opts$celltype.colors, drop=F, guide = 'none') +
  coord_polar() +
  guides(fill = guide_legend(ncol=1)) +
  theme_bw() +
  theme(
    legend.key.size = unit(0.25,"cm"),
    legend.position = "right",
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.title=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.line=element_blank(),
    axis.text.x = element_blank()
    # axis.text.x = element_text(angle= -76 - 360 / length(unique(to.plot$celltype)) * seq_along(to.plot$celltype))
  )

# for (i in unique(dt$celltype)) {
#   to.plot <- dt[celltype==i] %>% .[!is.na(sig)] 
#   p <- gg_volcano_plot(to.plot, top_genes=15)
#   
#   pdf(sprintf("%s/%s_vs_%s_%s_volcano.pdf",io$outdir,opts$groupA,opts$groupB,i), width=9, height=5, useDingbats = F)
#   print(p)
#   dev.off()
# }



# to.plot <- dt %>% copy %>%
#   .[!is.na(logFC) & sig==T] %>%
#   .[,sign:=c("up in Mutant", "up in WT")[as.numeric(logFC>0)+1]] %>%
#   .[,.N,by=c("sign","celltype","groupA","groupB")]