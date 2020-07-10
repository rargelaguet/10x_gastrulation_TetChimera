library(ggplot2)
linMap <- function(x, from, to) return( (x - min(x)) / max(x - min(x)) * (to - from) + from )

source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")

################
## Define I/O ##
################

io$outdir <- paste0(io$basedir,"/results/mapping/pdf")

####################
## Define options ##
####################

# Figure dimensions
# opts$width <- 3
# opts$height <- 4

opts$classes <- c(
  "E7.5_Host", 
  "E7.5_TET_TKO", 
  "E8.5_Host", 
  "E8.5_TET_TKO"
  # "E10.5_Host", 
  # "E10.5_TET_TKO"
)

opts$wt.classes <- c("E7.5_Host","E8.5_Host")
  
# opts$to.merge <- c(
#   "Erythroid3" = "Erythroid",
#   "Erythroid2" = "Erythroid",
#   "Erythroid1" = "Erythroid",
#   "Blood_progenitors_1" = "Blood_progenitors",
#   "Blood_progenitors_2" = "Blood_progenitors",
#   "Intermediate_mesoderm" = "Mixed_mesoderm",
#   "Paraxial_mesoderm" = "Mixed_mesoderm",
#   "Nascent_mesoderm" = "Mixed_mesoderm",
#   "Pharyngeal_mesoderm" = "Mixed_mesoderm",
#   "Visceral_endoderm" = "ExE_endoderm"
# )

# Update sample metadata
sample_metadata <- sample_metadata %>%
  .[class%in%opts$classes & !is.na(celltype.mapped)] %>%
  .[!celltype.mapped%in%c("ExE_ectoderm","ExE_endoderm","Parietal_endoderm","Visceral_endoderm")]# %>%
  # .[,celltype.mapped:=stringr::str_replace_all(celltype.mapped,opts$to.merge)]

# Calculate WT proportions
wt.dt <- sample_metadata[class%in%opts$wt.classes] %>%
  .[,ncells:=.N, by="class"] %>%
  # .[,.(proportion=.N/unique(ncells)),by=c("celltype.mapped","class")] %>%
  .[,.(proportion=.N/unique(ncells), N=.N),by=c("celltype.mapped","class")] %>%
  .[,.(proportion=mean(proportion), N=mean(N)),by="celltype.mapped"]
  # .[,c("class","class"):="WT"]


#############################################
## Plot number of cells for each cell type ##
#############################################

ko.dt <- sample_metadata[!class%in%opts$wt.classes] %>%
  .[,ncells:=.N, by="class"] %>%
  .[,.(proportion=.N/unique(ncells), N=.N),by=c("celltype.mapped","class","class")]
  

# Merge
dt <- merge(ko.dt, wt.dt, by="celltype.mapped", allow.cartesian=T, suffixes = c(".ko",".wt"))

# to.plot <- dt[,.(diff_proportion=log2(proportion.ko/proportion.wt), diff_N=abs(N.ko-N.wt)), by=c("class","class","celltype.mapped")]# %>%
to.plot <- dt[,.(diff_proportion=log2(proportion.ko/proportion.wt), diff_N=N.ko-N.wt), by=c("class","class","celltype.mapped")]# %>%
  # setnames(c("class","class","celltype.mapped","diff"))

to.plot[,diff_N_norm:=linMap(abs(diff_N), from=0.15, to=2.5)]

# sequence_length = length(unique(to.plot.test$celltype.mapped))
# first_sequence = c(1:(sequence_length%/%2)) 
# second_sequence = c((sequence_length%/%2+1):sequence_length) 
# first_angles = c(90 - 180/length(first_sequence) * first_sequence)
# second_angles = c(-90 - 180/length(second_sequence) * second_sequence)

to.plot.wt_line <- data.table(
  celltype.mapped = unique(dt$celltype.mapped),
  diff_proportion = log2(1),
  diff_N = 0 #+  0.01
)


ggplot(to.plot, aes(x=factor(celltype.mapped), y=diff_proportion, group=1)) +
  geom_point(aes(color = celltype.mapped, size = diff_N_norm), stat = 'identity') + 
  geom_polygon(color="black", fill=NA, alpha=0.5, linetype="dashed", data=to.plot.wt_line) +
  facet_wrap(~class, nrow=2) +
  scale_color_manual(values=opts$celltype.colors, drop=F) +
  coord_polar() +
  guides(colour = guide_legend(override.aes = list(size=2), ncol=1)) +
  scale_size(guide = 'none') +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.title=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.line=element_blank(),
    axis.text.x = element_blank()
    # axis.text.x = element_text(angle= -76 - 360 / length(unique(to.plot.test$celltype.mapped)) * seq_along(to.plot.test$celltype.mapped))
  )

# pdf(paste0(io$outdir,"/pdf/mapping_stats_day2.pdf"), width=opts$width, height=opts$height)
# print(p)
# dev.off()






########################################
## Plot difference in number of cells ##
########################################

# ggplot(to.plot, aes(x=factor(celltype.mapped), y=log2(diff_N+0.01), group=1)) +
ggplot(to.plot, aes(x=factor(celltype.mapped), y=diff_N, group=1)) +
  geom_point(aes(color = celltype.mapped, size = diff_N_norm), stat = 'identity') + 
  geom_polygon(color="black", fill=NA, alpha=0.5, linetype="dashed", data=to.plot.wt_line) +
  facet_wrap(~class, nrow=2) +
  scale_color_manual(values=opts$celltype.colors, drop=F) +
  coord_polar() +
  guides(colour = guide_legend(override.aes = list(size=2), ncol=1)) +
  scale_size(guide = 'none') +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.title=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.line=element_blank(),
    axis.text.x = element_blank(),
    # axis.text.x = element_text(angle= -76 - 360 / length(unique(to.plot.test$celltype.mapped)) * seq_along(to.plot.test$celltype.mapped))
  )
