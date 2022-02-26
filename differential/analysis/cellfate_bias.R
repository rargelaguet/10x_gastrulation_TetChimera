here::i_am("differential/analysis/cellfate_bias.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

# Load utils
source(here::here("differential/analysis/utils.R"))

##############
## Settings ##
##############

# I/O
io$indir <- file.path(io$basedir,"results_all/differential")
io$outdir <- file.path(io$basedir,"results_all/differential/pdf/cellfate_bias"); dir.create(io$outdir, showWarnings = F)

# Options
opts$min.cells <- 50

opts$remove_ExE <- TRUE

opts$rename_celltypes <- c(
  "Erythroid3" = "Erythroid",
  "Erythroid2" = "Erythroid",
  "Erythroid1" = "Erythroid",
  "Blood_progenitors_1" = "Blood_progenitors",
  "Blood_progenitors_2" = "Blood_progenitors"
)

###############
## Load data ##
###############

source(here::here("differential/analysis/load_data.R"))

if (opts$remove_ExE) {
  diff.dt <- diff.dt %>% .[!celltype%in%c("ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

stopifnot(opts$rename_celltypes%in%unique(diff.dt$celltype))

#######################
## Load marker genes ##
#######################

opts$min.marker.score <- 0.85

marker_genes.dt <- fread(io$atlas.marker_genes) %>%
  .[score>=opts$min.marker.score]

if (opts$remove_ExE) {
  marker_genes.dt <- marker_genes.dt %>% .[!celltype%in%c("ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
  
  stopifnot(opts$rename_celltypes%in%unique(diff.dt$celltype))
}

# Rename celltypes
marker_genes.dt <- marker_genes.dt %>% 
  .[,celltype:=stringr::str_replace_all(celltype,opts$rename_celltypes)] %>%
  .[,.(score=mean(score)), by=c("celltype","gene")]

####################
## Filter results ##
####################

# Filter
diff.dt <- diff.dt[groupA_N>=opts$min.cells & groupB_N>=opts$min.cells]

# Remove sex-specific differences
# dt.filt <- dt[gene!="Xist"]

# Remove Riks
# dt.filt <- dt.filt[!grep("Rik",gene)]

# Remove hemoglobins
diff.dt <- diff.dt[!grep("^Hb",gene)]

# Subset to lineage markers
diff_markers.dt <- diff.dt[gene%in%unique(marker_genes.dt$gene)]

# Filter out genes that are DE across all celltypes
diff_markers.dt <- diff_markers.dt[,N:=sum(sig),by=c("gene","class")] %>% .[N<=10] %>% .[,N:=NULL]
# diff_markers.dt[class=="E8.5_Dnmt1KO",sum(sig),by=c("gene","class")] %>% View
# diff_markers.dt[class=="E8.5_Dnmt1KO" & gene=="Hoxc8"]

######################################################
## Quantify cell fate disruption using the DE genes ##
######################################################

to.plot <- diff_markers.dt %>%
  merge(
    marker_genes.dt[,c("gene","celltype")] %>% setnames("celltype","celltype_marker"), by = "gene", allow.cartesian=TRUE
  ) %>% .[,sum(sig), by=c("celltype","celltype_marker","class","sign")]

to.plot[!celltype%in%c("Rostral_neurectoderm","Pharyngeal_mesoderm","")]
to.plot <- to.plot[sign=="Downregulated in TET_TKO"]
to.plot <- to.plot[V1>=5]

p <- ggplot(to.plot, aes(x=celltype, y=V1)) +
  geom_bar(aes(fill = celltype_marker), color="black", stat="identity") + 
  facet_wrap(~sign, scales="fixed") +
  scale_fill_manual(values=opts$celltype.colors[names(opts$celltype.colors)%in%unique(to.plot$celltype_marker)], drop=F) +
  labs(x="", y="Number of DE genes") +
  guides(x = guide_axis(angle = 90), fill=guide_legend(ncol=1)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    # axis.line = element_blank(),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black", size=rel(0.85))
  )

pdf(sprintf("%s/DE_barplots_marker_genes_fate_bias_legend.pdf",io$outdir), width=7, height=6)
print(p)
dev.off()

###################################################################
## Quantify cell fate disruption using the DE genes: polar plots ##
###################################################################

for (i in opts$ko.classes) {

  to.plot <- diff_markers.dt[class==i] %>%
    merge(
      marker_genes.dt[,c("gene","celltype")] %>% setnames("celltype","celltype_marker"), by = "gene", allow.cartesian=TRUE
    ) %>% .[,sum(sig), by=c("celltype","celltype_marker","class")]
  
  to.plot[V1>50,V1:=50]
  
  p <- ggplot(to.plot, aes(x=celltype_marker, y=V1)) +
    geom_bar(aes(fill = celltype_marker), color="black", stat = 'identity') + 
    # geom_polygon(color="black", fill=NA, alpha=0.5, linetype="dashed", data=foo) +
    facet_wrap(~celltype) +
    scale_fill_manual(values=opts$celltype.colors, drop=F) +
    coord_polar() + ylim(c(0,50)) +
    # guides(colour = guide_legend(override.aes = list(size=2), ncol=1)) +
    theme_bw() +
    theme(
      legend.position = "none",
      legend.text = element_text(size=rel(0.75)),
      legend.title = element_blank(),
      axis.title=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.line=element_blank(),
      axis.text.x = element_blank()
      # axis.text.x = element_text(angle= -76 - 360 / length(unique(to.plot$celltype)) * seq_along(to.plot$celltype))
    )
  
  pdf(file.path(io$outdir,sprintf("%s_DE_polar_plot_marker_genes_fate_bias.pdf",i)), width=11, height=7)
  print(p)
  dev.off()
  
}

