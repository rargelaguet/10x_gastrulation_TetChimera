here::i_am("processing/2_QC.R")

source(here::here("settings.R"))

#####################
## Define arguments ##
#####################

p <- ArgumentParser(description='')
p$add_argument('--metadata',       type="character",                    help='Metadata')
p$add_argument('--outdir',       type="character",                    help='Output directory')
p$add_argument('--min_nFeature_RNA',       type="integer",                    help='Minimum number of expressed genes')
p$add_argument('--max_nFeature_RNA',       type="integer",                    help='Maximum number of expressed genes')
p$add_argument('--mit_percent_RNA',       type="integer",                    help='Maximum percentage of mit reads')
p$add_argument('--rib_percent_RNA',       type="integer",                    help='Maximum percentage of rib reads')
# p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args <- list()
# # args$samples <- opts$samples
# args$metadata <- paste0(io$basedir,"/processed_all/metadata.txt.gz")
# args$min_nFeature_RNA <- 1500
# args$max_nFeature_RNA <- 10000
# args$mit_percent_RNA <- 30
# args$rib_percent_RNA <- 35
# args$outdir <- paste0(io$basedir,"/results_all/qc")
## END TEST ##

# Sanity checks
stopifnot(args$samples%in%opts$samples)

###############
## Load data ##
###############

metadata <- fread(args$metadata) %>% 
    # .[sample%in%args$samples] %>%
    .[,pass_rnaQC:=nFeature_RNA<=args$max_nFeature_RNA & nFeature_RNA>=args$min_nFeature_RNA & mit_percent_RNA<args$mit_percent_RNA & rib_percent_RNA<args$rib_percent_RNA]

table(metadata$sample)
table(metadata$stage)
table(metadata$pass_rnaQC)

#####################################
## Plot QC metrics before QC calls ##
#####################################

to.plot <- metadata %>%
    # .[,log_nFeature_RNA:=log10(nFeature_RNA)] %>%
    melt(id.vars=c("sample","cell","stage"), measure.vars=c("nFeature_RNA","mit_percent_RNA","rib_percent_RNA"))

facet.labels <- c("nFeature_RNA" = "Num. of genes", "mit_percent_RNA" = "Mitochondrial %", "rib_percent_RNA" = "Ribosomal %")

tmp <- data.table(
    variable = c("nFeature_RNA", "mit_percent_RNA", "rib_percent_RNA"),
    value = c(args$min_nFeature_RNA, args$mit_percent_RNA, args$rib_percent_RNA)
)

p <- gghistogram(to.plot, x="value", fill="sample", bins=50) +
    geom_vline(aes(xintercept=value), linetype="dashed", data=tmp) +
    facet_wrap(~variable, scales="free", nrow=1) +
    theme(
        axis.text =  element_text(size=rel(0.8)),
        axis.title.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=rel(0.75))
    )
    
pdf(file.path(args$outdir,"qc_metrics_histogram.pdf"), width=13, height=6)
print(p)
dev.off()

#####################################
## Plot QC metrics after QC calls ##
#####################################

to.plot <- metadata %>% .[pass_rnaQC==TRUE] %>%
    melt(id.vars=c("sample","cell","stage"), measure.vars=c("nFeature_RNA","mit_percent_RNA","rib_percent_RNA"))

p <- ggplot(to.plot, aes_string(x="sample", y="value")) +
    geom_boxplot(fill="gray70", outlier.shape=NA, coef=1) +
    facet_wrap(~variable, scales="free_y", nrow=1, labeller = as_labeller(facet.labels)) +
    # scale_fill_manual(values=opts$stage.colors) +
    theme_classic() +
    theme(
        axis.text.y = element_text(colour="black",size=rel(1)),
        axis.text.x = element_text(colour="black",size=rel(0.55), angle=20, hjust=1, vjust=1),
        axis.title.x = element_blank()
    )

pdf(file.path(args$outdir,"qc_metrics_boxplot.pdf"), width=12, height=5)
# pdf(sprintf("%s/qc_metrics_boxplot.pdf",args$outdir))
print(p)
dev.off()


#########################################################
## Plot fraction of cells that pass QC for each sample ##
#########################################################

to.plot <- metadata %>%
    .[,mean(pass_rnaQC,na.rm=T),by=c("sample","stage")]

p <- ggbarplot(to.plot, x="sample", y="V1", fill="gray70") +
    # scale_fill_manual(values=opts$stage.colors) +
    labs(x="", y="Fraction of cells that pass QC (RNA)") +
    # facet_wrap(~stage)
    theme(
        legend.position = "none",
        axis.text.y = element_text(colour="black",size=rel(0.8)),
        axis.text.x = element_text(colour="black",size=rel(0.50), angle=20, hjust=1, vjust=1),
    )

pdf(file.path(args$outdir,"qc_metrics_barplot.pdf"), width=8, height=6)
print(p)
dev.off()

##########
## Save ##
##########

fwrite(metadata, file.path(args$outdir,"sample_metadata_after_qc.txt.gz"), quote=F, na="NA", sep="\t")

