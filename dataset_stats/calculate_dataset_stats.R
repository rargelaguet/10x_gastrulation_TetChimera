here::i_am("processing/stats/calculate_dataset_stats.R")

source(here::here("settings.R"))

#####################
## Define settings ##
#####################

# I/O
# io$metadata <- file.path(io$basedir,"results/sex_assignment/sample_metadata_after_sex_assignment.txt.gz")
io$outdir <- file.path(io$basedir,"results/dataset_stats"); dir.create(io$outdir, showWarnings = F)

# Options
# opts$classes <- c("WT_tdTomato+", "WT_tdTomato-", "TET_TKO")
opts$stage_classes <- c("E7.5_WT", "E8.5_WT", "E7.5_TET_TKO","E8.5_TET_TKO")

###############
## Load data ##
###############

metadata <- fread(io$metadata) %>% 
    # .[,class:=sprintf("%s_%s",stage,class)] %>%
    .[,stage_class:=paste(stage,class2,sep="_")] %>%
    .[pass_rnaQC==TRUE & stage_class%in%opts$stage_classes] %>%
    .[,stage_class:=factor(stage_class, levels=opts$stage_classes)]

####################################
## Plot number of cells per class ##
####################################

to.plot <- metadata %>% .[,.N,by=c("stage_class")]

p <- ggbarplot(to.plot, x="stage_class", y="N", fill="stage_class", position=position_dodge(width = 0.75)) +
    labs(x="", y="Number of cells") +
    scale_fill_manual(values=opts$stage_class2_colors) +
    theme(
        legend.position = "none",
        legend.title = element_blank(),
        axis.text.y = element_text(colour="black",size=rel(0.8)),
        axis.text.x = element_text(colour="black",size=rel(0.80)),
    )

pdf(file.path(io$outdir,"ncells_per_class.pdf"), width=5, height=5.5)
print(p)
dev.off()

###################################################
## Plot number of samples per class and data set ##
###################################################

to.plot <- metadata %>% 
    .[,.(N=length(unique(alias))), by=c("stage_class")] %>%
    .[,N:=factor(N, levels=1:max(N))]

p <- ggbarplot(to.plot, x="stage_class", y="N", fill="stage_class", position=position_dodge(width = 0.75)) +
    labs(x="", y="Number of samples") +
    scale_fill_manual(values=opts$stage_class2_colors) +
    theme(
        legend.position = "none",
        legend.title = element_blank(),
        axis.text.y = element_text(colour="black",size=rel(0.8)),
        axis.text.x = element_text(colour="black",size=rel(0.80)),
    )

pdf(file.path(io$outdir,"nsamples_per_class.pdf"), width=4, height=5.5)
print(p)
dev.off()

