suppressMessages(library(RColorBrewer))
suppressMessages(library(MOFA2))

#####################
## Define settings ##
#####################

source("/Users/ricard/10x_gastrulation_TetChimera/settings.R")

################
## Load model ##
################

# file <- paste0(io$basedir,"/mofa/hdf5/model.hdf5")
# model <- load_model(file)

#########################
## Add sample metadata ##
#########################

cells <- as.character(unname(unlist(MOFA2::samples_names(model))))

sample_metadata.mofa <- copy(sample_metadata) %>%
  setnames("cell","sample") %>%
  setnames("batch","group") %>%
  .[sample%in%cells] %>% setkey(sample) %>% .[cells]
stopifnot(all(cells==sample_metadata.mofa$cell))

samples_metadata(model) <- sample_metadata.mofa

####################
## Subset factors ##
####################

# r2 <- model@cache$variance_explained$r2_per_factor
# factors <- sapply(r2, function(x) x[,"RNA"]>0.01)
# model <- subset_factors(model, which(apply(factors,1,sum)>=1))
# factors(model) <- paste("Factor",1:get_dimensions(model)[["K"]], sep=" ")

#############################
## Plot variance explained ##
#############################

plot_variance_explained(model, x="group", y="factor")

##################
## Plot factors ##
##################

plot_factor(model, factors = 1, groups="all", color_by = "celltype", group_by = "celltype", add_violin = T, dodge=T) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_factor(model, factors = 5, color_by = "BioClassification", group_by = "BioClassification", add_violin = T, dodge=T) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_factors(model, factors = c(1,2), color_by = "stage")



######################
## Sumarise factors ##
######################

# levels_df <- model@samples_metadata[,c("sample","BioClassification")] %>% setnames("BioClassification","leve2l")
# summarise_factors(model, levels_df, factors = 5, groups = "all", abs = F, return_data = F)

##########
## UMAP ##
##########

# model <- run_umap(model)
# plot_dimred(model, method="UMAP", color_by = "BioClassification")
# plot_dimred(model, method="UMAP", color_by = "TUBA1B")

##################
## Plot weights ##
##################

plot_weights(model, factor = 3, nfeatures = 20, scale = F)

###############
## Plot data ##
###############

plot_factor(model, factors = 5, color_by = "LYZ")

plot_factor(model, factors = 5, color_by = "IGLL1", group_by="BioClassification", dodge=TRUE) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_data_heatmap(model,
                  factor = 5,
                  features = 25,
                  view = 1,
                  denoise = FALSE,
                  legend = TRUE,
                  # min.value = 0, max.value = 6,
                  cluster_rows = T, cluster_cols = F,
                  show_colnames = F, show_rownames = T,
                  scale="row"
                  # annotation_samples = "Category",  annotation_colors = list("Category"=opts$colors), annotation_legend = F
)

plot_data_scatter(model, factor=1, view=1, color_by = "lab", features = 8, dot_size = 2)


