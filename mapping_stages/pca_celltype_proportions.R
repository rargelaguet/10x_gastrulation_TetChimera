here::i_am("mapping_stages/pca_celltype_proportions.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$cellype_proportions.atlas <- file.path(io$atlas.basedir,"results/celltype_proportions/celltype_proportions_noExE.txt.gz")
io$outdir <- file.path(io$basedir,"results_new/mapping_stages"); dir.create(io$outdir, showWarnings = F)

# Options
opts$remove_ExE_cells <- TRUE
opts$stages <- c("E7.5","E8.5")

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  setnames("celltype.mapped","celltype") %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & !is.na(celltype) & stage%in%opts$stages]

if (opts$remove_ExE_cells) {
  sample_metadata <- sample_metadata %>%
    .[!celltype%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

################
## Load query ##
################

dt.query <- sample_metadata %>% copy %>%
  .[,N:=.N,by="sample"] %>%
  .[,.(celltype_proportion=.N/unique(N)),by=c("sample","celltype")]

sample_metadata.query <- sample_metadata %>%
  .[,c("sample","stage","class")] %>% unique %>%
  .[,dataset:="Query"]

################
## Load atlas ##
################

sample_metadata.atlas <- fread(io$atlas.metadata) %>% 
  .[stage!="mixed_gastrulation",c("sample","stage")] %>% unique %>%
  .[,c("class","dataset"):="Atlas"]

# Load cell type proportions in the atlas
dt.atlas <- fread(io$cellype_proportions.atlas) %>%
  .[sample%in%sample_metadata.atlas$sample,c("sample","celltype","celltype_proportion")]

#################
## Concatenate ##
#################

sample_metadata <- rbind(sample_metadata.query, sample_metadata.atlas)
dt <- rbind(dt.query,dt.atlas)

##############################
## Dimensionality reduction ##
##############################

matrix <- dt %>% 
  dcast(sample~celltype, fill=0, value.var="celltype_proportion") %>% 
  matrix.please

# PCA
pca <- prcomp(matrix, rank.=5)
pca.var.explained <- 100*(pca$sdev**2 / sum(pca$sdev**2)) %>% round(4)

##########
## Plot ##
##########

to.plot <- pca$x %>% as.data.table %>% 
  .[,sample:=rownames(pca$x)] %>%
  merge(sample_metadata, by="sample")

# Define colors (stage)
stages <- unique(to.plot$stage) %>% sort
stage.colors <- viridis::viridis(n=length(stages))
names(stage.colors) <- rev(stages)
  
# Define shapes (conditions)
classes <- unique(to.plot$class)
shapes <- 21:(21+length(classes)-1)
names(shapes) <- classes

p <- ggplot(to.plot, aes(x=PC1, y=PC2, fill=stage, shape=class)) +
  geom_point(stroke=0.5, color="black", size=5) +
  scale_fill_manual(values=stage.colors) +
  scale_shape_manual(values=shapes) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  guides(shape=guide_legend(override.aes=list(fill="black"))) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text = element_text(color="black", size=rel(0.8))
  )

pdf(file.path(io$outdir,"pca_mapping_stages.pdf"), width=7, height=5)
print(p)
dev.off()
