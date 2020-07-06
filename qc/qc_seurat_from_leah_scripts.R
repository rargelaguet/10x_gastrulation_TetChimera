library(Seurat)
library(data.table)
library(purrr)

io <- list()
io$outdir <- "/Users/ricard/data/10x_gastrulation_TetChimera/results/first_batch"
seurat <- readRDS("/Users/ricard/data/10x_gastrulation_TetChimera/processed/first_batch/seurat.rds")
metadata <- fread("/Users/ricard/data/10x_gastrulation_TetChimera/processed/first_batch/sample_metadata.txt.gz")

batches <- unique(metadata$batch)

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "mt-")
# ribo.genes <- c(grep(pattern = "^Rpl", x = rownames(seurat), value = TRUE),grep(pattern = "^Rps", x = rownames(seurat), value = TRUE))
# seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, features = ribo.genes)

to.plot <- seurat@meta.data %>% as.data.table %>%
  melt(id.vars=c("cell","batch"), measure.vars=c("nCount_RNA","nFeature_RNA"))


p <- ggplot(to.plot, aes(x=value)) +
  facet_wrap(~batch+variable, nrow=length(batches), scales="free_y") +
  geom_histogram(position = 'identity', bins=250) +
  labs(x="") +
  # geom_vline(xintercept=io$min_nFeature_RNA[[opts$experiment]][[b]]) +
  # ggtitle(paste0(b, ", min=", io$min_nFeature_RNA[[opts$experiment]][[b]])) +
  theme_classic() +
  theme(
    axis.text = element_text(size = rel(0.75), color="black")
  )
  
pdf(paste0(io$outdir,"/qc_firstbatch.pdf"), width=10, height=8, useDingbats = F)
print(p)
dev.off()