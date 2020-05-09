library(SingleCellExperiment)
library(scran)
library(scater)

sce <- readRDS("/Users/ricard/data/10x_gastrulation_TetChimera/processed/sce.rds")

sf <- computeSumFactors(sce)

# sce <- calculateQCMetrics(sce)

# foo <- sce@assays$data$logcounts

# sce_filt <- sce

# # Temporarily remove the lowly expressed genes
# sce_filt <- sce_filt[!(rowMeans(counts(sce)) <= 1 | rowData(sce)$pct_dropout_by_counts > 90),]

# # Compute size factors without the lowly expressed genes
# sf = computeSumFactors(sce_filt, positive=TRUE, sf.out=T)

# # qplot(sf, sce_filt$total_counts, log="xy", ylab="Library size (mapped reads)", xlab="Size factor")
# ggplot(data.frame(sf=log(sf), counts=log(sce_filt$total_counts))) +
#   geom_point(aes(x=sf,y=counts)) +
#   labs(y="Library size (log)", x="Size factor (log)") +
#   theme_bw() +
#   theme(
#     axis.title = element_text(colour="black", size=15),
#     axis.text = element_text(colour="black", size=12)
#   )

# # Normalise and log transform with the lowly expressed genes
# sizeFactors(sce) <- sf; sce$sizeFactor <- sf
# sizeFactors(sce_filt) <- sf; sce_filt$sizeFactor <- sf
# sce <- normalize(sce, exprs_values="counts")
# sce_filt <- normalize(sce_filt, exprs_values="counts")

# # Update quality metrics
# sce = calculateQCMetrics(sce)

# foo <- data.frame(sd=apply(exprs(sce),1,sd), mean=apply(exprs(sce),1,mean))
# ggplot(foo, aes(x=mean, y=sd)) +
#   geom_point() +
#   stat_smooth() +
#   scale_color_manual(values=c("black","red")) +
#   xlab('Mean') + ylab('Standard deviation')

# # saveRDS(sce,io$out.file)

# #############################################
# ## Correlate Seurat vs scran normalisation ##
# #############################################

# bar <- sce@assays$data$logcounts

# foo[1:10,1:10]
# bar[1:10,1:10]

# saveRDS(sf, "/Users/ricard/data/10x_gastrulation_TetChimera/processed/size_factors.rds")










##########

# a <- readRDS("/Users/ricard/data/10x_gastrulation_TetChimera/processed/size_factors.rds")
# a@colData

df <- data.frame(
  library_size = colSums(counts(a)),
  size_factor = sizeFactors(a),
  cell = colnames(a)
)
fwrite(df, "/Users/ricard/data/10x_gastrulation_TetChimera/processed/size_factors.tab", sep="\t", col.names = T, quote=F)

plot(df$library_size,df$size_factor)
