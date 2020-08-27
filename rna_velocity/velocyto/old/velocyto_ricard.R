library(velocyto.R)
library(scater)
sce <- readRDS("/Users/ricard/data/gastrulation/rna/velocyto/sce/exon_spanning_sce.rds")

idx <- colData(sce)$stage%in%c("E5.5","E6.5","E7.5") & colData(sce)$lineage%in%c("EPI","Ectoderm","Mesoderm","Endoderm") & colData(sce)$KO_3b=="not"

sce <- sce[,idx]
smat <- as.matrix(readRDS("/Users/ricard/data/gastrulation/rna/velocyto/sce/exon_spanning_sce.rds")@assays$data[[1]])[,idx]
emat <- as.matrix(readRDS("/Users/ricard/data/gastrulation/rna/velocyto/sce/exons_sce.rds")@assays$data[[1]])[,idx]
nmat <- as.matrix(readRDS("/Users/ricard/data/gastrulation/rna/velocyto/sce/introns_sce.rds")@assays$data[[1]])[,idx]

sce$stage_lineage <- as.factor(paste(sce$stage,sce$lineage,sep="_"))

cell.colors <-  sce$stage_lineage
names(cell.colors) <- colnames(sce)

# filter expression matrices based on some minimum max-cluster averages
emat <- filter.genes.by.cluster.expression(emat, cell.colors, min.max.cluster.average = 5)
nmat <- filter.genes.by.cluster.expression(nmat, cell.colors, min.max.cluster.average = 1)
smat <- filter.genes.by.cluster.expression(smat, cell.colors, min.max.cluster.average = 0.5)

# look at the resulting gene set
str(intersect(intersect(rownames(emat),rownames(nmat)),rownames(smat)))

rvel.qf <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 5,fit.quantile = 0.02)
# rvel <- gene.relative.velocity.estimates(emat, nmat, smat=smat, deltaT=1, kCells = 5, min.nmat.emat.slope = 0.1, min.nmat.smat.correlation = 0.1)

palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#FFFFFF")
legend <- as.character(unique(cell.colors))

cell.colors.palette <- as.character(cell.colors)
for (i in seq_along(legend)){
    cell.colors.palette[cell.colors.palette==legend[i]] <- palette[i]
}
names(cell.colors.palette) <- names(cell.colors)


io <- list()
io$outdir <- "/Users/ricard/gastrulation/rna/velocyto/out"
pdf(file=paste0(io$outdir,"/velocyto2.pdf"), useDingbats=F)
pca.velocity.plot(rvel.qf, 
  nPcs=2, plot.cols=1, cell.colors=ac(cell.colors.palette,alpha=1),cex=1.2, pcount=0.1
)
legend("topleft", legend=legend, fill=palette[seq_along(legend)], cex=0.75)
dev.off()