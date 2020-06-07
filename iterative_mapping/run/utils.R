
recursive.fn <- function(sce.query, sce.atlas, dist) {
  
  celltypes.to.loop <- sce.query$celltype_mapped %>% unique %>% .[grep("%",.)]
  
  dt_list <- list()
  for (i in celltypes.to.loop) {
    ids <- which(sce.query$celltype_mapped==i)
    foo <- stringr::str_split(i,"%") %>% unlist
    
    # Subset
    dist.sub <- usedist::dist_subset(dist, foo)
    h.sub <- hclust(dist.sub, method="ward.D")
    sce.query.sub <- sce.query[,ids]
    sce.atlas.sub <- sce.atlas[,sce.atlas$celltype %in% foo]
    
    # Run cell type assignment prediction
    dt_list[[i]] <- mapping.fn(sce.query.sub, sce.atlas.sub, h.sub)
  }
  dt <- rbindlist(dt_list)
  return(dt)
    
}

mapping.fn <- function(sce.query, sce.atlas, h) {
  
  # Cut the tree into two groups
  cut <- cutree(h, k=2)
  
  # Define groups
  groupA <- names(which(cut==1))
  groupB <- names(which(cut==2))
  groups.parsed <- c(paste(groupB,collapse="%"),paste(groupA,collapse="%"))
  sce.atlas$celltype <- groups.parsed[as.numeric(sce.atlas$celltype%in%groupA)+1]
  # sce.atlas$celltype <- factor(sce.atlas$celltype,levels=groups.parsed)
  
  # Run MNN
  dt <- mnn.fn(sce.query, sce.atlas)
  
  return(dt)
    
}


mnn.fn <- function(sce.query, sce.atlas, npcs = 30, k = 25) {
  block <- c(rep("atlas",ncol(sce.atlas)),rep("query",ncol(sce.query))) %>% as.factor
  # block <- c(rep("atlas",ncol(sce.atlas)),sce.query$z)) %>% as.factor
  # block <- c(sce.atlas$sample,sce.query$z) %>% as.factor
  
  # Concatenate
  sce_all <- SingleCellExperiment(
    list(counts=Matrix::Matrix(cbind(counts(sce.atlas),counts(sce.query)),sparse=TRUE))
  )
  
  # Normalise
  # sce_all <- logNormCounts(sce_all)
  sce_all <- multiBatchNorm(sce_all, batch=as.factor(block))
  # sce_all <- multiBatchNorm(sce_all, batch=c(atlas_meta$sample, map_meta$batch))
  
  # Highly variable genes
  # hvgs <- getHVGs(sce_all, block=block)
  hvgs <- rownames(sce_all)
  
  # PCA
  # pca_all <- irlba::prcomp_irlba(t(assay(sce_all,"logcounts")), n = npcs)$x
  # rownames(pca_all) <- colnames(sce_all)
  
  pca_all <- multiBatchPCA(sce_all,
    batch = block,
    subset.row = hvgs,
    d = npcs,
    preserve.single = TRUE,
    assay.type = "logcounts"
  )[[1]]
  rownames(pca_all) <- colnames(sce_all)
  atlas_pca <- pca_all[1:ncol(sce.atlas),]
  query_pca <- pca_all[-(1:ncol(sce.atlas)),]

  # (IS THIS REQUIRED???) Batch effect correction for the atlas
  # order_df        <- atlas_meta[!duplicated(atlas_meta$sample), c("stage", "sample")]
  # order_df$ncells <- sapply(order_df$sample, function(x) sum(atlas_meta$sample == x))
  # order_df$stage  <- factor(order_df$stage, 
  #                           levels = rev(c("E8.5","E8.25","E8.0","E7.75","E7.5","E7.25","mixed_gastrulation","E7.0","E6.75","E6.5")))
  # order_df       <- order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
  # order_df$stage <- as.character(order_df$stage)
  # 
  # set.seed(42)
  # atlas_corrected <- doBatchCorrect(counts         = logcounts(atlas_sce[hvgs,]), 
  #                                   timepoints      = atlas_meta$stage, 
  #                                   samples         = atlas_meta$sample, 
  #                                   timepoint_order = order_df$stage, 
  #                                   sample_order    = order_df$sample, 
  #                                   pc_override     = atlas_pca,
  #                                   npc             = npcs)
  
  # MNN mapping
  # correct <- reducedMNN(rbind(atlas_pca, query_pca),
  correct <- reducedMNN(pca_all, batch = block)[["corrected"]]
  correct_atlas <- correct[1:nrow(atlas_pca),,drop=F]
  correct_query   <- correct[-(1:nrow(atlas_pca)),,drop=F]

  # get metadata
  meta_query <- as.data.frame(colData(sce.query)) %>% tibble::rownames_to_column("cell") %>% .[,c("cell"),drop=F]
  meta_atlas <- as.data.frame(colData(sce.atlas)) %>% tibble::rownames_to_column("cell") %>% .[,c("cell","celltype")]
  mapping <- get_meta(
    correct_atlas = correct_atlas,
    atlas_meta = meta_atlas,
    correct_query = correct_query,
    query_meta = meta_query,
    k_map = k
  )
  
  # Compute mapping scores
  out <- list()
  for (i in seq(from = 1, to = k)) {
    out$closest.cells[[i]]     <- sapply(mapping, function(x) x$cells.mapped[i])
    out$celltypes.mapped[[i]]  <- sapply(mapping, function(x) x$celltypes.mapped[i])
  }  
  multinomial.prob <- getMappingScore(out)
  

  #  Prepare output
  ct <- sapply(mapping, function(x) x$celltype.mapped); is.na(ct) <- lengths(ct) == 0
  cm <- sapply(mapping, function(x) x$cells.mapped[1]); is.na(cm) <- lengths(cm) == 0
  mapping.dt <- data.table(
    cell = names(mapping), 
    celltype_mapped = unlist(ct),
    mapping_score = multinomial.prob
  )  
  
  return(mapping.dt)
}


get_meta <- function(correct_atlas, atlas_meta, correct_query, query_meta, k_map = 10){
  knns <- BiocNeighbors::queryKNN(correct_atlas, correct_query, k = k_map, get.index = TRUE, get.distance = FALSE)
  
  #get closest k matching cells
  k.mapped  <- t(apply(knns$index, 1, function(x) atlas_meta$cell[x]))
  
  # get celltypes
  celltypes <- t(apply(k.mapped, 1, function(x) atlas_meta$celltype[match(x, atlas_meta$cell)]))
  celltype.mapped <- apply(celltypes, 1, function(x) getmode(x, 1:length(x)))
  
  out <- lapply(1:length(celltype.mapped), function(x){
    list(cells.mapped     = k.mapped[x,],
         celltype.mapped  = celltype.mapped[x],
         celltypes.mapped = celltypes[x,]
        )
  })
  names(out) <- query_meta$cell
  return(out)
}


getmode <- function(v, dist) {
  tab <- table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied <- names(tab)[tab == max(tab)]
    sub  <- dist[v %in% tied]
    names(sub) <- v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}

getHVGs <- function(sce, block, min.mean = 1e-3, p.value=0.01){
  decomp <- modelGeneVar(sce, block=block)
  decomp <- decomp[decomp$mean > min.mean,]
  decomp$FDR <- p.adjust(decomp$p.value, method = "fdr")
  return(rownames(decomp)[decomp$p.value < p.value])
}


getAreaFactorsDirectly = function(sce, transform = NULL) {
  # sce is a SingleCellExperiment object
  # transform is a function
  meta = colData(sce)
  counts = assay(sce, "counts")
  if (!is.null(transform)) {
    sizeFactors.area <- transform(meta$Area)/mean(transform(meta$Area))
  } else {
    sizeFactors.area <- meta$Area/mean(meta$Area)
  }
  logcounts.area <- log2(t(t(counts)/sizeFactors.area) + 1)
  return(list(sizeFactors = sizeFactors.area,
              logcounts.area = logcounts.area))
}


getMappingScore <- function(mapping){
  celltypes_accrossK <- matrix(unlist(mapping$celltypes.mapped),
                               nrow=length(mapping$celltypes.mapped[[1]]),
                               ncol=length(mapping$celltypes.mapped))
  celltype.score <- c()
  for (i in 1:nrow(celltypes_accrossK)){
    p <- max(table(celltypes_accrossK[i,]))
    index <- which(table(celltypes_accrossK[i,]) == p)
    p <- p/length(mapping$celltypes.mapped)
    celltype.score <- c(celltype.score,p)
  }
  return(celltype.score)  
}
