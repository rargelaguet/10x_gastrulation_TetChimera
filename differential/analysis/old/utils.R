
gg_volcano_plot <- function(to.plot, top_genes=10, xlim=NULL, ylim=NULL) {
  negative_hits <- to.plot[sig==TRUE & logFC<0,gene]
  positive_hits <- to.plot[sig==TRUE & logFC>0,gene]
  all <- nrow(to.plot)
  
  if (is.null(xlim))
    xlim <- max(abs(to.plot$logFC), na.rm=T)
  if (is.null(ylim))
    ylim <- max(-log10(to.plot$p.value+1e-100), na.rm=T)
  
  p <- ggplot(to.plot, aes(x=logFC, y=-log10(p.value+1e-100))) +
    labs(x="Log fold change", y=expression(paste("-log"[10],"(p.value)"))) +
    # geom_hline(yintercept = -log10(opts$threshold_fdr), color="blue") +
    geom_segment(aes(x=0, xend=0, y=0, yend=ylim-1), color="orange", size=0.5) +
    ggrastr::geom_point_rast(aes(color=sig, size=sig)) +
    scale_color_manual(values=c("black","red")) +
    scale_size_manual(values=c(0.75,1.25)) +
    scale_x_continuous(limits=c(-xlim-1.5,xlim+1.5)) +
    scale_y_continuous(limits=c(0,ylim+3)) +
    # annotate("text", x=0, y=ylim+3, size=4, label=sprintf("(%d)", all)) +
    # annotate("text", x=-xlim-0.5, y=ylim+3, size=4, label=sprintf("%d (-)",length(negative_hits))) +
    # annotate("text", x=xlim+0.5, y=ylim+3, size=4, label=sprintf("%d (+)",length(positive_hits))) +
    annotate("text", x=-xlim-0.5, y=0, size=3, label=sprintf("Up in Host (N=%s)",to.plot[1,paste0("N_",opts$groupA),with=F][[1]])) +
    annotate("text", x=xlim+0.5, y=0, size=3, label=sprintf("Up in TET TKO (N=%s)",to.plot[1,paste0("N_",opts$groupB),with=F][[1]])) +
    ggrepel::geom_text_repel(data=head(to.plot[sig==T],n=top_genes), aes(x=logFC, y=-log10(p.value+1e-100), label=gene), size=4) +
    theme_classic() +
    theme(
      axis.text = element_text(size=rel(0.75), color='black'),
      axis.title = element_text(size=rel(1.0), color='black'),
      legend.position="none"
    )
  return(p)
}