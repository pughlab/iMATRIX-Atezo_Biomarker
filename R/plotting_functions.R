# Function to align plots (from stackoverflow) 
align_plots1 <- function (...) {
  pl <- list(...)
  stopifnot(do.call(all, lapply(pl, inherits, "gg")))
  gl <- lapply(pl, ggplotGrob)
  bind2 <- function(x, y) gtable:::rbind_gtable(x, y, "first")
  combined <- Reduce(bind2, gl[-1], gl[[1]])
  wl <- lapply(gl, "[[", "widths")
  combined$widths <- do.call(grid::unit.pmax, wl)
  grid::grid.newpage()
  grid::grid.draw(combined)
}

gsea.fx <- function(resdf){
  #remove no genenames and dedup
  resdf <- resdf[!is.na(resdf$Gene_gencode),]
  resdf <- resdf[!duplicated(resdf$Gene_gencode),]
  ranks <- resdf$log2FoldChange
  names(ranks) <- resdf$Gene_gencode
  ranks <- ranks[!is.na(ranks)]
  ranks <- sort(ranks)
  rownames(Hs.H) <- Hs.H$V1
  Hs.H$V1 <- NULL
  Hs.H$V2 <- NULL
  Hs.H.list <- as.list(as.data.frame(t(Hs.H)))
  fgseaRes <- fgsea(Hs.H.list, ranks, minSize=5, maxSize = 500, eps = 0)
  p <- ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill= padj < 0.05)) + coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",title="Hallmark pathways NES from GSEA") + myplot + 
    theme(axis.text.x = element_text(size = 40, angle = 0, hjust = 0.5),
          axis.title.x = element_text(size = 40),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 25),
          plot.title = element_text(size = 40, hjust = 0.5)) + 
    labs(y = "Normalized Enrichment Score", title = "Differential pathway analysis (Hallmark gene sets)") 
  return(p)
}



volcano.fx <- function(resdf, fc, padj, ttl){
  
  # remove padj = NA
  
  resdf <- resdf[ !is.na(resdf$padj),]
  
  resdf$threshold <- NA
  resdf$threshold[ resdf$log2FoldChange > fc & resdf$padj < padj] <- "Up-regulated"
  resdf$threshold[ resdf$log2FoldChange < -fc & resdf$padj < padj] <- "Down-regulated"
  resdf$threshold[ is.na(resdf$threshold)] <- "not significant"
  
  res_upreg <- resdf[ resdf$threshold == "Up-regulated",]    
  res_upreg <- res_upreg[order(res_upreg$log2FoldChange, decreasing = T),]    
  res_downreg <- resdf[ resdf$threshold == "Down-regulated",]    
  res_downreg <- res_downreg[order(res_downreg$log2FoldChange, decreasing = F),]  
  
  if(nrow(res_upreg) < 10){
    resdf$genelabels[ rownames(resdf) %in% rownames(res_upreg)] <- "UP"}
  if(nrow(res_downreg) < 10){
    resdf$genelabels[rownames(resdf) %in% rownames(res_downreg)] <- "DOWN"}    
  
  if(nrow(res_upreg) >= 10){
    resdf$genelabels[rownames(resdf) %in% rownames(res_upreg)[1:10]] <- "UP"   }
  if(nrow(res_downreg) >= 10){
    resdf$genelabels[rownames(resdf) %in% rownames(res_downreg)[1:10]] <- "DOWN" } 
  
  p <- ggplot(resdf, aes(x=log2FoldChange, y=-log10(pvalue))) +
    geom_point(aes(color = threshold), size=2.5) +
    scale_colour_manual(values = c("Down-regulated"= "blue", "Up-regulated"="red",  "not significant"= "black")) +    
    geom_text_repel(data = subset(resdf, genelabels == "UP"),
                    label = subset(resdf, genelabels == "UP")$Gene_gencode,
                    size = 6, box.padding = 1, max.overlaps = Inf,  min.segment.length = 0,
                    direction = "both", nudge_x = 3, nudge_y = 0.5, vjust = 0.5, hjust = 0.5) + 
    geom_text_repel(data = subset(resdf, genelabels == "DOWN"),
                    label = subset(resdf, genelabels == "DOWN")$Gene_gencode,
                    size = 6, box.padding = 1, max.overlaps = Inf, direction = "both", 
                    nudge_x = -3, nudge_y = 0.5,
                    vjust = 0.5, hjust = 0.5, min.segment.length = 0) + 
    myplot + myaxis +
    theme(axis.text.x = element_text(size = 30, angle = 0, hjust = 0.5),
          axis.title = element_text(size = 30), axis.text.y = element_text(size = 30),
          plot.title = element_text(size = 30, hjust = 0.5), legend.position = "none") + 
    labs(x = "Fold change (Log2)" ,y = "p-value (-Log10)", title = ttl)     
  
  return(p)
}