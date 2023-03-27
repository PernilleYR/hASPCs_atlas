################################################################
#                                                              #
#                      UTILS DEGs SEURAT                       #
#                                                              #
################################################################



findMarkersOfClusts <- function(x, dataset = int, group.var = "batch", min.number.cells = 10, lgfc_treshold = 0.25){
  print(x)
  out <- FindConservedMarkers(dataset, ident.1 = x, grouping.var = group.var, 
                              verbose = T, min.cells.group = min.number.cells, 
                              logfc.threshold = lgfc_treshold)
  return(out)
}


add_avgFC <- function(x, DEGs = myDEGs_seurat){
  print(x)
  x <- DEGs[[x]]
  
  fc_col <- colnames(x)[grep("avg_log2FC", colnames(x))]
  
  dep <- unique(substr(fc_col, 1, 1))
  avg_logFC <- lapply(dep, function(i) return(rowMeans(x[,fc_col[startsWith(fc_col, i)], drop = F])))
  names(avg_logFC) <- paste0("avg_logFC_", dep)
  avg_logFC <- t(rlist::list.rbind(avg_logFC))
  x <- cbind(x, avg_logFC)
  x$avg_logFC_all <- rowMeans(x[,fc_col])
  
  log_avg_FC <- lapply(dep, function(i) return(log2(rowMeans(2^(x[,fc_col[startsWith(fc_col, i)], drop = F])))))
  names(log_avg_FC) <- paste0("log_avg_FC_", dep)
  log_avg_FC <- t(rlist::list.rbind(log_avg_FC))
  x <- cbind(x, log_avg_FC)
  x$log_avg_FC_all <- log2(rowMeans(2^(x[,fc_col])))
  
  x$geneID <- data.annot[rownames(x), "gene_short_name"]
  
  x <- x[order(x$avg_logFC_all, decreasing = T),]
  
  return(x)
}

plot_topN <- function(clust, topN = 20, colors =  mygradients_of_colors, data = myDEGs_seurat){
  if(clust %in% names(colors)){
    colors <- colors[[clust]]
  }else{
    colors <- colors[["other"]]
  }
  p <- lapply(rownames(data[[clust]])[1:topN], function(g)
    plot.colored_gene.continuous(int@reductions$tsne@cell.embeddings, 
                                 data.scale = GetAssayData(int, assay = "RNA"), 
                                 gene_name = g, data.annot = data.annot, gene.text = T, 
                                 col = colors))
  ggsave(gridExtra::marrangeGrob(p, nrow = 5, ncol = 4),width = 11.8, height = 11,
         filename = paste0("6.TopDEGs/plots/Top", topN,"_", clust, ".pdf" ))
}

