################################################################
#                                                              #
#                 UTILS - 7. comparison mASPCs                 #
#                                                              #
################################################################

mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
convert_mouse_to_human <- function(gene_list){
  
  output = c()
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  cat( length(output), "homolog genes found /", length(gene_list), "\n")
  return (output)
}

convert_human_to_mouse <- function(gene_list){
  
  output = c()
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="human"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      mouse_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
      for(mouse_gene in mouse_genes){
        output = append(output,mouse_gene)
      }
    }
  }
  cat( length(output), "homolog genes found /", length(gene_list), "\n")
  return (output)
}

get_ensembl_names <- function(gene_IDs, data.annot){
  cat( sum(gene_IDs %in% data.annot$gene_short_name), "found genes/", length(gene_IDs), "\n")
  out <- data.annot[data.annot$gene_short_name %in% gene_IDs, ]
  return(out)
}

#' calculate_score:
#' @author pernille
#' @description Calculate the score as the sum of the expression of a list of genes
#' @param exprs expression maytrix 
#' @param gene_list list of genes per clust 
#'          mylist$clust1$gene.Gene.stable.ID.1
#'                $clust2$gene.Gene.stable.ID.1
#' @param n_markers number of markers to consider from the list (if NULL all are considered)
#' @return a matrix as column the cells as row the scores
calculate_score <- function(exprs, gene_list, n_markers = NULL){
  
  n <- min(n_markers, length(gene_list$ens_id))
  genes_oi_all <- gene_list$ens_id[1:n]
  genes_oi <- genes_oi_all[genes_oi_all %in% rownames(exprs)]
  
  temp1=exprs[ genes_oi,]
  temp1[is.na(temp1)] <- 0
  
  print(paste0(length(genes_oi),"/",length(genes_oi_all), " markers in data"))
  
  out <- colSums(temp1)
  
  out <- data.frame(row.names = names(out), 
                    score = out)
  return(out)
}

calculate_score_scaled <- function(exprs, gene_list, n_markers = NULL){
  
  n <- min(n_markers, length(gene_list$ens_id))
  genes_oi_all <- gene_list$ens_id[1:n]
  genes_oi <- genes_oi_all[genes_oi_all %in% rownames(exprs)]
  
  temp1=exprs[ genes_oi,]
  temp1[is.na(temp1)] <- 0
  
  print(paste0(length(genes_oi),"/",length(genes_oi_all), " markers in data"))
  
  out <- colSums(temp1)
  
  out <- out/nrow(temp1)
  
  out <- data.frame(row.names = names(out), 
                    score = out)
  return(out)
}

#' anovaTukey
#' Performs ANOVA + Tukey test for the give expression table for ONE gene
#' @param inTab input table, with columns Sample, Gene, Expr, Cluster
#' @return data frame with columns comp, p adj, start, end, Expr, suitable for
#' plotting with facets in ggplot2
anovaTukey <- function(inTab) {
  # first - ANOVA, to see, if ANY groups differ
  anovaRes <- aov(Expr ~ Cluster, data = inTab)
  anovaPval <- summary(anovaRes)[[1]]$`Pr(>F)`[1]
  # if yes - do Tukey
  tukeyRes <- data.table(comp = list(), `p adj` = NA, start = NA, end = NA, 
                         Expr = NA)
  if (anovaPval < 0.05) {
    tukeyRes <- as.data.table(TukeyHSD(anovaRes)$Cluster, keep.rownames = T)
    tukeyRes <- tukeyRes[`p adj` < 0.05, ]
    if (nrow(tukeyRes) >= 1) {
      tukeyRes[ , rn := lapply(rn, function(x) sort(strsplit(x, '-')[[1]]))]
      setnames(tukeyRes, 'rn', 'comp')
      tukeyRes <- tukeyRes[, c('comp', 'p adj')]
      
      # for facets
      tukeyRes[, start := sapply(comp, function(x) x[[1]])]
      tukeyRes[, end := sapply(comp, function(x) x[[2]])]
      maxExpr <- max(inTab$Expr)
      tukeyRes[, Expr := maxExpr * seq(1.05, by = 0.05, 
                                       length.out = nrow(tukeyRes))]
    }
  }
  tukeyRes
}

#' boxplot_score_datatable:
#' @author pernille
#' @description Function to plot score results 
#' @param score_res the score to plot, matrix with scores as row, cells as column
#' @param clust clustering results, factor
#' @param clust_name clusters names 
#' @param col colors of the different clusters
#' @param stat_on If T stats will be displayed, default T
#' @return creates a pdf in pdf_path and return a list with the ggplots & results of pairwise ttest (corrected BH)
#' @example boxplot_score(Gs_scores, ...)
#' @return list plot: ggplot
boxplot_score_datatable <- function(score_res, clust, clust_name, col, stat_on = T){
  
  print(rownames(score_res))
  mashaGgplot2Theme <- list(
    theme_classic(base_size = fontSize) + 
      theme(text = element_text(size = fontSize)) +
      theme(axis.line.x = element_line(colour = 'black', size = 0.5,
                                       linetype = 'solid'),
            axis.line.y = element_line(colour = 'black', size=0.5,
                                       linetype ='solid'),
            panel.grid.minor = element_line(colour = "white", size = 0.5,
                                            linetype = 2))
  )
  # re-shape the table to the "long" format
  score_res <- data.table::as.data.table(score_res, keep.rownames = T)
  score_res <- reshape2::melt(score_res, id.vars = 'rn')
  score_res <- data.table::as.data.table(score_res)
  data.table::setnames(score_res, colnames(score_res), c( 'Gene','Sample', 'Expr'))
  # add cluster information to the exprMatrShort
  clust <- data.table::data.table(Sample = names(clust), Cluster = clust)
  data.table::setkey(score_res, Sample)
  data.table::setkey(clust, Sample)
  score_res <- merge(score_res, clust, all.x = T)
  
  # calculate statistical significance via ANOVA + Tukey
  statsSign <- score_res[, anovaTukey(.SD), by = Gene]
  # adjust for the number of tested genes
  statsSign[, `p adj` := p.adjust(`p adj`, method = 'BH')]
  statsSign <- statsSign[`p adj` < 0.05]
  statsSign[, `p adj` := format.pval(`p adj`, 1)]
  
  P <- list()
  for(group in names(table(score_res$Gene))){
    print(group)
    df <- score_res[score_res$Gene == group,]
    stats_subset <- statsSign[statsSign$Gene == group,]
    p <- ggplot(df, aes(x = Cluster, y = Expr)) + 
      ggrastr::rasterise(geom_point(alpha = 0.02, aes(color = Cluster), 
                 position = position_jitterdodge(dodge.width = 0.8,
                                                 jitter.width = 2.5), 
                 size = 2 , pch = 16), dpi = 650) + 
      geom_boxplot(aes(fill = Cluster), alpha = 0.9, outlier.color = NA) + 
      scale_discrete_manual(labels = clust_name , values = col, 
                            aesthetics = c("fill", "colour"))  +
      ylab('Score') +
      mashaGgplot2Theme +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggtitle(group)
    
    if(!stat_on){
      p1 <- p
    }else if(nrow(stats_subset) == 0){
      p1 <- p
    }else{
      p1 <- p +
        ggsignif::geom_signif(mapping = aes(xmin = start, xmax = end, annotations = `p adj`,
                                            y_position = Expr), 
                              data = stats_subset, textsize = 2, vjust = -0.2, manual = T)
    }
    
    P[[group]] <- p1 + theme(legend.position = "none")
  }
  # ggsave( marrangeGrob(pGrid, ncol = 1, nrow = 1),
  #         width = 6.3, height = 6.25,
  #         filename = pdf_path)
  if(stat_on){
    return(statsSign)
  }else{
    return(P)    
  }

}


plot_score <- function(s, score_data, integrated_data = int){
  integrated_data$score <- score_data[colnames(integrated_data), s]
  plot.tsne.colored_continuous(integrated_data@reductions$tsne@cell.embeddings, 
                               data.categories = integrated_data$score, axis = T, tittle = s,
                               col = rev(c('#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4')))
}
