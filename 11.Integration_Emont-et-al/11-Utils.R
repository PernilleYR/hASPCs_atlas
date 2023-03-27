myfunction_to_convert <- function(d){
  c <- GetAssayData(d, assay = "RNA", slot = "counts")
  c <- c[rownames(c) %in% data.annot$gene_short_name,]
  rownames(c) <- convert_geneID_to_data.annot(rownames(c), data.annot = data.annot)[, "ens_id"]
  return(c)
}

#' prepareData
#' @description filter the matrices from list_matrices (keep cells from the corresponding
#' matric in list_original_data) and create seurat object, normalized & set variable genes from the original datasets
#' @param list_matrices list of matrices to prepare 
prepareData <- function(list_matrices){
  for(d in names(list_matrices)){
    print(d)
    list_matrices[[d]] <- CreateSeuratObject(list_matrices[[d]])
    list_matrices[[d]] <- NormalizeData(list_matrices[[d]], normalization.method = "LogNormalize", scale.factor = 10000)
    list_matrices[[d]] <- FindVariableFeatures(list_matrices[[d]], nfeatures = 2000 )
    }
  return(list_matrices)
}

#' anovaTukey
#' Performs ANOVA + Tukey test for the give expression table for ONE gene
#' @param inTab input table, with columns Sample, Gene, Expr, Cluster
#' @return data frame with columns comp, p adj, start, end, Expr, suitable for
#' plotting with facets in ggplot2
anovaTukey <- function(inTab) {
  # first - ANOVA, to see, if ANY groups differ
  anovaRes <- aov(value ~ clust, data = inTab)
  anovaPval <- summary(anovaRes)[[1]]$`Pr(>F)`[1]
  # if yes - do Tukey
  tukeyRes <- data.table(comp = list(), 
                         `p adj` = NA, 
                         start = NA, 
                         end = NA, 
                         Expr = NA)
  if (anovaPval < 0.05) {
    tukeyRes <- as.data.table(TukeyHSD(anovaRes)$clust, keep.rownames = T)
    tukeyRes <- tukeyRes[`p adj` < 1, ]
    if (nrow(tukeyRes) >= 1) {
      tukeyRes[ , rn := lapply(rn, function(x) sort(strsplit(x, '-')[[1]]))]
      setnames(tukeyRes, 'rn', 'comp')
      tukeyRes <- tukeyRes[, c('comp', 'p adj')]
      
      # for facets
      tukeyRes[, start := sapply(comp, function(x) x[[1]])]
      tukeyRes[, end := sapply(comp, function(x) x[[2]])]
      tukeyRes[, max_score := sapply(comp, function(x) max(inTab[clust %in% x, "value"]))]
      tukeyRes[, mean_1 := sapply(start, function(x) mean(inTab[clust == x, "value"]$value))]
      tukeyRes[, mean_2 := sapply(end, function(x) mean(inTab[clust == x, "value"]$value))]
      
    }
  }
  tukeyRes
}

find_enriched_clust <- function(d, clustering){
  print(clustering)
  # Find data oi
  c <- colnames(d@meta.data)[grep("prediction", colnames(d@meta.data))]
  c <- c[!c %in% c("prediction.score.max")]
  df <- d@meta.data[, c]
  df <- cbind("sample" = colnames(d),
              "clust" = d@meta.data[, clustering],
              df)
  # Reshape
  df <- reshape2::melt(df, id.vars = c("sample", "clust"))
  df <- rename(df, "variable" = "pop", "value" = "score")
  
  # Calculate avg and median value per pop.score and clust
  df <- df %>% group_by(pop, clust) %>% 
    mutate(Avg_by_pop.clust = mean(score), Med_by_pop.clust = median(score))
  df <- as.data.table(df)
  
  # Find the max pop.score for each clust
  max_avg <- as.data.frame(df %>% group_by(clust) %>% 
                             dplyr::slice(which.max(Avg_by_pop.clust)))
  
  # if IGFBP2.score is the highest for a clust run t.test
  if("prediction.score.IGFBP2" %in% max_avg$pop){
    print("YES!")}
  
  comparisons <- combn(unique(paste(as.character(df$pop))), 2)

  # t test
  clust_oi <- (df %>% filter(pop == "prediction.score.IGFBP2") %>% dplyr::slice(which.max(Avg_by_pop.clust)))$clust
  
  tTestPvalsAll <- apply(comparisons, 2, function(x){
    testPval <- t.test(df[pop == x[1] & clust == as.character(clust_oi)]$score,
                       df[pop == x[2] & clust == as.character(clust_oi)]$score)$p.value
    return(c(x, as.character(clust_oi), testPval))
  })
  
  tTestPvalsAll <- as.data.frame(t(tTestPvalsAll), stringsAsFactors = F)
  colnames(tTestPvalsAll) <- c("pop_1", "pop_2", "clust",'pValue')
  
  # Filtering of comparisons with IGFBP2
  tTestPvalsAll <- tTestPvalsAll %>% filter(pop_1 == "prediction.score.IGFBP2" | pop_2 == "prediction.score.IGFBP2")
  
  # Adjust P value (FDR)
  tTestPvalsAll$p.adjusted <- p.adjust(tTestPvalsAll$pValue, method = 'fdr')
  
  # if(all(tTestPvalsAll$p.adjusted < 0.05)){
  #   print("IGFBP2 WON!")}
  
  return(list("data" = df,
              "max_avg" = max_avg,
              "ttest" = tTestPvalsAll))
}

plot_enriched_cluster <- function(d = NULL, clustering, c, enr_c = enriched_clust, seu = mylist){
  
  if(!is.null(d)){
    seu <- seu[[d]]
    enr_c <- enr_c[[d]]
  }
  
  dat <- enr_c[[clustering]]$data %>% filter(clust == c)
  dat$pop <- sapply(strsplit(as.character(dat$pop),"\\."), `[`, 3)
  dat$pop[dat$pop == "PR"] <- "PR specific"
  dat$pop[dat$pop == "res"] <- "res.0.2_8"
  dat$pop[dat$pop == "CHI3L1"] <- "CHI3L1-2"
  dat$pop <- factor(dat$pop, levels = names(myIntegratedColors))
  
  pb <- ggplot(dat, aes(y = score, x = clust, fill = pop)) + 
    geom_boxplot() + #geom_jitter() +
    scale_fill_manual(values = myIntegratedColors) + 
    mashaGgplot2Theme
  ps <- plot.tsne.colored_continuous(seu@reductions$umap@cell.embeddings, 
                                     axis = T, x.lab = "tSNE1", y.lab = "tSNE2",
                                     col = c("gray88",'#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b'),
                                     data.categories = seu$prediction.score.IGFBP2)
  seu$toplot <- "Rest"; seu$toplot[seu@meta.data[, clustering] == c] <- "IGFBP2"
  seu$toplot <- factor( seu$toplot, levels = c("Rest", "IGFBP2"))
  pc <- plot.tsne.colored_discrete_2(seu@reductions$umap@cell.embeddings,
                                     axis = T, x.lab = "tSNE1", y.lab = "tSNE2", 
                                     data.categories = seu$toplot,
                                     col = c("gray88","#1f78b4"))
  P <- list(wleg= grid.arrange(ps,pc,pb, ncol = 3),
            woleg = grid.arrange(ps+theme(legend.position = "none"),
                                 pc+theme(legend.position = "none"),
                                 pb+theme(legend.position = "none"), ncol = 3))
  ggsave(P$woleg, filename = paste0("11.Integration-emont/enriched_in_clust/Plot_", d, "_clustering.", clustering, "_clust", c, "_woleg.pdf"),
         height = 2.51, width = 7.68)
  ggsave(P$wleg, filename = paste0("11.Integration-emont/enriched_in_clust/Plot_", d, "_clustering.", clustering, "_clust", c, "_leg.pdf"),
         height = 2.51, width = 7.68)
  return(P)
}

