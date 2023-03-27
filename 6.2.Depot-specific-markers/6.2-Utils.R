FindMarkers_by_depots <- function(cell_pop, depot, all_dep_to_compare = c("SC", "EP", "MG", "PR")){
  cp_d_oi <- paste0(cell_pop, "_", depot)
  
  if(cp_d_oi %in% levels(int$myIntegratedCluster_depot)){
    ident.2 <- paste0(cell_pop, "_", all_dep_to_compare[all_dep_to_compare != depot])
    ident.2 <- ident.2[ident.2 %in% levels(int$myIntegratedCluster_depot)]
    
    out <- FindMarkers(object = int, ident.1 = cp_d_oi, ident.2 = ident.2)
    
    out <- out[order(out$avg_log2FC, decreasing = T),]
    
    out$signi <- T
    out$signi[out$p_val_adj > 0.05] <- F
    
    out$up <- F
    out$up[out$avg_log2FC > 0] <- T
    
    out$geneID <-  data.annot[rownames(out), "gene_short_name"]
    
  }else{
    out <- NULL
  }
  
  return(out)
}

FindMarkers_by_depotType <- function(cell_pop){
  cp_d_oi <- paste0(cell_pop, "_", c("SC", "PR"))
  
  if(sum(cp_d_oi %in% levels(int$myIntegratedCluster_depot))){
    ident.2 <- paste0(cell_pop, "_", c("OM", "MK"))
    ident.2 <- ident.2[ident.2 %in% levels(int$myIntegratedCluster_depot)]
    
    out <- FindMarkers(object = int, ident.1 = cp_d_oi, ident.2 = ident.2)
    
    out$geneID <-  data.annot[rownames(out), "gene_short_name"]
    
    out <- out[order(out$avg_log2FC, decreasing = T),]
    
    out$signi <- T
    out$signi[out$p_val_adj > 0.05] <- F
    
    out$up <- F
    out$up[out$avg_log2FC > 0] <- T
    
  }else{
    out <- NULL
  }
  
  return(out)
}

finding_specific_markers_v1 <- function(cp, dp_oi, all_depots_to_compare = c("SC", "EP", "PR", "MG")){
  # Check if the cell type is in depot of interest
  int <- SetIdent(int, value = "myIntegratedCluster_depot")
  if(paste0(cp, "_", dp_oi) %in% levels(int$myIntegratedCluster_depot)){
    dep_to_compare <- all_depots_to_compare[!all_depots_to_compare %in% dp_oi]
    # Run DE: for the depot where cell type exists
    DE_vsEachDepot <- lapply(dep_to_compare, function(d){
      if(paste0(cp, "_", d) %in% levels(int$myIntegratedCluster_depot)){
        out <- FindMarkers_by_depots(cp, dp_oi, all_dep_to_compare = d)
        colnames(out)[!colnames(out) %in% c("pct.1", "geneID")] <- paste0(colnames(out)[!colnames(out) %in% c("pct.1", "geneID")], ".", d)
        return(out)
      }else{
        return(data.frame("geneID" = NA))
      }
    })
    names(DE_vsEachDepot) <- dep_to_compare
    genes_to_keep <- DE_vsEachDepot[[1]][DE_vsEachDepot[[1]][,grep("signi", colnames(DE_vsEachDepot[[1]]))], "geneID"]
    for(d in names(DE_vsEachDepot)[-1]){
      genes_to_keep <- intersect(genes_to_keep,
                                 DE_vsEachDepot[[d]][DE_vsEachDepot[[d]][,grep("signi", colnames(DE_vsEachDepot[[d]]))], "geneID"])
    }
    DE_vsEachDepot_tobind <- lapply(DE_vsEachDepot, function(x){
      x$ens_id <- rownames(x)
      rownames(x) <- x$geneID
      x <- x[genes_to_keep,]
      return(x)
    })
    DE_res <- plyr::join_all(DE_vsEachDepot_tobind,  type='left')
    
    DE_res$avg_log2FC.all <- rowMeans(DE_res[, grep("avg_log2FC", colnames(DE_res))])
  }else{
    DE_res <- NULL
  }
  return(DE_res)
}

order_g <- function(d,m){
  dep_compared <- unique(sapply(strsplit(colnames(d),"\\."), `[`, 2))
  dep_compared <- dep_compared[dep_compared %in% c("SC", "PR", "EP", "MG")]
  
  # to_keep_signi <- rowSums(d[,paste0("signi.",dep_compared)]) >= length(dep_compared)
  # d <- d[to_keep_signi,]
  
  to_keep_up <- rowSums(d[,paste0("up.",dep_compared)]) >= length(dep_compared)
  d <- d[to_keep_up,]
  
  if(m == "avg"){
    d <- d[order(d$avg_log2FC.all, decreasing = T),]
  }
  
  if(m == "rank"){
    for(dep in dep_compared){
      d[, paste0("rank.", dep)] <- order(d[, paste0("avg_log2FC.", dep)], decreasing = T)
    }
    d$final.rank <- rowSums(d[, paste0("rank.", dep_compared)])
    d <- d[order(d$final.rank, decreasing = F), ]
  }
  return(d)
}

plot_heatmap_pop <- function(pop, n_top = 30, method = "avg", data){
  # Prepare avg expression data
  int <- SetIdent(int, value = "myIntegratedClustering")
  pop_sub <- subset(int, idents = pop)
  pop_sub <- SetIdent(pop_sub, value = "depot")
  avg.exp.cells <- log1p(AverageExpression(pop_sub, verbose = FALSE)$RNA)
  avg.exp.cells <- as.data.frame(avg.exp.cells)
  
  # Select genes
  data <- lapply(data, function(d) order_g(d = d, m = "rank"))
  to_plot <- unlist(lapply(data, function(d) return(d$ens_id[1:n_top])))
  
  # Subset data for genes oi
  to_plot <- avg.exp.cells[to_plot, ]
  rownames(to_plot) <- data.annot[rownames(to_plot), "gene_short_name"]
  
  p <- pheatmap::pheatmap(to_plot, cluster_rows = F, cluster_cols = F, fontsize_row = 6, border_color = NA)
  p_s <- pheatmap::pheatmap(to_plot, cluster_rows = F, cluster_cols = F, scale = "row", fontsize_row = 6, border_color = NA)
  
  return(list("p_non.scaled" = p,
              "p_scaled" = p_s ))
  
}
