
calculate_summary <- function(d_ref, d_projected, percent = T){
  res <- data.frame(row.names = rownames(colData(myscesets[[d_projected]])),
                    projected = as.character(colData(myscesets[[d_projected]])$cell_type1),
                    projection_on_ref = scmap_res_1000feat[[paste0(d_projected,"on", d_ref)]]$scmap_cluster_labs[, d_ref])
  
  mysummary <- table(res$projection_on_ref, res$projected)
  if(percent){
    mysummary <- round(t(mysummary)/colSums(mysummary)*100,2)
    
  }
  
  mysummary <- reshape2::melt(mysummary); colnames(mysummary) <- c("projected", "projection_on_ref", "value")
  
  mysummary$dataset_projected <- d_projected
  mysummary$dataset_ref <- d_ref
  
  return(mysummary)
}

scmap_res_table <- function(scmap_res, in_percent = T){
  
  output <- list()
  for(i in names(scmap_res_1000feat)){
    ref <- strsplit(i, "on")[[1]][2]
    projected <- strsplit(i, "on")[[1]][1]
    
    output[[i]] <- calculate_summary(d_ref = ref, d_projected = projected, percent = in_percent )
  }
  
  return(output)
}
    
    for(ii in names(list_markers)){
      #Scan the samples  
      if(i != ii){
        output[[paste0(i,"-",ii)]] <- data.frame(popRef = character(),
                                                 popToCompare = character(),
                                                 percent = numeric(), 
                                                 stringsAsFactors = F)
        for(j in names(list_markers[[i]])){
          for(jj in names(list_markers[[ii]])){
            #Scan the populations
            x <- rownames(list_markers[[i]][[j]])
            y <- rownames(list_markers[[ii]][[jj]])
            if(!is.null(n_top_genes)){
              #Select minimum between n_top_genes and the lenght of the gene list
              n_genes_x <- min(n_top_genes, nrow(list_markers[[i]][[j]]))
              if(n_genes_x < n_top_genes){cat(paste0(i, " ",j , " had only ", n_genes_x, "\n"))}
              n_genes_y <- min(n_top_genes, nrow(list_markers[[ii]][[jj]]))
              if(n_genes_y < n_top_genes){cat(paste0(ii, " ",jj , " had only ", n_genes_y, "\n"))}
              x <- rownames(list_markers[[i]][[j]])[1:n_genes_x]
              y <- rownames(list_markers[[ii]][[jj]])[1:n_genes_y]
            }
            #Calculate the % of shared markers 
            per <- round(sum(x %in% y)/length(x)*100, 2)
            output[[paste0(i,"-",ii)]][nrow(output[[paste0(i,"-",ii)]])+ 1,] <- list(popRef = j,popToCompare = jj, percent = per)
          }
        }
      }
    }
  }
  return(output)
}






calculate_label_ypos <- function(values, breaks, order, min_value = 0, shift_text = 10){
  value_out <- c()
  t <- table(values$ref)
  t0 <- which(t == 0) 
  val <- data.frame()
  for(i in 1:max(t)){
    val <- rbind(val, values[(i*sum(t!=0)-sum(t!=0)+1):(i*sum(t!=0)-sum(t!=0)+sum(t!=0)),])
    if(identical(names(t0), character(0))){
      for(ii in names(t0)){
        val <- rbind(val, list(0, ii))
      }
    }
    
  }
  
  for(i in seq(0, nrow(val)-breaks, by = breaks)){
    v <- values[(i+1):(i + breaks),]; rownames(v) <- v[,2]
    w <- v[order,1, drop = F]
    w_cumsum <- cumsum(w[,1]); names(w_cumsum) <- rownames(w)
    value_out <- c(value_out,w_cumsum[as.character(v[,2])] )
  }
  value_out[ values[,1] <= min_value ] <- NA
  return(value_out - min(min_value, shift_text))
}


barplot_shared_markers <- function(percent, col, min_value = 0, shift_text_yaxis = 10, n_top_genes = NULL){
  str_names <- strsplit(names(percent), "-")
  percent <- percent[[1]]
  
  mysummary$label_ypos <- calculate_label_ypos(values = mysummary[,c("value", "ref")], breaks = length(levels(mysummary$ref)),
                                              order = rev(levels(mysummary$ref)), min_value = 2, shift_text = 1)
  
  
  # which pops are there?
  pops <- names(col)[names(col) %in% levels(as.factor(percent$popToCompare))]
  percent$popToCompare <- factor(percent$popToCompare, levels = pops)
  col <- unname(col[pops])
  
  p <- ggplot(mysummary, aes(x = ref, y = percent, fill = projection)) +
    geom_bar(stat = "identity") +
    #xlab(str_names[[1]][1]) + ylab(paste0("% of shared top", n_top_genes," markers")) +
    #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_text(aes(y = label_ypos, label = percent), col = "white",fontface = "bold") + 
    labs(fill = str_names[[1]][2]) +
    scale_fill_manual(values = col)
  return(p)
}

