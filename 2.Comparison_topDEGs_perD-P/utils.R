################################################################
#                                                              #
#                       UTILS ---- PART 2                      #
#           Comparison barplot of top DEGs per D - P           #
#                                                              #
################################################################

##---------------------------------------------##
##------------------Functions------------------##
##---------------------------------------------##
calculate_label_ypos <- function(values, breaks, order, min_value = 0, shift_text = 10){
  value_out <- c()
  for(i in seq(0, (nrow(values)-breaks), by = breaks)){
    v <- values[(i+1):(i + breaks),]; rownames(v) <- v[,2]
    w <- v[order,1, drop = F]
    w_cumsum <- cumsum(w[,1]); names(w_cumsum) <- rownames(w)
    value_out <- c(value_out,w_cumsum[as.character(v[,2])] )
  }
  value_out[ values[,1] <= min_value ] <- NA
  return(value_out - shift_text)
}
calculate_percent_shared_markers <- function(list_markers, n_top_genes = NULL){
  ## list_markers: a list containing the top markers of the different sample in such form:
  # list_markers$sample1$pop1$gene
  #                     $pop2$gene
  #                     $[...]$gene
  #             $sample2$pop1$gene
  #                     $pop2$gene
  #                     $[...]$gene
  ##n_top_genes: the number of top genes to consider (by default all the genes of the list are considered)
  
  ## Calculate the percentage of markers between the populations of all the samples of a list 
  ## (!DO NOT CONSIDER MARKERS THAT ARE ASSIGNED TO MORE THAN ONE POP)
  
  output <- list()
  for(i in names(list_markers)){
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
calculate_percent_shared_markers_common <- function(list_markers, n_top_genes = NULL, list_common){
  ## list_markers: a list containing the top markers of the different sample in such form:
  # list_markers$sample1$pop1$gene
  #                     $pop2$gene
  #                     $[...]$gene
  #             $sample2$pop1$gene
  #                     $pop2$gene
  #                     $[...]$gene
  ##n_top_genes: the number of top genes to consider (by default all the genes of the list are considered)
  
  ## Calculate the percentage of markers between the populations of all the samples of a list 
  ## If a marker is assigned to more than one sample and is in shared with the markers of a population of another sample it goes in 
  ## the category "common" on the barplot 
  
  output <- list()
  for(i in names(list_markers)){
    for(ii in names(list_markers)){
      #Scan the samples  
      if(i != ii){
        output[[paste0(i,"-",ii)]] <- data.frame(popRef = character(),
                                                 popToCompare = character(),
                                                 percent = numeric(), 
                                                 stringsAsFactors = F)
        for(j in names(list_markers[[i]])){
          commons <-c()
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
            shared <- x[x %in% y]
            new_commons <- shared[shared %in% names(list_common[[ii]])]
            commons <- unique(c(commons, new_commons))
            
            per <- round((sum(x %in% y)-length(new_commons))/length(x)*100, 2)
            output[[paste0(i,"-",ii)]][nrow(output[[paste0(i,"-",ii)]])+ 1,] <- list(popRef = j,popToCompare = jj, percent = per)
          }
          output[[paste0(i,"-",ii)]][nrow(output[[paste0(i,"-",ii)]])+ 1,] <-list(popRef = j,popToCompare = "common", percent = round(length(commons)/length(x)*100,2))
        }
      }
    }
  }
  return(output)
}
calculate_percent_shared_markers_common_v2 <- function(list_markers, n_top_genes = NULL, list_common){
  ## list_markers: a list containing the top markers of the different sample in such form:
  # list_markers$sample1$pop1$gene
  #                     $pop2$gene
  #                     $[...]$gene
  #             $sample2$pop1$gene
  #                     $pop2$gene
  #                     $[...]$gene
  ##n_top_genes: the number of top genes to consider (by default all the genes of the list are considered)
  
  ## Calculate the percentage of markers between the populations of all the samples of a list 
  ## If a marker is a top marker of several populations we have two cases: 
  ##     - the marker is also a marker of the same population (for example SC0 darkgreen has a marker assigned to both Darkgreen and VSMP of SC1)
  ##       in that case it is considered as a marker of that population (darkgreen in the example)
  ##     - the marker is not a marker of the same population (for example SC0 darkgreen has a marker assigned to both LightGreen and VSMP of SC1)
  ##       in that case it is defined as "common" on the barplot
  
  output <- list()
  for(i in names(list_markers)){
    for(ii in names(list_markers)){
      #Scan the samples  
      if(i != ii){
        output[[paste0(i,"-",ii)]] <- data.frame(popRef = character(),
                                                 popToCompare = character(),
                                                 percent = numeric(), 
                                                 stringsAsFactors = F)
        for(j in names(list_markers[[i]])){
          commons <-c()
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
            shared <- x[x %in% y]
            new_commons <- shared[shared %in% names(list_common[[ii]])]
            commons <- unique(c(commons, new_commons))
            
            per <- round((sum(x %in% y)-length(new_commons))/length(x)*100, 2)
            output[[paste0(i,"-",ii)]][nrow(output[[paste0(i,"-",ii)]])+ 1,] <- list(popRef = j,popToCompare = jj, percent = per)
          }
          commons <- doubleAssigned[[ii]][commons]
          additional_counts <- 0
          for(k in names(commons)){
            pops <- strsplit(commons[k], ", ")
            if(j %in% pops[[1]]){
              additional_counts <- additional_counts + 1
            }
          }
          ind <- which(output[[paste0(i,"-",ii)]]$popRef == j & output[[paste0(i,"-",ii)]]$popToCompare == j)
          output[[paste0(i,"-",ii)]][ind,"percent"] <- output[[paste0(i,"-",ii)]][ind,"percent"] + round(additional_counts/length(x)*100,2)
          output[[paste0(i,"-",ii)]][nrow(output[[paste0(i,"-",ii)]])+ 1,] <-list(popRef = j,popToCompare = "common", percent = round((length(commons)-additional_counts)/length(x)*100,2))
        }
      }
    }
  }
  return(output)
}

barplot_shared_markers <- function(percent, col, min_value = 0, shift_text_yaxis = 10, n_top_genes = NULL){
  str_names <- strsplit(names(percent), "-")
  percent <- percent[[1]]
  
  percent$label_ypos <- calculate_label_ypos(values = percent[,c("percent", "popToCompare")], breaks = length(table(percent$popToCompare)),
                                             order = rev(levels(as.factor(percent$popToCompare))), 
                                             min_value = min_value, shift_text = shift_text_yaxis)
  
  # which pops are there?
  pops <- names(col)[names(col) %in% levels(as.factor(percent$popToCompare))]
  percent$popToCompare <- factor(percent$popToCompare, levels = pops)
  col <- unname(col[pops])
  
  p <- ggplot(percent, aes(x = popRef, y = percent, fill = popToCompare)) +
    geom_bar(stat = "identity") +
    xlab(str_names[[1]][1]) + ylab(paste0("% of shared top", n_top_genes," markers")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_text(aes(y = label_ypos, label = percent), col = "white",fontface = "bold") + 
    labs(fill = str_names[[1]][2]) +
    scale_fill_manual(values = col)
  return(p)
}

detect_doubleAssigned_markers <- function(list_markers){
  output <- c()
  for(i in 1:length(list_markers)){
    for(ii in 1:length(list_markers)){
      if(ii > i){
        shared_genes <- rownames(list_markers[[i]])[rownames(list_markers[[i]]) %in% rownames(list_markers[[ii]])]
        # Find genes already in output
        present_output <- shared_genes[shared_genes %in% names(output)]
        not_present_output <- shared_genes[!shared_genes %in% names(output)]
        output[present_output] <- paste0(output[present_output], ", ", names(list_markers)[i], ", ",names(list_markers)[ii])
        names_output <- names(output)
        output <- c(output, rep(paste0(names(list_markers)[i],", ",names(list_markers)[ii]), length(not_present_output)))
        names(output) <- c(names_output, not_present_output)
      }
    }
  }
  if(length(output) != 0){
    for(j in 1:length(output)){
      str_pops <- strsplit(output[j], ", ")
      output[j] <- paste0(unique(str_pops[[1]]), collapse = ", ")
    }
  }
  
  return(output)
}
