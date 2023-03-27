##---------------------------------------------##
##-----------------Functions-------------------##
##---------------------------------------------##

GOenrichment <- function(selectedGenes, allGenesList, ont = 'BP', 
                         topNodes = 100, nodeSize = 5, OrderBy = "classicFisher") {
  allGenesList_bg <- rep(1, length(allGenesList))
  names(allGenesList_bg) <- allGenesList
  allGenesList_bg[selectedGenes] <- 0.01
  
  tg.1 <- new("topGOdata", description = 'GO analysis',
              ontology =  ont, allGenes = allGenesList_bg,
              geneSel = function (x) {return (x < 0.05)},
              annot = annFUN.org ,
              nodeSize = nodeSize , # minimum number of genes in a GO categorie
              ID = "ENSEMBL", mapping = "org.Hs.eg.db")
  #GO.resKS <- runTest(tg.1, algorithm = "classic", statistic = "ks")
  GO.resF <- runTest(tg.1, algorithm = "classic", statistic = "fisher")
  #GO.resKS.elim <- runTest(tg.1, algorithm = "elim", statistic = "ks")
  GO.resF.elim <- runTest(tg.1, algorithm = "elim", statistic = "fisher")
  
  result <- GenTable(tg.1, classicFisher = GO.resF, eLimFisher = GO.resF.elim, orderBy = OrderBy,
                     ranksOf = "eLimFisher", topNodes = topNodes)
  return(result)
}
#' myEnrichRfunction
#' @description 
#' @param Markers list such as Markers_common (!! Please name your list by the cluster names)
#'                Markers[[1]]$geneID
#'                Markers[[2]]$geneID
#'                [...]            
#' @param libraryName databases from enrichR to test, default "Human_Gene_Atlas"
#'                  to check the different databases do:
#'                  require(enrichR)
#'                  dbs <- listEnrichrDbs()
#'                  if (is.null(dbs)) websiteLive <- FALSE else websiteLive <- TRUE
#'                  if (websiteLive) dbs$libraryName
#' @param n_genes number of genes send to enrichR (min between this number and number of Markers of the pop)
#' @param print if T print the results, default T
#' @param save boolean indicating if the data should be saved in .txt file
#' @param folder folder in which will save enrichR output
#' @columns_to_save the columns of enrichR to save (vector of numeric), default c(2,3,6)
#'                  1-"Index", 2-"Name", 3-"Adjusted_P-value", 4-"Z-score" 5-"Combined_Score", 6-"Genes", 7-"Overlap_P-value"
#' @return output of enrich & print them if print T
myEnrichRfunction <- function(Markers, libraryName = "Human_Gene_Atlas", 
                              n_genes = 50, print = T, save = T,
                              folder = NULL, columns_to_save = c(1,2,3,6)){
  
  # Test if Markers as names
  if(is.null(names(Markers))){
    stop("Please name your list of Markers")
  }
  
  # Load enrichR
  require(enrichR)
  # dbs <- listEnrichrDbs()
  # if (is.null(dbs)) websiteLive <- FALSE else websiteLive <- TRUE
  # if (websiteLive) dbs$libraryName
  
  # Perform enrichR for each clust
  enrichr_output <- list()
  for(i in names(Markers)){
    
    n <- min(length(Markers[[i]]$geneID), n_genes)
    
    if(n != n_genes){
      cat(paste0("!! Only ", n, " genes for clust ", i))
      cat("\n")
    }
    enrichr_output[[i]] <- enrichR::enrichr(as.character(Markers[[i]]$geneID[1:n]), 
                                            databases = libraryName)
  }
  
  # Print the output in the console
  if(print){
    for(i in names(Markers)){
      cat(paste("Clust",i,":", "\n"))
      print(enrichr_output[[i]][[1]][1:3, c("Term", "Adjusted.P.value")])
      cat("\n")
    }
  }
  
  # Save the ouput in .txt file 
  # if(save){
  #   # if folder not specified, saved in working directory
  #   if(is.null(folder)){
  #     folder <- getwd()
  #   }
  #   # check if folder ends with "/"
  #   if(!endsWith(folder, "/")){
  #     folder = paste0(folder,"/")
  #   }
  #   # save output
  #   sapply(names(enrichr_output),
  #          function(x) printEnrich(enrichr_output[[x]], 
  #                                  paste0(folder, "EnrichR_output_Clust", x,".txt" ), 
  #                                  sep = "\t", columns = columns_to_save))
  # }
  
  return(enrichr_output)
}
dotPlot_GO <- function(d, ref) {
  y <- d
  y$eLimFisher <-
    as.numeric(y$eLimFisher)
  y$classicFisher <- as.numeric(y$classicFisher)
  y <- y[order(y[, colnames(y) == ref], decreasing = F), ]
  y$Term <- factor(y$Term, levels = y$Term)
  c <- colnames(y)
  c[grep(ref, colnames(y))] <- "Pvalue"
  colnames(y) <- c
  
  p <- ggplot(y, aes(x = Pvalue, y = Term)) +
    geom_point(aes(size = Pvalue, color = "red")) +
    theme_bw(base_size = 14) + scale_size(trans = 'reverse')  + theme(legend.position = "none")
  
  return(p)
}