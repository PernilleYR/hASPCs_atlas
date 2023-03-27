################################################################
#                                                              #
#                sc map between all samples                    #
#                                                              #
################################################################


### Author: Pernille
### Date: 10.08.2022
### Datasets: scRNA-seq Depots: SC - EP - MK - MG - PR - GB 
###                     Patients: B0, B1, B7, L3, L11, L12 (B - Bariatric, L - Lean)
### Goal: Compute scmap between all samples

library(scran)
library(Seurat)
library(ggplot2)
library(scmap)

setwd("~/SVRAW1/prainer/hASPCs/PAPER/10X_scRNA-seq/")

source("3.scmap/utils.R")
source("Utility/General_utils.R")

##---------------------------------------------##
##----------------Loading data-----------------##
##---------------------------------------------##

myseu <- readRDS("0.data/List_seurat_objects.rds")
data.annot <- read.table("~/SVRAW1/prainer/Files/Human/data.annot/Homo_sapiens.GRCh38.92_data.annot.txt")

##---------------------------------------------##
##-----------------Loading data----------------##
##---------------------------------------------##

clust <- list()
for(i in names(clust_oi)){
  c <- Bar017_seurats[[i]]@meta.data[,clust_oi[i]]
  names(c) <- colnames(Bar017_seurats[[i]])
  clust[[i]] <- c 
}

##---------------------------------------------##
##---------SingleCellExperiment objects--------##
##---------------------------------------------##  

myscesets <- list()
for(d in myDepots){
  ### --- Create SingleCellExperiment Object
  myscesets[[d]] <- SingleCellExperiment(assays = list(counts = as.matrix(GetAssayData(myseu[[d]], slot = "counts"))))
  rowData(myscesets[[d]])$feature_symbol = data.annot[rownames(myscesets[[d]]),"gene_short_name"]
  ### --- Estimate normalisation factors 
  if(d != "GB7"){clusters <- quickCluster(myscesets[[d]])}
  ### --- Compute size factor based on clusters
  myscesets[[d]] <- computeSumFactors(myscesets[[d]])
  ### --- normalize(): Compute normalised expression values from an SCESet object using the size factors stored in the object. Return the object with the normalised expression values added.
  myscesets[[d]] <- logNormCounts(myscesets[[d]])
}

##---------------------------------------------##
##----------------Scmap Cluster----------------##
##---------------------------------------------##

### --- Prepare the data
# - Add clustering to SingleCellExperiments
for(d in c("PR11", "PR12")){
  myscesets[[d]]$cell_type1 <- myseu[[d]]$mySelectedClustering[colnames(myscesets[[d]])]
}
rm(d)

# - Feature selection
# Select 500 features by default in a simplified M3Drop 
#Let's select 1000
for(d in myDepots){
  myscesets[[d]] <- selectFeatures(myscesets[[d]], suppress_plot = FALSE, n_features = 1000)
}

# - Calculate the number of shared features
selectedFeatures <- list()
for(d in myDepots){
  selectedFeatures[[d]] <- as.character(rowData(myscesets[[d]])[["feature_symbol"]][rowData(myscesets[[d]])[["scmap_features"]]])
}

common_feat <- c(); names_feat <- c(); pops <- myDepots
for(i in 1:length(myDepots)){
  for(ii in 1:length(myDepots)){
    if( i < ii ){
      common_feat <- c(common_feat, sum(selectedFeatures[[i]] %in% selectedFeatures[[ii]]))
      names_feat <- c(names_feat, paste0(pops[i], "-", pops[[ii]]))
      cat(tail(names_feat, 1), "\n")
    }
  }
}
names(common_feat) <- names_feat

rm(selectedFeatures, names_feat, pops, i, ii)
common_feat

#1000 features:
# EP0-EP1   EP0-EP7   EP0-SC0   EP0-SC1   EP0-SC7   EP0-MG7   EP0-MK7   EP0-GB7   EP0-PR3  EP0-PR11  EP0-PR12   EP1-EP7   EP1-SC0   EP1-SC1 
# 786       426       680       535       376       376       391       286       268       320       341       512       680       625 
# EP1-SC7   EP1-MG7   EP1-MK7   EP1-GB7   EP1-PR3  EP1-PR11  EP1-PR12   EP7-SC0   EP7-SC1   EP7-SC7   EP7-MG7   EP7-MK7   EP7-GB7   EP7-PR3 
# 469       464       474       332       340       401       431       487       530       715       786       757       489       509 
# EP7-PR11  EP7-PR12   SC0-SC1   SC0-SC7   SC0-MG7   SC0-MK7   SC0-GB7   SC0-PR3  SC0-PR11  SC0-PR12   SC1-SC7   SC1-MG7   SC1-MK7   SC1-GB7 
# 639       648       586       520       484       517       353       344       416       446       592       529       559       373 
# SC1-PR3  SC1-PR11  SC1-PR12   SC7-MG7   SC7-MK7   SC7-GB7   SC7-PR3  SC7-PR11  SC7-PR12   MG7-MK7   MG7-GB7   MG7-PR3  MG7-PR11  MG7-PR12 
# 419       510       523       766       737       481       546       662       692       821       510       564       686       681 
# MK7-GB7   MK7-PR3  MK7-PR11  MK7-PR12   GB7-PR3  GB7-PR11  GB7-PR12  PR3-PR11  PR3-PR12 PR11-PR12 
# 503       555       668       667       439       463       470       687       612       774 

#700 features:
# EP0-EP1 EP0-EP7 EP0-SC0 EP0-SC1 EP0-SC7 EP0-MG7 EP0-MK7 EP0-GB7 EP1-EP7 EP1-SC0 EP1-SC1 EP1-SC7 EP1-MG7 EP1-MK7 EP1-GB7 EP7-SC0 
# 537     319     436     378     275     291     298     210     378     456     451     342     352     356     258     385 
# EP7-SC1 EP7-SC7 EP7-MG7 EP7-MK7 EP7-GB7 SC0-SC1 SC0-SC7 SC0-MG7 SC0-MK7 SC0-GB7 SC1-SC7 SC1-MG7 SC1-MK7 SC1-GB7 SC7-MG7 SC7-MK7 
# 373     486     538     521     359     437     414     390     411     273     424     386     405     272     521     507 
# SC7-GB7 MG7-MK7 MG7-GB7 MK7-GB7 
# 343     563     378     371

# - Index the reference dataset
# Find the median gene expression for each cluster (use cell_type1 by default)
for(d in myDepots){
  myscesets[[d]] <- indexCluster(myscesets[[d]], cluster_col = "cell_type1")
}
rm(d)

# - Heatmap cluster index
# pheatmap::pheatmap(metadata(myscesets[["SC0"]])$scmap_cluster_index)
# pheatmap::pheatmap(metadata(myscesets[["GB7"]])$scmap_cluster_index)
# pheatmap::pheatmap(metadata(myscesets[["MK7"]])$scmap_cluster_index)
# pheatmap::pheatmap(metadata(myscesets[["EP7"]])$scmap_cluster_index)


# - Projection 
scmap_res_1000feat <- list()
n = 1
names_res <- c()
for(i in myDepots){
  for(ii in myDepots){
    if(i != ii){
      
      ind_list <- list(metadata(myscesets[[ii]])$scmap_cluster_index)
      names(ind_list) <- ii
      
      scmap_res_1000feat[[n]] <- scmapCluster( projection = myscesets[[i]],
                                              index_list = ind_list )
      
      names_res[n] <- paste0(i, "on", ii)
      cat(n, ": ",tail(names_res,1), "\n")
      
      rm(ind_list)
      n = n + 1
    }
  }
}
names(scmap_res_1000feat) <- names_res

rm(i, ii, n, names_res)

##---------------------------------------------##
##-----------------Sankey plot-----------------##
##---------------------------------------------##
# - Sankey plot
for( i in myDepots){
  for(ii in myDepots){
    if(i != ii){
      plot(getSankey( colData(myscesets[[i]])$cell_type1,
                      scmap_res_1000feat[[paste0(i, "on", ii)]]$scmap_cluster_labs[, ii]))
      #colors = col[[i]]))
    }
  }
}
rm(i,ii)

##---------------------------------------------##
##-----------Calculate data summary------------##
##---------------------------------------------##

projection_results_percent <- scmap_res_table(scmap_res_1000feat, in_percent = T)
projection_results_ncells <- scmap_res_table(scmap_res_1000feat, in_percent = F)

saveRDS(scmap_res_1000feat, "output_table/scmap_res.rds")
saveRDS(projection_results_percent, "output_table/projection_results_in-percent.rds")
saveRDS(projection_results_ncells, "output_table/projection_results_in-number-of-cells.rds")
