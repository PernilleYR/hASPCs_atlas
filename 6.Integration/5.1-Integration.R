################################################################
#                                                              #
#                         Integration                          #
#                                                              #
################################################################

### Author: Pernille
### Date: 10.08.2022
### Datasets: scRNA-seq Depots: SC - EP - MK - MG - PR - GB 
###                     Patients: B0, B1, B7, L3, L11, L12 (B - Bariatric, L - Lean)
### Goal: Integration

library(ggplot2); library(data.table); library(harmony); library(Seurat); library(stringr)

setwd("~/SVRAW1/prainer/hASPCs/PAPER/10X_scRNA-seq/")
source("Utility/General_utils.R")
source("5.Integration/5-utils.R")

##---------------------------------------------##
##-----------------Loading data----------------##
##---------------------------------------------##

## Seurat objs
myseu <- readRDS("0.data/List_seurat_objects.rds")
myseu <- lapply(myDepots, function(d) 
  FindVariableFeatures(myseu[[d]], 
                       selection.method = "vst",
                       nfeatures = 2000))
names(myseu) <- myDepots

## Raw data
setwd("~/SVRAW1/prainer/hASPCs/10X/")
SC0_r <- Read10X_h5("10X_Bar0_141118/cellRangerOutput/output_Bar0_SC0_141118/outs/filtered_gene_bc_matrices_h5.h5", use.names = F)
SC1_r <- Read10X_h5("10X_Bar1_201118/cellRangerOutput/output_Bart1_SC1_201118/outs/filtered_gene_bc_matrices_h5.h5", use.names = F)
SC7_r <- Read10X_h5("10X_Bar7_SC_190613/output_Bar7_SC/outs/filtered_feature_bc_matrix.h5", use.names = F)
EP0_r <- Read10X_h5("10X_Bar0_141118/cellRangerOutput/output_Bar0_EP0_141118/outs/filtered_gene_bc_matrices_h5.h5", use.names = F)
EP1_r <- Read10X_h5("10X_Bar1_201118/cellRangerOutput/output_Bart1_EP1_201118/outs/filtered_gene_bc_matrices_h5.h5", use.names = F)
EP7_r <- Read10X_h5("10X_Bar7_EP_190613/output_Bar7_EP/outs/filtered_feature_bc_matrix.h5", use.names = F)
MG7_r <- Read10X_h5("10X_Bar7_MG_190613/output_Bar7_MG/outs/filtered_feature_bc_matrix.h5", use.names = F)
MK7_r <- Read10X_h5("10X_Bar7_MK_190613/output_Bar7_MK/outs/filtered_feature_bc_matrix.h5", use.names = F)
#GB7_r <- Read10X_h5("10X_Bar7_GB_190613/output_Bar7_GB/outs/filtered_feature_bc_matrix.h5", use.names = F)
PR3_12 <- Read10X_h5("PR_LP12_03_frozen/RScriptData/filtered_feature_bc_matrix.h5", use.names = F)
PR11_r <- Read10X_h5("PR_LP11/filtered_feature_bc_matrix.h5", use.names = F)

myseu_raw <- list("EP0" = EP0_r, "EP1" = EP1_r, "EP7" = EP7_r,
                  "SC0" = SC0_r, "SC1" = SC1_r, "SC7" = SC7_r,
                  "MG7" = MG7_r, "MK7" = MK7_r, #"GB7" = GB7_r,
                  "PR3" = PR3_12, "PR11" = PR11_r, "PR12" = PR3_12)

rm(EP0_r, EP1_r, EP7_r, 
   SC0_r, SC1_r, SC7_r,
   MK7_r, MG7_r, GB7_r,
   PR11_r, PR3_12)

##---------------------------------------------##
##------------Create new list objs-------------##
##---------------------------------------------##

myseu_raw <- prepareData(list_matrices = myseu_raw, list_original_data = myseu)

##---------------------------------------------##
##-----------Extract clustering info-----------##
##---------------------------------------------##

myDepots <- myDepots[-which(myDepots == "GB7")]

recover_clustering <- function(d){
  myclust <- myseu[[d]]$mySelectedClustering
  output <- data.frame(integration_names = paste0(d, "_", colnames(myseu[[d]])),
                       original_names = colnames(myseu[[d]]),
                       mySelectedClustering = myseu[[d]]$mySelectedClustering,
                       mySelectedClustering_perSample = paste0(d, "_", myseu[[d]]$mySelectedClustering))
}
clustering_metadata <- as.data.frame(rbindlist(lapply(myDepots, recover_clustering)))
rownames(clustering_metadata) <- clustering_metadata$integration_names
clustering_metadata$mySelectedClustering_perSample <- factor(clustering_metadata$mySelectedClustering_perSample,
                                                             levels = c(paste0(myDepots, "_", "ASCs"),
                                                                        paste0(myDepots, "_", "PreAs"),
                                                                        paste0(myDepots[startsWith(myDepots, "EP")], "_", "IGFBP2"),
                                                                        paste0(myDepots[startsWith(myDepots, "EP")], "_", "Meso"),
                                                                        paste0(myDepots, "_", "VSMPs"),
                                                                        paste0(myDepots[grep("7", myDepots)], "_", "Endo"),
                                                                        paste0(myDepots[grep("7", myDepots)], "_", "Immune"),
                                                                        "PR11_Unknown_VSMPs"))
clustering_metadata$mySelectedClustering_perSample_simplified <- clustering_metadata$mySelectedClustering_perSample
levels(clustering_metadata$mySelectedClustering_perSample_simplified) <- c(levels(clustering_metadata$mySelectedClustering_perSample_simplified)[1:39], #rep("VSMPs", length(myDepots)), 
                                                                           rep("Endo", length(myDepots[grep("7", myDepots)])),
                                                                           rep("Immune", length(myDepots[grep("7", myDepots)])), 
                                                                           "PR11_VSMPs")

#####
# ##---------------------------------------------##
# ##----------Create obj for integration---------##
# ##---------------------------------------------##
# 
# # Find features
# myDepots <- myDepots[-which(myDepots == "GB7")]
# myanchorfeatures <- c()
# for(d in myDepots){
#   myanchorfeatures <- unique(c(myanchorfeatures, VariableFeatures(myseu[[d]])))
# }
# 
# # Bind counts
# bind_counts <- GetAssayData(myseu_raw$EP0, slot = "counts")
# for(d in myDepots[-1]){
#   bind_counts <- cbind(bind_counts, GetAssayData(myseu_raw[[d]], slot = "counts"))
# }
# rm(d)
# 
# # Create obj
# harmony_int <- CreateSeuratObject(counts = bind_counts, min.cells = 3) 
# harmony_int <- FindVariableFeatures(harmony_int)
# 
# # Add metadata
# harmony_int$dataset <- factor(sapply(strsplit(colnames(harmony_int),"_"), `[`, 1), levels = myDepots)
# 
# # Add variable features (i.e., in features)
# #VariableFeatures(harmony_int) <- myanchorfeatures
# 
# # Clean space 
# rm(myseu, myseu_raw, bind_counts, myanchorfeatures)
# 
# # Normalize, scale and PCA
# harmony_int <- harmony_int %>% NormalizeData(verbose = FALSE)  %>% 
#   ScaleData(verbose = FALSE) %>% 
#   RunPCA(pc.genes = VariableFeatures(harmony_int), npcs = 60, verbose = FALSE)
# 
# ##---------------------------------------------##
# ##--------------Harmony & DimRed---------------##
# ##---------------------------------------------##
# 
# #Harmony
# harmony_int <- harmony_int %>% RunHarmony(group.by.vars = "dataset", reduction = "pca", 
#                                           dims.use = 1:60, plot_convergence = TRUE)
# 
# # Run Umap
# harmony_int <- RunUMAP(harmony_int, reduction = "harmony", dims = 1:60)
# # Run tsne
# harmony_int <- RunTSNE(harmony_int, reduction = "harmony", dims = 1:60)
# 
# FeaturePlot(harmony_int, features = "ENSG00000115457", reduction = "umap", slot = "data")
# 
# ##---------------------------------------------##
# ##-------Add metadata clustering per D-P-------##
# ##---------------------------------------------##
# 
# harmony_int$mySelectedClustering <- clustering_metadata[colnames(harmony_int), "mySelectedClustering"]
#   
# DimPlot(harmony_int, group.by = "mySelectedClustering", reduction ="umap", cols = myColors)
# 
# dir.create("5.Integration/output/")
# saveRDS(harmony_int, "5.Integration/output/Harmony.Rds")
#####

##---------------------------------------------##
##-------------Seurat - integration------------##
##---------------------------------------------##

# Anchors & integration
rm(myseu)
anchors <- FindIntegrationAnchors(object.list = myseu_raw, dims = 1:100)
int <- IntegrateData(anchorset = anchors, dims = 1:100) 
saveRDS(anchors, "5.Integration/output/anchors_seurat_100dims.Rds")
rm(anchors)

# Pipeline 
int <- ScaleData(int, verbose = FALSE)
int <- RunPCA(int , npcs = 60, verbose = FALSE)
int <- RunUMAP(int , reduction = "pca", dims = 1:60)
int_test <- RunTSNE(int , reduction = "pca", dims = 1:30, perplexity = 15)

# Add metadata
int$batch <- factor(sapply(strsplit(colnames(int),"_"), `[`, 1), levels = myDepots)
int$mySelectedClustering <- clustering_metadata[colnames(int), "mySelectedClustering"]
int$mySelectedClustering_perSample <- clustering_metadata[colnames(int), "mySelectedClustering_perSample"]
int$mySelectedClustering_perSample_simplified <- clustering_metadata[colnames(int), "mySelectedClustering_perSample_simplified"]

setwd("~/SVRAW1/prainer/hASPCs/PAPER/10X_scRNA-seq/")
saveRDS(object = int, file = "5.Integration/output/Seurat_2000HVGs.Rds")



all_genes <- c()
for(i in names(reorganised_list_orthologs)){
   for(ii in names(reorganised_list_orthologs[[i]])){
      all_genes <- c(all_genes, reorganised_list_orthologs[[i]][[ii]]$gene.Gene.stable.ID.1)
   }
}
all_genes <- unique(all_genes)
genes_to_int <- all_genes[all_genes %in% rownames(anchors@object.list[[1]])]
integrated_feat <- IntegrateData(anchorset = anchors, dims = 1:40, features.to.integrate = genes_to_int)

