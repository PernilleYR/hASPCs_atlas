################################################################
#                                                              #
#               Seurat objects & Clustering                    #
#                                                              #
################################################################


### Author: Pernille
### Date: 10.08.2022
### Datasets: scRNA-seq Depots: SC - EP - MK - MG - PR - GB 
###                     Patients: B0, B1, B7, L3, L11, L12 (B - Bariatric, L - Lean)
### Goal: Collect Seurat object and select final clustering

library(Seurat)

setwd("~/SVRAW1/prainer/hASPCs/")

##---------------------------------------------##
##----------------Loading data-----------------##
##---------------------------------------------##

EP0 <- readRDS("10X/10X_Bar0_141118/analysis/EP0/RScriptData/EP0_Seurat_2.Rds")
EP1 <- readRDS("10X/10X_Bar1_201118/analysis/EP1/RScriptData/EP1_Seurat_2.Rds")
EP7 <- readRDS("10X/10X_Bar7_EP_190613/Analysis_noDBL_mergedASPCsresults/RScriptData/EP7_noDBL_seurat.Rds")

SC0 <- readRDS("10X/10X_Bar0_141118/analysis/SC0/RScriptData/SC0_Seurat_2.Rds")
SC1 <- readRDS("10X/10X_Bar1_201118/analysis/SC1/RScriptData/SC1_Seurat_2.Rds")
SC7 <- readRDS("10X/10X_Bar7_SC_190613/Analysis_strictFiltBarcodes/RScriptData/SC7_seurat.Rds")

MG7 <- readRDS("10X/10X_Bar7_MG_190613/Analysis_strictFiltBarcodes/RScriptData/MG7_seurat.Rds")
MK7 <- readRDS("10X/10X_Bar7_MK_190613/Analysis_strictFiltBarcodes/RScriptData/MK7_sce_seurat.Rds")

GB7 <- readRDS("10X/10X_Bar7_GB_190613/Analysis_ASPCs/RScriptData/GB7_ASPCs_seurat.Rds")

PR3 <- readRDS("10X/PR_LP12_03_frozen/9.Identify_Cell_Patients/output/LP03_seurat.Rds")
PR11 <- readRDS("10X/PR_LP11/RScriptData/PR_Seurat.Rds")
  PR11 <- RenameCells(PR11, new.names = sapply(strsplit(colnames(PR11),"_"), `[`, 3))
PR12 <- readRDS("10X/PR_LP12_03_frozen/9.Identify_Cell_Patients/output/LP12_seurat.Rds")

##---------------------------------------------##
##-----------------Clustering------------------##
##---------------------------------------------##

EP0$mySelectedClustering <- as.factor(EP0$RNA_snn_res.0.3)
  levels(EP0$mySelectedClustering) <- c("Meso", "ASCs", "IGFBP2", "PreAs", "VSMPs")
EP1$mySelectedClustering <- EP1$RNA_snn_res.0.3
  levels(EP1$mySelectedClustering) <- c("PreAs", "Meso", "ASCs", "VSMPs", "IGFBP2")
EP7$mySelectedClustering <- EP7$RNA_snn_res.0.2
  levels(EP7$mySelectedClustering) <- c("PreAs", "Meso", "ASCs", "-", "Immune", "VSMPs")
  EP7$mySelectedClustering <- as.character(EP7$mySelectedClustering)
  EP7$mySelectedClustering[EP7$RNA_snn_res.0.9 == 7] <- "IGFBP2"
  EP7$mySelectedClustering[EP7$RNA_snn_res.0.4 == 8] <- "VSMPs"
  EP7$mySelectedClustering[EP7$RNA_snn_res.0.4 %in% c(4,7)] <- "Endo"
  EP7$mySelectedClustering <- as.factor(EP7$mySelectedClustering)

SC0$mySelectedClustering <- SC0$RNA_snn_res.0.3
  levels(SC0$mySelectedClustering) <- c("PreAs", "ASCs", "VSMPs")
SC1$mySelectedClustering <- SC1$RNA_snn_res.0.1
  levels(SC1$mySelectedClustering) <- c("PreAs", "ASCs", "VSMPs")
SC7$mySelectedClustering <- SC7$RNA_snn_res.0.1
  levels(SC7$mySelectedClustering) <- c("PreAs", "ASCs", "VSMPs", "Endo", "Immune")

GB7$mySelectedClustering <- GB7$RNA_snn_res.0.5
  levels(GB7$mySelectedClustering) <- c("ASCs", "PreAs")
  
MG7$mySelectedClustering <- MG7$RNA_snn_res.0.1
  levels(MG7$mySelectedClustering) <- c("PreAs", "ASCs", "-", "Immune") 
  MG7$mySelectedClustering <- as.character(MG7$mySelectedClustering)
  MG7$mySelectedClustering[MG7$RNA_snn_res.0.2 == 5] = "VSMPs"
  MG7$mySelectedClustering[MG7$RNA_snn_res.0.2 == 3] = "Endo"
  MG7$mySelectedClustering <- as.factor(MG7$mySelectedClustering)
MK7$mySelectedClustering <- MK7$RNA_snn_res.0.3
  levels(MK7$mySelectedClustering) <- c("0","0", 2:5) #For simplicity merge two clusters of preAs
  levels(MK7$mySelectedClustering) <- c("PreAs", "ASCs", "Endo", "VSMPs", "Immune")

PR3$mySelectedClustering <- PR3$RENAMED_RNA_snn_res.0.2
levels(PR3$mySelectedClustering) <- c("ASCs", "PreAs", "VSMPs")
PR11$mySelectedClustering <- PR11$RNA_snn_res.0.1
  levels(PR11$mySelectedClustering) <- c("PreAs", "ASCs", "VSMPs", "Unknown_VSMPs")
PR12$mySelectedClustering <- PR12$RENAMED_RNA_snn_res.0.2
  levels(PR12$mySelectedClustering) <- c("ASCs", "PreAs", "VSMPs")
  
##---------------------------------------------##
##--------------------tSNE---------------------##
##---------------------------------------------##

EP0 <- RunTSNE(EP0, dims = 1:29)
EP1 <- RunTSNE(EP1, dims = 1:31)
EP7 <- RunTSNE(EP7, dims = 1:32)

SC0 <- RunTSNE(SC0, dims = 1:33)
SC1 <- RunTSNE(SC1, dims = 1:22)
SC7 <- RunTSNE(SC7, dims = 1:36)

MG7 <- RunTSNE(MG7, dims = 1:45)
MK7 <- RunTSNE(MK7, dims = 1:39)

PR3 <- RunTSNE(PR3, dims = 1:15)
PR11 <- RunTSNE(PR11, dims = 1:45)
PR12 <- RunTSNE(PR12, dims = 1:16)

##---------------------------------------------##
##--------------------SAVE---------------------##
##---------------------------------------------##

myseu <- list("EP0"=EP0, "EP1"=EP1,"EP7"=EP7, "SC0"=SC0, "SC1"=SC1, "SC7"=SC7, "MG7"=MG7, "MK7"=MK7, "PR3"=PR3, "PR11"=PR11, "PR12"=PR12, "GB7" = GB7)

saveRDS(myseu, file = "PAPER/10X_scRNA-seq/0.data/Seurats_objects.rds")

