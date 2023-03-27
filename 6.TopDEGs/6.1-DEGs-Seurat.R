################################################################
#                                                              #
#                  Investigate top DE markers                  #
#                                                              #
################################################################

### Author: Pernille
### Date: 17.08.2022 - adapted from old script
### Datasets: scRNA-seq Depots: SC - EP - MK - MG - PR - GB 
###                     Patients: B0, B1, B7, L3, L11, L12 (B - Bariatric, L - Lean)
### Goal: Compare markers between human and mouse ASPCs 

library(ggplot2); library(data.table); library(Seurat); library(dplyr)

setwd("~/SVRAW1/prainer/hASPCs/PAPER/10X_scRNA-seq/")
dir.create("6.TopDEGs/")

source("Utility/General_utils.R")
source("6.TopDEGs/6-utils.R")

##---------------------------------------------##
##-----------------Loading data----------------##
##---------------------------------------------##

int <- readRDS("5.Integration/output/Seurat_2000HVGs.Rds")

##---------------------------------------------##
##-------------------Set data------------------##
##---------------------------------------------##

DefaultAssay(int) <- "RNA"
int <- SetIdent(int, value = int$myIntegratedClustering)
myclust <- levels(int$myIntegratedClustering)

##---------------------------------------------##
##---------------------DEGs--------------------##
##---------------------------------------------##

sink("6.TopDEGs/console_output.txt", type = c("output", "message"))
myDEGs_seurat <- lapply(myclust, findMarkersOfClusts)
sink()

names(myDEGs_seurat) <- myclust

myDEGs_seurat <- lapply(myclust, add_avgFC)
names(myDEGs_seurat) <- myclust

##---------------------------------------------##
##---------------------DEGs--------------------##
##---------------------------------------------##
DefaultAssay(int) <- "RNA"
IGFBP2_vs_ASPCs <- FindConservedMarkers(int, ident.1 = "IGFBP2", 
                                        ident.2 = c("ASCs", "PreAs", "CHI3L1-2", "IFIT", "HHIP"), 
                                        grouping.var = "batch", verbose = T, min.cells.group = 10)
IGFBP2_vs_ASPCs$geneID <- data.annot[rownames(IGFBP2_vs_ASPCs), "gene_short_name"]
IGFBP2_vs_Meso <- FindConservedMarkers(int, ident.1 = "IGFBP2", ident.2 = "Meso", 
                                       grouping.var = "batch",
                                       verbose = T, min.cells.group = 10)
IGFBP2_vs_Meso$geneID <- data.annot[rownames(IGFBP2_vs_Meso), "gene_short_name"]

save(IGFBP2_vs_ASPCs, IGFBP2_vs_Meso, file = "6.TopDEGs/DEGs_IGFBP2-vs-Meso-or-ASPCs.RData")

##---------------------------------------------##
##---------------------save--------------------##
##---------------------------------------------##

saveRDS(myDEGs_seurat, "6.TopDEGs/DEGs.Rds")
