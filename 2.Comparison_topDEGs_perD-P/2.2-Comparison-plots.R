################################################################
#                                                              #
#          Comparison of shared markers between pops           #
#                          for Bar017                          #
#                                                              #
################################################################

### Author: Pernille
### Date: 10.08.2022
### Datasets: scRNA-seq Depots: SC - EP - MK - MG - PR - GB 
###                     Patients: B0, B1, B7, L3, L11, L12 (B - Bariatric, L - Lean)
### Goal: Compare percent of shared markers between sub pops across Depots and Patients

setwd("~/SVRAW1/prainer/hASPCs/PAPER/10X_scRNA-seq/2.Comparison_topDEGs_perD-P/")
source("utils.R")

##---------------------------------------------##
##-----------------Loading data----------------##
##---------------------------------------------##
myDepots <- c(paste0("EP", c(0,1,7)),
              paste0("SC", c(0,1,7)),
              "MG7", "MK7",
              paste0("PR", c(3,11,12)))

topMarkers <- list()
for(d in myDepots){
  topMarkers[[d]] <- readRDS(paste0("DEGs/", d, "/Markers_", d, ".rds"))
}

myColors <- c("ASCs" = "#33A02C", "PreAs" = "#E31A1C", "IGFBP2" = "#284724", "Meso" = "#810F7C", 
              "VSMPs" = "#FC8D62", "Endo" = "darkgoldenrod1", "Immune" = "#E889BD", "Unknown_VSMPs" = "gray")

##---------------------------------------------##
##----------------Create Plots-----------------##
##---------------------------------------------##

## -- Print the number of markers per sample and cell types:
for(i in names(topMarkers)){
  cat(i,":", "\n")
  for(ii in names(topMarkers[[i]])){
    cat(paste0(ii,": ",nrow(topMarkers[[i]][[ii]]), ", "))
  }
  cat("\n")
}
# EP0 : 
#   Meso: 428, ASCs: 461, IGFBP2: 206, PreAs: 335, VSMPs: 310, 
# EP1 : 
#   PreAs: 478, Meso: 499, ASCs: 433, VSMPs: 504, IGFBP2: 44, 
# EP7 : 
#   ASCs: 525, Endo: 1367, IGFBP2: 63, Immune: 1008, Meso: 704, PreAs: 506, VSMPs: 275, 
# SC0 : 
#   PreAs: 198, ASCs: 229, VSMPs: 402, 
# SC1 : 
#   PreAs: 287, ASCs: 385, VSMPs: 454, 
# SC7 : 
#   PreAs: 339, ASCs: 511, VSMPs: 860, Endo: 651, Immune: 547, 
# MG7 : 
#   SCs: 522, Endo: 1229, Immune: 737, PreAs: 391, VSMPs: 637, 
# MK7 : 
#   PreAs: 340, ASCs: 393, Endo: 895, VSMPs: 297, Immune: 173, 
# PR3 : 
#   ASCs: 589, PreAs: 348, VSMPs: 192, 
# PR11 : 
#   PreAs: 409, ASCs: 693, VSMPs: 895, Unknown_VSMPs: 176, 
# PR12 : 
#   ASCs: 371, PreAs: 319, VSMPs: 648, 

## -- Detect markers assigned to more than one pop
doubleAssigned <- list()
for(i in names(topMarkers)){
  doubleAssigned[[i]] <- detect_doubleAssigned_markers(topMarkers[[i]])
}

## -- Calculate the % of shared genes
percent_markers <- calculate_percent_shared_markers(topMarkers, n_top_genes = 100 )
# EP1 IGFBP2 had only 44
# EP7 IGFBP2 had only 63

percent_markers_common <- calculate_percent_shared_markers_common_v2(topMarkers, n_top_genes = 100,  doubleAssigned)
# EP1 IGFBP2 had only 44
# EP7 IGFBP2 had only 63

## -- Create the barplots 
p <- list()
for(i in names(percent_markers)){
  n <- strsplit(i, "-")
  p[[i]] <- barplot_shared_markers(percent_markers[i], myColors, 
                                   shift_text_yaxis = 1.5, min_value = 7, n_top_genes = 100 ) + 
    theme_bw()
  
}

p2 <- list()
cols <- c(myColors, "common" = "gray")
for(i in names(percent_markers_common)){
  n <- strsplit(i, "-")
  p2[[i]] <- barplot_shared_markers(percent_markers_common[i], cols, shift_text_yaxis = 2.5, min_value = 7, n_top_genes = 100 ) +
    theme_bw()
}

##---------------------------------------------##
##-----------------SAVE PLOTS------------------##
##---------------------------------------------##
dir.create("Barplot_perComparison")
dir.create("Barplot_perComparison/CommonNotHighlighted")
dir.create("Barplot_perComparison/Common")

for(i in names(p)){
  ggsave(p[[i]], width = 8.35, height = 5.58,
         file = paste0("Barplot_perComparison/CommonNotHighlighted/", "Barplot_", i, "_top100Markers.pdf"))
}
for(i in names(p2)){
  ggsave(p2[[i]], width = 8.35, height = 5.58,
         file = paste0("Barplot_perComparison/Common/Barplot_", i, "_top100Markers.pdf"))
}

saveRDS(percent_markers_common, file = "Barplot_perComparison/Common/percent_markers_common.Rds")
saveRDS(percent_markers, file = "Barplot_perComparison/CommonNotHighlighted/percent_markers.Rds")
saveRDS(doubleAssigned, file = "Barplot_perComparison/Common/Markers_Assigned_to_more_than_onePopulation.Rds")
