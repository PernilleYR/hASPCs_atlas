################################################################
#                                                              #
#               Seurat objects & Clustering                    #
#                                                              #
################################################################


### Author: Pernille
### Date: 10.08.2022
### Datasets: scRNA-seq Depots: SC - EP - MK - MG - PR - GB 
###                     Patients: B0, B1, B7, L3, L11, L12 (B - Bariatric, L - Lean)
### Goal: Plot tSNEs

library(Seurat); library(ggrastr)

setwd("~/SVRAW1/prainer/hASPCs/PAPER/10X_scRNA-seq/")

##---------------------------------------------##
##----------------Loading data-----------------##
##---------------------------------------------##

myseu <- readRDS("0.data/List_seurat_objects.rds")

##---------------------------------------------##
##-----------Colors & plot function------------##
##---------------------------------------------##

myColors <- c("ASCs" = "#33A02C", "PreAs" = "#E31A1C", "IGFBP2" = "#1f78b4", "Meso" = "#810F7C", 
              "VSMPs" = "#FC8D62", "Endo" = "darkgoldenrod1", "Immune" = "#E889BD", "Unknown_VSMPs" = "gray")

plot_tsne <- function(data, name.category.plot = "", 
                      axis = T, x.lab = "", y.lab = "", labels.category = NULL, 
                      size.point = 0.5, size.text = 3, 
                      col = myColors,
                      legend = "none", leg.point.size = 3){
  
  data.plot = data.frame(data@reductions$tsne@cell.embeddings[,1:2])
  colnames(data.plot) = c("x", "y")
  data.plot$category = data$mySelectedClustering[rownames(data.plot)]
  data.plot$cell = rownames(data.plot)
  col <- myColors[levels(data$mySelectedClustering)]
  p <- ggplot(data.plot, aes(x=x, y=y, col=category)) + 
    rasterise(geom_point(shape=19, size = size.point), dpi = 650) +
    scale_color_manual(values = unname(col), name = NULL) +
    xlab(label = "tsne1") + ylab(label = "tsne2") +
    theme(axis.line = element_line(colour = "black"),legend.position = legend, panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), axis.text=element_blank(), 
          axis.ticks=element_blank(), plot.background=element_blank(), panel.background=element_blank()) +
    theme(legend.key=element_blank()) + 
    guides(colour = guide_legend(override.aes = list(size=3)))
  if(axis == F){ p <- p + theme(axis.line = element_blank())+ xlab(label = "") + ylab(label = "")}
  
  return(p)
}

##---------------------------------------------##
##-----------Colors & plot function------------##
##---------------------------------------------##

p <- lapply(myseu, plot_tsne)
names(p) <- names(myseu)

ggsave(marrangeGrob(list(p$SC0, p$SC1, p$SC7, p$EP0, p$EP1, p$EP7,p$PR11, p$PR3, p$PR12, p$MG7, p$MK7, p$MG7), 
                    ncol = 3, nrow = 4,layout_matrix = matrix(1:12,  nrow = 4, ncol=3, byrow=TRUE)), 
      width = 7.2, height = 8.8, filename = "1.tSNEs_perD-P/tSNEs_blue_woLeg.pdf")
