################################################################
#                                                              #
#                    Clustering and plots                      #
#                       of integration                         #
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
int <- readRDS("5.Integration/output/Seurat_2000HVGs.Rds")

##---------------------------------------------##
##-----------------Plot batches----------------##
##---------------------------------------------##

reord <- sample(ncol(int))
p <- plot.tsne.colored_discrete_2(int@reductions$tsne@cell.embeddings[reord,], 
                             data.categories = int$batch[reord],
                             size.point = 0.8, axis = T, x.lab = "tSNE1", y.lab = "tSNE2",
                             col = unname(myBatchColors[levels(int$batch)]))

dir.create("5.Integration/plots/")
ggsave(p, 
       filename = "5.Integration/plots/tSNE_integration_colored_batch.pdf")
#Saving 9.77 x 5.85 in image
ggsave(p + theme(legend.position = "none"), width = 6.36, height = 6, 
       filename = "5.Integration/plots/tSNE_integration_colored_batch_woLeg.pdf")

##---------------------------------------------##
##----------Plot individual clustering---------##
##---------------------------------------------##

myGradientColors <- c(colorRampPalette(c("#D4FF88","#41ab5d","#004529" ))(11), #ASCs
                      colorRampPalette(c("#fee0d2","#ef3b2c","#67000d" ))(11), #PreAs
                      colorRampPalette(c("#bdd7e7", "#2171b5"))(3),#IGFBP2
                      colorRampPalette(c("#bcbddc", "#6E1779"))(3), #Meso
                      colorRampPalette(c("#fee6ce", "darkorange1", "#a63603"))(11), #VSMPs 
                      "darkgoldenrod1", #Endo        
                      "#E889BD" #Immune
)
names(myGradientColors) <- levels(clustering_metadata$mySelectedClustering_perSample_simplified)

reord <- sample(ncol(int))
p <- plot.tsne.colored_discrete_2(int@reductions$umap@cell.embeddings[reord,], 
                                  data.categories = int$mySelectedClustering_perSample_simplified[reord],
                                  size.point = 0.8, axis = T, x.lab = "tSNE1", y.lab = "tSNE2",
                                  col = unname(myGradientColors))

ggsave(p, 
       filename = "5.Integration/plots/tSNE_integration_colored_mySelectedClustering_perSample_simplified.pdf")
#Saving 9.77 x 5.85 in image
ggsave(p + theme(legend.position = "none"), width = 6.36, height = 6,
       filename = "5.Integration/plots/tSNE_integration_colored_mySelectedClustering_perSample_simplified_woLeg.pdf")

##---------------------------------------------##
##------------------Clustering-----------------##
##---------------------------------------------##

int <- FindNeighbors(object = int, dims = 1:60, reduction = "pca")

## -- Clustering for res from 0.2 to 1.5
for( i in seq( 0.1, 1.5, 0.1)){
  int <- FindClusters(object = int, reduction.type = "pca", resolution = i, print.output = T)}

saveRDS(int, "5.Integration/output/Seurat_2000HVGs.Rds")

##---------------------------------------------##
##---------------Plot clustering---------------##
##---------------------------------------------##

m <- grep("integrated_snn_res", colnames(int@meta.data))
p <- lapply( m, function (i) plot.tsne.colored_discrete(int@reductions$tsne@cell.embeddings, 
                                                        int@meta.data[i], 
                                                        size.point = 1, size.text = 4, 
                                                        col = c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
                                                                '#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
                                                                '#cab2d6','#6a3d9a','#ffff99',"#FFFF33",
                                                                "#E5C494", '#b15928',"#fccde5","#ccebc5",
                                                                "black", "gray", "orange")))

dir.create("5.Integration/plots/Clustering/")
pdf( "5.Integration/plots/Clustering/tSNE_integration_Clustering_res0.1-1.5.pdf", paper = "a4r")
  marrangeGrob(p, nrow = 2, ncol = 1)
dev.off()

clustree::clustree(int, prefix = "integrated_snn_res.")

##---------------------------------------------##
##--------------Select clustering--------------##
##---------------------------------------------##
int$myIntegratedClustering <- int$integrated_snn_res.0.1
levels(int$myIntegratedClustering) <- c("PreAs", "ASCs", "Meso", "VSMPs", "4", "5", "PR specific")

int$myIntegratedClustering <- as.character(int$myIntegratedClustering)

int$myIntegratedClustering[int$integrated_snn_res.0.2 == 8] <- "res.0.2_8"
int$myIntegratedClustering[int$integrated_snn_res.0.2 == 6] <- "res.0.2_6"

int$myIntegratedClustering[int$integrated_snn_res.0.3 == 7] <- "Endo"
int$myIntegratedClustering[int$integrated_snn_res.0.3 == 8] <- "Immune"

int$myIntegratedClustering[int$myIntegratedClustering == "res.0.2_6" & int$integrated_snn_res.0.7 == 16] <- "IFIT"

int$myIntegratedClustering[int$myIntegratedClustering == "res.0.2_6" & int$integrated_snn_res.0.7 == 14] <- "IGFBP2"
int$myIntegratedClustering[int$myIntegratedClustering == "res.0.2_6"] <- "CHI3L1-2"

int$myIntegratedClustering[int$integrated_snn_res.0.8 == 13] <- "HHIP"
int$myIntegratedClustering[int$integrated_snn_res.0.8 == 7] <- "CILP"

int$myIntegratedClustering[int$myIntegratedClustering == 5] <- "PreAs"

#int$myIntegratedClustering <- as.factor(int$myIntegratedClustering)

int$myIntegratedClustering <- factor(as.character(int$myIntegratedClustering), 
                                     levels = names(myIntegratedColors))
table(int$myIntegratedClustering)

p <- plot.tsne.colored_discrete_2(int@reductions$tsne@cell.embeddings, 
                                  data.categories = int$myIntegratedClustering,
                                  size.point = 0.8, axis = T, x.lab = "tSNE1", y.lab = "tSNE2",
                                  col = unname(myIntegratedColors))

ggsave(p, height = 3.92, width = 5.48,
       filename = "5.Integration/plots/Clustering/tSNE_integration_colored_myIntegrationClustering.pdf")
#Saving 5.48 x 3.92 in image
ggsave(p + theme(legend.position = "none"), width = 6.36, height = 6,
       filename = "5.Integration/plots/Clustering/tSNE_integration_colored_myIntegrationClustering_woLeg.pdf")


##---------------------------------------------##
##----------Plot Lib Size & nFeatures----------##
##---------------------------------------------##

VlnPlot(int, features = c("nCount_RNA", "nFeature_RNA"), cols = myIntegratedColors)

df <- int@meta.data[, c("nCount_RNA", "nFeature_RNA", "myIntegratedClustering")]
df <- df %>% filter(myIntegratedClustering != "res.0.2_8")

p_counts <- ggplot(df, aes(x = myIntegratedClustering, y = nCount_RNA, fill = myIntegratedClustering)) + 
  scale_fill_manual(values = myIntegratedColors) + 
  ggrastr::rasterise(geom_jitter(alpha = 0.04)) + 
  geom_violin(scale = "width") + 
  mashaGgplot2Theme + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 60, hjust=1))

p_Feat <- ggplot(df, aes(x = myIntegratedClustering, y = nFeature_RNA, fill = myIntegratedClustering)) + 
  scale_fill_manual(values = myIntegratedColors) + 
  ggrastr::rasterise(geom_jitter(alpha = 0.04)) + 
  geom_violin(scale = "width") + 
  mashaGgplot2Theme + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 60,  hjust=1))

ggsave(grid.arrange(p_counts, p_Feat, ncol = 2), 
       filename = "5.Integration/plots/Violin_nCounts-nFeat.pdf",
       width = 12, height = 4.9)


