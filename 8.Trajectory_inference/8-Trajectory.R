

library(dyno)
library(tidyverse)
library(DoubletFinder)
setwd("~/SVRAW1/prainer/hASPCs/PAPER/10X_scRNA-seq/")

##---------------------------------------------##
##-----------------Loading data----------------##
##---------------------------------------------##

myseu <- readRDS("0.data/List_seurat_objects.rds")

int <- readRDS("5.Integration/output/Seurat_2000HVGs.Rds")

##---------------------------------------------##
##------Detect Doublets in individual EPs------##
##---------------------------------------------##

# EP0
bcmvn <- pKidentification(myseu$EP0) 
pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_low <- 0.16
nExps <- homotypicDblPropEstimate(myseu$EP0, clust = "RNA_snn_res.0.2")
myseu$EP0 <- RunDblFinder(myseu$EP0, pK.selected = c(pK_low, pK$pK), nExps.poi = nExps)
myseu$EP0 <- prepare(myseu$EP0)
DimPlot(myseu$EP0, reduction = "tsne", group.by ="DF_hi.lo_0.07", cols = c('#e41a1c','#377eb8', "gray")) #SELECTED
DimPlot(myseu$EP0, reduction = "tsne", group.by ="DF_hi.lo_0.16", cols = c('#e41a1c','#377eb8', "gray"))
DimPlot(myseu$EP0, reduction = "tsne", group.by ="DF_hi.lo_0.27", cols = c('#e41a1c','#377eb8', "gray"))

# EP1
bcmvn <- pKidentification(myseu$EP1) 
pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
nExps <- homotypicDblPropEstimate(myseu$EP1, clust = "mySelectedClustering")
myseu$EP1 <- RunDblFinder(myseu$EP1, pK.selected = c(0.07, 0.1, pK$pK), nExps.poi = nExps)
myseu$EP1 <- prepare(myseu$EP1)
DimPlot(myseu$EP1, reduction = "tsne", group.by ="DF_hi.lo_0.07", cols = c('#e41a1c','#377eb8', "gray"))
DimPlot(myseu$EP1, reduction = "tsne", group.by ="DF_hi.lo_0.1", cols = c('#e41a1c','#377eb8', "gray")) #SELECTED
DimPlot(myseu$EP1, reduction = "tsne", group.by ="DF_hi.lo_0.29", cols = c('#e41a1c','#377eb8', "gray"))

# EP7
bcmvn <- pKidentification(myseu$EP7) 
pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
nExps <- homotypicDblPropEstimate(myseu$EP7, clust = "mySelectedClustering")
myseu$EP7 <- RunDblFinder(myseu$EP7, pK.selected = c(pK$pK), nExps.poi = nExps)
myseu$EP7 <- prepare(myseu$EP7)
DimPlot(myseu$EP7, reduction = "tsne", group.by ="DF_hi.lo_0.07") #SELECTED

##---------------------------------------------##
##-----------Select cells of interest----------##
##---------------------------------------------##
int$DblFinderEPs <- NA
int$DblFinderEPs[paste0("EP7_",colnames(myseu$EP7))] <- myseu$EP7$DF_hi.lo_0.07
int$DblFinderEPs[paste0("EP1_",colnames(myseu$EP1))] <- myseu$EP1$DF_hi.lo_0.1
int$DblFinderEPs[paste0("EP0_",colnames(myseu$EP0))] <- myseu$EP0$DF_hi.lo_0.07

DimPlot(int, reduction = "tsne", group.by ="DblFinderEPs", cols = c('#e41a1c','#377eb8', "gray"), na.value = "gray95")

#Makes sense
cells_oi <- int@meta.data %>% filter(depot == "EP" & DblFinderEPs == "Singlet" &
                                       myIntegratedClustering %in% c("ASCs", "PreAs", "IGFBP2", "Meso", "VSMPs")) 

##---------------------------------------------##
##------------------Save data------------------##
##---------------------------------------------##

# Int 
for_PAGA <- subset(int, cells = rownames(cells_oi))
DefaultAssay(for_PAGA) <- "integrated"
write.csv(t(GetAssayData(for_PAGA, slot = "data")), 
            file = "8.Trajectory_inference/On-integration/data/Int-EPs_normalizedData_batchCorrected-Seurat_PreAs-ASCs-IGFBP2-Meso-VSMPs_Singlets.csv")
write.csv(for_PAGA@meta.data[,c("myIntegratedClustering", "batch")], 
            file = "8.Trajectory_inference/On-integration/data/Int-EPs_meta_PreA-ASCs-IGFBP2-Meso-VSMPs_Singlets.csv")

# EP0

# EP1

# EP7
# myseu$EP7$myIntegratedClustering <- int$myIntegratedClustering[paste0("EP7_", colnames(myseu$EP7))]
# cells_oi_EP7 <- cells_oi[grep("EP7", rownames(cells_oi)),]
# cells_oi_EP7 <- sapply(strsplit(rownames(cells_oi_EP7),"_"), `[`, 2)
# for_PAGA <- subset(myseu$EP7, cells = cells_oi_EP7)
# write.csv(t(GetAssayData(for_PAGA, slot = "data")), 
#           file = "8.Trajectory_inference/On-each-dataset/data/EP7_normalized_ASCs-PreAs-IGFBP2-Meso-VSMPs_Singlets.csv")
# write.csv(for_PAGA@meta.data[,c("myIntegratedClustering"), drop = F], 
#           file = "8.Trajectory_inference/On-each-dataset/data/EP7_meta_ASCs-PreAs-IGFBP2-Meso-VSMPs_Singlets.csv")

##---------------------------------------------##
##------------------PYTHON---------------------##
##---------------------------------------------##


##---------------------------------------------##
##---Create dynverse data & import trajectory--##
##---------------------------------------------##

data_int <- subset(x = int, cells = rownames(cells_oi))
DefaultAssay(data_int) <- "integrated"
dataset_int <- wrap_expression(
  counts = t(GetAssayData(data_int, assay = "RNA", slot = "counts")[rownames(data_int),]),
  expression = t(GetAssayData(data_int, slot = "data")), 
  CellType = data_int$myIntegratedClustering
)


##---------------------------------------------##
##---Identify genes changing along pseudotime--##
##---------------------------------------------##

### Identify genes 
# traj <- readRDS("8.Trajectory_inference/On-integration/int-EPs_PreAs-ASCs-Meso-VSMPs_Singlets_PAGA-model.rds")
# branch_feature_importance <- calculate_branch_feature_importance(test, expression_source = dataset_int)
# branch_feature_importance <- as.data.frame(branch_feature_importance)
# branch_feature_importance$feature_id <- as.character(branch_feature_importance$feature_id)
# branch_feature_importance$gene_id <- data.annot[branch_feature_importance$feature_id, "gene_short_name"]
# saveRDS(branch_feature_importance, file = "8.Trajectory_inference/On-integration/branch_feat_importance.Rds")


branch_feature_importance <- readRDS("8.Trajectory_inference/On-integration/branch_feat_importance.Rds")
# 1 ASCs, 4 PreAs, 2 IGFBP2, 3 Begging Meso, 6 Meso

features_toMeso <- branch_feature_importance %>% 
  filter(to == 6) %>% arrange(order(importance, decreasing = T))
features_toEarlyMeso <- branch_feature_importance %>% 
  filter(to == 3)  %>% arrange(order(importance, decreasing = T))
features_toIGFBP2 <- branch_feature_importance %>% 
  filter(to == 2)  %>% arrange(order(importance, decreasing = T))

##---------------------------------------------##
##--------Heatmap genes along pseudotime-------##
##---------------------------------------------##

## -- Load pseudotime
pseudotime <- read.csv("8.Trajectory_inference/On-integration/output_Python/pseudotime.csv", row.names = 1)
pseudotime <- pseudotime %>% arrange(dpt_pseudotime)
pseudotime$clust <- data_int$myIntegratedClustering[rownames(pseudotime)]

## -- Annot of heatmap
annot_col <- data.frame(row.names = rownames(pseudotime),
                        pop = as.factor(as.character(data_int$myIntegratedClustering[rownames(pseudotime)])),
                        pseudotime = pseudotime$dpt_pseudotime)
my_colour = list(
  pop = c("ASCs" = "#33A02C", "PreAs" = "#E31A1C", "IGFBP2" = "#1f78b4", "Meso" = "#810F7C", "VSMPs" =  "#FC8D62"),
  pseudotime = c('#fff7f3','#fde0dd','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a')
)

## -- Gene selection
#####
# d_heat <- GetAssayData(d, slot = "data")[features_toIGFBP2$feature_id[1:100], rownames(pseudotime)]
# rownames(d_heat) <- data.annot[rownames(d_heat), "gene_short_name"]
# d_heat <- dynutils::scale_quantile(t(d_heat))
# 
# p <- pheatmap::pheatmap(t(d_heat), 
#                         annotation_col = annot_col[rownames(d_heat),],
#                         clustering_method = "ward.D2",
#                         cluster_rows = T, cluster_cols = F,
#                         show_colnames = F, annotation_colors = my_colour,
#                         border_color = "NA")
# 
# g_oi <- p$tree_row$labels[p$tree_row$order][c(1:49,68:100)]
# p <- pheatmap::pheatmap(t(d_heat)[g_oi,], 
#                         annotation_col = annot_col[rownames(d_heat),],
#                         clustering_method = "ward.D2",
#                         cluster_rows = T, cluster_cols = F,
#                         show_colnames = F, annotation_colors = my_colour,
#                         border_color = "NA")
# g_oi <- p$tree_row$labels[p$tree_row$order][c(1:40, 42:43, 47:48, 66:79, 82)]
# p <- pheatmap::pheatmap(t(d_heat)[g_oi,], 
#                         annotation_col = annot_col[rownames(d_heat),],
#                         clustering_method = "ward.D2",
#                         cluster_rows = T, cluster_cols = F,
#                         show_colnames = F, annotation_colors = my_colour,
#                         border_color = "NA")
# g_oi <- g_oi[!g_oi %in% c("GAS6", "CXCL14", "FST", "EDNRB")]
# 
# 
# g <- features_toMeso[1:100,]
# g <- g[!g$gene_id %in% features_toIGFBP2$gene_id[1:100],"feature_id"]
# d_heat <- GetAssayData(d, slot = "data")[g, rownames(pseudotime)]
# rownames(d_heat) <- data.annot[rownames(d_heat), "gene_short_name"]
# d_heat <- dynutils::scale_quantile(t(d_heat))
# 
# p <- pheatmap::pheatmap(t(d_heat), 
#                         annotation_col = annot_col[rownames(d_heat),],
#                         clustering_method = "ward.D2",
#                         cluster_rows = T, cluster_cols = F,
#                         show_colnames = F, annotation_colors = my_colour,
#                         border_color = "NA")
# g <- p$tree_row$labels[p$tree_row$order][c(1:7,24:25,29:31,45)]
# 
# g_oi <- c(g_oi, g)
# g_oi <- g_oi[g_oi != "C7"]
# 
# g <- features_toEarlyMeso[1:100,]
# g <- g[!g$gene_id %in% features_toIGFBP2$gene_id[1:100],]
# g <- g[!g$gene_id %in% features_toMeso$gene_id[1:100],"feature_id"]
# 
# d_heat <- GetAssayData(d, slot = "data")[g, rownames(pseudotime)]
# rownames(d_heat) <- data.annot[rownames(d_heat), "gene_short_name"]
# d_heat <- dynutils::scale_quantile(t(d_heat))
# 
# p <- pheatmap::pheatmap(t(d_heat), 
#                         annotation_col = annot_col[rownames(d_heat),],
#                         clustering_method = "ward.D2",
#                         cluster_rows = T, cluster_cols = F,
#                         show_colnames = F, annotation_colors = my_colour,
#                         border_color = "NA")
# g <-  p$tree_row$labels[p$tree_row$order][66:90]
# p <- pheatmap::pheatmap(t(d_heat)[g,], 
#                         annotation_col = annot_col[rownames(d_heat),],
#                         clustering_method = "ward.D2",
#                         cluster_rows = T, cluster_cols = F,
#                         show_colnames = F, annotation_colors = my_colour,
#                         border_color = "NA")
# 
#####
write.table(g_oi, "8.Trajectory_inference/On-integration/output_forHeatmap/genes_forheatmap.txt", sep="\t")
g_oi <- read.table("8.Trajectory_inference/On-integration/output_forHeatmap/genes_forheatmap.txt", sep = "\t")
# g <- c("SPARCL1", "COL3A1", "APOE", "FGL2", "IGFBP2", "RBP1", "APOC1", "KIF21A", "G0S2", "APOA1", "PENK",
#        "WNT4", "WNT6", "SNAI2", "PLCXD3", "PTGFR", "PHLDA1", "DDIT4", "C7", "COL1A1","CRABP2", "PDMIL3", 
#        "TMEM100", "PLAT", "CYP1B1", "MDK", "NBL1", "LINC01697", "NPW", "EMILIN1")
g_specific <- c("SPARCL1", "WNT4", "PLCXD3", "RBP1", "IGFBP2", "G0S2", "COL1A1", "PLAT","TMEM100", "APOE", "APOC1", "APOA1")
g_oi_specific<- data.annot[data.annot$gene_short_name %in% g_specific, ]

d_heat <- GetAssayData(data_int, slot = "data")[, rownames(pseudotime)]
d_heat <- d_heat[rownames(d_heat) %in% g_oi$ens_id,]
rownames(d_heat) <- data.annot[rownames(d_heat), "gene_short_name"]
d_heat <- dynutils::scale_quantile(t(d_heat))

## -- pheatmap plot
big_pheat <- pheatmap::pheatmap(t(d_heat), 
                   annotation_col = annot_col,
                   clustering_method = "ward.D2",
                   cluster_rows = T, cluster_cols = F,
                   show_colnames = F, annotation_colors = my_colour,
                   border_color = "NA")
ggsave(big_pheat, filename = "8.Trajectory_inference/On-integration/figures/heatmap.pdf",
       height = 7.25, width = 11.34)
ggsave(big_pheat, filename = "8.Trajectory_inference/On-integration/figures/heatmap.png",
       height = 7.25, width = 11.34)

saveRDS(big_pheat$tree_row, file = "8.Trajectory_inference/On-integration/output_forHeatmap/dendrogram.Rds")

## -- ggplot heatmap
g_ord <- big_pheat$tree_row$labels[big_pheat$tree_row$order]

df <- as.data.frame(t(d_heat)[rev(g_ord),])
df$gene_id <- rownames(df)
df_reshaped <- reshape2::melt(df)
df_reshaped$gene_id <- factor(df_reshaped$gene_id,
                              levels = df$gene_id)
df_reshaped$variable <- factor(df_reshaped$variable, levels = colnames(df))

  # heatmap
p_heat <- ggplot(df_reshaped, aes(x = variable, y = gene_id, fill = value)) + 
  ggrastr::rasterise(geom_tile(), dpi = 650) + xlab("") + theme(axis.text.x = element_blank()) + 
  scale_fill_distiller(palette="RdYlBu", type = "div")

  # bar pseudotime
bar_dpt <- ggplot(annot_col, aes(x = 1:nrow(annot_col), y = 1, fill = pseudotime)) + 
  geom_bar(stat="identity") + theme_void() + theme(legend.position = "top") + 
  scale_fill_gradientn(colors = rev(c('#a50026','#d73027','#f46d43','#fdae61',
                                      '#fee08b','#ffffbf','#d9ef8b','#a6d96a',
                                      '#66bd63','#1a9850','#006837')))
  # bar cell type 
bar_ct <- ggplot(annot_col, aes(x = 1:nrow(annot_col), y = 1, fill = pop)) + 
  geom_bar(stat = "identity") + theme_void() + theme(legend.position = "bottom") +
  scale_fill_manual(values = myIntegratedColors)

  # Dendrogram
p_dendro <- plot(big_pheat$tree_row, hang = -1)
  
ggsave(p_heat, filename = "8.Trajectory_inference/On-integration/figures/Heatmap_ggplot.pdf",
       width = 12.46, height = 9.84)
# ggsave(p_heat, filename = "8.Trajectory_inference/On-integration/figures/Heatmap_ggplot_IGFBP2-specific.pdf",
#        width = 12.46, height = 1.47)
ggsave(gridExtra::grid.arrange(bar_dpt, bar_ct, ncol = 1),
       filename = "8.Trajectory_inference/On-integration/figures/Heatmap_barLegend.pdf",
       width = 12.46 , height = 2.4 )
#p_dendro, "8.Trajectory_inference/On-integration/figures/Heatmap_gene_dendro.pdf"
#       width = 12.46, height = 


## ZOOM 

d_heat <- GetAssayData(data_int, slot = "data")[c(g_oi$ens_id, g_oi_specific$ens_id), rownames(pseudotime)]
rownames(d_heat) <- data.annot[rownames(d_heat), "gene_short_name"]
d_heat <- dynutils::scale_quantile(t(d_heat))

cat <- cutree(big_pheat$tree_row, k = 2)
dat <- as.data.frame(t(d_heat)) 

d <- data.frame( cells = colnames(dat),
                 d_up = colMeans(dat[names(cat[cat == 1]),]),
                 d_down = colMeans(dat[names(cat[cat == 2]),]),
                 d_specific = colMeans(dat[g_oi_specific$gene_short_name,]))
d$cell_type <- int$myIntegratedClustering[rownames(d)]
which(d$cell_type == "IGFBP2")

to_select <- 4250:4700
d_filt <- d[to_select, ]
df <- d_filt[, c("cells", "d_up", "d_down")]
df <- reshape2::melt(df)
df$cells <- factor(df$cells, levels = rownames(d_filt))
p_updown <- ggplot(df, aes(x = cells, y = value, col = variable)) + 
  geom_point(aes(fill = variable), alpha = 0.7, shape = 16, size = 1.5) +
  scale_color_manual(values = c("#E31A1C",  "#810F7C")) + 
  scale_fill_manual(values = c( "#E31A1C", "#810F7C")) + 
  geom_smooth(aes(group = variable)) +
  mashaGgplot2Theme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  xlab("pseudotime") + ylab("normalized mean expr.")

df <- d_filt[, c("cells", "d_specific")]
df <- reshape2::melt(df)
df$cells <- factor(df$cells, levels = rownames(d_filt))
p_specific <- ggplot(df, aes(x = cells, y = value, col = variable)) + 
  geom_point(aes(fill = variable), alpha = 0.7, shape = 16, size = 1.5) +
  scale_color_manual(values = c( "#1f78b4")) + 
  scale_fill_manual(values = c( "#1f78b4")) + 
  geom_smooth(aes(group = variable), method = "gam") +
  mashaGgplot2Theme +
  scale_y_continuous(breaks=c(0.2,0.4,0.6,0.8)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  xlab("pseudotime") + ylab("normalized mean expr.")

ggsave(gridExtra::grid.arrange(p_updown, p_specific, ncol = 1), 
       filename = "8.Trajectory_inference/On-integration/figures/Scatter-ZOOM_UpDownSpecific.pdf",
       width = 5, height = 8.5)

bar_dpt <- ggplot(annot_col[to_select,], aes(x = 1:length(to_select), y = 1, fill = pseudotime)) + 
  ggrastr::rasterise(geom_bar(stat="identity")) + theme_void() + theme(legend.position = "top") + 
  scale_fill_gradientn(colors = rev(c('#d9ef8b','#a6d96a',
                                      '#66bd63','#1a9850','#006837')))
# bar cell type 
bar_ct <- ggplot(annot_col[to_select,], aes(x = 1:length(to_select), y = 1, fill = pop)) + 
  ggrastr::rasterise(geom_bar(stat = "identity")) + theme_void() + theme(legend.position = "bottom") +
  scale_fill_manual(values = myIntegratedColors)

ggsave(gridExtra::grid.arrange(bar_dpt , bar_ct + theme(legend.position = "none"), ncol = 1),
       filename = "8.Trajectory_inference/On-integration/figures/Scatter-ZOOM_barLegend.pdf",
       width = 5.11, height = 2)
       
ggsave(p, filename = "8.Trajectory_inference/On-integration/figures/Scatter-ZOOM.pdf",
       width = 5.11, height = 4)


library(dyno)
paga_traj <- readRDS("8.Trajectory_inference/On-integration/output_Python/int-EPs_PreAs-ASCs-Meso-VSMPs_Singlets_PAGA-model.rds")
pseudotime <- read.csv("8.Trajectory_inference/On-integration/output_Python/pseudotime.csv")
rownames(pseudotime) <- pseudotime$cell_id

df <- as.data.frame(paga_traj$dimred)
df$pseudotime <- pseudotime[rownames(df),"dpt_pseudotime"]
p <- ggplot(df,aes(comp_1, comp_2, col = pseudotime)) + geom_point(size = 0.8) + 
  scale_color_gradientn(colors = rev(c('#a50026','#d73027','#f46d43','#fdae61',
                                      '#fee08b','#ffffbf','#d9ef8b','#a6d96a',
                                      '#66bd63','#1a9850','#006837'))) + 
  theme_void() +
  theme(axis.line = element_blank())
ggsave(p, filename = "8.Trajectory_inference/On-integration/figures/FA2_coloredPseudotime.pdf",
       width = 5.02, height = 2.85)


## SOME EMT genes 
g <- c("WT1", "VEGFA", "TWIST1", "TWIST2", "ALDH1A2", "SNAI2",
       "ZFPM2", "TGFB1", "TGFB3", "FERMT2", "GLIPR2", "DCLK1")
g <- data.annot[data.annot$gene_short_name %in% g,]  
