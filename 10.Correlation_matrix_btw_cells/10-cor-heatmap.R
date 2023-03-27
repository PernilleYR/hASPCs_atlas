

int <- readRDS(file = "5.Integration/output/Seurat_2000HVGs.Rds")
int_harm <- readRDS("5.Integration/output/Harmony.Rds")

##---------------------------------------------##
##------Randomly sample fraction of cells------##
##---------------------------------------------##
# int$depot <- substr(int$batch, 1,2)
# int$depot[int$depot %in% c("MG", "MK")] <- "ME"

#set.seed(1721)
set.seed(1234)
cells_for_plot <- c()
n = 1000
c_oi <- c("ASCs", "PreAs")
for(l in c_oi){
  cells_oi <- names(int$myIntegratedClustering[int$myIntegratedClustering == l])
  if(length(cells_oi) < 35){
    cells_for_plot <- cells_for_plot
  }else if(length(cells_oi) < n){
    cells_for_plot <- c(cells_for_plot, cells_oi)
  }else{
    cells_for_plot <- c(cells_for_plot,
                        sample(cells_oi, size = n, replace = F))
  }
}

cells_for_plot <- unique(cells_for_plot)
table(int$depot[cells_for_plot])

##---------------------------------------------##
##---------CALCULATE COR IN PCA SPACE----------##
##---------------------------------------------##
# 
# n_d = 60; m_hc = "ward.D2"; n_cells = 1000; suffix = ""
# plot_myCorrHeatmap(reduc = int@reductions$pca@cell.embeddings, n_dim = n_d,  
#                    cells_to_plot = cells_for_plot, method_hclust = m_hc, 
#                    path = paste0("10.Correlation_matrix_btw_cells/heatmap_PCA_",n_d, "dim_", m_hc, "_", n_cells, "cells_", suffix, ".pdf"))
# DefaultAssay(int) <- "integrated"
# test <- subset(int, cells = cells_for_plot)
# test <-  ScaleData(test, vars.to.regress = c("nCount_RNA", "nFeature_RNA"))
# test <- RunPCA(test)

d <- int@reductions$pca@cell.embeddings[names(x[x %in% c(1,3)]), ]
d <- cor(t(d))
  
  annot_c <- data.frame(row.names = colnames(d), 
                        cell_type= int$myIntegratedClustering[colnames(d)],
                        depot = int$depot[colnames(d)])
  annot_c$depot <- factor(annot_c$depot, levels = c("EP", "MG", "SC", "PR"))
  levels(annot_c$depot) <- c("EP", "MC", "SC", "PR")
  annot_c <- annot_c[ order(annot_c$depot),]
  annot_c$cell_type <- factor(as.character(annot_c$cell_type),
                              levels = c("ASCs", "PreAs"))
  annot_c <- annot_c[ order(annot_c$cell_type),]
  my_cols <-  list(cell_type = myIntegratedColors[c("ASCs", "PreAs")],
                    depot = myDepotsColors)
  
  d <- d[rownames(annot_c), rownames(annot_c)]
  d <- d[names(x[ x %in% c(1,3)]), names(x[ x %in% c(1,3)])]
  big_pheat <- pheatmap::pheatmap(d,
                                  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                                  annotation_col = annot_c,
                                  annotation_row = annot_c,
                                  clustering_method = "ward.D2",
                                  cluster_rows =T, cluster_cols = T,
                                  show_colnames = F, scale = "none",
                                  annotation_colors = my_cols,
                                  border_color = "NA", show_rownames = F )
ggsave(big_pheat, "10.Correlation_matrix_btw_cells/heatmap_ASCs-PreAs_PCAofint_2000cells.pdf")
  
  
  
df <- reshape2::melt(d)
df$catVar1 <- annot_c[as.character(df$Var1) , "cell_type"]
df$catVar2 <- annot_c[as.character(df$Var2) , "cell_type"]

df$depotVar1 <- annot_c[as.character(df$Var1) , "depot"]
df$depotVar2 <- annot_c[as.character(df$Var2) , "depot"]

df$Var1 <- factor(as.character(df$Var1), levels = rownames(annot_c))
df$Var2 <- factor(as.character(df$Var2), levels = rownames(annot_c))

p <- ggplot(df, aes(x = Var1, y = Var2, fill = value)) + 
  ggrastr::rasterise(geom_tile(), dpi = 200) + 
  theme(axis.text = element_blank(), axis.title = element_blank()) + 
  scale_fill_distiller(palette="RdYlBu", type = "div")
ggsave(p, filename = "10.Correlation_matrix_btw_cells/heatmap_PCA-50_ASCs-PreAs.pdf",
       width = 2.45, height = 1.62)
#Saving 2.45 x 1.62 in image

t <- df %>% filter(catVar1 == catVar2)
p <- ggplot(t, aes(x = catVar1, y = value, fill = catVar1)) + geom_boxplot() + 
  mashaGgplot2Theme + scale_fill_manual(values = myIntegratedColors) + 
  theme(legend.position = "none") + ylab("Correlation in PCA space") + xlab("")
ggsave(p, filename = "10.Correlation_matrix_btw_cells/Boxplot_Corr-ASCs-PreAs_PC50.pdf",
       width = 2.72, height = 3.34)
#Saving 2.72 x 3.34 in image
