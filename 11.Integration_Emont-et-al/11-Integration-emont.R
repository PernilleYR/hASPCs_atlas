################################################################
#                                                              #
#                 Transfer emont to our data                   #
#                                                              #
################################################################

### Author: Pernille
### Date: 04.09.2022 - adapted from old script
### Datasets: scRNA-seq Depots: SC - EP - MK - MG - PR - GB 
###                     Patients: B0, B1, B7, L3, L11, L12 (B - Bariatric, L - Lean)
###           Emont scRNAseq published in Nature
### Goal: Compare markers between human and mouse ASPCs 

library(ggplot2); library(data.table); library(Seurat); library(biomaRt); library(dplyr)

setwd("~/SVRAW1/prainer/hASPCs/PAPER/10X_scRNA-seq/")
source("Utility/General_utils.R")

library(gridExtra)

##---------------------------------------------##
##-----------------Loading data----------------##
##---------------------------------------------##

int <- readRDS("~/SVRAW1/prainer/hASPCs/PAPER/10X_scRNA-seq/5.Integration/output/Seurat_2000HVGs.Rds")
DefaultAssay(int) <- "integrated"

#myseu <- readRDS("0.data/List_seurat_objects.rds")
emont <- readRDS("~/SVRAW1/prainer/hASPCs/PAPER/10X_scRNA-seq/11.Integration-emont/data/human_all.rds")

emont_meta <- emont@meta.data

##---------------------------------------------##
##------------Extract Meso & ASPCs-------------##
##---------------------------------------------##
dir.create("11.Integration-emont/emont_int_ASPCs-Meso/")

# - Subset Meso & ASPCs
sub_emont <- subset(emont, cell_type2 %in% c("mesothelium", "ASPC"))
rm(emont)

# - Run dim red
DefaultAssay(sub_emont) <- "integrated"
sub_emont <- RunPCA(sub_emont, npcs = 50)
sub_emont <- RunUMAP(sub_emont, dims = 1:50) 
sub_emont <- RunTSNE(sub_emont, dims = 1:50)

# - Clustering
DefaultAssay(sub_emont) <- "integrated"
sub_emont <- FindNeighbors(sub_emont, dims = 1:50)
for(i in seq( 0.1, 2, 0.1)){
  sub_emont <- FindClusters(object = sub_emont, resolution = i, print.output = T)
}

saveRDS(sub_emont, file = "11.Integration-emont/emont_int_ASPCs-Meso/emont_ASPCs-Meso_seurat.Rds")

# - Plot clustering
m <- colnames(sub_emont@meta.data)[grep("integrated_snn_res", colnames(sub_emont@meta.data))]
p <- lapply( m, function(i) DimPlot(sub_emont, group.by = i, reduction = "tsne"))

pdf("11.Integration-emont/emont_int_ASPCs-Meso/tsne_int-emont-ASPCs-meso_colored_clustering_res.0.1-2.pdf")
marrangeGrob(p, ncol = 2, nrow = 2)
dev.off()

##---------------------------------------------##
##--------------Split by samples---------------##
##---------------------------------------------##

# - Split by Samples
sub_emont_split <- SplitObject(sub_emont, split.by = "sample")

sub_emont_split <- sub_emont_split[c("Hs_SAT_253","Hs_SAT_09", "Hs_SAT_10" )]

mylist <- readRDS("11.Integration-emont/output/datasets_withPredictions.Rds")

# - Extract counts to rename from symbol to ens_id
n <- names(sub_emont_split)
sub_emont_split <- lapply(sub_emont_split, myfunction_to_convert)
names(sub_emont_split) <- n

# - Create seurat obj and normalize
sub_emont_split <- prepareData(sub_emont_split)

##---------------------------------------------##
##-----------Preprocess each sample------------##
##---------------------------------------------##

# Find integration feat
features <- SelectIntegrationFeatures(object.list = append(sub_emont_split, myseu))

# - Scale and PCA
sub_emont_split <- lapply(X = sub_emont_split, FUN = function(x) {
  x <- FindVariableFeatures(x, nFeatures = 2000)
  x <- ScaleData(x, verbose = FALSE, features = features)
  x <- RunPCA(x, verbose = FALSE, npcs = 50, features = features)
})

# Generate tSNE
doTSNE <- function(d){
  print(substr(colnames(d), 1, 10)[1])
  d <- JackStraw(d, num.replicate = 50, dims = 25, reduction = "pca")
  d <- ScoreJackStraw(d, dims = 1:25)
  tokeep <- rle(d@reductions$pca@jackstraw@overall.p.values[,"Score"] < 0.05)$lengths[1]
  print(tokeep)
  d <- RunTSNE(d, dims = 1:tokeep)
  return(d)
}
sub_emont_split <- lapply(sub_emont_split, doTSNE)  

# Clustering
perform_clustering <- function(d){
  d <- FindNeighbors(d, reduction = "pca", dims = d@commands$RunTSNE$dims)
  for( i in seq( 0.1, 3, 0.1)){
    d <- FindClusters(object = d, reduction.type = "pca", resolution = i, print.output = T)}
  return(d)
}
n <- names(sub_emont_split)
sub_emont_split <- lapply(sub_emont_split, perform_clustering)  
names(sub_emont_split) <- n

# Plot clustering
for(n in names(sub_emont_split)){
  m <- colnames(sub_emont_split[[n]]@meta.data)[grep("RNA_snn_res", colnames(sub_emont_split[[n]]@meta.data))]
  p <- lapply( m, function(i) DimPlot(sub_emont_split[[n]], group.by = i))
  pdf( paste0("11.Integration-emont/clustering/tSNE_", n, "_Clustering_res0.1-3.pdf"))
    marrangeGrob(p, nrow = 2, ncol = 2)
  dev.off()  
}

# Add emont metadata
add_meta <- function(d){
  d@meta.data <- cbind(d@meta.data, emont_meta[colnames(d), c("cell_type", "cell_type2")])
  return(d)
}
n <- names(mylist)
t <- lapply(mylist, add_meta)
names(t) <- n

s <- mylist

mylist <- readRDS("11.Integration-emont/output/datasets_withPredictions.Rds")

mylist[["Hs_SAT_253"]] <- sub_emont_split$Hs_SAT_253
mylist[["Hs_SAT_09"]] <- sub_emont_split$Hs_SAT_09
mylist[["Hs_SAT_10"]] <- sub_emont_split$Hs_SAT_10

##---------------------------------------------##
##---------------Transfer on int---------------##
##---------------------------------------------##
for(i in names(mylist)){
  mylist[[i]]@meta.data[grep("predict", colnames(mylist[[i]]@meta.data))] <- NULL
}
myTransferFunction <- function(mydata.query, mydata.ref = int, ref.data.clust = int$myIntegratedClustering){
  anchors.prediction <- FindTransferAnchors(reference = mydata.ref, query = mydata.query,
                                            dims = 1:30, reference.reduction = "pca")
  predictions <- TransferData(anchorset = anchors.prediction, 
                              refdata = ref.data.clust,
                              dims = 1:30)
  mydata.query <- AddMetaData(mydata.query, metadata = predictions)
  return(mydata.query)
}
n <- names(mylist)
mylist <- lapply(mylist, myTransferFunction)
names(mylist) <- n

dir.create("11.Integration-emont/output")
saveRDS(mylist, "11.Integration-emont/output/datasets_withPredictions.Rds")

##---------------------------------------------##
##--------------Extract prediction-------------##
##---------------------------------------------##

library(tidyr); library(tibble)
#####
# predic <- grep("prediction", colnames(data[[d]]@meta.data))
# df <- data[[d]]@meta.data[, predic, drop = F]
# df <- df[, -grep("res.0.2_6", colnames(df))]
# df <- df[, -grep("5", colnames(df))]
# #df <- df[, -grep("max", colnames(df))]
# 
# colnames(df) <- sapply(strsplit(colnames(df),"\\."), `[`, 3)
# 
# df <- df %>% 
#   rownames_to_column('id') %>%  # creates an ID number
#   gather("celltype", "score.prediction", 2:14) %>% 
#   # http://statseducation.com/Introduction-to-R/modules/tidy%20data/gather/ kind of a melt
#   group_by(id) %>% 
#   dplyr::slice(which.max(score.prediction)) 
# df <- as.data.frame(df)
#####

# Plot predicted.id
p <- lapply(names(mylist), function(d)
  DimPlot(mylist[[d]], group.by = "predicted.id", cols = myIntegratedColors) + 
    ggtitle(d))
ggsave(filename = paste0("11.Integration-emont/prediction/tSNE_predictedID.pdf"),
       plot = marrangeGrob(p, nrow = 2, ncol = 2))

# Extract prediction
prediction <- lapply(names(mylist), function(d){
  c <- colnames(mylist[[d]]@meta.data)[grep("prediction", colnames(mylist[[d]]@meta.data))]
  out <- mylist[[d]]@meta.data[, c("predicted.id", "cell_type", "cell_type2",c)]
  out$cellid <- rownames(out)
  out$sample <- d
  return(out)})
prediction <- as.data.frame(rbindlist(prediction))
prediction$Depot <- sapply(strsplit(prediction$sample,"_"), `[`, 2)

# Convert to % of emont cell_type2
  # All
prediction_percent <- prediction %>% 
  group_by(sample, cell_type2) %>% 
  dplyr::count(predicted.id) %>% 
  dplyr::mutate(Freq = round(n/sum(n)*100, 2))
prediction_percent$Depot <- sapply(strsplit(prediction_percent$sample,"_"), `[`, 2)


# OATs
prediction_percent_OAT <- prediction %>% 
  filter(Depot == "OAT") %>% 
  group_by(sample, cell_type2) %>% 
  dplyr::count(predicted.id) %>% 
  dplyr::mutate(Freq = round(n/sum(n)*100, 2))
  
  # SATs
prediction_percent_SAT <- prediction %>% 
  filter(Depot == "SAT") %>% 
  group_by(sample, cell_type2) %>% 
  dplyr::count(predicted.id) %>% 
  dplyr::mutate(Freq = round(n/sum(n)*100, 2))

# prediction_percent_predicted.id <- prediction %>% 
#   group_by(sample, predicted.id) %>% 
#   dplyr::count(cell_type2) %>% 
#   dplyr::mutate(Freq = round(n/sum(n)*100, 2))


##---------------------------------------------##
##----Show that Meso/ASPCs match btw studies---##
##---------------------------------------------##

predic_percent <- prediction_percent

# Calculate mean and std
predic_percent <- as.data.table(predic_percent)
dataSum <- predic_percent[, .(M = mean(Freq, na.rm = T), 
                                  S = sd(Freq, na.rm = T)),
                              by = .(cell_type2, predicted.id)]

# Order
predic_percent$cell_type2 <- factor(as.character(predic_percent$cell_type2), 
                                        levels = c("ASPC", "mesothelium"))
predic_percent$predicted.id <- factor(as.character(predic_percent$predicted.id), 
                                          levels = names(myIntegratedColors))
dataSum$cell_type2 <- factor(as.character(dataSum$cell_type2), 
                             levels = c("ASPC", "mesothelium"))
dataSum$predicted.id <- factor(as.character(dataSum$predicted.id), 
                               levels = names(myIntegratedColors))

# Plot
thePlot <- ggplot(dataSum, aes(x = cell_type2, y = M, 
                               fill = predicted.id)) + 
  geom_bar(stat = "identity", color = "black", 
           position = position_dodge(.9), alpha = 0.6) +
  geom_errorbar(aes(ymin = M, ymax = M + S), width = .2,
                position = position_dodge(.9)) +
  geom_jitter(data = predic_percent, 
              mapping = aes(x = cell_type2, 
                            y = Freq,
                            col = predicted.id, 
                            group = predicted.id,
                            shape = Depot),
              size = 2, 
              position = position_jitterdodge(dodge.width = 0.9, 
                                              jitter.width = 0.1)) +
  scale_shape_manual(values = c(1,19)) + 
  scale_fill_manual(values = myIntegratedColors) +
  scale_color_manual(values = myIntegratedColors) + 
  mashaGgplot2Theme + theme(legend.position = "bottom") +
  ylab(label = "%") + xlab("") 
ggsave(plot = thePlot, "11.Integration-emont/prediction/Barplot_predicted.id_per_emont-celltype2.pdf")
#Saving 7.71 x 6.84 in image


##---------------------------------------------##
##--Plot individ. tsne IGFBP2 prediction score-##
##---------------------------------------------##

## Plot tSNEs
P <- list()
n = names(mylist); n <- n[grep("OAT",n)]
for(i in n){
  P[[i]] <- plot.tsne.colored_continuous(mylist[[i]]@reductions$tsne@cell.embeddings, 
                                         axis = T, x.lab = "tSNE1", y.lab = "tSNE2",
                                         col = c('#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b'),
                                         data.categories = mylist[[i]]$prediction.score.IGFBP2) + 
    ggtitle(i)
}
ggsave(gridExtra::grid.arrange(P$Hs_OAT_01, P$Hs_OAT_09, P$Hs_OAT_10, P$Hs_OAT_12, P$Hs_OAT_13,
                        P$Hs_OAT_253, P$Hs_OAT_254, P$Hs_OAT_255, P$Hs_OAT_256, P$Hs_OAT_266, 
                        nrow = 2, ncol = 5), 
       filename = "11.Integration-emont/prediction/tSNEs_colored_prediction.score.IGFBP2.leg.pdf",
       width = 12.85, height = 3.83)

t <- theme(legend.position = "none")
ggsave(gridExtra::grid.arrange(P$Hs_OAT_01 + t, P$Hs_OAT_09 + t, P$Hs_OAT_10 + t, P$Hs_OAT_12 + t, P$Hs_OAT_13 + t,
                               P$Hs_OAT_253 + t, P$Hs_OAT_254 + t, P$Hs_OAT_255 + t, P$Hs_OAT_256 + t, P$Hs_OAT_266 + t, 
                               nrow = 2, ncol = 5), 
       filename = "11.Integration-emont/prediction/tSNEs_colored_prediction.score.IGFBP2.woleg.pdf",
       width = 9.56, height = 3.83)
       

##---------------------------------------------##
##----A cluster enriched for IGFBP2 score?-----##
##---------------------------------------------##

enriched_clust <- lapply(names(mylist), FUN = function(n){
  print(n)
  out <- lapply(seq(0.1,3,0.1), function(c) find_enriched_clust(mylist[[n]], 
                                                                  clustering = paste0("RNA_snn_res.", c)))
  names(out) <- paste0("res.", seq(0.1,3,0.1))
  return(out)
})
names(enriched_clust) <- names(mylist)

saveRDS(enriched_clust, "11.Integration-emont/output/enriched_clust.Rds")

#####
# [1] "Hs_OAT_256"
# [1] "RNA_snn_res.0.1"
# [1] "RNA_snn_res.0.2"
# [1] "RNA_snn_res.0.3"
# [1] "RNA_snn_res.0.4"
# [1] "RNA_snn_res.0.5"
# [1] "RNA_snn_res.0.6"
# [1] "RNA_snn_res.0.7"
# [1] "RNA_snn_res.0.8"
# [1] "RNA_snn_res.0.9"
# [1] "RNA_snn_res.1"
# [1] "RNA_snn_res.1.1"
# [1] "RNA_snn_res.1.2"
# [1] "RNA_snn_res.1.3"
# [1] "RNA_snn_res.1.4"
# [1] "RNA_snn_res.1.5"
# [1] "RNA_snn_res.1.6"
# [1] "RNA_snn_res.1.7"
# [1] "RNA_snn_res.1.8"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.1.9"
# [1] "RNA_snn_res.2"
# [1] "RNA_snn_res.2.1"
# [1] "RNA_snn_res.2.2"
# [1] "RNA_snn_res.2.3"
# [1] "RNA_snn_res.2.4"
# [1] "RNA_snn_res.2.5"
# [1] "RNA_snn_res.2.6"
# [1] "RNA_snn_res.2.7"
# [1] "RNA_snn_res.2.8"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.2.9"
# [1] "RNA_snn_res.3"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# 
# [1] "Hs_OAT_01"
# [1] "RNA_snn_res.0.1"
# [1] "RNA_snn_res.0.2"
# [1] "RNA_snn_res.0.3"
# [1] "RNA_snn_res.0.4"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.0.5"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.0.6"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.0.7"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.0.8"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.0.9"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.1"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.1.1"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.1.2"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.1.3"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.1.4"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.1.5"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.1.6"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.1.7"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.1.8"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.1.9"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.2"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.2.1"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.2.2"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.2.3"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.2.4"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.2.5"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.2.6"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.2.7"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.2.8"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.2.9"
# [1] "YES!"
# [1] "IGFBP2 WON!"
# [1] "RNA_snn_res.3"
# [1] "YES!"
# [1] "IGFBP2 WON!"
#####




dir.create("11.Integration-emont/enriched_in_clust/")
# Signi
P <- plot_enriched_cluster("Hs_OAT_256", clustering = "RNA_snn_res.1.8", c = 15)
P <- plot_enriched_cluster("Hs_OAT_01", clustering = "RNA_snn_res.0.4", c = 6)
# not signi
P <- plot_enriched_cluster("Hs_OAT_09", clustering = "RNA_snn_res.0.8", c = 6)
P <- plot_enriched_cluster("Hs_OAT_10", clustering = "RNA_snn_res.1.2", c = 8)
P <- plot_enriched_cluster("Hs_OAT_13", clustering = "RNA_snn_res.0.3", c = 4)
P <- plot_enriched_cluster("Hs_OAT_255", clustering = "RNA_snn_res.1.5", c = 12)
P <- plot_enriched_cluster("Hs_OAT_266", clustering = "RNA_snn_res.0.4", c = 5)

## plot clustering
l <- list("Hs_OAT_01" = c("res.0.4", 6), #Signi
     "Hs_OAT_09" = c("res.0.8", 4),
     "Hs_OAT_10" = c("res.1.2", 8),
     "Hs_OAT_12" = NULL, 
     "Hs_OAT_13" = c("res.0.3", 4),
     "Hs_OAT_253" = NULL,
     "Hs_OAT_254" = c("res.2.7", 19),
     "Hs_OAT_255" = c("res.1.5", 12),
     "Hs_OAT_256" = c("res.1.8", 15), #Signi, 
     "Hs_OAT_266" = c("res.0.4", 5))
P <- list()
for(n in names(mylist)){
  print(n)
  if(!is.null(l[[n]])){
    mylist[[n]]$toplot <- "Rest"; mylist[[n]]$toplot[ mylist[[n]]@meta.data[, paste0("RNA_snn_", l[[n]][1]) ] == l[[n]][2]] <- "IGFBP2"
    mylist[[n]]$toplot <- factor(  mylist[[n]]$toplot, levels = c("Rest", "IGFBP2"))
    pc <- plot.tsne.colored_discrete_2( mylist[[n]]@reductions$tsne@cell.embeddings,
                                        axis = T, x.lab = "tSNE1", y.lab = "tSNE2", 
                                        data.categories =  mylist[[n]]$toplot,
                                        col = c("gray88","#1f78b4")) + 
      ggtitle(n)
    P[[n]] <- pc
  }else{
    P[[n]] <- ggplot()
  }

  
}
t <- theme(legend.position = "none")
ggsave(gridExtra::grid.arrange(P$Hs_OAT_01 + t, P$Hs_OAT_09 + t, P$Hs_OAT_10 + t, P$Hs_OAT_12 + t, P$Hs_OAT_13 + t,
                               P$Hs_OAT_253 + t, P$Hs_OAT_254 + t, P$Hs_OAT_255 + t, P$Hs_OAT_256 + t, P$Hs_OAT_266 + t, 
                               nrow = 2, ncol = 5), 
       filename = "11.Integration-emont/prediction/tSNEs_colored_potential.IGFBP2.clust.woleg.pdf",
       width = 9.56, height = 3.83)


## plot boxplot
P <- PB <- list()
for(n in names(l)){
  print(n)
  if(!is.null(l[[n]])){
    t <- enriched_clust[[n]][[l[[n]][1]]]$data
    t$clust_plot <- t$clust
    levels(t$clust_plot)[levels(t$clust_plot) != l[[n]][2]] = 0
    t$pop <- sapply(strsplit(as.character(t$pop),"\\."), `[`, 3)
    t$pop[t$pop == "PR"] <- "PR specific"
    t$pop[t$pop == "res"] <- "res.0.2_8"
    t$pop[t$pop == "CHI3L1"] <- "CHI3L1-2"
    
    print(max(as.numeric(t$clust)))
    steps = round(length(8:90)/max(as.numeric(t$clust)))
    mygrays <- paste0("gray", seq(90, 8, -steps))
    p <- ggplot(t %>% filter(pop == "IGFBP2"),
                aes(x = clust_plot, y = score, fill = clust_plot)) + 
      geom_jitter(aes(color = clust), width = 0.2, size = 0.5) +
      geom_boxplot(alpha = 0.8) + 
      scale_color_manual(values = c(mygrays[1:as.numeric(l[[n]][2])],
                                    "#1f78b4",
                                    mygrays[(as.numeric(l[[n]][2])+1):length(mygrays)]) ) +
      scale_fill_manual(values = c("gray", "#1f78b4"))+
      mashaGgplot2Theme+ 
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.y = element_text(size = 5))
    
    pb <- ggplot(t %>% filter(clust == l[[n]][2]), aes(y = score, x = pop, fill = pop)) + 
      #geom_jitter(aes(group = pop, col = pop), width = 0.3) + 
      geom_boxplot(alpha = 0.8) + 
      scale_fill_manual(values = myIntegratedColors) + 
      scale_color_manual(values = myIntegratedColors) + 
      mashaGgplot2Theme + 
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.y = element_text(size = 5),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    
    P[[n]] <- p
    PB[[n]] <- pb
  }else{
    P[[n]] <- ggplot()
    PB[[n]] <- ggplot()
  }
  
}
ggsave(gridExtra::grid.arrange(P$Hs_OAT_01, P$Hs_OAT_09, P$Hs_OAT_10, P$Hs_OAT_12, P$Hs_OAT_13,
                               P$Hs_OAT_253, P$Hs_OAT_254, P$Hs_OAT_255, P$Hs_OAT_256, P$Hs_OAT_266, 
                               nrow = 2, ncol = 5), 
       filename = "11.Integration-emont/prediction/boxplot_score.IGFBP2_selectedClust.pdf",
       width = 9.56, height = 3.83)

ggsave(gridExtra::grid.arrange(PB$Hs_OAT_01, PB$Hs_OAT_09, PB$Hs_OAT_10, PB$Hs_OAT_12, PB$Hs_OAT_13,
                               PB$Hs_OAT_253, PB$Hs_OAT_254, PB$Hs_OAT_255, PB$Hs_OAT_256, PB$Hs_OAT_266, 
                               nrow = 2, ncol = 5), 
       filename = "11.Integration-emont/prediction/boxplot_scores.pops_selectedCluster.pdf",
       width = 9.56, height = 3.83)

#.../11.Integration-emont/prediction/tSNEs_colored_prediction.score.IGFBP2.pdf

## Plot prediction IGFBP2 per depot
p <- ggplot(prediction, aes(y = prediction.score.IGFBP2, 
                            x = Depot, 
                            group = sample)) + 
  #geom_point(position=position_jitterdodge(jitter.width = 0.05), alpha = 0.3)+
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.2, fill = "#6F3996") + 
  
  mashaGgplot2Theme + 
  theme(legend.position = "none") 
ggsave(p, filename = "11.Integration-emont/prediction/Boxplot_Score_byDepot.pdf")  
#Saving 3.91 x 4.98 in image

## Plot on integrated data

##---------------------------------------------##
##-----Enriched clust in integration emont-----##
##---------------------------------------------##

sub_emont@meta.data <- cbind(sub_emont@meta.data, prediction[colnames(sub_emont),])

out <- lapply(seq(0.1,2,0.1), function(c) find_enriched_clust(sub_emont, 
                                                              clustering = paste0("integrated_snn_res.", c)))
names(out) <- paste0("res.", seq(0.1,2,0.1))

names(out) <- paste0("integrated_snn_", names(out))
P <- plot_enriched_cluster(enr_c = out, seu = sub_emont, clustering = "integrated_snn_res.1", c = 19)

t <- out$integrated_snn_res.0.6$data
t$clust_plot <- t$clust
levels(t$clust_plot)[levels(t$clust_plot) != 14] = 0
t$pop <- sapply(strsplit(as.character(t$pop),"\\."), `[`, 3)
t$pop[t$pop == "PR"] <- "PR specific"
t$pop[t$pop == "res"] <- "res.0.2_8"
t$pop[t$pop == "CHI3L1"] <- "CHI3L1-2"
t$clust <- as.factor(as.numeric(as.character(t$clust)))
t$cell_type2 <- sub_emont$cell_type2[t$sample]
t$cell_type2[t$clust == 14] <- "IGFBP2"
print(max(as.numeric(t$clust)))
steps = round(length(8:90)/max(as.numeric(t$clust)))
mygrays <- paste0("gray", seq(90, 8, -steps))
p <- ggplot(t %>% filter(pop == "IGFBP2"),
            aes(x = clust_plot, y = score, fill = clust_plot)) + 
  geom_jitter(aes(color = clust), width = 0.2, size = 0.5) +
  geom_boxplot(alpha = 0.8) + 
  scale_color_manual(values = c(mygrays[1:14],
                                "#1f78b4",
                                mygrays[(15):length(mygrays)]) ) +
  scale_fill_manual(values = c("gray", "#1f78b4"))+
  mashaGgplot2Theme+ 
  theme(legend.position = "none",
        axis.title = element_blank())

pb <- ggplot(t %>% filter(clust == 14), aes(y = score, x = pop, fill = pop)) + 
  #geom_jitter(aes(group = pop, col = pop), width = 0.3) + 
  geom_boxplot(alpha = 0.8) + 
  scale_fill_manual(values = myIntegratedColors) + 
  scale_color_manual(values = myIntegratedColors) + 
  mashaGgplot2Theme + 
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p1 <- ggplot(t %>% filter(pop == "PreAs"),
            aes(x = clust, y = score)) + 
  geom_jitter(aes(color = cell_type2), width = 0.2, size = 0.5) +
  geom_boxplot(alpha = 0.8) + 
  scale_color_manual(values = c("#E31A1C", "#1f78b4", "#810F7C" ))+
  mashaGgplot2Theme+ 
  theme(legend.position = "none",
        axis.title = element_blank()) + 
  ggtitle("score PreAs")
p2 <- ggplot(t %>% filter(pop == "IGFBP2"),
             aes(x = clust, y = score)) + 
  geom_jitter(aes(color = cell_type2), width = 0.2, size = 0.5) +
  geom_boxplot(alpha = 0.8) + 
  scale_color_manual(values = c("#E31A1C", "#1f78b4", "#810F7C" ))+
  mashaGgplot2Theme+ 
  theme(legend.position = "none",
        axis.title = element_blank()) + 
  ggtitle("Score IGFBP2")

ggsave(grid.arrange(p,pb, ncol = 2),
       width = 7.47, height = 5, 
       filename = "11.Integration-emont/emont_int_ASPCs-Meso/boxplot_scores_selectedClust.pdf")

ggsave(grid.arrange(p1,p2, ncol = 2),
       width = 7.47, height = 5, 
       filename = "11.Integration-emont/emont_int_ASPCs-Meso/boxplot_scoresPreAs-IGFBP2_acrossClust.pdf")

##---------------------------------------------##
##-----Enriched clust in integration emont-----##
##---------------------------------------------##

ggsave(plot = DimPlot(sub_emont, group.by = "predicted.id", reduction = "umap",
        order = c("IGFBP2", "CILP", "HHIP", "CHI3L1-2", "IFIT", "ASCs", "Meso"),
        cols = myIntegratedColors), width = 6.83, height = 4.98,
       "11.Integration-emont/emont_int_ASPCs-Meso/Umap_colored_predicted.id.pdf")

##---------------------------------------------##
##--------Potential stats-------##
##---------------------------------------------##

meta <- sub_emont@meta.data[, c("predicted.id", "prediction.score.IGFBP2", "sample", "individual", "compartment",
                        "bmi.range", "sex", "race", "ethnicity", "age", "age.range", "bmi")]
meta_f <- meta %>% filter(prediction.score.IGFBP2 > 0)
ggplot(meta_f, aes(y = `prediction.score.IGFBP2`, x = age, col = sample)) + geom_boxplot()
ggplot(meta_f, aes(y = `prediction.score.IGFBP2`, x = bmi, col = sample)) + geom_boxplot()
ggplot(meta_f, aes(y = `prediction.score.IGFBP2`, x = bmi.range, col = sample)) + geom_boxplot()
ggplot(meta_f, aes(y = `prediction.score.IGFBP2`, x = sex, col = sample)) + geom_boxplot()


## extract cells of potential IGFBP2 cluster
cells_to_keep <- c()
for(i in names(l)){
  print(i)
  if(!is.null(l[[i]])){
    x <- mylist[[i]]@meta.data[, paste0("RNA_snn_", l[[i]][1]), drop = F]
    x <- x[as.character(x[,1]) == l[[i]][2],, drop = F]
    cells_to_keep <- c(cells_to_keep,
                       rownames(x))
  }
}

meta_f_test <- meta_f[rownames(meta_f) %in% cells_to_keep,]
meta_test <- meta[cells_to_keep,]
ggplot(meta_f_test, aes(y = `prediction.score.IGFBP2`, x = bmi.range, col = sample)) + geom_violin()


##---------------------------------------------##
##----A cluster enriched for IGFBP2 score?-----##
##---------------------------------------------##

myDEGs <- readRDS(file = "6.TopDEGs/DEGs.Rds")
IGFBP2_g <- myDEGs$IGFBP2 %>%  filter(avg_logFC_all > 0)
IGFBP2_g$geneID

IGFBP2_EdgeR <- read.table("~/SVRAW1/prainer/hASPCs/10X/Joint_Analysis/Seurat/SC0_SC1_SC7_EP0_EP1_EP7_MG7_MK7/DE_EdgeR/DEEdgeR_filtered_IGFBP2.txt")
IGFBP2_EdgeR <- IGFBP2_EdgeR %>% filter(padj < 0.05) %>% filter(logFC > 0.25)

sub_emont <- AddModuleScore(sub_emont, features = list(IGFBP2_EdgeR$gene_id,
                                                       IGFBP2_EdgeR$gene_id[1:100],
                                                       IGFBP2_g$geneID), 
                            name = c("IGFBP2oldFull_score", "IGFBP2old_score", "IGFBP2_new"))
int <- AddModuleScore(int, features = list(rownames(IGFBP2_EdgeR),
                                           rownames(IGFBP2_EdgeR)[1:100],
                                           rownames(IGFBP2_g)), 
                            name = c("IGFBP2oldFull_score", "IGFBP2old_score", "IGFBP2_new"))

p <- lapply(c("IGFBP2oldFull_score1", "prediction.score.IGFBP2"), 
       function(x) plot.tsne.colored_continuous_2(sub_emont@reductions$umap@cell.embeddings, 
                                                  data.categories = sub_emont@meta.data[, x, drop = F], label.legend = "",
                                                  col = c("#878787","#BABABA","#E0E0E0","#F7F7F7",'#D1E5F0','#92C5DE','#4393C3','#2166AC', "#053061")) + ggtitle(x))
ggsave(gridExtra::grid.arrange(p[[1]], p[[2]], ncol = 2),
       filename = "11.Integration-emont/emont_int_ASPCs-Meso/UMAP_colored_Module-Prediction-score-IGFBP2.pdf",
       height = 5.24/2, width = 6.83)

p1 <- DimPlot(sub_emont, group.by = "predicted.id", order = c("IGFBP2", "CILP", "HHIP", "IFIT", "Meso", "ASCs"), cols = myIntegratedColors, 
              raster = T, raster.dpi = c(600, 600) )
p2 <- DimPlot(sub_emont, group.by = "cell_type", order = c(paste0("hASPC", c(6,5,4,3,2)), "hMes3"), raster = T,raster.dpi = c(600, 600),
              cols = c("hASPC1" = '#cb181d', "hASPC2" = '#238b45', "hASPC3" = '#fb6a4a',"hASPC4" = "#8da0cb", 
                       "hASPC5" = '#74c476', "hASPC6" = '#1f78b4',"hMes1" ='#6a51a3',  "hMes2" ='#9e9ac8', "hMes3" ='#cbc9e2'))
p3 <- DimPlot(sub_emont, group.by = "depot", cols = unname(myDepotsColors[c("SC", "EP")]), shuffle = T, raster = T, raster.dpi = c(600, 600))

ggsave(gridExtra::grid.arrange(p3+theme(legend.position = "none"),
                               p2+theme(legend.position = "none"),
                               p1+theme(legend.position = "none"), ncol = 3),
       filename = "11.Integration-emont/emont_int_ASPCs-Meso/UMAP_colored_Depot-CT-PredictID.pdf",
       height = 5.24/2, width =8.196)


# ggsave(FeaturePlot(sub_emont, feature = c("IGFBP2old_score2", "prediction.score.IGFBP2"),
#             reduction = "tsne", order = T, raster = T, raster.dpi = c(600,600),
#             cols = c("#BABABA","#BABABA","#BABABA","#E0E0E0","#F7F7F7",'#D1E5F0','#92C5DE','#4393C3','#2166AC', "#053061")),
#        filename = "11.Integration-emont/emont_int_ASPCs-Meso/tSNE_emont-Meso-ASPCs_colored_Prediction-ModuleScore-IGFBP2_RASTER.pdf", height = 3.15, width = 7.8)
#tSNE_emont-Meso-ASPCs_colored_Prediction-ModuleScore-IGFBP2.pdf
#7.80x3.15

p <- lapply(c("IGFBP2", "KIF21A", "WNT4", "PLAT"), 
            function(x) plot.colored_gene.continuous(sub_emont@reductions$umap@cell.embeddings, 
                                                     GetAssayData(sub_emont, assay = "RNA"), data.annot = data.annot, 
                                                     col = c("#878787","#BABABA","#E0E0E0","#F7F7F7",'#D1E5F0','#92C5DE','#4393C3','#2166AC', "#053061"),
                                                     x, gene.text = F) + ggtitle(x))
ggsave(gridExtra::grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], ncol = 2),
       filename = "11.Integration-emont/emont_int_ASPCs-Meso/Umap_colored-IGFBP2-markers.pdf",
       height = 5.24, width = 6.83)

p <- lapply(convert_geneID_to_data.annot(c("IGFBP2", "KIF21A", "WNT4", "PLAT"), data.annot = data.annot)[, "ens_id"], 
            function(x) plot.colored_gene.continuous(int@reductions$tsne@cell.embeddings, 
                                                     GetAssayData(int, assay = "RNA"), data.annot = data.annot, 
                                                     col = c("#878787","#BABABA","#E0E0E0","#F7F7F7",'#D1E5F0','#92C5DE','#4393C3','#2166AC', "#053061"),
                                                     x, gene.text = F) + ggtitle(x))
ggsave(gridExtra::grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], ncol = 2),
       filename = "11.Integration-emont/emont_int_ASPCs-Meso/tSNE-OURINT_colored-IGFBP2-markers.pdf",
       height = 5.24, width = 6.83)


p <- lapply(c("COL1A2", "C7", "RARRES1", "DMKN"), 
            function(x) plot.colored_gene.continuous(sub_emont@reductions$umap@cell.embeddings, 
                                                     GetAssayData(sub_emont, assay = "RNA"), data.annot = data.annot, 
                                                     col = c("#878787","#BABABA","#E0E0E0","#F7F7F7",'#D1E5F0','#92C5DE','#4393C3','#2166AC', "#053061"),
                                                     x, gene.text = F) + ggtitle(x))
ggsave(gridExtra::grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], ncol = 2),
       filename = "11.Integration-emont/emont_int_ASPCs-Meso/Umap_colored-Meso-ASPCs-markers.pdf",
       height = 5.24, width = 6.83)

##---------------------------------------------##
##-----Enriched clust in integration emont-----##
##---------------------------------------------##
emont <- readRDS("~/SVRAW1/prainer/hASPCs/PAPER/10X_scRNA-seq/11.Integration-emont/emont_int_ASPCs-Meso/emont_ASPCs-Meso_seurat.Rds")

table(emont$predicted.id, emont$Depot)
#             OAT   SAT
# IGFBP2      481     0

df <- table(emont$predicted.id[emont$Depot == "OAT"], 
            emont$sample[emont$Depot == "OAT"])

#df <- t(t(df)/colSums(df)*100)
df <- df/rowSums(df)*100
df <- data.frame(row.names = colnames(df),
                 "IGFBP2" = df["IGFBP2",])

BMI <- emont@meta.data[, c("bmi", "sample", "age", "sex")]
BMI <- BMI[!duplicated(BMI$sample),]
rownames(BMI) <- BMI$sample
BMI <- BMI[rownames(df),]
df <- cbind(df, BMI)
df$sex <- as.factor(df$sex)


ggplot(df, aes(x = bmi, y = IGFBP2)) + 
  geom_point(aes(color = sex)) + geom_smooth(method = "lm") + 
  mashaGgplot2Theme

summary(lm(IGFBP2~bmi,df))

##---------------------------------------------##
##-----------------------------------------##
##---------------------------------------------##
FeaturePlot(emont, features = c(#"TIMP1", "WT1", "RARRES1",
                                "DDIT4", "COL1A2", "GPX3", "APOD", "IGFBP7",
                                "MFAP4", "FBLN1"), order = F)
FeaturePlot(emont, features = c("TIMP1", "WT1", "RARRES1", ""
  "DDIT4", "COL1A2", "GPX3", "APOD", "IGFBP7"), order = T)


##---------------------------------------------##
##--------------------HHIP---------------------##
##---------------------------------------------##

df <- sub_emont@meta.data[, c("prediction.score.HHIP", "cell_type", "cell_type2")]
df <- df %>% filter(cell_type2 == "ASPC")

p_v <- ggplot(df, aes(x = cell_type, y = prediction.score.HHIP)) + 
  ggrastr::rasterise(geom_jitter(alpha = 0.1, aes(col = cell_type), width = 0.2), dpi = 550) + 
  geom_violin(scale = "width", fill = NA) + mashaGgplot2Theme +
  scale_color_manual(values = c("hASPC1" = '#cb181d', "hASPC2" = '#238b45', 
                                "hASPC3" = '#fb6a4a',"hASPC4" = '#8da0cb', 
                                "hASPC5" = '#74c476', "hASPC6" = '#1f78b4')) + 
  theme(legend.position = "none")
p_b <- ggplot(df, aes(x = cell_type, y = prediction.score.HHIP)) + 
  ggrastr::rasterise(geom_jitter(alpha = 0.1, aes(col = cell_type), width = 0.2), dpi = 550) + 
  geom_boxplot(outlier.shape = NA, fill = NA) + mashaGgplot2Theme +
  scale_color_manual(values = c("hASPC1" = '#cb181d', "hASPC2" = '#238b45', 
                               "hASPC3" = '#fb6a4a',"hASPC4" = '#8da0cb', 
                               "hASPC5" = '#74c476', "hASPC6" = '#1f78b4')) + 
  theme(legend.position = "none")

ggsave(plot = grid.arrange(p_v, p_b, ncol = 2), 
       "11.Integration-emont/prediction/ViolinPlot_prediction.score.HHIP_inASPCs.pdf",
       width = 5.44, height = 3.66)
df$cell_type
comparisons <- combn(levels(as.factor(df$cell_type)), 2)
comparisons <- t(comparisons)

WilcoxPvalsAll <- apply(comparisons, 1, function(x){
  # pair1 <- strsplit(x[1], '//')[[1]]
  # pair2 <- strsplit(x[2], '//')[[1]]
  testPval <- wilcox.test(df[df$cell_type == x[1], "prediction.score.HHIP"],
                          df[df$cell_type == x[2], "prediction.score.HHIP"])$p.value
  return(c(x, testPval))
})
WilcoxPvalsAll <- as.data.frame(t(WilcoxPvalsAll), stringsAsFactors = F)
colnames(WilcoxPvalsAll) <- c("cat_1", "cat_2",'pValue')

WilcoxPvalsAll <- WilcoxPvalsAll %>% filter(cat_1 =="hASPC4" | cat_2 =="hASPC4" )
WilcoxPvalsAll$p.adjusted <- p.adjust(WilcoxPvalsAll$pValue, method = 'fdr')

f_addstars <- function(d){
  d$stars <- ""
  d$stars[d$p.adjusted < 0.05] <- "*"
  d$stars[d$p.adjusted < 0.01] <- "**"
  d$stars[d$p.adjusted < 0.001] <- "***"
  d$stars[d$p.adjusted < 0.0001] <- "****"
  return(d)
}
WilcoxPvalsAll <- f_addstars(WilcoxPvalsAll)
WilcoxPvalsAll
# cat_1  cat_2 pValue p.adjusted stars
#   1 hASPC1 hASPC4      0          0  ****
#   2 hASPC2 hASPC4      0          0  ****
#   3 hASPC3 hASPC4      0          0  ****
#   4 hASPC4 hASPC5      0          0  ****
#   5 hASPC4 hASPC6      0          0  ****

##---------------------------------------------##
##----------------CILP or SFRP4----------------##
##---------------------------------------------##

p <- FeaturePlot(sub_emont, features = "prediction.score.CILP", order = T, 
            cols = c("gray", "#a6d854" ,"#41ab5d"))
ggsave(plot = p, "11.Integration-emont/emont_int_ASPCs-Meso/Umap_colored_Prediction.score.CILP.pdf",
       height = 4.63, width = 6.4)

DefaultAssay(sub_emont) <- "RNA"
FeaturePlot(sub_emont, features = c("SFRP4", "SFRP2"), order = T)

sub_emont <- AddModuleScore(sub_emont, features = list("CILP" = myDEGs_seurat_f$CILP$geneID), name = "CILP.Module.Score")
p <- FeaturePlot(sub_emont, features = c("CILP.Module.Score1"), order = F, 
            cols = c("gray", "#a6d854" ,"#41ab5d"))
ggsave(plot = p, "11.Integration-emont/emont_int_ASPCs-Meso/Umap_colored_Module.Score.CILP.pdf",
       height = 4.63, width = 6.4)

m <- sub_emont@meta.data
m <- m %>% filter(cell_type2 == "ASPC")

tab <- table(m$predicted.id, m$sample)
tab <- as.data.frame(t(tab)); colnames(tab) <- c("sample", "celltype", "number")

n_cells <- table(m$sample)
tab$ncells <- n_cells[tab$sample]

tab$percent <-  tab$number/tab$ncells*100

bmi_meta <- sub_emont@meta.data[, c("sample", "bmi")]
bmi_meta <- bmi_meta[!duplicated(bmi_meta$sample),]; rownames(bmi_meta) <- bmi_meta$sample

tab$bmi <- bmi_meta[tab$sample, "bmi"]

ggplot(tab, aes(x = bmi, y = percent, col = celltype)) + geom_point() + geom_smooth()

#### 


t <- table(sub_emont$cell_type, sub_emont$sample)

t_igfbp2 <- t["hASPC6",]
#t <- t[t > 0]
df_igfbp2 <- data.frame(sample = names(t_igfbp2), nIGFBP2 = t_igfbp2)
df_igfbp2 <- left_join(bmi_meta, df_igfbp2)
#df_igfbp2 <- df_igfbp2[grep("OAT", rownames(df_igfbp2)),]

df_igfbp2$nCells <- as.vector(table(sub_emont$sample)[df_igfbp2$sample])
df_igfbp2$freq <- df_igfbp2$nIGFBP2/df_igfbp2$nCells*100
df_igfbp2$depot <- sapply(strsplit(df_igfbp2$sample, "_"), `[`, 2)

p <- ggplot(df_igfbp2, aes(x = bmi, 
                      y = freq, 
                      col = depot)) + 
  geom_smooth(method ="lm", fill = "#cab2d6") +
  geom_point(size = 2) +
  mashaGgplot2Theme + 
  xlab("BMI") + ylab("% of IGFBP2+ cells") +
  scale_color_manual(values = c("#6a3d9a", "#FEC010"))
ggsave(plot = p, filename = "11.Integration-emont/emont_int_ASPCs-Meso/BMI-vs-percentIGFBP2.pdf",
       width = 4.15, height = 2.85)
summary(lm(freq~bmi, df_igfbp2 %>% filter(depot == "OAT")))
# Call:
#   lm(formula = freq ~ bmi, data = df_igfbp2 %>% filter(depot == "OAT"))
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.61441 -0.40240  0.08863  0.29194  0.53713 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.56661    0.54878  -2.855   0.0213 *  
#   bmi          0.12813    0.01479   8.663 2.45e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4479 on 8 degrees of freedom
# Multiple R-squared:  0.9037,	Adjusted R-squared:  0.8916 
# F-statistic: 75.04 on 1 and 8 DF,  p-value: 2.452e-05

t <- table(sub_emont$predicted.id, sub_emont$cell_type); t <- t["IGFBP2",]
t <- as.data.frame(t); t$x <-  as.factor(1); t$clust <- rownames(t)
p <- ggplot(t, aes(x = x,y = t, fill = clust)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = c("hASPC1" = '#cb181d', "hASPC2" = '#238b45', "hASPC3" = '#fb6a4a',"hASPC4" = "#8da0cb", 
                               "hASPC5" = '#74c476', "hASPC6" = '#1f78b4',"hMes1" ='#6a51a3',  "hMes2" ='#9e9ac8', "hMes3" ='#cbc9e2')) +
  xlab("IGFBP2") + ylab("# of cells") + mashaGgplot2Theme
ggsave(p,filename = "11.Integration-emont/emont_int_ASPCs-Meso/DistributionPredictedIGFBP2_per_Clust.pdf",
       width = 2.82, height = 3.34)




t_meso <- t[c("hMes1", "hMes2", "hMes3"),]
t_meso <- colSums(t_meso)
#t <- t[t > 0]
df_meso <- data.frame(sample = names(t_meso ), nMeso = t_meso )
df_meso  <- left_join(bmi_meta, df_meso)
#df_igfbp2 <- df_igfbp2[grep("OAT", rownames(df_igfbp2)),]

df_meso$nCells <- as.vector(table(sub_emont$sample)[df_meso$sample])
df_meso$freq <- df_meso$nMeso/df_meso$nCells*100
df_meso$depot <- sapply(strsplit(df_meso$sample, "_"), `[`, 2)

p <- ggplot(df_meso, aes(x = bmi, 
                           y = freq, 
                           col = depot)) + 
  geom_smooth(method ="lm", fill = "#cab2d6") +
  geom_point(size = 2) +
  mashaGgplot2Theme + 
  xlab("BMI") + ylab("% of Mesothelial cells") +
  scale_color_manual(values = c("#6a3d9a", "#FEC010"))
ggsave(plot = p, filename = "11.Integration-emont/emont_int_ASPCs-Meso/BMI-vs-percentIGFBP2.pdf",
       width = 4.15, height = 2.85)
summary(lm(freq~bmi, df_meso %>% filter(depot == "OAT")))
