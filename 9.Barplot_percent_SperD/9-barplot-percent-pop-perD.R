################################################################
#                                                              #
#      Barplot of % of cells per depots in each cluster        #
#                                                              #
################################################################

### Author: Pernille
### Date: 17.08.2022 - adapted from old script
### Datasets: scRNA-seq Depots: SC - EP - MK - MG - PR - GB 
###                     Patients: B0, B1, B7, L3, L11, L12 (B - Bariatric, L - Lean)
### Goal: Compare markers between human and mouse ASPCs 

library(ggplot2); library(data.table); library(Seurat); library(biomaRt)

setwd("~/SVRAW1/prainer/hASPCs/PAPER/10X_scRNA-seq/")
source("Utility/General_utils.R")
source("6.TopDEGs/6-utils.R")

##---------------------------------------------##
##-----------------Loading data----------------##
##---------------------------------------------##
int <- readRDS("5.Integration/output/Seurat_2000HVGs.Rds")

myDepots <- myDepots[-which(myDepots == "GB7")]

##---------------------------------------------##
##--------------Calculate percent--------------##
##---------------------------------------------##
df <- table(int$myIntegratedClustering,int$batch)
df_percent <- reshape2::melt(t(df/rowSums(df)*100))
colnames(df_percent) <- c("sample", "cluster", "percent")

df_percent$percent <- round(df_percent$percent, 2)
df_percent$sample <- factor(as.character(df_percent$sample),
                            levels = myDepots)
df_percent$cluster <- factor(as.character(df_percent$cluster),
                             levels = levels(int$myIntegratedClustering))

##---------------------------------------------##
##-------------------Plot %--------------------##
##---------------------------------------------##
calculate_label_ypos <- function(df, xaxis = "cluster", fill = "sample", values = "percent",
                                 min_value = 6, shift_text = 2){
  
  #rename
  df <- rename(df, "xaxis" = xaxis, "fill" = fill, "values" = values)
  df <- df %>% group_by(xaxis) %>% mutate(y_pos = rev(values)) %>% mutate(y_pos = rev(cumsum(y_pos)))
  
  df$y_pos <- df$y_pos - shift_text
  df$y_pos[df$values < min_value | df$y_pos < 0] <- NA
  
  return(df$y_pos)  
}

df_percent$label_ypos <- calculate_label_ypos(df_percent)

p <- ggplot(df_percent, aes(x = cluster, y = percent, fill = sample)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = unname(myBatchColors[levels(df_percent$sample)])) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_text(aes(y = label_ypos, label = percent), col = "white",fontface = "bold") + 
  mashaGgplot2Theme 

ggsave(p, filename = "9.Barplot_percent_SperD/Barplot_brown.pdf", height = 6.75, width = 9.24)

##---------------------------------------------##
##-----------------Plot n cells----------------##
##---------------------------------------------##
df <- as.data.frame(table(int$myIntegratedClustering))
colnames(df) <- c("pop", "n_cells")
df$pop <- as.character(df$pop)
df <- df[order(df$n_cells, decreasing = T),]
df$pop <- factor(as.character(df$pop), levels = rev(unique(df$pop)))
p <- ggplot(df, aes(y = pop, x = n_cells, fill = pop)) + 
  geom_bar(stat = "identity") + mashaGgplot2Theme + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = myIntegratedColors) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "b")  
ggsave(plot = p, filename = "9.Barplot_percent_SperD/Barplot_numberCells-per-Clust.pdf",
       width = 2.75, height = 4.11)

##---------------------------------------------##
##-----------------Plot n cells----------------##
##---------------------------------------------##
df <- table(int$myIntegratedClustering, int$batch)

df <- reshape2::melt(df)
colnames(df) <- c("pop", "batch", "n_cells")
df$depot <- str_sub(df$batch, 1, 2)
df$depot[df$depot == "MK"] <- "MG"
df <- df %>% group_by(pop) %>% mutate(sum_pop = sum(n_cells))
df <- df[order(df$sum_pop, decreasing = F),]
df$pop <- factor(as.character(df$pop), levels = unique(df$pop))
P <- lapply(c("SC", "EP", "MG", "PR"), function(d){
  data <- df %>% filter(depot == d)
  data <- data %>% group_by(pop) %>% mutate(sum_pop_depot = sum(n_cells))
  data_plot1 <- data[!duplicated(data$pop),]
  p_log <- ggplot(data_plot1, aes(y = pop, x = sum_pop_depot)) + 
    geom_bar(stat = "identity") + mashaGgplot2Theme + 
    theme(legend.position = "none") + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks(sides = "b") +
    ylab("") + xlab("") + theme(axis.text.y = element_blank())
  
  p_forProp <- ggplot(data, aes(y = pop, x = n_cells, fill = batch)) + 
    geom_bar(stat = "identity") + mashaGgplot2Theme + 
    theme(legend.position = "none") + 
    scale_fill_manual(values = myBatchColors) +
    ylab("") + xlab("") + theme(axis.text.y = element_blank())
  
  return(list(p_log = p_log, p_forProp = p_forProp))
})
names(P) <- c("SC", "EP", "MG", "PR")

ggsave(plot = marrangeGrob(grobs = list(p, P$SC$p_log, P$PR$p_log, P$EP$p_log, P$MG$p_log), 
                           layout_matrix = matrix(rep(c(1,1,1,2,2,3,3,4,4,5,5), 6), ncol = 11, byrow = T)),
      "9.Barplot_percent_SperD/Barplot_numberCells-per-Clust_dividedPerDepots.pdf", 
      width = 9.14, height = 4.11)

ggsave(plot = marrangeGrob(grobs = list(P$SC$p_forProp, P$PR$p_forProp, P$EP$p_forProp, P$MG$p_forProp), ncol = 4, nrow = 1),
       "9.Barplot_percent_SperD/Barplot_nCells-per-Clust_dividedPerDepots_FOR-PROP.pdf", 
       width = 6.64, height = 4.11)

#---------------------------------------------##
##----------Dot plot N cells per pops----------##
##---------------------------------------------##

t <- table(int$myIntegratedClustering)
t <- t[names(t) != "res.0.2_8"]
t <- as.data.frame(t)
t$Var1 <- factor(as.character(t$Var1),
                 levels = c("ASCs", "PreAs", "HHIP", "IFIT", "CILP", "CHI3L1-2", "PR specific", "IGFBP2",
                            "Meso", "VSMPs", "Endo", "Immune"))
t <- t[order(t$Var1)]
t$perc <- t$Freq/sum(t$Freq)
t$x <- as.factor("x")
t$y <- paste(t$Var1, "-", t$Freq)
t$y <- factor(t$y, levels = rev(t$y))
p <- ggplot(t, aes(x = x, y = y, size = perc, color = Var1 )) + 
  geom_point() + 
  mashaGgplot2Theme + 
  #scale_size(range = c(3,13)) +
  scale_color_manual(values = myIntegratedColors)
p
ggsave(p, filename = "9.Barplot_percent_SperD/DotPlot_Percent-of-cells.pdf",
       width = 3.54, height = 3)


##-----------------------------------------------------##
##-% cells corrected for # cells/batch & # batch/depot-##
##-----------------------------------------------------##

df <- table(int$myIntegratedClustering, int$batch)

# ## VERSION 1 - ALL
# df <- df[rownames(df) != "res.0.2_8",]
# ## VERSION 2 - without Immune, Endo
df <- df[!rownames(df) %in% c("Endo", "Immune","res.0.2_8"),]
# ## VERSION 3 - without Immune, Endo, VSMPs
# df <- df[!rownames(df) %in% c("Endo", "Immune", "VSMPs","res.0.2_8"),]
# ## VERSION 4 - only ASPCs, and IGFBP2
# df <- df[c("ASCs", "PreAs", "CILP", "IFIT", "HHIP", "CHI3L1-2", "PR specific", "IGFBP2"),]
## VERSION 5 - onlyASPCs
df <- df[c("ASCs", "PreAs", "CILP", "IFIT", "HHIP", "CHI3L1-2", "PR specific", "IGFBP2"),]

# 1. Correct by the # of cells in a batch
df <- t(t(df)/colSums(df))

# 2. Reshape and add depot info
df_percent <- reshape2::melt(df)
colnames(df_percent) <- c("cluster", "sample", "freq")
df_percent$depot <- stringr::str_sub(df_percent$sample, 1,2)
df_percent$depot[df_percent$depot %in% c("MK", "MG")] <- "MC"

# 3. Correct by the number of samples per depot (3 for SC, OM, PR and 2 for MC)
df_percent$freq_corrected <- df_percent$freq
df_percent$freq_corrected[df_percent$depot %in% c("EP", "SC", "PR")] <- df_percent$freq_corrected[df_percent$depot %in% c("EP", "SC", "PR")] *11/3
df_percent$freq_corrected[df_percent$depot %in% c("MC")] <- df_percent$freq_corrected[df_percent$depot %in% c("MC")] *11/2

# 4. Convert in % per pop 
df_percent <- df_percent %>% group_by(cluster) %>% mutate(sum_freq_corr = sum(freq_corrected)) 
df_percent$percent_corrected <- round(df_percent$freq_corrected/df_percent$sum_freq_corr*100,2)

# 5. Prepare for plot
df_percent$sample <- factor(as.character(df_percent$sample),
                            levels = rev(c("SC0", "SC1", "SC7",
                                           "PR3", "PR11", "PR12", 
                                           "EP0", "EP1", "EP7", 
                                           "MG7", "MK7")))
df_percent$cluster <- factor(as.character(df_percent$cluster),
                             levels = c("ASCs", "PreAs", "HHIP", "IFIT", "CILP", "CHI3L1-2", "PR specific", "IGFBP2"))
                                       # "Meso", "VSMPs"))
df_percent <- df_percent[order(df_percent$sample),]
df_percent$label_ypos_batch <- calculate_label_ypos(df_percent, values = "percent_corrected")


# 6. Plot
p <- ggplot(df_percent, aes(x = cluster, y = percent_corrected, fill = sample)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = unname(myBatchColors[levels(df_percent$sample)])) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_text(aes(y = label_ypos_batch, label = percent_corrected), col = "white",fontface = "bold") + 
  mashaGgplot2Theme 
p
ggsave(plot = p, "9.Barplot_percent_SperD/Barplot_PercentCellsPerPop_correctedFornCellsPerDatasets_normalizedFornBatch_leg.pdf",
       width = 5.72, height = 4.71)
ggsave(plot = p + theme(legend.position = "none"), "9.Barplot_percent_SperD/Barplot_PercentCellsPerPop_correctedFornCellsPerDatasets_normalizedFornBatch_NOleg.pdf",
       width = 4.99, height = 4.71)

# 7. Calculate % for each depot
df_percent <- df_percent %>% group_by(depot, cluster) %>% mutate(percent_corrected_depot = sum(percent_corrected))
df_percent_for_depot <- df_percent %>% filter(sample %in% c("SC1", "EP1", "PR11", "MG7"))
df_percent_for_depot$depot <- factor(as.character(df_percent_for_depot$depot), levels = rev(c("SC", "PR", "EP","MC")))
df_percent_for_depot <- df_percent_for_depot[order(df_percent_for_depot$depot),]
df_percent_for_depot$label_ypos_depot <- calculate_label_ypos(df_percent_for_depot, values = "percent_corrected_depot", fill = "depot")

names(myDepotsColors) <- c("EP", "MC", "SC", "PR")
p <- ggplot(df_percent_for_depot, aes(x = cluster, y = percent_corrected_depot, fill = depot)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = unname(myDepotsColors[levels(df_percent_for_depot$depot)])) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_text(aes(y = label_ypos_depot, label = percent_corrected_depot), col = "white",fontface = "bold") + 
  mashaGgplot2Theme 
p
ggsave(plot = p, "9.Barplot_percent_SperD/Barplot_PercentCellsPerPop_correctedFornCellsPerDatasets_normalizedFornBatch_colDep_leg.pdf",
       width = 5.72, height = 4.71)
ggsave(plot = p + theme(legend.position = "none"), "9.Barplot_percent_SperD/Barplot_PercentCellsPerPop_correctedFornCellsPerDatasets_normalizedFornBatch_colDep_NOleg.pdf",
       width = 4.99, height = 4.71)

ggsave(plot = p, "9.Barplot_percent_SperD/Barplot_PercentCellsPerPop_correctedFornCellsPerDatasets_onlyASPCs.pdf",
       width = 7, height = 5.3)
