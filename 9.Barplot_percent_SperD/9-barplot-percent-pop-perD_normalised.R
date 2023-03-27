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
