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
### Goal: summary barplot of scmap results

library(ggplot2); library(data.table)

setwd("~/SVRAW1/prainer/hASPCs/PAPER/10X_scRNA-seq/3.scmap/")

##---------------------------------------------##
##-----------------Loading data----------------##
##---------------------------------------------##
myColors <- c("ASCs" = "#33A02C", "PreAs" = "#E31A1C", "IGFBP2" = "#284724", "Meso" = "#810F7C", 
              "VSMPs" = "#FC8D62", "Endo" = "darkgoldenrod1", "Immune" = "#E889BD")

data <- readRDS("output_table/projection_results_in-percent.rds")

##---------------------------------------------##
##-----------------Prepare data----------------##
##---------------------------------------------##

# unlist
data <- data.table::rbindlist(data)

# Does the pop exist in the ref
data$exist <- T 
data <- data %>% mutate(exist=replace(exist,
                                      (projected %in% c("Meso", "IGFBP2")) & (!dataset_ref %in% c("EP0", "EP1", "EP7")),
                                      values = F))
data <- data %>% mutate(exist=replace(exist,
                                      (projected == "Endo") & (!dataset_ref %in% c("EP7", "SC7", "MK7", "MG7")),
                                      values = F))
data <- data %>% mutate(exist=replace(exist,
                                      (projected == "Immune") & (!dataset_ref %in% c("EP7", "SC7", "MG7")),
                                      values = F))
data <- data %>% mutate(exist=replace(exist,
                                      (projected %in% c("IGFBP2", "Immune",
                                                        "Meso", "VSMPs", 
                                                        "Endo", "Immune")) & (dataset_ref == c("GB7")),
                                      values = F))

# Calculate mean and std
datasum <- data %>% filter(exist == T)
dataSum <- datasum[, .(M = mean(value, na.rm = T), S = sd(value, na.rm = T)),
                   by = .(projection_on_ref, projected)]
to_add <- as.data.table(data.frame(projection_on_ref = rep(c("Endo", "Immune"),2),
                                   projected = rep(c("IGFBP2", "Meso"), each = 2), 
                                   M = rep(0,4), S = rep(0,4)))
dataSum <- rbindlist(list(dataSum, to_add))

# Filter out common & unknown VSMPs 
data <- data %>% filter(projected != "Unknown_VSMPs") %>% filter(projection_on_ref != "Unknown_VSMPs")
dataSum <- dataSum %>% filter(projected != "Unknown_VSMPs") %>% filter(projection_on_ref != "Unknown_VSMPs")

# Order
data$projected <- factor(as.character(data$projected), levels = names(myColors))
data$projection_on_ref <- factor(as.character(data$projection_on_ref), levels = c(names(myColors), "unassigned"))
dataSum$projected <- factor(as.character(dataSum$projected), levels = names(myColors))
dataSum$projection_on_ref <- factor(as.character(dataSum$projection_on_ref), levels = c(names(myColors), "unassigned"))

##---------------------------------------------##
##---------------------Plot--------------------##
##---------------------------------------------##
myColors <- c(myColors, c("unassigned" = "gray"))

thePlot <- ggplot(dataSum, aes(x = projected, y = M, fill = projection_on_ref)) + 
  geom_bar(stat = "identity", color = "black", 
           position = position_dodge(.9), alpha = 0.6) +
  geom_errorbar(aes(ymin = M, ymax = M + S), width = .2,
                position = position_dodge(.9)) +
  geom_jitter(data = data, 
              mapping = aes(x = projected, 
                            y = value,
                            col = projection_on_ref, 
                            group = projection_on_ref,
                            shape = exist,
                            alpha = exist),
              size = 0.8, 
              position = position_jitterdodge(dodge.width = 0.9, 
                                              jitter.width = 0.3)) +
  scale_alpha_manual(values = c(0.3, 1))+
  scale_shape_manual(values = c(1,19))+
  scale_fill_manual(values = myColors) +
  scale_color_manual(values = myColors) + mashaGgplot2Theme +
  ylab(label = "%") + xlab("") 
# 
ggsave(plot = thePlot, filename = "Summary_Barplot_percent_projected_cells.pdf")
#Saving 11.2 x 5.76 in image

ggsave(thePlot + theme(legend.position = "none"),
       filename = "Summary_Barplot_percent_projected_cells_woLeg.pdf",
       width = 11.25, height = 4.11)
