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
### Goal: 

library(ggplot2); library(data.table)

##---------------------------------------------##
##-----------------Loading data----------------##
##---------------------------------------------##
myColors <- c("ASCs" = "#33A02C", "PreAs" = "#E31A1C", "IGFBP2" = "#284724", "Meso" = "#810F7C", 
              "VSMPs" = "#FC8D62", "Endo" = "darkgoldenrod1", "Immune" = "#E889BD")

data <- readRDS("2.Comparison_topDEGs_perD-P/Barplot_perComparison/Common/percent_markers_common.Rds")

##---------------------------------------------##
##-----------------Prepare data----------------##
##---------------------------------------------##

# add comparison info
for(c in names(data)){
  data[[c]]$comparison <- c
}

# unlist
data <- data.table::rbindlist(data)

# Calculate mean and std
dataSum <- data[, .(M = mean(percent, na.rm = T), S = sd(percent, na.rm = T)),
                by = .(popRef, popToCompare)]

# Filter out common & unknown VSMPs 
data <- data %>% filter(!popRef %in% c( "Unknown_VSMPs")) %>% filter(!popToCompare %in% c("common", "Unknown_VSMPs"))
dataSum <- dataSum %>% filter(!popRef %in% c( "Unknown_VSMPs")) %>% filter(!popToCompare %in% c("common", "Unknown_VSMPs"))

# Order
data$popRef <- factor(data$popRef, levels = names(myColors))
data$popToCompare <- factor(data$popToCompare, levels = names(myColors))
dataSum$popToCompare <- factor(dataSum$popToCompare, levels = names(myColors))
dataSum$popRef <- factor(dataSum$popRef, levels = names(myColors))

##---------------------------------------------##
##---------------------Plot--------------------##
##---------------------------------------------##
myColors <- myColors[levels(data$popToCompare)]

thePlot <- ggplot(dataSum, aes(x = popRef, y = M, fill = popToCompare)) + 
  geom_bar(stat = "identity", color = "black", 
           position = position_dodge(.9), alpha = 0.6) +
  geom_errorbar(aes(ymin = M, ymax = M + S), width = .2,
                position = position_dodge(.9)) +
  geom_jitter(data = data, 
              mapping = aes(x = popRef, 
                            y = percent,
                            col = popToCompare, 
                            group = popToCompare),
              size = 0.8, 
              position = position_jitterdodge(dodge.width = 0.9, 
                                              jitter.width = 0.3)) +
  scale_fill_manual(values = myColors) +
  scale_color_manual(values = myColors) + mashaGgplot2Theme +
  ylab(label = "% shared top 100 markers") + xlab("") 
  
ggsave(thePlot, filename = "Summary_Barplot_percent-top-100-shared-markers.pdf")
#Saving 11.2 x 5.76 in image

ggsave(thePlot + theme(legend.position = "none"), 
       filename = "Summary_Barplot_percent-top-100-shared-markers_woLeg.pdf", 
       width = 11.25, height = 4.11)

##---------------------------------------------##
##---------------------Plot--------------------##
##---------------------------------------------##

head(data)
data$s1 <- substr(sapply(strsplit(data$comparison,"-"), `[`, 1), start = 1, stop = 2); data$s1[data$s1 == "MK"] <- "MG"
data$s2 <- substr(sapply(strsplit(data$comparison,"-"), `[`, 2), start = 1, stop = 2); data$s2[data$s2 == "MK"] <- "MG"

f_comp_type <- function(s1,s2){
  if(s1 == s2){
    return("same")
  }else{
    if(s1 %in% c("EP", "SC") & s2 %in% c("EP", "SC")){ return("EP-SC")}
    else if(s1 %in% c("EP", "MG") & s2 %in% c("EP", "MG")){ return("EP-MG")}
    else if(s1 %in% c("EP", "PR") & s2 %in% c("EP", "PR")){ return("EP-PR")}
    else if(s1 %in% c("SC", "MG") & s2 %in% c("SC", "MG")){ return("MG-SC")}
    else if(s1 %in% c("SC", "PR") & s2 %in% c("SC", "PR")){ return("PR-SC")}
    else if(s1 %in% c("MG", "PR") & s2 %in% c("MG", "PR")){ return("PR-MG")}
    else{return("AW")}
  }
}

f_comp_type <- function(s1,s2){
  if(s1 %in% c("EP") & s2 %in% c("EP")){ return("EP")}
    if(s1 %in% c("SC") & s2 %in% c("SC")){ return("SC")}
    else if(s1 %in% c("MG") & s2 %in% c("MG")){ return("MG")}
    else if(s1 %in% c("PR") & s2 %in% c("PR")){ return("PR")}
    else if(s1 %in% c("EP", "SC") & s2 %in% c("EP", "SC")){ return("EP-SC")}
    else if(s1 %in% c("EP", "MG") & s2 %in% c("EP", "MG")){ return("EP-MG")}
    else if(s1 %in% c("EP", "PR") & s2 %in% c("EP", "PR")){ return("EP-PR")}
    else if(s1 %in% c("SC", "MG") & s2 %in% c("SC", "MG")){ return("MG-SC")}
    else if(s1 %in% c("SC", "PR") & s2 %in% c("SC", "PR")){ return("PR-SC")}
    else if(s1 %in% c("MG", "PR") & s2 %in% c("MG", "PR")){ return("PR-MG")}
    else{return("AW")}
  }

t <- unlist(apply(data, 1, function(x) f_comp_type(x["s1"], x["s2"])))
data$comp_type_all <- t

data$ind1 <- substr(sapply(strsplit(data$comparison,"-"), `[`, 1), start = 3, stop = 10000)
data$ind2 <- substr(sapply(strsplit(data$comparison,"-"), `[`, 2), start = 3, stop = 10000)
f_comp_type <- function(ind1,ind2){
  if(ind1 == ind2){
    return("same")
  }else{
    return("not")
  }
}
t <- unlist(apply(data, 1, function(x) f_comp_type(x["ind1"], x["ind2"])))
data$ind_comp <- t


d <- data %>% filter(popRef == popToCompare)
d <- d %>% filter(popRef %in% c("ASCs", "PreAs", "VSMPs"))
d$comp_type <- factor(d$comp_type, levels = c("same", "PR-MG", "EP-MG", "MG-SC", 
                                              "EP-PR", "PR-SC", "EP-SC"))
d$comp_type_all <- factor(d$comp_type_all, levels = c("EP", "SC", "MG", "PR", "PR-MG", "EP-MG", "MG-SC", 
                                              "EP-PR", "PR-SC", "EP-SC"))

p <- ggplot(d, aes(y = percent, x = popRef, fill = comp_type_all)) + 
  geom_boxplot(alpha = 0.5) +
  geom_jitter(aes( col = comp_type_all),
              position = position_jitterdodge(dodge.width = 0.75, 
                                              jitter.width = 0.1)) + 
  scale_fill_manual(values = c('#bae4b3','#74c476','#31a354','#006d2c','#ffffb2','#fed976','#feb24c','#fd8d3c','#f03b20','#bd0026')) + 
  scale_color_manual(values = c('#bae4b3','#74c476','#31a354','#006d2c','#ffffb2','#fed976','#feb24c','#fd8d3c','#f03b20','#bd0026')) + 
  mashaGgplot2Theme

ggsave(p, filename = "2.Comparison_topDEGs_perD-P/Percent_comparison_type.pdf")
#Saving 12 x 3.64 in image

p <- ggplot(d, aes(y = percent, x = popRef, fill = ind_comp)) + 
  geom_boxplot(alpha = 0.5) +
  geom_jitter(aes( col = ind_comp),
              position = position_jitterdodge(dodge.width = 0.75, 
                                              jitter.width = 0.1)) + 
  scale_fill_manual(values = c('#006d2c','#bd0026')) + 
  scale_color_manual(values = c('#006d2c','#bd0026')) + 
  mashaGgplot2Theme
ggsave(p, filename = "2.Comparison_topDEGs_perD-P/Percent_comparison_ind.pdf")
#Saving 7.32 x 3.64 in image


