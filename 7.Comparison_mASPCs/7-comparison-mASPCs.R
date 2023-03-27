################################################################
#                                                              #
#          Comparison of shared markers between pops           #
#     mouse and human integration of Subcutaneous&Visceral     #
#                                                              #
################################################################

### Author: Pernille
### Date: 17.08.2022 - adapted from old script
### Datasets: scRNA-seq Depots: SC - EP - MK - MG - PR - GB 
###                     Patients: B0, B1, B7, L3, L11, L12 (B - Bariatric, L - Lean)
### Goal: Compare markers between human and mouse ASPCs 

library(ggplot2); library(data.table); library(Seurat); library(biomaRt); library(dplyr)

setwd("~/SVRAW1/prainer/hASPCs/PAPER/10X_scRNA-seq/")
source("Utility/General_utils.R")
source("7.Comparison_mASPCs/7-Utils.R")

##---------------------------------------------##
##-----------------Loading data----------------##
##---------------------------------------------##

# integration
int <- readRDS("5.Integration/output/Seurat_2000HVGs.Rds")

# mouse markers
mouseMarkers_sc <- readRDS("~/SVRAW1/prainer/Integration_SingleCell_Adipo_dataset/Subcutaneous_Mouse/Seurat/dim30_SelectedOne/TopMarkers_fisherMethod/AggregatedFisherPval_PerCluster_orderFC.Rds")
names(mouseMarkers_sc) <- c("ASCs", "PreAs", "Aregs", "Endo", "Immune 1", "Immune 2", "Cilp+", "Immune 3")

mouseMarkers <- readRDS("~/SVRAW1/prainer/Integration_Sc_datasets/SC_and_Visc_Mouse/dim60/TopMarkers_fisherMethod/AggregatedFisherPval_PerCluster_orderFC.Rds")

myDEGs_mouse <- readRDS("7.Comparison_mASPCs/ReRun_mASPCs/DEGs_mouse_SC-OM.Rds")

##---------------------------------------------##
##-------------Filter mouse data---------------##
##---------------------------------------------##

## -- Filter on logFC and keep only ASPCs and meso
pop_to_keep <- c("ASCs", "PreAs", "Aregs", "Ifit+", "Cilp+", 
                 "Meso", "VSMPs")

# markers based on sc & om
# mouseMarkers_filter <- lapply(names(mouseMarkers), 
#                               FUN = function(x) return(mouseMarkers[[x]][mouseMarkers[[x]]$avg_overallFC > 0 &
#                                                                          mouseMarkers[[x]]$`Fisher Aggregated p.adj` < 0.05, ]))
# names(mouseMarkers_filter) <- names(mouseMarkers)
# mouseMarkers_filter <- mouseMarkers_filter[pop_to_keep]; 
# rm(mouseMarkers)

myDEGs_mouse_filter <- lapply(names(myDEGs_mouse), FUN = function(x) return(myDEGs_mouse[[x]][myDEGs_mouse[[x]]$avg_logFC_all >0.25 &
                                                                                              myDEGs_mouse[[x]]$minimump_p_val < 0.05, ]))
names(myDEGs_mouse_filter) <- names(myDEGs_mouse)
myDEGs_mouse_filter <- myDEGs_mouse_filter[pop_to_keep]; 
rm(myDEGs_mouse)

# # markers based on sc
# mouseMarkers_sc_filter <- lapply(names(mouseMarkers_sc), FUN = function(x) return(mouseMarkers_sc[[x]][mouseMarkers_sc[[x]]$avg_overallFC >0 &
#                                                                                                        mouseMarkers_sc[[x]]$`Fisher Aggregated p.adj` < 0.05, ]))
# names(mouseMarkers_sc_filter) <- names(mouseMarkers_sc)
# mouseMarkers_sc_filter <- mouseMarkers_sc_filter[pop_to_keep]; 
# rm(mouseMarkers_sc)

##---------------------------------------------##
##---------Convert to human orthologs----------##
##---------------------------------------------##

mouseMarkers_filter_ortho <- lapply(names(myDEGs_mouse_filter), 
                                    function(p) convert_mouse_to_human(myDEGs_mouse_filter[[p]]$geneID))
mouseMarkers_filter_ortho_2 <- lapply(mouseMarkers_filter_ortho, function(d) get_ensembl_names(d, data.annot = data.annot))
names(mouseMarkers_filter_ortho_2) <- names(myDEGs_mouse_filter)
# 182 homolog genes found / 187 
# 93 homolog genes found / 82 
# 60 homolog genes found / 62 
# 120 homolog genes found / 104 
# 55 homolog genes found / 58 
# 805 homolog genes found / 801 
# 735 homolog genes found / 762 
# -----------------------------
# 181 found genes/ 182 
# 93 found genes/ 93 
# 60 found genes/ 60 
# 120 found genes/ 120 
# 54 found genes/ 55 
# 801 found genes/ 805 
# 722 found genes/ 735 

# mouseMarkers_sc_filter_ortho <- lapply(names(mouseMarkers_sc_filter), 
#                                        function(p) convert_mouse_to_human(mouseMarkers_sc_filter[[p]]$geneID))
# mouseMarkers_sc_filter_ortho <- lapply(mouseMarkers_sc_filter_ortho, function(d) get_ensembl_names(d, data.annot = data.annot))
# names(mouseMarkers_sc_filter_ortho) <- names(mouseMarkers_sc_filter)
# mouseMarkers_sc_filter_ortho <- mouseMarkers_sc_filter_ortho[!is.na(names(mouseMarkers_sc_filter_ortho))]
# # 182 found genes/ 183 
# # 83 found genes/ 83 
# # 85 found genes/ 87 
# # 72 found genes/ 73 


##---------------------------------------------##
##---------Calculate % shared markers----------##
##---------------------------------------------##
to_compare <- list(mouse_sc = mouseMarkers_sc_filter_ortho,
                   mouse_sc_om = mouseMarkers_filter_ortho,
                   human = myDEGs_seurat)

## -- Calculate the % of shared genes
percent_markers <- calculate_percent_shared_markers(to_compare, n_top_genes = 100 )
# mouse_sc_om PreAs had only 38
# mouse_sc_om Aregs had only 37
# mouse_sc_om Ifit+ had only 53
# mouse_sc_om Cilp+ had only 37
# 
# mouse_sc PreAs had only 83
# mouse_sc Aregs had only 85
# 
# human Endo had only 98
# human PreAs had only 92
# human HHIP had only 46
# human CHI3L1-2 had only 53
# human res.0.7_16 had only 53


## -- Create the barplots 
p <- list()
for(i in names(percent_markers)){
  n <- strsplit(i, "-")
  p[[i]] <- barplot_shared_markers(percent_markers[i], myColors, 
                                   shift_text_yaxis = 1.5, min_value = 7, n_top_genes = 100 ) + 
    theme_bw()
}


##---------------------------------------------##
##-------Calculate score mouse sc & om---------##
##---------------------------------------------##
n_top = 1000

## Select genes 
all_genes <- c()
for(i in names(mouseMarkers_filter_ortho)){
  n <- min(n_top, nrow(mouseMarkers_filter_ortho[[i]]))
  all_genes <- c(all_genes, mouseMarkers_filter_ortho[[i]]$ens_id[1:n])
}
all_genes <- unique(all_genes)
length(all_genes)
genes_to_int <- all_genes[all_genes %in% rownames(GetAssayData(int, assay = "RNA"))]
length(genes_to_int)

## Scale data
data_scaled <- GetAssayData(int, assay = "RNA")[genes_to_int, ]
data_scaled <- apply(data_scaled, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X)))
data_scaled <- t(data_scaled)

## Calculate score
Scores <- lapply(names(mouseMarkers_filter_ortho), 
                 function(x) calculate_score(data_scaled, 
                                             gene_list = mouseMarkers_filter_ortho[[x]],
                                             n_markers = n_top))
names(Scores) <- names(mouseMarkers_filter_ortho)
Scores_scaled <- lapply(names(mouseMarkers_filter_ortho), 
                 function(x) calculate_score_scaled(data_scaled, 
                                                    gene_list = mouseMarkers_filter_ortho[[x]],
                                                    n_markers = n_top))
names(Scores_scaled) <- names(mouseMarkers_filter_ortho)
# [1] "181/181 markers in data"
# [1] "90/90 markers in data"
# [1] "60/60 markers in data"
# [1] "103/103 markers in data"
# [1] "54/54 markers in data"
# [1] "783/783 markers in data"
# [1] "724/724 markers in data"

## Reshape
Scores_sc_om <- rlist::list.cbind(Scores); colnames(Scores_sc_om) <- names(Scores)
Scores_sc_om_scaled <- rlist::list.cbind(Scores_scaled); colnames(Scores_sc_om_scaled) <- names(Scores_scaled)

write.table(Scores_sc_om, "7.Comparison_mASPCs/Scores_sc-om_top1000genes.txt")
write.table(Scores_sc_om_scaled, "7.Comparison_mASPCs/Scores_sc-om_top1000genes_scaledByNumberOfGenes.txt")

##---------------------------------------------##
##--------------Plot scores tSNE---------------##
##---------------------------------------------##

P_sc_om <- lapply(colnames(Scores_sc_om), function(s) plot_score(s, Scores_sc_om, integrated_data = int )); names(P_sc_om) <- colnames(Scores_sc_om)
P_sc_om_scaled <- lapply(colnames(Scores_sc_om_scaled), 
                         function(s) plot_score(s, Scores_sc_om_scaled, integrated_data = int ))
names(P_sc_om_scaled) <- colnames(Scores_sc_om_scaled)

gridExtra::grid.arrange(P_sc_om$ASCs, P_sc_om$PreAs, P_sc_om$Aregs, P_sc_om$`Ifit+`, P_sc_om$`Cilp+`, P_sc_om$Meso, P_sc_om$VSMPs, nrow = 1, ncol = 7)
gridExtra::grid.arrange(P_sc_om_scaled$ASCs, P_sc_om_scaled$PreAs, P_sc_om_scaled$Aregs, P_sc_om_scaled$`Ifit+`, P_sc_om_scaled$`Cilp+`, 
                        P_sc_om_scaled$Meso, P_sc_om_scaled$VSMPs, nrow = 1, ncol = 7)
"7.Comparison_mASPCs/tSNE_colored_mouse-SC-OM-pops-score.pdf"


##---------------------------------------------##
##----------------Boxplot scores---------------##
##---------------------------------------------##

clust <- int$myIntegratedClustering
clust <- factor(clust, levels = c("ASCs", "PreAs", "IGFBP2", 
                                  "HHIP", "CHI3L1-2", "IFIT", "CILP", "PR specific", 
                                  "Meso", "VSMPs", "Endo", "Immune", 
                                  "res.0.2_8", "res.0.2_6", "5"))
clust <- clust[!clust %in% c("res.0.2_8", "res.0.2_6", "5")]
levels(clust) <- c("ASCs", "PreAs", "IGFBP2", "HHIP", "CHI3L12", "IFIT", "CILP", "PRspecific",
                   "Meso", "VSMPs", "Endo", "Immune", NA, NA, NA)

cols <- myIntegratedColors; names(cols) <- gsub("-", "", names(cols)); names(cols) <- gsub(" ", "", names(cols))
p <- boxplot_score_datatable(score_res = t(Scores_sc_om_scaled[names(clust), ]), clust = clust,
                        clust_name = levels(clust), col = cols, stat_on = F)
stats <- boxplot_score_datatable(score_res = t(Scores_sc_om_scaled[names(clust), ]), clust = clust,
                                      clust_name = levels(clust), col = cols, stat_on = T)

ggsave(gridExtra::grid.arrange(p$ASCs, p$PreAs, p$Aregs, p$`Cilp+`, p$`Ifit+`, p$Meso, p$VSMPs, ncol = 3 ),
      filename = "7.Comparison_mASPCs/Boxplot_scores_stricterFiltering.pdf", height = 8.80, width = 7.28)

stats %>% filter(Gene == "ASCs") %>% filter(start == "ASCs" | end == "ASCs")
# Gene            comp  p adj start        end      Expr
# 1: ASCs      ASCs,PreAs <2e-16  ASCs      PreAs 0.3514231
# 2: ASCs     ASCs,IGFBP2 <2e-16  ASCs     IGFBP2 0.3681575
# 3: ASCs       ASCs,HHIP <2e-16  ASCs       HHIP 0.3848919
# 4: ASCs    ASCs,CHI3L12 <2e-16  ASCs    CHI3L12 0.4016263
# 5: ASCs       ASCs,IFIT <2e-16  ASCs       IFIT 0.4183608
# 6: ASCs       ASCs,CILP <2e-16  ASCs       CILP 0.4350952
# 7: ASCs ASCs,PRspecific <2e-16  ASCs PRspecific 0.4518296
# 8: ASCs       ASCs,Meso <2e-16  ASCs       Meso 0.4685641
# 9: ASCs      ASCs,VSMPs <2e-16  ASCs      VSMPs 0.4852985
# 10: ASCs       ASCs,Endo <2e-16  ASCs       Endo 0.5020329
# 11: ASCs     ASCs,Immune <2e-16  ASCs     Immune 0.5187674

stats %>% filter(Gene == "PreAs") %>% filter(start == "PreAs" | end == "PreAs")
# 1: PreAs    ASCs,PreAs <2e-16    ASCs PreAs 0.3400850
# 2: PreAs  IGFBP2,PreAs <2e-16  IGFBP2 PreAs 0.5182248
# 3: PreAs CHI3L12,PreAs  1e-09 CHI3L12 PreAs 0.5344193
# 4: PreAs    IFIT,PreAs  5e-06    IFIT PreAs 0.5506139
# 5: PreAs    CILP,PreAs <2e-16    CILP PreAs 0.5668084
# 6: PreAs    Meso,PreAs <2e-16    Meso PreAs 0.5830029
# 7: PreAs   PreAs,VSMPs <2e-16   PreAs VSMPs 0.5991974
# 8: PreAs    Endo,PreAs <2e-16    Endo PreAs 0.6153920
# 9: PreAs  Immune,PreAs <2e-16  Immune PreAs 0.6315865

stats %>% filter(Gene == "Aregs") %>% filter(start == "HHIP" | end == "HHIP")
# 1: Aregs       ASCs,HHIP <2e-16    ASCs       HHIP 0.3670814
# 2: Aregs      HHIP,PreAs <2e-16    HHIP      PreAs 0.5107219
# 3: Aregs     HHIP,IGFBP2 <2e-16    HHIP     IGFBP2 0.6543625
# 4: Aregs    CHI3L12,HHIP <2e-16 CHI3L12       HHIP 0.7660829
# 5: Aregs       HHIP,IFIT <2e-16    HHIP       IFIT 0.7820429
# 6: Aregs       CILP,HHIP <2e-16    CILP       HHIP 0.7980030
# 7: Aregs HHIP,PRspecific <2e-16    HHIP PRspecific 0.8139631
# 8: Aregs       HHIP,Meso <2e-16    HHIP       Meso 0.8299231
# 9: Aregs      HHIP,VSMPs <2e-16    HHIP      VSMPs 0.8458832
# 10: Aregs       Endo,HHIP <2e-16    Endo       HHIP 0.8618432
# 11: Aregs     HHIP,Immune <2e-16    HHIP     Immune 0.8778033

stats %>% filter(Gene == "Ifit+") %>% filter(start == "IFIT" | end == "IFIT")
# 1: Ifit+       ASCs,IFIT <2e-16    ASCs       IFIT 0.4373585
# 2: Ifit+      IFIT,PreAs <2e-16    IFIT      PreAs 0.5831447
# 3: Ifit+     IFIT,IGFBP2 <2e-16    IFIT     IGFBP2 0.7107076
# 4: Ifit+       HHIP,IFIT <2e-16    HHIP       IFIT 0.8200472
# 5: Ifit+    CHI3L12,IFIT <2e-16 CHI3L12       IFIT 0.9111636
# 6: Ifit+       CILP,IFIT <2e-16    CILP       IFIT 1.0387265
# 7: Ifit+ IFIT,PRspecific <2e-16    IFIT PRspecific 1.0569498
# 8: Ifit+       IFIT,Meso <2e-16    IFIT       Meso 1.0751730
# 9: Ifit+      IFIT,VSMPs <2e-16    IFIT      VSMPs 1.0933963
# 10: Ifit+       Endo,IFIT <2e-16    Endo       IFIT 1.1116196
# 11: Ifit+     IFIT,Immune <2e-16    IFIT     Immune 1.1298429

stats %>% filter(Gene == "Cilp+") %>% filter(start == "CILP" | end == "CILP")
# 1: Cilp+       ASCs,CILP <2e-16    ASCs       CILP 0.5236909
# 2: Cilp+      CILP,PreAs <2e-16    CILP      PreAs 0.7049685
# 3: Cilp+     CILP,IGFBP2 <2e-16    CILP     IGFBP2 0.8862462
# 4: Cilp+       CILP,HHIP  1e-10    CILP       HHIP 1.0473818
# 5: Cilp+    CHI3L12,CILP <2e-16 CHI3L12       CILP 1.1883755
# 6: Cilp+       CILP,IFIT <2e-16    CILP       IFIT 1.2890853
# 7: Cilp+ CILP,PRspecific <2e-16    CILP PRspecific 1.4099371
# 8: Cilp+       CILP,Meso <2e-16    CILP       Meso 1.4300790
# 9: Cilp+      CILP,VSMPs <2e-16    CILP      VSMPs 1.4502210
# 10: Cilp+       CILP,Endo <2e-16    CILP       Endo 1.4703630
# 11: Cilp+     CILP,Immune <2e-16    CILP     Immune 1.4905049

stats %>% filter(Gene == "Meso") %>% filter(start == "Meso" | end == "Meso")
# 1: Meso       ASCs,Meso <2e-16    ASCs       Meso 0.2019438
# 2: Meso      Meso,PreAs <2e-16    Meso      PreAs 0.2860870
# 3: Meso     IGFBP2,Meso <2e-16  IGFBP2       Meso 0.3197443
# 4: Meso       HHIP,Meso <2e-16    HHIP       Meso 0.3786446
# 5: Meso    CHI3L12,Meso <2e-16 CHI3L12       Meso 0.4207162
# 6: Meso       IFIT,Meso <2e-16    IFIT       Meso 0.4627878
# 7: Meso       CILP,Meso <2e-16    CILP       Meso 0.5048594
# 8: Meso Meso,PRspecific <2e-16    Meso PRspecific 0.5385167
# 9: Meso      Meso,VSMPs <2e-16    Meso      VSMPs 0.5637597
# 10: Meso       Endo,Meso <2e-16    Endo       Meso 0.5721740
# 11: Meso     Immune,Meso <2e-16  Immune       Meso 0.5805884

stats %>% filter(Gene == "VSMPs") %>% filter(start == "VSMPs" | end == "VSMPs")
# 1: VSMPs       ASCs,VSMPs <2e-16       ASCs VSMPs 0.3395169
# 2: VSMPs      PreAs,VSMPs <2e-16      PreAs VSMPs 0.4526892
# 3: VSMPs     IGFBP2,VSMPs <2e-16     IGFBP2 VSMPs 0.5532868
# 4: VSMPs       HHIP,VSMPs <2e-16       HHIP VSMPs 0.6161603
# 5: VSMPs    CHI3L12,VSMPs <2e-16    CHI3L12 VSMPs 0.7041832
# 6: VSMPs       IFIT,VSMPs <2e-16       IFIT VSMPs 0.7544820
# 7: VSMPs       CILP,VSMPs <2e-16       CILP VSMPs 0.8173555
# 8: VSMPs PRspecific,VSMPs <2e-16 PRspecific VSMPs 0.8676543
# 9: VSMPs       Meso,VSMPs <2e-16       Meso VSMPs 0.9053784
# 10: VSMPs       Endo,VSMPs <2e-16       Endo VSMPs 0.9431024
# 11: VSMPs     Immune,VSMPs <2e-16     Immune VSMPs 0.9556771


##---------------------------------------------##
##----------------Score  human to mouse---------------##
##---------------------------------------------##






