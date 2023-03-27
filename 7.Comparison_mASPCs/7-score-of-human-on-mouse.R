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
# my human markers
myDEGs_human <- readRDS("6.TopDEGs/DEGs.Rds")
myDEGs_human_filter <- lapply(names(myDEGs_human), FUN = function(x) return(myDEGs_human[[x]][myDEGs_human[[x]]$avg_logFC_all >0.25 &
                                                                                                   myDEGs_human[[x]]$minimump_p_val < 0.05, ]))
names(myDEGs_human_filter) <- names(myDEGs_human)
rm(myDEGs_human)

# data.annot_mouse
data.annot_mouse <- read.table("~/SVRAW1/prainer/Files/Mouse/data.annot/Mus_musculus.GRCm38.100_data.annot.txt", header =T)
rownames(data.annot_mouse) <- data.annot_mouse$ens_id <- data.annot_mouse$Ensembl
data.annot_mouse <- data.annot_mouse[, c("ens_id", "Name")]
colnames(data.annot_mouse) <- c("ens_id", "gene_short_name")

# integration mouse SC OM 
int_mASPCs <- readRDS("~/SVRAW1/prainer/Integration_SingleCell_Adipo_dataset/SC_and_Visc_Mouse/Seurat/Sc_Visc_integrated_StdWorkflow_seurat.Rds")
int_mASPCs <- int_mASPCs$`60`

##---------------------------------------------##
##-----------------Loading data----------------##
##---------------------------------------------##

humanMarkers_ortho <- lapply(names(myDEGs_human_filter), 
                             function(p) convert_human_to_mouse(myDEGs_human_filter[[p]]$geneID))
humanMarkers_ortho <- lapply(humanMarkers_ortho, function(d) get_ensembl_names(d, data.annot = data.annot_mouse))
names(humanMarkers_ortho) <- names(myDEGs_human_filter)
# 144 homolog genes found / 148 
# 40 homolog genes found / 40 
# 100 homolog genes found / 101 
# 418 homolog genes found / 475 
# 215 homolog genes found / 227 
# 11 homolog genes found / 17 
# 11 homolog genes found / 12 
# 29 homolog genes found / 35 
# 37 homolog genes found / 41 
# 77 homolog genes found / 91 
# 622 homolog genes found / 671 
# 34 homolog genes found / 42 
# 31 homolog genes found / 32 

# 142 found genes/ 144 
# 40 found genes/ 40 
# 99 found genes/ 100 
# 418 found genes/ 418 
# 214 found genes/ 215 
# 11 found genes/ 11 
# 11 found genes/ 11 
# 29 found genes/ 29 
# 37 found genes/ 37 
# 76 found genes/ 77 
# 621 found genes/ 622 
# 34 found genes/ 34 
# 31 found genes/ 31 

##---------------------------------------------##
##----------Calculate score mouse sc-----------##
##---------------------------------------------##
n_top <- 1000

## Select genes 
all_genes <- c()
for(i in names(humanMarkers_ortho)){
  n <- min(n_top, nrow(humanMarkers_ortho[[i]]))
  all_genes <- c(all_genes, humanMarkers_ortho[[i]]$ens_id[1:n])
}
all_genes <- unique(all_genes)
length(all_genes)
genes_to_int <- all_genes[all_genes %in% rownames(GetAssayData(int_mASPCs, assay = "RNA"))]
length(genes_to_int)

## Scale data
data_scaled <- GetAssayData(int_mASPCs, assay = "RNA")[genes_to_int,]
data_scaled <- apply(data_scaled, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X)))
data_scaled <- t(data_scaled)

## Calculate score
Scores <- lapply(names(humanMarkers_ortho), 
                 function(x) calculate_score(data_scaled, 
                                             gene_list = humanMarkers_ortho[[x]],
                                             n_markers = n_top))
names(Scores) <- names(humanMarkers_ortho)
Scores_scaled <- lapply(names(humanMarkers_ortho), 
                        function(x) calculate_score_scaled(data_scaled, 
                                                           gene_list = humanMarkers_ortho[[x]],
                                                           n_markers = n_top))
names(Scores_scaled) <- names(humanMarkers_ortho)
# [1] "141/142 markers in data"
# [1] "39/40 markers in data"
# [1] "98/99 markers in data"
# [1] "408/416 markers in data"
# [1] "214/214 markers in data"
# [1] "10/11 markers in data"
# [1] "11/11 markers in data"
# [1] "29/29 markers in data"
# [1] "37/37 markers in data"
# [1] "75/76 markers in data"
# [1] "615/623 markers in data"
# [1] "34/34 markers in data"
# [1] "31/31 markers in data"

## Reshape
Scores_mouse <- rlist::list.cbind(Scores); colnames(Scores_mouse) <- names(Scores)
Scores_scaled_mouse <- rlist::list.cbind(Scores_scaled); colnames(Scores_scaled_mouse) <- names(Scores_scaled)

saveRDS(Scores_mouse, "7.Comparison_mASPCs/score_ofHuman_in_mouse.Rds")
saveRDS(Scores_scaled_mouse, "7.Comparison_mASPCs/score_ofHuman_in_mouse_scaledByNumberOfGenes.Rds")

##---------------------------------------------##
##--------------Plot scores tSNE---------------##
##---------------------------------------------##

P <- lapply(colnames(Scores_mouse), function(s) plot_score(s, Scores_mouse, integrated_data = int_mASPCs)); names(P) <- colnames(Scores_mouse)
P_scaled <- lapply(colnames(Scores_scaled_mouse), 
                         function(s) plot_score(s, Scores_scaled_mouse, integrated_data = int_mASPCs))
names(P_scaled) <- colnames(Scores_scaled_mouse)

ggsave(gridExtra::grid.arrange(P_scaled$ASCs, P_scaled$PreAs, P_scaled$IGFBP2, 
                        P_scaled$HHIP, P_scaled$`CHI3L1-2`, P_scaled$IFIT, 
                        P_scaled$CILP, P_scaled$Meso, P_scaled$VSMPs,
                         ncol = 3),
       filename = "7.Comparison_mASPCs/tSNE_mouseIntegration_coloredHumanBasedScores.pdf",
       height = 7.18, width = 9.43)

##---------------------------------------------##
##----------------Boxplot scores---------------##
##---------------------------------------------##

clust <- int_mASPCs$myIntegratedClustering
clust <- factor(clust, levels = c("ASCs", "PreAs", "Aregs", "Ifit+", "Cilp+", 
                                 "Meso", "VSMPs", "Endo", "Immune"))
cols <- c(myIntegratedColors[c("ASCs", "PreAs")], "Aregs" = "#377EB8", 
          "Ifit+" = unname(myIntegratedColors["IFIT"]),
          "Cilp+" = unname(myIntegratedColors["CILP"]),
          myIntegratedColors[c("Meso", "VSMPs", "Endo", "Immune")])
p <- boxplot_score_datatable(score_res = t(Scores_scaled_mouse), clust = clust,
                             clust_name = levels(clust), col = cols, stat_on = F)
stats <- boxplot_score_datatable(score_res = t(Scores_scaled_mouse[names(clust), ]), clust = clust,
                                 clust_name = levels(clust), col = cols, stat_on = T)

ggsave(gridExtra::grid.arrange(p$ASCs, p$PreAs,p$IGFBP2, p$HHIP, p$CILP, p$IFIT, 
                               p$`PR specific`,p$`CHI3L1-2`, p$Meso, p$VSMPs,  ncol = 5),
       filename = "7.Comparison_mASPCs/Boxplot_scores_human_in_mouse_2.pdf", 
       height = 5.8, width = 15)

stats %>% filter(Gene == "ASCs") %>% filter(start == "ASCs" | end == "ASCs")
# Gene        comp  p adj start    end      Expr
# 1: ASCs  ASCs,PreAs <2e-16  ASCs  PreAs 0.4027382
# 2: ASCs  Aregs,ASCs <2e-16 Aregs   ASCs 0.4219163
# 3: ASCs  ASCs,Ifit+ <2e-16  ASCs  Ifit+ 0.4410943
# 4: ASCs  ASCs,Cilp+ <2e-16  ASCs  Cilp+ 0.4602723
# 5: ASCs   ASCs,Meso <2e-16  ASCs   Meso 0.4794503
# 6: ASCs  ASCs,VSMPs <2e-16  ASCs  VSMPs 0.4986283
# 7: ASCs   ASCs,Endo <2e-16  ASCs   Endo 0.5178063
# 8: ASCs ASCs,Immune <2e-16  ASCs Immune 0.5369843

stats %>% filter(Gene == "IGFBP2") %>% filter(start == "Meso" | end == "Meso")
# Gene        comp  p adj  start   end      Expr
# 1: IGFBP2   ASCs,Meso <2e-16   ASCs  Meso 0.3132937
# 2: IGFBP2  Meso,PreAs <2e-16   Meso PreAs 0.4010159
# 3: IGFBP2  Aregs,Meso <2e-16  Aregs  Meso 0.4511429
# 4: IGFBP2  Ifit+,Meso  2e-06  Ifit+  Meso 0.5012698
# 5: IGFBP2  Cilp+,Meso <2e-16  Cilp+  Meso 0.5513968
# 6: IGFBP2  Meso,VSMPs <2e-16   Meso VSMPs 0.6015238
# 7: IGFBP2   Endo,Meso <2e-16   Endo  Meso 0.6140556
# 8: IGFBP2 Immune,Meso <2e-16 Immune  Meso 0.6265873

stats %>% filter(Gene == "Meso") %>% filter(start == "Meso" | end == "Meso")
# Gene        comp  p adj  start   end      Expr
# 1: Meso   ASCs,Meso <2e-16   ASCs  Meso 0.3353750
# 2: Meso  Meso,PreAs <2e-16   Meso PreAs 0.4331928
# 3: Meso  Aregs,Meso <2e-16  Aregs  Meso 0.4890886
# 4: Meso  Ifit+,Meso <2e-16  Ifit+  Meso 0.5589584
# 5: Meso  Cilp+,Meso <2e-16  Cilp+  Meso 0.6148542
# 6: Meso  Meso,VSMPs <2e-16   Meso VSMPs 0.6567761
# 7: Meso   Endo,Meso <2e-16   Endo  Meso 0.6707501
# 8: Meso Immune,Meso <2e-16 Immune  Meso 0.6847240

stats %>% filter(Gene == "HHIP") %>% filter(start %in% c("Aregs", "Cilp+") | 
                                              end %in% c("Aregs", "Cilp+"))
# Gene         comp  p adj start    end      Expr
# 1: HHIP   Aregs,ASCs <2e-16 Aregs   ASCs 0.4019751
# 2: HHIP  Aregs,PreAs <2e-16 Aregs  PreAs 0.5298762
# 3: HHIP  Aregs,Ifit+ <2e-16 Aregs  Ifit+ 0.6395058
# 4: HHIP  Aregs,Cilp+  0.003 Aregs  Cilp+ 0.6577774
# 5: HHIP   Aregs,Meso <2e-16 Aregs   Meso 0.6760490
# 6: HHIP  Aregs,VSMPs <2e-16 Aregs  VSMPs 0.6943206
# 7: HHIP   Aregs,Endo <2e-16 Aregs   Endo 0.7125922
# 8: HHIP Aregs,Immune <2e-16 Aregs Immune 0.7308638

stats %>% filter(Gene == "CILP") %>% filter(start == "Cilp+" | end == "Cilp+")
# Gene         comp  p adj start    end      Expr
# 1: CILP   ASCs,Cilp+ <2e-16  ASCs  Cilp+ 0.6893715
# 2: CILP  Cilp+,PreAs <2e-16 Cilp+  PreAs 0.8617143
# 3: CILP  Aregs,Cilp+ <2e-16 Aregs  Cilp+ 1.0053334
# 4: CILP  Cilp+,Ifit+ <2e-16 Cilp+  Ifit+ 1.1489525
# 5: CILP   Cilp+,Meso <2e-16 Cilp+   Meso 1.2925715
# 6: CILP  Cilp+,VSMPs <2e-16 Cilp+  VSMPs 1.3212953
# 7: CILP   Cilp+,Endo <2e-16 Cilp+   Endo 1.3500191
# 8: CILP Cilp+,Immune <2e-16 Cilp+ Immune 1.3787430

stats %>% filter(Gene == "VSMPs") %>% filter(start == "VSMPs" | end == "VSMPs")
Gene         comp  p adj  start   end      Expr
1: VSMPs   ASCs,VSMPs <2e-16   ASCs VSMPs 0.4637592
2: VSMPs  PreAs,VSMPs <2e-16  PreAs VSMPs 0.5529437
3: VSMPs  Aregs,VSMPs <2e-16  Aregs VSMPs 0.6242913
4: VSMPs  Ifit+,VSMPs <2e-16  Ifit+ VSMPs 0.6956389
5: VSMPs  Cilp+,VSMPs <2e-16  Cilp+ VSMPs 0.7669864
6: VSMPs   Meso,VSMPs <2e-16   Meso VSMPs 0.8204971
7: VSMPs   Endo,VSMPs <2e-16   Endo VSMPs 0.8740078
8: VSMPs Immune,VSMPs <2e-16 Immune VSMPs 0.8918447

stats %>% filter(Gene == "PreAs") %>% filter(start == "PreAs" | end == "PreAs")
Gene         comp  p adj  start   end      Expr
1: PreAs   ASCs,PreAs <2e-16   ASCs PreAs 0.4134248
2: PreAs  Aregs,PreAs <2e-16  Aregs PreAs 0.5709200
3: PreAs  Cilp+,PreAs <2e-16  Cilp+ PreAs 0.5906069
4: PreAs   Meso,PreAs <2e-16   Meso PreAs 0.6102938
5: PreAs  PreAs,VSMPs <2e-16  PreAs VSMPs 0.6299807
6: PreAs   Endo,PreAs <2e-16   Endo PreAs 0.6496676
7: PreAs Immune,PreAs <2e-16 Immune PreAs 0.6693545

stats %>% filter(Gene == "IFIT") %>% filter(start == "Ifit+" | end == "Ifit+")
Gene         comp  p adj start    end      Expr
1: IFIT   ASCs,Ifit+ <2e-16  ASCs  Ifit+ 0.5262537
2: IFIT  Ifit+,PreAs <2e-16 Ifit+  PreAs 0.6864179
3: IFIT  Aregs,Ifit+ <2e-16 Aregs  Ifit+ 0.8237014
4: IFIT  Cilp+,Ifit+ <2e-16 Cilp+  Ifit+ 0.9609850
5: IFIT   Ifit+,Meso <2e-16 Ifit+   Meso 0.9838656
6: IFIT  Ifit+,VSMPs <2e-16 Ifit+  VSMPs 1.0067462
7: IFIT   Endo,Ifit+ <2e-16  Endo  Ifit+ 1.0296268
8: IFIT Ifit+,Immune <2e-16 Ifit+ Immune 1.0525074

stats %>% filter(Gene == "CHI3L1-2") %>% filter(start == "PreAs" | end == "PreAs")
Gene         comp  p adj  start   end      Expr
1: CHI3L1-2   ASCs,PreAs <2e-16   ASCs PreAs 0.5546067
2: CHI3L1-2  Aregs,PreAs <2e-16  Aregs PreAs 0.7394756
3: CHI3L1-2  Cilp+,PreAs <2e-16  Cilp+ PreAs 0.7658854
4: CHI3L1-2   Meso,PreAs <2e-16   Meso PreAs 0.7922953
5: CHI3L1-2  PreAs,VSMPs <2e-16  PreAs VSMPs 0.8187051
6: CHI3L1-2   Endo,PreAs <2e-16   Endo PreAs 0.8451150
7: CHI3L1-2 Immune,PreAs <2e-16 Immune PreAs 0.8715248

stats %>% filter(Gene == "PreAs") %>% filter(start %in% c("PreAs","Aregs","Ifit+","Clip+") | 
                                               end %in% c("PreAs","Aregs","Ifit+","Clip+"))
Gene         comp  p adj  start    end      Expr
1: PreAs   ASCs,PreAs <2e-16   ASCs  PreAs 0.4134248
2: PreAs   Aregs,ASCs <2e-16  Aregs   ASCs 0.4331117
3: PreAs   ASCs,Ifit+ <2e-16   ASCs  Ifit+ 0.4527986
4: PreAs  Aregs,PreAs <2e-16  Aregs  PreAs 0.5709200
5: PreAs  Cilp+,PreAs <2e-16  Cilp+  PreAs 0.5906069
6: PreAs   Meso,PreAs <2e-16   Meso  PreAs 0.6102938
7: PreAs  PreAs,VSMPs <2e-16  PreAs  VSMPs 0.6299807
8: PreAs   Endo,PreAs <2e-16   Endo  PreAs 0.6496676
9: PreAs Immune,PreAs <2e-16 Immune  PreAs 0.6693545
10: PreAs  Aregs,Ifit+  5e-08  Aregs  Ifit+ 0.6890414
11: PreAs   Aregs,Meso <2e-16  Aregs   Meso 0.7087283
12: PreAs  Aregs,VSMPs <2e-16  Aregs  VSMPs 0.7284152
13: PreAs   Aregs,Endo <2e-16  Aregs   Endo 0.7481021
14: PreAs Aregs,Immune <2e-16  Aregs Immune 0.7677890
15: PreAs  Cilp+,Ifit+  3e-05  Cilp+  Ifit+ 0.7874759
16: PreAs   Ifit+,Meso <2e-16  Ifit+   Meso 0.8071628
17: PreAs  Ifit+,VSMPs <2e-16  Ifit+  VSMPs 0.8268497
18: PreAs   Endo,Ifit+ <2e-16   Endo  Ifit+ 0.8465366
19: PreAs Ifit+,Immune <2e-16  Ifit+ Immune 0.8662235

p_100 <- boxplot_score_datatable(score_res = t(Scores_scaled_mouse_100), clust = int_mASPCs$myIntegratedClustering,
                             clust_name = levels(int_mASPCs$myIntegratedClustering), col = cols, stat_on = F)

