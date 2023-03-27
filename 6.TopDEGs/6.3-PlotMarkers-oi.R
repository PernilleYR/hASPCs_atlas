################################################################
#                                                              #
#                    Plot markers oi                        #
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

##---------------------------------------------##
##-----------------Loading data----------------##
##---------------------------------------------##
int <- readRDS("5.Integration/output/Seurat_2000HVGs.Rds")

##---------------------------------------------##
##--------------IGFBP2-TM4SF1-MSLN-------------##
##---------------------------------------------##

DefaultAssay(int) <- "RNA"
df <- as.data.frame(GetAssayData(int)[data.annot[data.annot$gene_short_name %in% c( "TM4SF1", "MSLN"), "ens_id"],])
df$gene_id <- data.annot[rownames(df), "gene_short_name"]
df <- reshape2::melt(as.data.frame(df))
df$variable <- as.character(df$variable)
df$clust <- int$myIntegratedClustering[df$variable]
df$depot <- int$depot[df$variable]
df <- df %>% filter(clust %in% c("ASCs", "PreAs", "IGFBP2", "CILP", "IFIT", "HHIP", "PR specific", "Meso", "VSMPs"))
df <- df %>% filter(depot == "EP")

p <- ggplot(df, aes(x = clust, y = value, fill = gene_id)) + 
  geom_boxplot(outlier.colour =  "gray80", outlier.alpha = 0.5) + 
  mashaGgplot2Theme + ylab("log norm. gene expr.") + xlab("") + 
  scale_fill_manual(values = c("#810F7C", "#CAB2D6"))
ggsave(p, filename = "5.Integration/plots/Markers_oi/Boxplot_IGFBP2-TM4SF1-MSLN.pdf",width = 6.67, height = 3)

ggsave(p, filename = "5.Integration/plots/Markers_oi/Boxplot_TM4SF1-MSLN_onlyEPs.pdf",width = 6.17, height = 3)



##---------------------------------------------##
##-----------QC-markers of cell types----------##
##---------------------------------------------##
g <- c("THY1", "PDGFRA",
       "ACTA2", "TAGLN", "MYH11", 
       "MSLN", "UPK3B", 
       "PECAM1", "PTPRC")
d <- GetAssayData(int, assay = "RNA")[data.annot[data.annot$gene_short_name %in% g,"ens_id"],]

rownames(d) <- data.annot[rownames(d), "gene_short_name"]
d <- as.data.frame(t(d))
d$batch <- int$batch[rownames(d)]
d$pop <- int$mySelectedClustering[rownames(d)]
d <- reshape2::melt(d)
d$variable <- factor(d$variable, levels = g)
d <- d %>% filter(pop != "Unknown_VSMPs")
d$pop <- factor(as.character(d$pop), levels = c("ASCs", "PreAs", "IGFBP2", "VSMPs", "Meso", "Endo", "Immune"))
p <- ggplot(d, aes(y = value, x = variable, fill = pop)) + 
  geom_boxplot(outlier.size = 0.5) + 
  mashaGgplot2Theme + 
  scale_fill_manual(values = myIntegratedColors)

ggsave(plot = p + theme(legend.position = "none"), "5.Integration/plots/Markers_oi/BOXPLOT_QCMarkers_mySelectedClustering.pdf",
       width = 12, height = 4.4)

##---------------------------------------------##
##-------PDGFRA & Immune for IFIT SuppFig------##
##---------------------------------------------##
g <- c( "PDGFRA", "PTPRC", "ITGA4", "NCAM1")
d <- GetAssayData(int, assay = "RNA")[data.annot[data.annot$gene_short_name %in% g,"ens_id"],]

rownames(d) <- data.annot[rownames(d), "gene_short_name"]
d <- as.data.frame(t(d))

d$batch <- int$batch[rownames(d)]
d$pop <- int$myIntegratedClustering[rownames(d)]
d <- reshape2::melt(d)
d$variable <- factor(d$variable, levels = g)
# d <- d %>% filter(pop != "Unknown_VSMPs")
# d$pop <- factor(as.character(d$pop), levels = c("ASCs", "PreAs", "IGFBP2", "VSMPs", "Meso", "Endo", "Immune"))
p <- ggplot(d, aes(y = value, x = variable, fill = pop)) + 
  geom_boxplot(outlier.size = 0.5) + 
  mashaGgplot2Theme + 
  scale_fill_manual(values = myIntegratedColors)

ggsave(plot = p + theme(legend.position = "none"), 
       "6.TopDEGs/Boxplot/Boxplot_PDGFRA-PTPRC_myIntegratedClustering.pdf", width = 5.5, height = 2.43)

##---------------------------------------------##
##-------------WT1, ALDH1A2, CD200-------------##
##---------------------------------------------##
g <- c( "WT1", "ALDH1A2", "CD200")
d <- GetAssayData(int, assay = "RNA")[data.annot[data.annot$gene_short_name %in% g,"ens_id"],]

rownames(d) <- data.annot[rownames(d), "gene_short_name"]
d <- as.data.frame(t(d))

d$depot <- int$depot[rownames(d)]
d$pop <- int$myIntegratedClustering[rownames(d)]
d <- reshape2::melt(d)
d$variable <- factor(d$variable, levels = g)
# d <- d %>% filter(pop != "Unknown_VSMPs")
# d$pop <- factor(as.character(d$pop), levels = c("ASCs", "PreAs", "IGFBP2", "VSMPs", "Meso", "Endo", "Immune"))
p <- ggplot(d, aes(y = value, x = variable, fill = pop)) + 
  geom_boxplot(outlier.size = 0) + 
  #geom_point(aes(shape = depot, group = pop), position = position_jitterdodge()) + scale_shape_manual(values = c(0,5,16,6)) +
  mashaGgplot2Theme + 
  scale_fill_manual(values = myIntegratedColors) 
p
ggsave(plot = p , 
       "6.TopDEGs/Boxplot/Boxplot_WT1-ALDH1A2-CD200_myIntegratedClustering.pdf", width = 10.24, height = 3.24)

##---------------------------------------------##
##----------------SFRP2 & SFRP4----------------##
##---------------------------------------------##
g <- c( "SFRP2", "SFRP4")
d <- GetAssayData(int, assay = "RNA")[data.annot[data.annot$gene_short_name %in% g,"ens_id"],]

rownames(d) <- data.annot[rownames(d), "gene_short_name"]
d <- as.data.frame(t(d))
d$cell <- rownames(d)
d$depot <- int$depot[rownames(d)]
d$pop <- int$myIntegratedClustering[rownames(d)]
d$cell_type <- int$cell_type[rownames(d)]
d <- reshape2::melt(d)
d$variable <- factor(d$variable, levels = g)
# d <- d %>% filter(pop != "Unknown_VSMPs")
# d$pop <- factor(as.character(d$pop), levels = c("ASCs", "PreAs", "IGFBP2", "VSMPs", "Meso", "Endo", "Immune"))
#d <- d %>% filter(cell_type %in% c( "ASPCs"))
p_sfrp2 <- ggplot(d %>% filter(variable == "SFRP2"), 
            aes(y = value, x = cell_type, fill = depot)) + 
  geom_boxplot(outlier.size = 0) + 
  #geom_point(aes(shape = depot, group = pop), position = position_jitterdodge()) + scale_shape_manual(values = c(0,5,16,6)) +
  mashaGgplot2Theme + ylab("SFRP2")+
  scale_fill_manual(values = myDepotsColors ) 

ggsave(plot = grid.arrange(p_sfrp2, p_sfrp4, ncol = 1) , 
       "6.TopDEGs/Boxplot/Boxplot_SFRP2-4_CT-byDep.pdf", width = 5.43, height = 5.52)

d <- d %>% filter(cell_type %in% c( "ASPCs"))
p_sfrp4 <- ggplot(d %>% filter(variable == "SFRP4"), 
            aes(y = value, x = pop, fill = depot)) + 
  geom_boxplot(outlier.size = 0) + 
  #geom_point(aes(shape = depot, group = pop), position = position_jitterdodge()) + scale_shape_manual(values = c(0,5,16,6)) +
  mashaGgplot2Theme + ylab("SFRP4")+
  scale_fill_manual(values = myDepotsColors ) 
ggsave(plot = grid.arrange(p_sfrp2, p_sfrp4, ncol = 1) , 
       "6.TopDEGs/Boxplot/Boxplot_SFRP2-4_PopOfASPCs-byDep.pdf", width = 7.78, height = 5.52)

##---------------------------------------------##
##-------------MSLN, UPK3B, LRRN4--------------##
##---------------------------------------------##
g <- c("MSLN", "UPK3B", "LRRN4")
d <- GetAssayData(int, assay = "RNA")[get_ensid(g),]

rownames(d) <- data.annot[rownames(d), "gene_short_name"]
d <- as.data.frame(t(d))
d$cell <- rownames(d)
d$depot <- int$depot[rownames(d)]
d$pop <- int$myIntegratedClustering[rownames(d)]
d$cell_type2 <- int$cell_type2[rownames(d)]
d$cell_type2 <- factor(d$cell_type2, 
                       levels = c("ASPCs", "IGFBP2", "Meso",
                                  "VSMPs", "Endo", "Immune"))
d <- reshape2::melt(d)
d$variable <- factor(d$variable, levels = g)
d <- d %>% filter(pop != "res.0.2_8")
#d <- d %>% filter(cell_type2 %in% c("ASPCs", "IGFBP2", "Meso", "VSMPs"))
# d$pop <- factor(as.character(d$pop), levels = c("ASCs", "PreAs", "IGFBP2", "VSMPs", "Meso", "Endo", "Immune"))
#d <- d %>% filter(cell_type %in% c( "ASPCs"))
p_upk3b <- ggplot(d %>% filter(variable == "CDH1"), 
                  aes(y = value, x = depot, 
                      fill = cell_type2)) + 
  ggrastr::rasterise(geom_point(aes(color = cell_type2), 
                                pch = 16, alpha = 0.1, size = 0.1, #color = "gray80", 
                                position = position_jitterdodge(jitter.width = 0.4)), dpi = 600) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.4) + 
  mashaGgplot2Theme + ylab("UPK3B")+
  scale_fill_manual(values = unname(c("#E31A1C",
                               myIntegratedColors[levels(d$cell_type2)[-1]]))) + 
  scale_color_manual(values = unname(c("#E31A1C",
                                      myIntegratedColors[levels(d$cell_type2)[-1]])))
p_msln
ggsave(plot = p_msln , 
       "6.TopDEGs/Boxplot/Boxplot_MSLN_CT-byDep.pdf", 
       width = 5.96, height = 4)

P <- grid.arrange(p_msln + theme(axis.title.x=element_blank()), 
                  p_upk3b + theme(axis.title.x=element_blank()), 
                  ncol = 1)
ggsave(P, filename = "6.TopDEGs/Boxplot/Boxplot_MSLN-UPK3B_CT-byDep_withoutEndoImmune.pdf",
       width = 5.24, height = 5.6)


##---------------------------------------------##
##-------------MSLN, UPK3B, LRRN4--------------##
##---------------------------------------------##
emt <- read.table("Utility/EMT-genes.txt", sep = "\t")
myemt <- read.gmt("Utility/GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION.v2022.1.Hs.gmt")
myemt[myemt$gene == "TASOR", "gene"] <- "FAM208A"
table(myemt$gene %in% IGFBP2_EdgeR$gene_id)
g_emt <- unique(c(emt$gene_short_name)#,
                "SNAI1", 
                 "MDK", "PTN", 
                "WNT4", "WNT6", 
                "UPK1B", 
                "ZFPM2",
                "ZEB1", "ZEB2", 
                "FGF1", "FGF2",
                "DCLK1", "KRAS",  
                "AIFM2", "COL1A1", "ACTA2", "DES", "VIM", "DDR2", "MMP9",
                "MMP2", "MMP19", "MMP23B"))
d <- GetAssayData(int, assay = "RNA")[rownames(g_emt),]

rownames(d) <- data.annot[rownames(d), "gene_short_name"]
d <- as.data.frame(t(d))
d$cell <- rownames(d)
d$depot <- int$depot[rownames(d)]
d$pop <- int$myIntegratedClustering[rownames(d)]
d$cell_type2 <- int$cell_type2[rownames(d)]
d$cell_type2 <- factor(d$cell_type2, 
                       levels = c("ASPCs", "IGFBP2", "Meso",
                                  "VSMPs", "Endo", "Immune"))
d <- reshape2::melt(d)
#d$variable <- factor(d$variable, levels = g_emt$)
d <- d %>% filter(pop != "res.0.2_8")
#d <- d %>% filter(cell_type2 %in% c("ASPCs", "IGFBP2", "Meso", "VSMPs"))
# d$pop <- factor(as.character(d$pop), levels = c("ASCs", "PreAs", "IGFBP2", "VSMPs", "Meso", "Endo", "Immune"))
#d <- d %>% filter(cell_type %in% c( "ASPCs"))
p <- ggplot(d %>% filter(variable == "FN1"), 
                  aes(y = value, x = depot, 
                      fill = cell_type2)) + 
  # ggrastr::rasterise(geom_point(aes(color = cell_type2), 
  #                               pch = 16, alpha = 0.1, size = 0.1, #color = "gray80", 
  #                               position = position_jitterdodge(jitter.width = 0.4)), dpi = 600) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.4) + 
  mashaGgplot2Theme + ylab("UPK3B")+
  scale_fill_manual(values = unname(c("#E31A1C",
                                      myIntegratedColors[levels(d$cell_type2)[-1]]))) + 
  scale_color_manual(values = unname(c("#E31A1C",
                                       myIntegratedColors[levels(d$cell_type2)[-1]])))
p
ggsave(plot = p_msln , 
       "6.TopDEGs/Boxplot/Boxplot_MSLN_CT-byDep.pdf", 
       width = 5.96, height = 4)

P <- grid.arrange(p_msln + theme(axis.title.x=element_blank()), 
                  p_upk3b + theme(axis.title.x=element_blank()), 
                  ncol = 1)
ggsave(P, filename = "6.TopDEGs/Boxplot/Boxplot_MSLN-UPK3B_CT-byDep_withoutEndoImmune.pdf",
       width = 5.24, height = 5.6)



[1] "UPP1"         "MLLT11"       "PITPNC1"      "CDCP1"        "NFKB2"        "PLAU"         "S100A2"      
[8] "TIMP3"        "INHBA"        "CUL4B"        "HMOX1"        "SPRY4"        "ARL4A"        "CXCL3"       
[15] "TFPI2"        "DUSP5"        "CCL26"        "GDAP2"        "CXCL8"        "BDKRB1"       "TIPIN"       
[22] "ITGA2"        "MMP3"         "RPSAP52"      "CXCL1"        "CBWD1"        "IL33"         "LIF"         
[29] "DIRAS3"       "CXCL6"        "SIGLEC15"     "PRR9"         "AJAP1"        "IVL"          "SERPINB2"    
[36] "CXCL5"        "BMP2"         "MMP9"         "CSF3"         "BCL2A1"       "PLXNA2"       "AREG"        
[43] "PI3"          "CD24"         "CORO7"        "EFNB2"        "ZNF549"       "LCEP3"        "ATP6V0E2-AS1"
[50] "PPBP"         "PCYT1B"       "RAP1GAP"      "AC207130.1"   "ARHGAP39"     "SMAD5-AS1"   


int <- SetIdent(int, value = "depot")
FeaturePlot(subset(int, depot %in% c("EP")), features = get_ensid(c("WT1", "ALDH1A2", "CD200", "MME", "DKK2")), 
            reduction ="tsne", order =T)

FeaturePlot(subset(int, depot %in% c("EP")), features = get_ensid(c("CXCL1", "CXCL3", "CXCL6", "CXCL8")), 
            reduction ="tsne", order =T)


FeaturePlot(subset(int, depot %in% c("SC")), features = get_ensid(c("IL33")), 
            reduction ="tsne", order =T)




aspcs <- subset(int, myIntegratedClustering %in% c("PreAs", "ASCs"))
aspcs <- RunPCA(aspcs, npcs = 1:30)



d <- dist(int@reductions$pca@cell.embeddings[int$myIntegratedClustering %in% c("PreAs", "ASCs"), 1:30])
df <- reshape2::melt(as.matrix(d))
df$Var1 <- as.character(df$Var1)
df$Var2 <- as.character(df$Var2)
df$cat_var1 <- int$myIntegratedClustering[df$Var1]
df$cat_var2 <- int$myIntegratedClustering[df$Var2]

t <- df[c(1:100, 2000:2010,]




