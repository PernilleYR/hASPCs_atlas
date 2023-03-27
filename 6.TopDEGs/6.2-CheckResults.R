################################################################
#                                                              #
#                  Investigate top DE markers                  #
#                                                              #
################################################################

### Author: Pernille
### Date: 17.08.2022 - adapted from old script
### Datasets: scRNA-seq Depots: SC - EP - MK - MG - PR - GB 
###                     Patients: B0, B1, B7, L3, L11, L12 (B - Bariatric, L - Lean)
### Goal: Compare markers between human and mouse ASPCs 

library(ggplot2); library(data.table); library(Seurat); library(dplyr)

setwd("~/SVRAW1/prainer/hASPCs/PAPER/10X_scRNA-seq/")
source("Utility/General_utils.R")
source("6.TopDEGs/6-utils.R")

##---------------------------------------------##
##-----------------Loading data----------------##
##---------------------------------------------##

myDEGs_seurat <- readRDS("6.TopDEGs/DEGs.Rds")
int <- readRDS("5.Integration/output/Seurat_2000HVGs.Rds")

##---------------------------------------------##
##-----------------tSNE top 20-----------------##
##---------------------------------------------##

dir.create("6.TopDEGs/plots")
mygradients_of_colors <- list("ASCs" = c('#f7fcf5','#e5f5e0','#c7e9c0','#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b'),
                              "PreAs" = c('#fff5f0','#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d'),
                              "IGFBP2" = c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b'),
                              "Meso" = c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d'),
                              "VSMPs" = c('#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506'),
                              "HHIP" = c('#fff7fb','#ece2f0','#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a','#016c59','#014636'),
                              "PR specific" = c('#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506'),
                              "other" = c('#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b'))

lapply(names(myDEGs_seurat), plot_topN)


##---------------------------------------------##
##--------------Dot plot ALL POPs--------------##
##---------------------------------------------##
dir.create("6.TopDEGs/Dotplots")

DefaultAssay(int) <- "RNA"
int$order_dot_plot <- as.character(int$myIntegratedClustering)
int$order_dot_plot <- factor(int$order_dot_plot,
                             levels = rev(c("ASCs", "PreAs", 
                                            "HHIP", "IFIT", "CILP",
                                            "CHI3L1-2", "PR specific",
                                            "IGFBP2",
                                            "Meso", 
                                            "VSMPs", "Endo", "Immune", "res.0.2_8")))
int <- SetIdent(int, value = "order_dot_plot")


##---------------##
##----MINI V1----##
##---------------##

mydotgenes <- c("CD55", "MFAP5", #ASCs
                "CXCL14", "GPC3", #PreAs
                "HHIP", "GDF10", "NOV", #Aregs
                "ISG15", "IFIT3", "STAT1", "MX1", #IFIT
                "SFRP2", #CLIP
                "CHI3L1", "CHI3L2", "RBP5", #CHI3L1
                "GPX3", #PR
                "IGFBP2", "RBP1", #IGFBP2
                
                
                "MSLN", "TM4SF1", #Meso
                
                "ACTA2", "MYH11", #VSMPs
                
                "PECAM1", #Endo
                
                "PTPRC" #Immune
                )
mydotgenes <- convert_geneID_to_data.annot(mydotgenes, data.annot = data.annot)

p <- DotPlot(int, features = rownames(mydotgenes), 
        cols = c("white", "red")) + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(labels= mydotgenes$gene_short_name)
ggsave(p, filename = "6.TopDEGs/Dotplots/Dotplot_long.pdf",
       height = 4.18, width = 10)


##---------------##
##----MINI V2----##
##---------------##
mydotgenes <- c("MFAP5", #ASCs
                "CXCL14", #PreAs
                "IGFBP2", #IGFBP2
                "GDF10", #Aregs
                "RBP5", #CHI3L1
                "IFIT3", #IFIT
                "SFRP2", "CILP",#CILP
                "GPX3", #PR
                
                "MSLN", #Meso
                
                "ACTA2", #VSMPs
                
                "PECAM1", #Endo
                
                "PTPRC" #Immune
)
mydotgenes <- convert_geneID_to_data.annot(mydotgenes, data.annot = data.annot)


p <- DotPlot(int, features = mydotgenes$ens_id, 
        cols = c("white", "red")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_x_discrete(labels= mydotgenes$gene_short_name)
  
ggsave(p, filename = "6.TopDEGs/Dotplots/Dotplot_simple.pdf",
       height = 4.18, width = 7.63)



##---------------##
##----BIG TEST---##
##---------------##

for(n in names(myDEGs_seurat)){
  #myDEGs_seurat[[n]] <- myDEGs_seurat[[n]] %>% filter(avg_logFC_all > 0 & minimump_p_val < 0.05)
  print(n)
  print(nrow(myDEGs_seurat[[n]]))
}
mydotgenes <-unique(c(rownames(myDEGs_seurat$ASCs)[1:30],
                rownames(myDEGs_seurat$PreAs)[1:30],
                rownames(myDEGs_seurat$HHIP)[1:30],
                rownames(myDEGs_seurat$IFIT)[1:30],
                rownames(myDEGs_seurat$CILP)[1:30],
                rownames(myDEGs_seurat$`CHI3L1-2`)[1:30],
                rownames(myDEGs_seurat$`PR specific`)[1:30],
                rownames(myDEGs_seurat$IGFBP2)[1:30],
                
                rownames(myDEGs_seurat$Meso)[1:30],
                
                rownames(myDEGs_seurat$VSMPs)[1:30],
                
                data.annot[data.annot$gene_short_name == "PECAM1", "ens_id"],
                
                data.annot[data.annot$gene_short_name == "PTPRC", "ens_id"]))

mydotgenes <- data.annot[rownames(myDEGs_seurat$ASCs)[1:40],]
#mydotgenes <- data.annot[mydotgenes, ]

##---------------##
##-----FINAL-----##
##---------------##

mydotgenes <- c("MFAP5", "PRG4", "CD55", "PI16", "TPPP3", "SEMA3C", "ACKR3", "PCOLCE2", "C1QTNF3", "HTRA3", #ASCs
                "CXCL14", "CXCL12", "APOD", "ADM", "IGFBP3", "GPC3", "SRPX", "GPX3", "GAS6", "TMEM176A", #PreAs
                "IGFBP7", "NOV", "TSPAN8", "ALDH1A1", "OMD", "HHIP", "GDF10", "CFH", "FAM180B", "F3", #HHIP
                "ISG15", "IFI6", "IFI27", "MX1", "IFI44L", "IFIT3", "STAT1", "MX2", "XAF1", "RSAD2", #IFIT
                "SFRP4", "SFRP2", "CTHRC1", "COL14A1", "WISP2", "ELN", "ASPN", "LUM", "CTSK", "LTBP2", #CILP
                "TIMP1", "COL3A1", "COL6A3", "TYMP", "RBP5", "PSME2", "COL6A1", "COL4A2", "SERPINE1", "PKM", #CHI3L
                "FMO2", "LMO3", "MT-RNR2", "MT-RNR1", "ABCA6", "SVEP1", "SFRP1", "DEPP1", "ABCA10", "MT-ATP6", #PR
                "IGFBP2", "G0S2", "C7", "APOE", "RBP1", "SNAI2", "DDIT4", "PHLDA1","PTN", "MARCKSL1", #IGFBP2)
                "ITLN1", "KRT19", "TM4SF1", "MSLN", "KRT8", "UPK3B", "TFPI2", "SLPI", "BCHE", "C19orf33",#Meso
                "ACTA2", "TAGLN", "MYH11", "TPM2", "MYL9", "RGS5", "PLN", "RERGL", "SORBS2", "MCAM",
                "PECAM1", 
                "PTPRC")
mydotgenes <- convert_geneID_to_data.annot(mydotgenes, data.annot = data.annot)
p <- DotPlot(int, features = mydotgenes$ens_id, col.max = 3) + 
        #scale_color_distiller(palette = "YlGnBu", direction = 1) + 
        scale_color_distiller(palette = "YlOrRd", direction = 1) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(labels= mydotgenes$gene_short_name)
ggsave(p, filename = "6.TopDEGs/Dotplots/DotPlot_10MarkersPerClust_red.pdf",
       height = 6.2, width = 22.16)


##---------------------------------------------##
##---------------Dot plot HHIP+----------------##
##---------------------------------------------##
murine_Aregs_markers <- readRDS("~/SVRAW1/prainer/Integration_sc_datasets/SC_and_Visc_Mouse/dim60/TopMarkers_fisherMethod/AggregatedFisherPval_PerCluster_orderFC.Rds")
murine_sc_Aregs <- read.table("~/SVRAW1/prainer/hASPCs/PAPER/Files/MyMurineAregsMarkers_Integration_Top20BRBseq_EMBOJournal.txt")
HHIP <- myDEGs_seurat$HHIP %>% filter(avg_logFC_all > 0 & minimump_p_val < 0.05)
HHIP$geneID[HHIP$geneID %in% toupper(murine_sc_Aregs$geneID)]
#"IGFBP7"  "APOD"    "CFH"     "BGN"     "HHIP"    "GDF10"   "MGP"     "INMT"    "F3"      "CLEC11A" "FGF7"  

mydotgenes <- c("IGFBP7", "MYOC", "NOV", "TSPAN8", "APOD","CFH","BGN","HHIP","GDF10","MGP","INMT","F3","CLEC11A", "FGF7", "EPHA3")
mydotgenes <- convert_geneID_to_data.annot(mydotgenes, data.annot = data.annot)
int <- SetIdent(int, value = "myIntegratedClustering")
p <- DotPlot(int, features = mydotgenes$ens_id) +
  scale_color_distiller(palette = "YlOrRd", direction = 1) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_x_discrete(labels= mydotgenes$gene_short_name)
ggsave(plot = p, filename = "6.TopDEGs/Dotplots/Dotplot_HHIP_shared-with-murineAregs_v2.pdf",
       height = 3.97, width = 6.26)
#Saving 6.26 x 3.97 in image

mydotgenes <- myDEGs$HHIP[1:34,]
int <- SetIdent(int, value = "myIntegratedClustering")
p <- DotPlot(int, features = c(rownames(mydotgenes),"ENSG00000044524") , 
             cols = c("white", "red")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_x_discrete(labels= c(mydotgenes$geneID, "EPHA3"))
ggsave(plot = p, filename = "6.TopDEGs/Dotplots/Dotplot_TOP35-HHIPmarkers.pdf",
       width = 12.75, height = 3.97)


##---------------------------------------------##
##---------------Dot plot IFIT+----------------##
##---------------------------------------------##
DefaultAssay(int) <- "RNA"
mydotgenes <- myDEGs_seurat$IFIT[1:20,]
int <- SetIdent(int, value = "myIntegratedClustering")
p <- DotPlot(int, features = rownames(mydotgenes)) +
  scale_color_distiller(palette = "YlOrRd", direction = 1) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_x_discrete(labels= mydotgenes$geneID)
ggsave(plot = p, filename = "6.TopDEGs/Dotplots/Dotplot_IFIT.pdf",
       height = 3.97, width = 8)
#Saving 6.26 x 3.97 in image

##---------------------------------------------##
##-----------Dot plot CILP+ or SRFP4+----------##
##---------------------------------------------##
DefaultAssay(int) <- "RNA"

mydotgenes <- myDEGs_seurat$CILP[1:20,]
int <- SetIdent(int, value = "myIntegratedClustering")
p <- DotPlot(int, features = rownames(mydotgenes)) +
  scale_color_distiller(palette = "YlGn", direction = 1) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_x_discrete(labels= mydotgenes$geneID)
ggsave(plot = p, filename = "6.TopDEGs/Dotplots/Dotplot_CILP.pdf",
       height = 4.17, width = 8.05)


##---------------------------------------------##
##--------------Dot plot IGFBP2+---------------##
##---------------------------------------------##
IGFBP2_EdgeR <- readRds("~/SVRAW1/prainer/hASPCs/")

##---------------------------------------------##
##-------------Dot plot emt genes--------------##
##---------------------------------------------##
int_eps <- subset(int, depot == "EP" )
int_eps <- subset(int_eps, cell_type %in% c("ASPCs","Meso") )
int_eps <- subset(int_eps, myIntegratedClustering %in% c("IGFBP2","Meso", "PreAs", "ASCs") )

emt <- read.table("Utility/EMT-genes.txt", sep = "\t")
g_emt <- unique(c(emt$gene_short_name)#,
                  "WT1", "ALDH1A2", "SNAI1", "SNAI2", "TWIST1", "TWIST2",
                  "GATA4", "TBX18", "MDK", "PTN", "WNT4", "WNT6", "UPK1B", 
                  "VEGFA", "VIM", "WNT5A", "SF1", "TCF21", "ZFPM2", "MUC16",
                  "LGR5", "ZEB1", "ZEB2", "DCLK1", "KRAS", "FGF1", "FGF2", 
                  "AIFM2", "COL1A1", "ACTA2", "DES", "VIM", "DDR2", "MMP9",
                  "MMP2", "MMP19", "MMP23B"))
g_emt <- data.annot[data.annot$gene_short_name %in% g_emt,]
DefaultAssay(int_eps) <- "RNA"
p <- DotPlot(int_eps, features = rownames(g_emt), dot.scale = 10, 
             scale = T) +
  scale_color_distiller(palette = "RdYlBu", direction = -1) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_x_discrete(labels= g_emt$gene_short_name)
p
ggsave(plot = p, filename = "6.TopDEGs/Dotplots/Dotplot_CILP.pdf",
       height = 4.17, width = 8.05)


##---------------------------------------------##
##------------------Heatmap--------------------##
##---------------------------------------------##
int <- FindVariableFeatures(int)
int <- ScaleData(int, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), 
                 features = VariableFeatures(int))
pdf("test.pdf")
DoHeatmap(int, features = unique(c(rownames(myDEGs$ASCs)[1:10],
                            rownames(myDEGs$PreAs)[1:10],
                            rownames(IGFBP2_EdgeR)[1:10],
                            rownames(myDEGs$CILP)[1:10],
                            rownames(myDEGs$HHIP)[1:10],
                            rownames(myDEGs$`CHI3L1-2`)[1:10],
                            rownames(myDEGs$Meso)[1:10])))
dev.off()

##---------------------------------------------##
##------------------Heatmap--------------------##
##---------------------------------------------##
my_genes <- unlist(lapply(myDEGs, function(x) return(rownames(x[1:20,]))))
my_genes <- unique(my_genes)

DefaultAssay(int) <- "RNA"
DoHeatmap(int, my_genes, slot = "data")

d <- GetAssayData(int, assay = "RNA")[my_genes,]
rownames(d) <- data.annot[rownames(d), "gene_short_name"]
my_annot_c <- data.frame(row.names = colnames(d),
                         cellType = int$myIntegratedClustering[colnames(d)])
my_colour = list(
  cat = myIntegratedColors
)

#sample randomly 200 cells from each pop
out_final <- list()
meta <- int@meta.data[, c("myIntegratedClustering", "depot")]
for(l in levels(int$myIntegratedClustering)){
  m <- meta %>% filter(myIntegratedClustering == l)
  to_grab <- min(200, nrow(m))
  n_dep <- table(meta$depot); n_dep <- names(n_dep[n_dep > 6])
  to_grab_perD <- round(to_grab/length(n_dep))
  m$cell_id <- rownames(m)
  out <- c()
  for(d in m$depot){
    grab <- min(to_grab_perD, nrow(m %>% filter(depot == d)))
    out <- c(out,sample(m[m$depot == d, "cell_id"], size = grab))
  }
  if(length(out) < to_grab){
    out <- c(out,sample(m[!m$cell_id %in% out, "cell_id"], size = (to_grab - length(out))))
  }
  out_final[[l]] <- out
}

pheatmap::pheatmap(d, scale = "row", clustering_method = "ward.D2",
                   cluster_rows = F,
                   annotation_colors = my_colour,
                   #color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100),
                   annotation = my_annot_c,
                   annotation_row = my_annot_r,
                   border_color = "NA",
                   cutree_col = 2,
                   fontsize_row = 7,
                   cellwidth = 6, cellheight = 8.5)

DE_in_all <- function(degs){
  lgfc <- colnames(degs)[grep("avg_log2FC", colnames(degs))]
  degs <- degs[rowSums(degs[, lgfc] > 0.25) == length(lgfc),]
  
  pval <- colnames(degs)[grep("_p_val_adj", colnames(degs))]
  degs <- degs[rowSums(degs[, pval] < 0.1) == length(pval),]
  return(degs)
  }

myDEGs_f <- lapply(myDEGs, DE_in_all)
names(myDEGs_f) <- names(myDEGs)

DefaultAssay(int) <- "RNA"


mygenes <- lapply(myDEGs_f, function(x) return(rownames(x[1:min(10, nrow(x)),])))
data.annot[unlist(mygenes),]

dir.create("6.TopDEGs/Heatmap/")

mygenes <- unique(c(rownames(myDEGs_f$ASCs)[1:20],
                    rownames(myDEGs_f$PreAs)[1:20]))
mycells <- int$myIntegratedCluster_depot[int$myIntegratedClustering %in% c("ASCs", "PreAs")]
mycells <- sort(mycells)
d <- GetAssayData(int, assay = "RNA")[unlist(mygenes),names(mycells)]
rownames(d) <- data.annot[rownames(d), "gene_short_name"]
d <- dynutils::scale_quantile(t(d))


d <- as.data.frame(t(d))
d$gene_id <- rownames(d)
d <- reshape2::melt(d)
d$variable <- factor(d$variable, levels = names(mycells))
d$gene_id <- factor(d$gene_id, levels = rev(data.annot[mygenes, "gene_short_name"]))
p_heat <- ggplot(d, aes(x = variable, y = gene_id, fill = value)) + 
  ggrastr::rasterise(geom_tile(), dpi = 650) + xlab("") + theme(axis.text.x = element_blank()) + 
  scale_fill_distiller(palette="RdYlBu", type = "div")

ggsave(p_heat, filename = "6.TopDEGs/Heatmap/Heatmap.pdf", width = 9.57, height = 9.18)


meta_heat <- int@meta.data[levels(d$variable), c("myIntegratedClustering", "depot")]
meta_heat$cell_id <- rownames(meta_heat)
meta_heat$cell_id <- factor(meta_heat$cell_id, levels = meta_heat$cell_id)

bar_clust <- ggplot(meta_heat, aes(x = cell_id, y = 1, fill = myIntegratedClustering)) + 
  ggrastr::rasterise(geom_bar(stat = "identity")) + theme_void() + 
  theme(legend.position = "bottom", axis.text.x = element_blank()) + 
  scale_fill_manual(values = myIntegratedColors)

bar_dep <- ggplot(meta_heat, aes(x = cell_id, y = 1, fill = depot)) + 
  ggrastr::rasterise(geom_bar(stat = "identity")) + theme_void() + 
  theme(legend.position = "bottom", axis.text.x = element_blank()) +
  scale_fill_manual(values = c("#FEC010", "#945200","#6F3996", "#0096FF"))

ggsave(gridExtra::grid.arrange(bar_dep, bar_clust, ncol = 1),
       filename = "6.TopDEGs/Heatmap/legend_bar.pdf")
#Saving 12.5 x 2.72 in image



