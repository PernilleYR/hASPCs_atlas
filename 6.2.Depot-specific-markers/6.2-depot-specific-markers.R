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
dir.create("6.TopDEGs/")

source("Utility/General_utils.R")
source("6.TopDEGs/6-utils.R")

dir.create("6.2.Depot-specific-markers")

##---------------------------------------------##
##-----------------Loading data----------------##
##---------------------------------------------##

int <- readRDS("5.Integration/output/Seurat_2000HVGs.Rds")

##---------------------------------------------##
##-------------------Set data------------------##
##---------------------------------------------##

DefaultAssay(int) <- "RNA"
int <- SetIdent(int, value = int$myIntegratedClustering)
myclust <- levels(int$myIntegratedClustering)

##---------------------------------------------##
##----------DEGs within pops 1 vs All----------##
##---------------------------------------------##

int <- SetIdent(int, value = "myIntegratedClustering")
int$myIntegratedCluster_depot <- paste0(int$myIntegratedClustering, "_", int$depot)
int$myIntegratedCluster_depot <- as.factor(int$myIntegratedCluster_depot)
int <- SetIdent(int, value = "myIntegratedCluster_depot")

cell_pop_oi <- c("ASCs", "PreAs", "HHIP", "IFIT", "CHI3L1-2", "CILP")
mydepot <- c("SC", "EP", "PR", "MG")
DE_res <- lapply(cell_pop_oi, function(cp){
  print(cp)
  res <- lapply(mydepot, function(d){
    print(d)
    out <- FindMarkers_by_depots(cell_pop = cp, depot = d)
    return(out)
  })
  names(res) <- mydepot
  return(res)
})
names(DE_res) <- cell_pop_oi

# DE_PR.SC_vs_OM.MK <- lapply(cell_pop_oi, FindMarkers_by_depotType)
# names(DE_PR.SC_vs_OM.MK) <- cell_pop_oi

##---------------------------------------------##
##------------DEGs specific to a pop-----------##
##---------------------------------------------##

# Within Pop perform 1 vs 1 DE, and identify genes DE in all case 
# (i.e. fully specific to this depot)

ASCs_SC <- finding_specific_markers_v1("ASCs", "SC")
ASCs_EP <- finding_specific_markers_v1("ASCs", "EP")
ASCs_PR <- finding_specific_markers_v1("ASCs", "PR")
ASCs_MG <- finding_specific_markers_v1("ASCs", "MG")
save(ASCs_SC, ASCs_EP, ASCs_PR, ASCs_MG, 
     file = "/Volumes/UPDEPLA/prainer/hASPCs/PAPER/10X_scRNA-seq/6.2.Depot-specific-markers/ASCs.Rdata")

PreAs_SC <- finding_specific_markers_v1("PreAs", "SC") #APOE, CD36, APOC1, KLF4, KLF6, WISP2, PDGFRB
PreAs_EP <- finding_specific_markers_v1("PreAs", "EP") #PTN, RARRES1, IGFBP3, IGFBP7,IGFBP6
PreAs_PR <- finding_specific_markers_v1("PreAs", "PR")
PreAs_MG <- finding_specific_markers_v1("PreAs", "MG")
save(PreAs_SC, PreAs_EP, PreAs_PR, PreAs_MG, 
     file = "/Volumes/UPDEPLA/prainer/hASPCs/PAPER/10X_scRNA-seq/6.2.Depot-specific-markers/PreAs.Rdata")

CILP_SC <- finding_specific_markers_v1("CILP", "SC")
CILP_EP <- finding_specific_markers_v1("CILP", "EP")
CILP_PR <- finding_specific_markers_v1("CILP", "PR")
CILP_MG <- finding_specific_markers_v1("CILP", "MG")
save(CILP_SC, CILP_EP, CILP_PR, CILP_MG, 
     file = "/Volumes/UPDEPLA/prainer/hASPCs/PAPER/10X_scRNA-seq/6.2.Depot-specific-markers/CILP.Rdata")

IFIT_SC <- finding_specific_markers_v1("IFIT", "SC")
IFIT_EP <- finding_specific_markers_v1("IFIT", "EP")
IFIT_PR <- finding_specific_markers_v1("IFIT", "PR")
IFIT_MG <- finding_specific_markers_v1("IFIT", "MG")
save(IFIT_SC, IFIT_EP, IFIT_PR, IFIT_MG, 
     file = "/Volumes/UPDEPLA/prainer/hASPCs/PAPER/10X_scRNA-seq/6.2.Depot-specific-markers/IFIT.Rdata")

CHI3L1.2_SC <- finding_specific_markers_v1("CHI3L1-2", "SC")
CHI3L1.2_EP <- finding_specific_markers_v1("CHI3L1-2", "EP")
CHI3L1.2_PR <- finding_specific_markers_v1("CHI3L1-2", "PR")
CHI3L1.2_MG <- finding_specific_markers_v1("CHI3L1-2", "MG")
save(CHI3L1.2_SC, CHI3L1.2_EP, CHI3L1.2_PR, CHI3L1.2_MG, 
     file = "/Volumes/UPDEPLA/prainer/hASPCs/PAPER/10X_scRNA-seq/6.2.Depot-specific-markers/CHI3L1.2.Rdata")

HHIP_SC <- finding_specific_markers_v1("HHIP", "SC")
HHIP_EP <- finding_specific_markers_v1("HHIP", "EP")
HHIP_PR <- finding_specific_markers_v1("HHIP", "PR")
HHIP_MG <- finding_specific_markers_v1("HHIP", "MG")
save(HHIP_SC, HHIP_EP, HHIP_PR, HHIP_MG, 
     file = "/Volumes/UPDEPLA/prainer/hASPCs/PAPER/10X_scRNA-seq/6.2.Depot-specific-markers/HHIP.Rdata")

##---------------------------------------------##
##-----------------Plot heatmap----------------##
##---------------------------------------------##
P_ascs <- plot_heatmap_pop(pop = "ASCs", 
                           data = list("SC" = ASCs_SC,
                                       "PR" = ASCs_PR,
                                       "EP" = ASCs_EP,
                                       "MG" = ASCs_MG),
                           method = "avg")
ggsave(P_ascs$p_non.scaled,
       filename = "6.2.Depot-specific-markers/Heatmap_ASCs_non-scaled-2.pdf",
       width = 1.9, height = 8.47)
ggsave(P_ascs$p_scaled,
       filename = "6.2.Depot-specific-markers/Heatmap_ASCs_scaled-2.pdf",
       width = 1.9, height = 8.47)

P_preas <- plot_heatmap_pop(pop = "PreAs", 
                            data = list("SC" = PreAs_SC,
                                        "PR" = PreAs_PR,
                                        "EP" = PreAs_EP,
                                        "MG" = PreAs_MG),
                           method = "avg")
ggsave(P_preas$p_non.scaled,
       filename = "6.2.Depot-specific-markers/Heatmap_PreAs_non-scaled-2.pdf",
       width = 1.9, height = 8.47)
ggsave(P_preas$p_scaled,
       filename = "6.2.Depot-specific-markers/Heatmap_PreAs_scaled-2.pdf",
       width = 1.9, height = 8.47)


##---------------------------------------------##
##----------------------GO---------------------##
##---------------------------------------------##

GO_res_CLASSIC_PreAs_perdepot <- lapply(list("SC" = PreAs_SC, "EP" = PreAs_EP,
                           "PR" = PreAs_PR, "MG" = PreAs_MG), 
                      function(m){
                        allG <- GetAssayData(int, assay = "RNA")[, int$myIntegratedClustering == "PreAs"]
                        allG <- rownames(allG)[rowSums(allG) > 0]
                      
                        # Find genes specific, ordered by avg exp 
                        m <- order_g(d = m, m = "avg")
                        out <- GOenrichment(m$ens_id, allGenesList = allG)
                        return(out)
                        })
saveRDS(GO_res_elim_PrAs_perdepot, "/Volumes/UPDEPLA/prainer/hASPCs/PAPER/10X_scRNA-seq/6.2.Depot-specific-markers/GO_res_elim_PreAs.Rds")

GO_res_elim_PrAs_perdepot <- readRDS("6.2.Depot-specific-markers/GO_res_elim_PreAs.Rds")

GO_res_elim_PrAs_perdepot$SC[grep("lip", GO_res_elim_PrAs_perdepot$SC$Term),]
GO.ID                                        Term Annotated Significant Expected classicFisher eLimFisher
10 GO:0045834 positive regulation of lipid metabolic p...       145           7     1.03       7.8e-05    7.8e-05
31 GO:0006869                             lipid transport       375          13     2.66       2.6e-06    0.00048
50 GO:0032368               regulation of lipid transport       109           5     0.77       0.00107    0.00107
69 GO:0016042                     lipid catabolic process       317           8     2.25       0.00190    0.00190
84 GO:0042157               lipoprotein metabolic process       136           5     0.96       0.00284    0.00284
GO_res_elim_PrAs_perdepot$SC[grep("diff", GO_res_elim_PrAs_perdepot$SC$Term),]
9  GO:0045597 positive regulation of cell differentiat...       821          17     5.82       6.7e-05    6.7e-05
71 GO:0045667    regulation of osteoblast differentiation       125           5     0.89       0.00197    0.00197


GO_res_elim_PrAs_perdepot$EP[grep("diff", GO_res_elim_PrAs_perdepot$EP$Term),]


GO_res_elim_PrAs_perdepot$EP[grep("diff", GO_res_elim_PrAs_perdepot$EP$Term),]

adipo_genes <- data.annot[data.annot$gene_short_name %in% c("FABP4", "APOE", "APOC1", "PDGFRA", 
                                                            "PPARG", "CEBPD", "CEBPE", "CEBPA",
                                                            "FAS", "LPL", "DLK1"),]

d <- as.data.frame(t(GetAssayData(int, assay = "RNA")[adipo_genes$ens_id,]))
colnames(d) <- data.annot[colnames(d), "gene_short_name"]
d$depot <- int$depot[rownames(d)]
d$clust <- int$myIntegratedClustering[rownames(d)]
d <- reshape2::melt(d)
ggplot(d, aes(x = clust, fill = depot, y = value)) + geom_boxplot() + 
  facet_wrap(~ variable, scales = "free_y", nrow = 10) + 
  mashaGgplot2Theme
  


int <- SetIdent(int, value = "myIntegratedClustering")
VlnPlot(int, features = adipo_genes$ens_id, idents = c("ASCs", "PreAs"), 
        group.by = "myIntegratedClustering", split.by = "depot") 


EnrichR_PreAs_perdepot <- lapply(list("SC" = PreAs_SC, "EP" = PreAs_EP,
                                             "PR" = PreAs_PR, "MG" = PreAs_MG), 
                                        function(m){
                                          allG <- GetAssayData(int, assay = "RNA")[, int$myIntegratedClustering == "PreAs"]
                                          allG <- rownames(allG)[rowSums(allG) > 0]
                                          m <- order_g(d = m, m = "avg")
                                          out <- myEnrichRfunction(list("aw" = m), 
                                                                   libraryName = c("GO_Biological_Process_2018",
                                                                                   "GO_Biological_Process_2021",
                                                                                   "WikiPathways_2016"))
                                          return(out)
                                        })
saveRDS(EnrichR_PreAs_perdepot, file = "6.2.Depot-specific-markers/EnrichR_preAs_perdepot.Rds")



## PR
EnrichR_PreAs_perdepot$PR$aw$GO_Biological_Process_2021[grep("fat",EnrichR_PreAs_perdepot$PR$aw$GO_Biological_Process_2021$Term),] 

cell morphogenesis involved in differentiation (GO:0000904)
positive regulation of cell differentiation (GO:0045597)
positive regulation of glucose import (GO:0046326)
positive regulation of glucose transmembrane transport (GO:0010828)
inositol lipid-mediated signaling (GO:0048017)

EnrichR_PreAs_perdepot$PR$aw$GO_Biological_Process_2021[grep("thermo",EnrichR_PreAs_perdepot$PR$aw$GO_Biological_Process_2021$Term),] 
negative regulation of cold-induced thermogenesis (GO:0120163)    ADAMTS5;LAMA4

EnrichR_PreAs_perdepot$PR$aw$GO_Biological_Process_2021[grep("thermo",EnrichR_PreAs_perdepot$PR$aw$GO_Biological_Process_2021$Term),] 

EnrichR_PreAs_perdepot$PR$aw$GO_Biological_Process_2021[grep("resp",EnrichR_PreAs_perdepot$PR$aw$GO_Biological_Process_2021$Term),] 
response to growth factor (GO:0070848)

EnrichR_PreAs_perdepot$PR$aw$GO_Biological_Process_2021[grep("BMP",EnrichR_PreAs_perdepot$PR$aw$GO_Biological_Process_2021$Term),] 
negative regulation of BMP signaling pathway (GO:0030514) 0.1839703

EnrichR_PreAs_perdepot$PR$aw$GO_Biological_Process_2021[grep("hormon",EnrichR_PreAs_perdepot$PR$aw$GO_Biological_Process_2021$Term),] 
cellular response to growth hormone stimulus (GO:0071378)

EnrichR_PreAs_perdepot$PR$aw$GO_Biological_Process_2021[grep("grow",EnrichR_PreAs_perdepot$PR$aw$GO_Biological_Process_2021$Term),] 
insulin-like growth factor receptor signaling pathway (GO:0048009)

EnrichR_PreAs_perdepot$PR$aw$GO_Biological_Process_2021[grep("reti",EnrichR_PreAs_perdepot$PR$aw$GO_Biological_Process_2021$Term),] 
retinal metabolic process (GO:0042574)
retinol metabolic process (GO:0042572)

## SC
EnrichR_PreAs_perdepot$SC$aw$GO_Biological_Process_2021[grep("morpho",EnrichR_PreAs_perdepot$SC$aw$GO_Biological_Process_2021$Term),] 
long-chain fatty acid transport (GO:0015909)
positive regulation of lipid biosynthetic process (GO:0046889)
fatty acid transport (GO:0015908)
triglyceride metabolic process (GO:0006641)
cellular response to low-density lipoprotein particle stimulus (GO:0071404)
regulation of cell morphogenesis (GO:0022604)
positive regulation of cell morphogenesis involved in differentiation (GO:0010770)
response to glucocorticoid (GO:0051384)
negative regulation of chemokine production (GO:0032682)
positive regulation of extracellular matrix organization (GO:1903055)
extracellular matrix organization (GO:0030198)
positive regulation of cell differentiation (GO:0045597)
fat cell differentiation (GO:0045444)
regulation of fat cell differentiation (GO:0045598)
liver development (GO:0001889)
neuron development (GO:0048666)

## OM
EnrichR_PreAs_perdepot$EP$aw$GO_Biological_Process_2021[grep("imm",EnrichR_PreAs_perdepot$EP$aw$GO_Biological_Process_2021$Term),] 

regulation of insulin-like growth factor receptor signaling pathway (GO:0043567)
leukocyte chemotaxis involved in inflammatory response (GO:0002232)
regulation of chemokine production (GO:0032642)
regulation of humoral immune response (GO:0002920) 
positive regulation of interleukin-8 production (GO:0032757)
integrin-mediated signaling pathway (GO:0007229)
regulation of type 2 immune response (GO:0002828)
negative regulation of canonical Wnt signaling pathway (GO:0090090)
negative regulation of Wnt signaling pathway (GO:0030178)
immunoglobulin mediated immune response (GO:0016064)

## MC
EnrichR_PreAs_perdepot$MG$aw$GO_Biological_Process_2021[grep("epith",EnrichR_PreAs_perdepot$MG$aw$GO_Biological_Process_2021$Term),] 
response to unfolded protein (GO:0006986)
chaperone cofactor-dependent protein refolding (GO:0051085)
cellular response to topologically incorrect protein (GO:0035967)
regulation of transforming growth factor beta receptor signaling pathway (GO:0017015)



to_plot <- c( "positive regulation of extracellular matrix organization (GO:1903055)",
              "positive regulation of lipid biosynthetic process (GO:0046889)",
              "fatty acid transport (GO:0015908)",
              "triglyceride metabolic process (GO:0006641)",
              
              "regulation of cell morphogenesis (GO:0022604)",
              "positive regulation of cell morphogenesis involved in differentiation (GO:0010770)",
              "liver development (GO:0001889)",
              "neuron development (GO:0048666)",
              
              "positive regulation of cell differentiation (GO:0045597)",
              "fat cell differentiation (GO:0045444)", ## ------- SC
              
              
              "positive regulation of glucose import (GO:0046326)",
              "inositol lipid-mediated signaling (GO:0048017)",
              "response to growth factor (GO:0070848)",
              "cellular response to growth hormone stimulus (GO:0071378)",
              "insulin-like growth factor receptor signaling pathway (GO:0048009)",
              "retinal metabolic process (GO:0042574)",
              "retinol metabolic process (GO:0042572)",  ## ------ PR
              
              "regulation of insulin-like growth factor receptor signaling pathway (GO:0043567)",
              "leukocyte chemotaxis involved in inflammatory response (GO:0002232)",
              "regulation of chemokine production (GO:0032642)",
              "regulation of humoral immune response (GO:0002920)", 
              "regulation of type 2 immune response (GO:0002828)", 
              "negative regulation of Wnt signaling pathway (GO:0030178)",
              "immunoglobulin mediated immune response (GO:0016064)", ## ------ OM
              
              "response to unfolded protein (GO:0006986)",
              "chaperone cofactor-dependent protein refolding (GO:0051085)",
              "cellular response to topologically incorrect protein (GO:0035967)",
              "regulation of transforming growth factor beta receptor signaling pathway (GO:0017015)" ## ----- MC
)

out <- data.frame()
for(dep in c("SC", "PR", "EP", "MG")){
  o <- data.frame()
  for(t in to_plot){
    print(t)
    g <- which(EnrichR_PreAs_perdepot[[dep]]$aw$GO_Biological_Process_2021$Term == t)
    o <- rbind(o,EnrichR_PreAs_perdepot[[dep]]$aw$GO_Biological_Process_2021[g,])
  }
  o$Depot <- dep
  out <- rbind(out, o)
}

out$GOTerm <- paste0("(",sapply(strsplit(out$Term,"\\("), `[`, 2))
out$Term <- factor(out$Term, levels = unique(to_plot))
out <- out[order(out$Term),]
out$GR <- as.numeric(sapply(strsplit(out$Overlap,"/"), `[`, 1))/as.numeric(sapply(strsplit(out$Overlap,"/"), `[`, 2))
out$Depot <- factor(out$Depot, levels = c("SC", "PR", "EP", "MG"))
out$Adjusted.P.value <- as.numeric(out$Adjusted.P.value)
p <- ggplot(out, aes(x= Depot, y = Term, 
                     size = GR, col = Adjusted.P.value)) + 
  geom_point() + 
  scale_color_distiller(palette="RdYlBu", direction = 1) + 
  mashaGgplot2Theme
p
ggsave(p, filename = "6.2.Depot-specific-markers/DotPlot_GO_PreAs_depotSpecific_vnew.pdf",
       width = 8.93, height = 6.60)











EnrichR_ASCs_perdepot <- lapply(list("SC" = ASCs_SC, "EP" = ASCs_EP,
                                      "PR" = ASCs_PR, "MG" = ASCs_MG), 
                                 function(m){
                                   allG <- GetAssayData(int, assay = "RNA")[, int$myIntegratedClustering == "PreAs"]
                                   allG <- rownames(allG)[rowSums(allG) > 0]
                                   m <- order_g(d = m, m = "avg")
                                   out <- myEnrichRfunction(list("aw" = m), 
                                                            libraryName = c("GO_Biological_Process_2018",
                                                                            "GO_Biological_Process_2021",
                                                                            "WikiPathways_2016"))
                                   return(out)
                                 })


EnrichR_ASCs_perdepot$SC$aw$GO_Biological_Process_2021


### heatmap

ASCs_SC_f <- ASCs_SC %>% filter(signi.MG & signi.EP & signi.PR) %>% filter(up.MG & up.EP & up.PR)
ASCs_SC_f <- ASCs_SC_f[order(ASCs_SC_f$avg_log2FC.all, decreasing = T),]

PreAs_MG_f$rank.SC <- order(PreAs_MG_f$avg_log2FC.SC, decreasing = T)
PreAs_MG_f$rank.EP <- order(PreAs_MG_f$avg_log2FC.EP, decreasing = T)
PreAs_MG_f$rank.PR <- order(PreAs_MG_f$avg_log2FC.PR, decreasing = T)

PreAs_MG_f$final.rank <- rowSums(PreAs_MG_f[, c("rank.SC", "rank.EP", "rank.PR")])
PreAs_MG_f <- PreAs_MG_f[order(PreAs_MG_f$final.rank, decreasing = F), ]


int <- SetIdent(int, value = "myIntegratedClustering")
PreAs <- subset(int, idents = "PreAs")
PreAs <- SetIdent(PreAs, value = "depot")
avg.PreAs.cells <- log1p(AverageExpression(PreAs, verbose = FALSE)$RNA)
avg.PreAs.cells <- as.data.frame(avg.PreAs.cells)
avg.PreAs.cells$geneID <- data.annot[rownames(avg.PreAs.cells),"gene_short_name"]

to_plot <- c(PreAs_EP_f$ens_id[1:30], PreAs_SC_f$ens_id[1:30], PreAs_PR_f$ens_id[1:30], PreAs_MG_f$ens_id[1:30])
to_plot <- avg.PreAs.cells[to_plot, ]
rownames(to_plot) <- to_plot$geneID; to_plot <- to_plot[, c("EP","SC", "PR", "MG")]
pheatmap::pheatmap(to_plot, cluster_rows = F, cluster_cols = F)
pheatmap::pheatmap(to_plot, cluster_rows = F, cluster_cols = F, scale = "row", fontsize_row = 6)

to_plot <- c(PreAs_EP_f$ens_id[1:30], PreAs_SC_f$ens_id[1:30], PreAs_PR_f$ens_id[1:30], PreAs_MG_f$ens_id[1:30])
to_plot <- avg.PreAs.cells[to_plot, ]
rownames(to_plot) <- to_plot$geneID; to_plot <- to_plot[, c("EP","SC", "PR", "MG")]
pheatmap::pheatmap(to_plot, cluster_rows = F, cluster_cols = F)
pheatmap::pheatmap(to_plot, cluster_rows = F, cluster_cols = F, scale = "row", fontsize_row = 6)



avg.ASCs.cells$DE_all <- F
avg.ASCs.cells$DE_all[rownames(avg.ASCs.cells) %in% rownames(myDEGs$ASCs %>% filter(avg_logFC_all > 0))] <- T

avg.ASCs.cells$DE_SC <- avg.ASCs.cells$DE_EP <- avg.ASCs.cells$DE_PR <- avg.ASCs.cells$DE_MG <- F
avg.ASCs.cells$DE_SC[rownames(avg.ASCs.cells) %in% rownames(DE_res$ASCs$SC %>% filter(signi &  up))] <- T
avg.ASCs.cells$DE_EP[rownames(avg.ASCs.cells) %in% rownames(DE_res$ASCs$EP %>% filter(signi &  up))] <- T
avg.ASCs.cells$DE_PR[rownames(avg.ASCs.cells) %in% rownames(DE_res$ASCs$PR %>% filter(signi &  up))] <- T
avg.ASCs.cells$DE_MG[rownames(avg.ASCs.cells) %in% rownames(DE_res$ASCs$MG %>% filter(signi &  up))] <- T

avg.ASCs.cells_topPlot <- avg.ASCs.cells[rowSums(avg.ASCs.cells[, c("DE", "DE_SC", "DE_EP", "DE_PR", "DE_MG")]) > 0,]
pheatmap::pheatmap(avg.ASCs.cells[, c("DE", "DE_SC", "DE_EP", "DE_PR", "DE_MG")])

ggplot(avg.ASCs.cells, aes(x=EP,y = SC, col = DE)) + geom_point()
ggplot(avg.ASCs.cells, aes(x=EP,y = PR, , col = DE)) + geom_point()
ggplot(avg.ASCs.cells, aes(x=EP,y = MG)) + geom_point()
ggplot(avg.ASCs.cells, aes(x=SC,y = PR)) + geom_point()
ggplot(avg.ASCs.cells, aes(x=SC,y = MG)) + geom_point()
ggplot(avg.ASCs.cells, aes(x=PR,y = MG)) + geom_point()

avg.t.cells$gene <- rownames(avg.t.cells)

cd14.mono <- subset(immune.combined, idents = "CD14 Mono")
Idents(cd14.mono) <- "stim"
avg.cd14.mono <- as.data.frame(log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA))
avg.cd14.mono$gene <- rownames(avg.cd14.mono)

genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)

##---------------------------------------------##
##---------------------DEGs--------------------##
##---------------------------------------------##
DefaultAssay(int) <- "RNA"
IGFBP2_vs_ASPCs <- FindConservedMarkers(int, ident.1 = "IGFBP2", 
                                        ident.2 = c("ASCs", "PreAs", "CHI3L1-2", "IFIT", "HHIP"), 
                                        grouping.var = "batch", verbose = T, min.cells.group = 10)
IGFBP2_vs_ASPCs$geneID <- data.annot[rownames(IGFBP2_vs_ASPCs), "gene_short_name"]
IGFBP2_vs_Meso <- FindConservedMarkers(int, ident.1 = "IGFBP2", ident.2 = "Meso", 
                                       grouping.var = "batch",
                                       verbose = T, min.cells.group = 10)
IGFBP2_vs_Meso$geneID <- data.annot[rownames(IGFBP2_vs_Meso), "gene_short_name"]

save(IGFBP2_vs_ASPCs, IGFBP2_vs_Meso, file = "6.TopDEGs/DEGs_IGFBP2-vs-Meso-or-ASPCs.RData")

##---------------------------------------------##
##---------------------save--------------------##
##---------------------------------------------##

saveRDS(myDEGs_seurat, "6.TopDEGs/DEGs.Rds")










SC.vs.EP <- FindMarkers_by_depots("ASCs", "SC", all_dep_to_compare = "EP")
SC.vs.PR <- FindMarkers_by_depots("ASCs", "SC", all_dep_to_compare = "PR")
SC.vs.MG <- FindMarkers_by_depots("ASCs", "SC", all_dep_to_compare = "MG")


ASCs_EP.vs.PR <- FindMarkers_by_depots("ASCs", "SC", all_dep_to_compare = "EP")


PreAs_SC.vs.EP <- FindMarkers_by_depots("PreAs", "SC", all_dep_to_compare = "EP")
PreAs_SC.vs.PR <- FindMarkers_by_depots("PreAs", "SC", all_dep_to_compare = "PR")
PreAs_SC.vs.MG <- FindMarkers_by_depots("PreAs", "SC", all_dep_to_compare = "MG")


common <- PreAs_SC.vs.EP[PreAs_SC.vs.EP$signi,"geneID"][PreAs_SC.vs.EP[PreAs_SC.vs.EP$signi ,"geneID"] %in% PreAs_SC.vs.PR[PreAs_SC.vs.PR$signi,"geneID"]]
common <- common[common %in% PreAs_SC.vs.MG[PreAs_SC.vs.MG$signi,"geneID"]]

PreAs_SC.vs.EP[PreAs_SC.vs.EP$geneID %in% common,]
library(ggplot2)
library(cowplot)


Idents(cd14.mono) <- "stim"
avg.cd14.mono <- as.data.frame(log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA))


