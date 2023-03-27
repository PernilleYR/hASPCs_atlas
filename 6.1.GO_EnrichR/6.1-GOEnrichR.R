#################################################
#                                               #
#         Gene Set Enrichment Analysis          #
#                                               #
#################################################

########### FILE DESCRIPTION

# NAME: GSEA_script.R
# DATE: 5 december 2020
# AUTHOR: Pernille
#
# DESCRIPTION: Perform GSEA for Manuscript Aregs


library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)

##---------------------------------------------##
##------------Load and prepare data------------##
##---------------------------------------------##
data.annot <- read.table("/Volumes/UPDEPLA/prainer/Files/Human/data.annot/Homo_sapiens.GRCh38.92_data.annot.txt")

myDEGs <- readRDS("/Volumes/UPDEPLA/prainer/hASPCs/PAPER/10X_scRNA-seq/6.TopDEGs/DEGs.Rds")
for(n in names(myDEGs)){
  myDEGs[[n]] <- myDEGs[[n]] %>% filter(minimump_p_val < 0.05) %>% filter(avg_logFC_all > 0)
}

library(dplyr)
int <- readRDS("/Volumes/UPDEPLA/prainer/hASPCs/PAPER/10X_scRNA-seq/5.Integration/output/Seurat_2000HVGs.Rds")


GO_res_2 <- lapply(myDEGs, function(m){
  t <- m %>%  filter(minimump_p_val < 0.05) %>% filter(avg_logFC_all > 0)
  out <- GOenrichment(rownames(t), allGenesList = rownames(GetAssayData(int, assay = "RNA")))
  return(out)
})

##---------------------------------------------##
##------------EnrichR------------##
##---------------------------------------------##
enrichR_res_test <- myEnrichRfunction(myDEGs_seurat_f, n_genes = 1000, libraryName = c("GO_Biological_Process_2021"))
enrichR_res_tabula_sapiens <- myEnrichRfunction(myDEGs, n_genes = 1000, libraryName = "WikiPathways_2019_Human")

enrichR_res_IGFBP2 <- myEnrichRfunction(list(IGFBP2 = IGFBP2_EdgeR), n_genes = 1000, libraryName = c("GO_Biological_Process_2021"))

saveRDS(enrichR_res_IGFBP2, "6.1.GO_EnirhcR/EnrichR_res_IGFBP2-EdgeR.Rds")

##---------------------------------------------##
##-------------------Explore-------------------##
##---------------------------------------------##
#### ASCs
EnrichR_res$ASCs$GO_Biological_Process_2018
3                                                        extracellular matrix organization (GO:0030198)   9/230 5.522733e-05       0.01809753           0                    0   5.751457       56.38758
6                                                            extracellular matrix assembly (GO:0085029)    3/16 2.071823e-04       0.03753452           0                    0  31.574005      267.80792
14                                                                   bone cell development (GO:0098751)     2/7 1.114725e-03       0.07573163           0                    0  54.375342      369.70597
20                   negative regulation of cell morphogenesis involved in differentiation (GO:0010771)    3/31 1.533243e-03       0.08190597           0                    0  14.648276       94.92625
27                                                  positive regulation of lipase activity (GO:0060193)    2/10 2.354117e-03       0.08190597           0                    0  33.979452      205.62969
25                                                                        bone development (GO:0060348)    3/35 2.184800e-03       0.08190597           0                    0  12.814655       78.50553

EnrichR_res$ASCs$GO_Biological_Process_2018[grep("deve", EnrichR_res$ASCs$GO_Biological_Process_2018$Term),]
52                                                            substantia nigra development (GO:0021762)    3/44 0.004211774       0.08804227           0                    0  9.9971405    54.68307290
57                                                                myeloid cell development (GO:0061515)    2/15 0.005361340       0.09752873           0                    0 20.9051633   109.30351151

EnrichR_res$ASCs$GO_Biological_Process_2018[grep("stem", EnrichR_res$ASCs$GO_Biological_Process_2018$Term),]
62                                 regulation of stem cell proliferation (GO:0072091)    2/16 0.006097651        0.1069056           0                    0 19.4109589    98.99301104
EnrichR_res$ASCs$GO_Biological_Process_2018[grep("diff", EnrichR_res$ASCs$GO_Biological_Process_2018$Term),]
20                   negative regulation of cell morphogenesis involved in differentiation (GO:0010771)    3/31 0.001533243       0.08190597           0                    0 14.6482759     94.9262488
127                                        negative regulation of fat cell differentiation (GO:0045599)    2/33 0.024720527       0.20615021           0                    0  8.7587274     32.4083540
EnrichR_res$ASCs$GO_Biological_Process_2018[grep("beta", EnrichR_res$ASCs$GO_Biological_Process_2018$Term),]
65                         cellular response to transforming growth factor beta stimulus (GO:0071560)   4/101 0.006750704        0.1128925           0                    0   5.657216       28.27538
74     negative regulation of transforming growth factor beta receptor signaling pathway (GO:0030512)    3/56 0.008265456        0.1214129           0                    0   7.728953       37.06551
115                           transforming growth factor beta receptor signaling pathway (GO:0007179)    3/77 0.019499614        0.1827248           0                    0   5.529730       21.77254
119             regulation of transforming growth factor beta receptor signaling pathway (GO:0017015)    3/80 0.021555150        0.1968945           0                    0   5.313480       20.38857
167                    negative regulation of transforming growth factor beta production (GO:0071635)     1/6 0.043591786        0.2332262           0                    0  27.002721       84.59646
252                    positive regulation of transforming growth factor beta production (GO:0071636)     1/8 0.057698788        0.2387369           0                    0  19.285714       55.01287
264 negative regulation of cellular response to transforming growth factor beta stimulus (GO:1903845)    2/54 0.060623662        0.2387369           0                    0   5.216017       14.62086
EnrichR_res$ASCs$GO_Biological_Process_2018[grep("Wnt", EnrichR_res$ASCs$GO_Biological_Process_2018$Term),]
213                                   Wnt signaling pathway (GO:0016055)   3/110 0.04832862        0.2343818           0                    0  3.8179181     11.5672661 WNT10B;CPE;CD24
859            negative regulation of Wnt signaling pathway (GO:0030178)   2/175 0.37224085        0.4696455           0                    0  1.5582390      1.5398738  TAX1BP3;IGFBP6
907                         canonical Wnt signaling pathway (GO:0060070)    1/74 0.42342186        0.5057797           0                    0  1.8431647      1.5839905          WNT10B
952     Wnt signaling pathway, planar cell polarity pathway (GO:0060071)    1/89 0.48444519        0.5531428           0                    0  1.5278293      1.1072958            RHOA
989                     regulation of Wnt signaling pathway (GO:0030111)   1/110 0.55924581        0.6146615           0                    0  1.2321663      0.7160934         TAX1BP3
1008                    non-canonical Wnt signaling pathway (GO:0035567)   1/128 0.61471142        0.6622312           0                    0  1.0565644      0.5141267            RHOA
1031 negative regulation of canonical Wnt signaling pathway (GO:0090090)   1/149 0.67071371        0.7071443           0                    0  0.9056812      0.3617407          IGFBP6
1059          regulation of canonical Wnt signaling pathway (GO:0060828)   1/214 0.79770729        0.8187987           0                    0  0.6272237      0.1417610          IGFBP6

EnrichR_res$ASCs$GO_Biological_Process_2018[grep("morpho", EnrichR_res$ASCs$GO_Biological_Process_2018$Term),]
24                                          embryonic organ morphogenesis (GO:0048562)    3/34 0.002008340       0.08190597           0                    0  13.228699      82.156126
34                                            embryonic eye morphogenesis (GO:0048048)    2/11 0.002863318       0.08190597           0                    0  30.202435     176.858640

EnrichR_res$ASCs$GO_Biological_Process_2018[grep("prolif", EnrichR_res$ASCs$GO_Biological_Process_2018$Term),]
61                                   negative regulation of cell proliferation (GO:0008285)   8/364 0.005793947        0.1032462           0                    0   3.129374      16.119222
62                                       regulation of stem cell proliferation (GO:0072091)    2/16 0.006097651        0.1069056           0                    0  19.410959      98.993011
114                            positive regulation of fibroblast proliferation (GO:0048146)    2/29 0.019376884        0.1827248           0                    0  10.058346      39.666842
207                                     regulation of fibroblast proliferation (GO:0048145)    2/47 0.047341215        0.2343818           0                    0   6.029528      18.392316

EnrichR_res$ASCs$GO_Biological_Process_2018[grep("migratio", EnrichR_res$ASCs$GO_Biological_Process_2018$Term),]
29                                              regulation of cell migration (GO:0030334)   8/317 0.002528636       0.08190597           0                    0   3.614055      21.612319


### PreAs enrichR BP
EnrichR_res$PreAs$GO_Biological_Process_2018
1                                         positive regulation of fibroblast proliferation (GO:0048146)    3/29 3.028987e-05       0.01039960           0                    0  58.970414      613.56930
5                                                   positive regulation of cell migration (GO:0030335)   5/222 9.805774e-05       0.01189573           0                    0  12.293561      113.46900
15                                        positive regulation of myoblast differentiation (GO:0045663)    2/14 3.855442e-04       0.01634707           0                    0  83.108333      653.30253
39                                        positive regulation of fat cell differentiation (GO:0045600)    2/40 3.192378e-03       0.05075880           0                    0  26.210526      150.63162
62                                    regulation of mesenchymal stem cell differentiation (GO:2000739)     1/6 1.253548e-02       0.08976420           0                    0  97.331707      426.23429

EnrichR_res$PreAs$GO_Biological_Process_2018[grep("fat", EnrichR_res$PreAs$GO_Biological_Process_2018$Term),]
39               positive regulation of fat cell differentiation (GO:0045600)    2/40 0.003192378        0.0507588           0                    0  26.210526     150.631615 CEBPB;ZFP36
58                        regulation of fat cell differentiation (GO:0045598)    2/73 0.010295787        0.0897642           0                    0  14.004930      64.086845 CEBPB;ZFP36
92  regulation of transcription involved in cell fate commitment (GO:0060850)     1/7 0.014609764        0.0897642           0                    0  81.105691     342.757942       CEBPB
543                                 fatty acid metabolic process (GO:0006631)   1/107 0.201910526        0.2358593           0                    0   4.567879       7.308289        GGT5

EnrichR_res$PreAs$GO_Biological_Process_2018[grep("lipid", EnrichR_res$PreAs$GO_Biological_Process_2018$Term),]
40                      inositol lipid-mediated signaling (GO:0048017)    2/40 0.003192378       0.05075880           0                    0
97                                        lipid transport (GO:0006869)    2/92 0.015993376       0.08976420           0                    0

EnrichR_res$PreAs$GO_Biological_Process_2018[grep("insulin", EnrichR_res$PreAs$GO_Biological_Process_2018$Term),]
19           "regulation of insulin-like growth factor receptor signaling pathway" (GO:0043567)    2/17 0.0005739015       0.01895259           0                    0  66.476667      496.11887

EnrichR_res$PreAs$GO_Biological_Process_2018[grep("diff", EnrichR_res$PreAs$GO_Biological_Process_2018$Term),]
15                        positive regulation of myoblast differentiation (GO:0045663)    2/14 0.0003855442       0.01634707           0                    0  83.108333     653.302530
23                            positive regulation of cell differentiation (GO:0045597)   4/195 0.0007340396       0.02029779           0                    0  10.893910      78.620779
29                                 regulation of myoblast differentiation (GO:0045661)    2/30 0.0018041721       0.03923812           0                    0  35.589286     224.840775
35                     positive regulation of muscle cell differentiation (GO:0051149)    2/35 0.0024514342       0.04454606           0                    0  30.189394     181.470923
39                        positive regulation of fat cell differentiation (GO:0045600)    2/40 0.0031923775       0.05075880           0                    0  26.210526     150.631615
58                                 regulation of fat cell differentiation (GO:0045598)    2/73 0.0102957865       0.08976420           0                    0  14.004930      64.086845
62                    regulation of mesenchymal stem cell differentiation (GO:2000739)     1/6 0.0125354766       0.08976420           0                    0  97.331707     426.234286


9                          positive regulation of phosphatidylinositol 3-kinase signaling (GO:0014068)    3/54 1.982441e-04       0.01392963
16       negative regulation of platelet-derived growth factor receptor signaling pathway (GO:0010642)    2/15 4.442669e-04       0.01765961
19                    regulation of insulin-like growth factor receptor signaling pathway (GO:0043567)    2/17 5.739015e-04       0.01895259
22                                  regulation of phosphatidylinositol 3-kinase signaling (GO:0014066)    3/81 6.555927e-04       0.01895259
23                                            positive regulation of cell differentiation (GO:0045597)   4/195 7.340396e-04       0.02029779
39                                        positive regulation of fat cell differentiation (GO:0045600)    2/40 3.192378e-03       0.05075880
40                                                      inositol lipid-mediated signaling (GO:0048017)    2/40 3.192378e-03       0.05075880
58                                                 regulation of fat cell differentiation (GO:0045598)    2/73 1.029579e-02       0.08976420
62                                    regulation of mesenchymal stem cell differentiation (GO:2000739)     1/6 1.253548e-02       0.08976420
97                                                                        lipid transport (GO:0006869)    2/92 1.599338e-02       0.08976420

enrichr$PreAs$GO_Biological_Process_2018[grep("stem",enrichr$PreAs$GO_Biological_Process_2018$Term ),]
519                 regulation of stem cell differentiation (GO:2000736)    1/89 0.17097859       0.20952290           0  0

## IGFBP2
enrichR_res_IGFBP2$IGFBP2$GO_Biological_Process_2021[grep("(GO:0033629)",
                                                          enrichR_res_IGFBP2$IGFBP2$GO_Biological_Process_2021$Term), 
                                                     "Genes"]
# positive regulation of epithelial cell proliferation (GO:0050679) 
# "MDK;OSR1;FZD7;VASH2;IGF1;PTN;CLDN1;VEGFA"
# epithelium development (GO:0060429) 
# "WT1;OSR1;ALDH1A2;GATA6;SNAI2;WNT4;VEGFA"
# positive regulation of epithelial to mesenchymal transition (GO:0010718) 
# "CRB2;COL1A1;MDK;TWIST1"
# mesenchymal to epithelial transition (GO:0060231) 
# "WT1;FZD7"
# regulation of epithelial to mesenchymal transition (GO:0010717) 
# "CRB2;COL1A1;MDK;TWIST1"
# positive regulation of morphogenesis of an epithelium (GO:1905332) 
# "MDK;VEGFA"
# positive regulation of cell motility (GO:2000147) 
# "PDGFRA;SEMA3A;SEMA3B;TWIST1;IGF1;SOD2;CLDN1;VEGFA;COL1A1;CXCL12;MDK;CCL3;SNAI2;HAS2;ROR2"
# regulation of immune effector process (GO:0002697) 
# "C3;C6;CFH;C7;CFI;HLA-DRA;CLU;HLA-DRB1;C2"
# positive regulation of cell population proliferation (GO:0008284) 
# "PTGFR;CD74;PDGFRA;IL11RA;OSR1;CDCA7L;CLEC11A;IGF1;PTN;MEIS2;VEGFA;MEIS1;MDK;S100A13;S100A6;HAS2;TIMP1;ZFPM2;NES;NKX3-1"
# positive regulation of cell migration (GO:0030335) 
# "PDGFRA;SEMA3A;SEMA3B;IGF1;SOD2;CLDN1;VEGFA;COL1A1;CXCL12;MDK;CCL3;SNAI2;HAS2;ROR2"
# regulation of insulin-like growth factor receptor signaling pathway (GO:0043567) 
# "IGFBP4;IGFBP2;IGF1;NKX3-1"
# regulation of angiogenesis (GO:0045765) 
# "SPARC;CXCL8;PTGIS;GATA6;VASH2;CYP1B1;TWIST1;EMILIN1;GLUL;DCN;VEGFA"
# cellular response to chemokine (GO:1990869)
# "CXCL8;CCL8;CXCL12;DUSP1;CCL3;CXCL2"
# regulation of neutrophil chemotaxis (GO:0090022)
# "CD74;CXCL8;MDK;RAC2"
# positive regulation of vasculature development (GO:1904018)
# "CXCL8;PTGIS;VASH2;GATA6;CYP1B1;TWIST1;VEGFA"
# regulation of inflammatory response (GO:0050727)
# "NFKBIA;PTGIS;MDK;NFKBIZ;PLA2G2A;CCL3;APOA1;APOE;PTGS2"
# integrin-mediated signaling pathway (GO:0007229)
# "COL3A1;APOA1;FBLN1;PTN"
# negative regulation of cell adhesion mediated by integrin (GO:0033629)
# "SNAI2;CYP1B1"

##---------------------------------------------##
##-----------Dot Plot PreAs and ASCs-----------##
##---------------------------------------------##
selected_go_terms <- c("extracellular matrix organization (GO:0030198)",
                       "extracellular matrix assembly (GO:0085029)",
                       "bone cell development (GO:0098751)",
                       "bone development (GO:0060348)",
                       "myeloid cell development (GO:0061515)",
                       "substantia nigra development (GO:0021762)",
                       "negative regulation of cell morphogenesis involved in differentiation (GO:0010771)" ,
                       "embryonic organ morphogenesis (GO:0048562)",
                       "regulation of cell migration (GO:0030334)",
                       "positive regulation of lipase activity (GO:0060193)",
                       "regulation of stem cell proliferation (GO:0072091)",
                       "cellular response to transforming growth factor beta stimulus (GO:0071560)",
                       "negative regulation of fat cell differentiation (GO:0045599)",

                       "positive regulation of fibroblast proliferation (GO:0048146)",
                       "positive regulation of cell migration (GO:0030335)",
                       "positive regulation of myoblast differentiation (GO:0045663)",
                       "positive regulation of fat cell differentiation (GO:0045600)",
                       "regulation of mesenchymal stem cell differentiation (GO:2000739)",
                       "regulation of transcription involved in cell fate commitment (GO:0060850)",
                       "inositol lipid-mediated signaling (GO:0048017)",
                       "lipid transport (GO:0006869)",
                       "regulation of insulin-like growth factor receptor signaling pathway (GO:0043567)",
                       "positive regulation of phosphatidylinositol 3-kinase signaling (GO:0014068)"
)

list_GO_res <- list("ASCs" = EnrichR_res$ASCs$GO_Biological_Process_2018,
                    "PreAs" = EnrichR_res$PreAs$GO_Biological_Process_2018)

out <- data.frame()
for(res in names(list_GO_res)){
 g <- list_GO_res[[res]]$Term %in% selected_go_terms
 o <- list_GO_res[[res]][g,]
 o$comp <- res
 out <- rbind(out, o)
}

out$Term <- factor(out$Term, levels = selected_go_terms)
out <- out[order(out$Term),]

out$GR <- as.numeric(sapply(strsplit(out$Overlap,"/"), `[`, 1))/as.numeric(sapply(strsplit(out$Overlap,"/"), `[`, 2))
out$comp <- factor(out$comp, levels = names(list_GO_res))
out$Adjusted.P.value <- as.numeric(out$Adjusted.P.value)


myPallette <-
  c(rev(brewer.pal(5, "YlOrRd"))
    , "white"
    , brewer.pal(9, "Blues"))

p <- ggplot(out, aes(x= comp, y = Term, 
                     size = GR, col = Adjusted.P.value)) + 
  geom_point() + 
  #scale_color_manual(values = myPallette, drop = F) +
  scale_color_distiller(palette="RdYlBu", direction = 1, limits = c(0,0.22)) + 
  mashaGgplot2Theme + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(p, filename = "6.1.GO_EnirhcR/DotPlot_EnrichR_GO_2018_ASCs-PreAs.pdf",
       width = 8.6, height = 7.2)


### IGFBP2 enrichR BP

313                  retinal metabolic process (GO:0042574)    1/11 0.05887746        0.1946377           0                    0  18.238532
765         cellular response to retinoic acid (GO:0071300)    1/43 0.21132603        0.2986189           0                    0   4.335518

##---------------------------------------------##
##----------------Dot Plot IFIT----------------##
##---------------------------------------------##
enrichr$IFIT$GO_Biological_Process_2018 

df <- enrichr$IFIT$GO_Biological_Process_2018[, c( "Term", "Adjusted.P.value", "Old.Adjusted.P.value", "Overlap", "Combined.Score")]
df$overlap_TOP <- as.numeric(sapply(strsplit(df$Overlap, "/"), `[`, 1))
df$overlap_BOTTOM <- as.numeric(sapply(strsplit(df$Overlap, "/"), `[`, 2))
df$overlap_numeric <- df$overlap_TOP/df$overlap_BOTTOM
df$log_comb.score <- log(df$Combined.Score)
term_to_keep <- c(1,2,5,6,8,9,11,13,17,31)
df <- df[term_to_keep,]
df <- df[order(df$overlap_numeric),]
df$Term <- factor(df$Term, levels = df$Term)
p <- ggplot(df, 
            aes(x = overlap_numeric, y = Term)) + 
  geom_point(aes(size = log_comb.score, color = Adjusted.P.value)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient( low="red") +
  ylab(NULL) +
  ggtitle("IFIT")
ggsave(p, filename = "6.1.GO_EnirhcR/DotPlot_BP-enrichR_IFIT.pdf")
#Saving 11.7 x 3.91 in image

##---------------------------------------------##
##-----------------CILP/SFRP4------------------##
##---------------------------------------------##
EnrichR_res$CILP$GO_Biological_Process_2018[grep("Wnt", EnrichR_res$CILP$GO_Biological_Process_2018$Term),]
16  positive regulation of canonical Wnt signaling pathway (GO:0090263)   3/114 0.0007930178       0.01655425           0
22            positive regulation of Wnt signaling pathway (GO:0030177)   3/141 0.0014646670       0.02162518           0
34           regulation of canonical Wnt signaling pathway (GO:0060828)   3/214 0.0047667097       0.04630251           0
42                         canonical Wnt signaling pathway (GO:0060070)    2/74 0.0062345168       0.04957925           0
65                     regulation of Wnt signaling pathway (GO:0030111)   2/110 0.0133521253       0.06656134           0

4                                                        sulfur compound catabolic process (GO:0044273)    3/36 2.562681e-05     2.139838e-03
5                                                      glycosaminoglycan catabolic process (GO:0006027)    3/57 1.026515e-04     6.857119e-03
6                                                        keratan sulfate catabolic process (GO:0042340)    2/12 1.620532e-04     8.751450e-03

#SFRP4 and SFRP2:
23                                                 negative regulation of cellular process (GO:0048523)   5/535 1.489159e-03     2.162518e-02

ggplot(df[1:10,], aes(x = Combined.Score, y = Term, fill = Adjusted.P.value)) + geom_dotplot()

GO_res <- readRDS("6.1.GO_EnirhcR/GOEnrichment_results_elimFisher.Rds")
GO_res$IFIT

##---------------------------------------------##
##---------------Dot Plot IGFBP2---------------##
##---------------------------------------------##

## IGFBP2
# positive regulation of epithelial cell proliferation (GO:0050679) 
# "MDK;OSR1;FZD7;VASH2;IGF1;PTN;CLDN1;VEGFA"
# epithelium development (GO:0060429) 
# "WT1;OSR1;ALDH1A2;GATA6;SNAI2;WNT4;VEGFA"
# positive regulation of epithelial to mesenchymal transition (GO:0010718) 
# "CRB2;COL1A1;MDK;TWIST1"
# mesenchymal to epithelial transition (GO:0060231) 
# "WT1;FZD7"
# regulation of epithelial to mesenchymal transition (GO:0010717) 
# "CRB2;COL1A1;MDK;TWIST1"
# positive regulation of morphogenesis of an epithelium (GO:1905332) 
# "MDK;VEGFA"
# positive regulation of cell motility (GO:2000147) 
# "PDGFRA;SEMA3A;SEMA3B;TWIST1;IGF1;SOD2;CLDN1;VEGFA;COL1A1;CXCL12;MDK;CCL3;SNAI2;HAS2;ROR2"
# regulation of immune effector process (GO:0002697) 
# "C3;C6;CFH;C7;CFI;HLA-DRA;CLU;HLA-DRB1;C2"
# positive regulation of cell population proliferation (GO:0008284) 
# "PTGFR;CD74;PDGFRA;IL11RA;OSR1;CDCA7L;CLEC11A;IGF1;PTN;MEIS2;VEGFA;MEIS1;MDK;S100A13;S100A6;HAS2;TIMP1;ZFPM2;NES;NKX3-1"
# positive regulation of cell migration (GO:0030335) 
# "PDGFRA;SEMA3A;SEMA3B;IGF1;SOD2;CLDN1;VEGFA;COL1A1;CXCL12;MDK;CCL3;SNAI2;HAS2;ROR2"
# regulation of insulin-like growth factor receptor signaling pathway (GO:0043567) 
# "IGFBP4;IGFBP2;IGF1;NKX3-1"
# regulation of angiogenesis (GO:0045765) 
# "SPARC;CXCL8;PTGIS;GATA6;VASH2;CYP1B1;TWIST1;EMILIN1;GLUL;DCN;VEGFA"
# cellular response to chemokine (GO:1990869)
# "CXCL8;CCL8;CXCL12;DUSP1;CCL3;CXCL2"
# regulation of neutrophil chemotaxis (GO:0090022)
# "CD74;CXCL8;MDK;RAC2"
# positive regulation of vasculature development (GO:1904018)
# "CXCL8;PTGIS;VASH2;GATA6;CYP1B1;TWIST1;VEGFA"
# regulation of inflammatory response (GO:0050727)
# "NFKBIA;PTGIS;MDK;NFKBIZ;PLA2G2A;CCL3;APOA1;APOE;PTGS2"
# integrin-mediated signaling pathway (GO:0007229)
# "COL3A1;APOA1;FBLN1;PTN"
# negative regulation of cell adhesion mediated by integrin (GO:0033629)
# "SNAI2;CYP1B1"
# negative regulation of cell adhesion (GO:0007162)    6/73 0.0008506121       0.02667614
# cellular response to hypoxia (GO:0071456)   8/131 0.0008940171       0.02691588
# positive regulation of transcription from RNA polymerase II promoter in response to hypoxia (GO:0061419)     2/5 0.0022341012
# regulation of inflammatory response (GO:0050727)   9/206 0.004389029       0.05832799
# leukocyte chemotaxis involved in inflammatory response (GO:0002232)     2/6 0.003317546       0.05061499 -> MDK, PTN
# regulation of angiogenesis (GO:0045765)  11/203 0.0002899366       0.01582553
# positive regulation of angiogenesis (GO:0045766)   8/116 0.0003991764       0.01802680
# cellular response to transforming growth factor beta stimulus (GO:0071560)   6/114 0.007960878       0.07645637

selected_terms <- c("regulation of insulin-like growth factor receptor signaling pathway (GO:0043567)",
                    
                    "positive regulation of epithelial cell proliferation (GO:0050679)",
                    "epithelium development (GO:0060429)",
                    
                    "positive regulation of cell motility (GO:2000147)",
                    "positive regulation of cell population proliferation (GO:0008284)",
                    "positive regulation of cell migration (GO:0030335)",
                    
                    "negative regulation of cell adhesion (GO:0007162)",
                    
                    "positive regulation of epithelial to mesenchymal transition (GO:0010718)",
                    "mesenchymal to epithelial transition (GO:0060231)",
                    
                    "regulation of angiogenesis (GO:0045765)",
                    "positive regulation of transcription from RNA polymerase II promoter in response to hypoxia (GO:0061419)",
                    
                    "cellular response to hypoxia (GO:0071456)",
                    
                    "cellular response to transforming growth factor beta stimulus (GO:0071560)",
                    
                    "cellular response to chemokine (GO:1990869)",
                    "regulation of inflammatory response (GO:0050727)"
                    )


df <- enrichR_res_IGFBP2$IGFBP2$GO_Biological_Process_2021[, c( "Term", "P.value", "Adjusted.P.value", "Old.Adjusted.P.value", "Overlap", "Combined.Score")]
df <- df[df$Term %in% selected_terms,]

df$overlap_TOP <- as.numeric(sapply(strsplit(df$Overlap, "/"), `[`, 1))
df$overlap_BOTTOM <- as.numeric(sapply(strsplit(df$Overlap, "/"), `[`, 2))
df$overlap_numeric <- df$overlap_TOP/df$overlap_BOTTOM
#df$Adjusted.P.value <- round(df$Adjusted.P.value, 3)
# df <- df[order(df$overlap_numeric),]
# df$Term <- paste(sapply(strsplit(as.character(df$Term),"\\("), `[`, 1),
#                  "-", round(df$Adjusted.P.value,2))
df$Term <- factor(df$Term, levels = selected_terms)
p <- ggplot(df, aes(x= overlap_numeric, y = Term, 
                    size = Combined.Score, 
                    col = Adjusted.P.value)) + 
  geom_point() + 
  # scale_color_gradient2(low = "#d7191c",
  #                       mid = "#ffffbf",
  #                       high = "#2c7bb6",
  #   midpoint = 0.06) +
#  scale_color_gradientn(colors = myPallette) + 
scale_color_distiller(palette = "RdYlBu", direction = 1, limits = c(0,0.1)) + 
  mashaGgplot2Theme
p
ggsave(p, filename = "6.1.GO_EnirhcR/DotPlot_BP-enrichR_IGFBP2.pdf",
       width = 8.75, height = 3.52)
#Saving 11.7 x 3.91 in imag



myPallette <-
  c(rev(brewer.pal(9, "YlOrRd"))
    , brewer.pal(6, "Blues")[2:6])

p <- ggplot(df, aes(x= Combined.Score, y = Term, 
                    size = overlap_numeric, 
                    col = Adjusted.P.value)) + 
  geom_point() + 
  scale_color_manual(values = myPallette, drop = F) +
  #scale_color_distiller(palette="RdYlBu", direction = 1, limits = c(0,0.22)) + 
  mashaGgplot2Theme 
p
