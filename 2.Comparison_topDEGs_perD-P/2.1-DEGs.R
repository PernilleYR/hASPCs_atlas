################################################################
#                                                              #
#              Differential Expression Analysis                #
#                                                              #
################################################################


### Author: Pernille
### Date: 10.08.2022
### Datasets: scRNA-seq Depots: SC - EP - MK - MG - PR - GB 
###                     Patients: B0, B1, B7, L3, L11, L12 (B - Bariatric, L - Lean)
### Goal: Compute DEGs

library(Seurat); library(ggrastr)

setwd("~/SVRAW1/prainer/hASPCs/PAPER/10X_scRNA-seq/")
out_dir = "~/SVRAW1/prainer/hASPCs/PAPER/10X_scRNA-seq/2.Comparison_topDEGs_perD-P/DEGs/"
dir.create(out_dir)

##---------------------------------------------##
##----------------Loading data-----------------##
##---------------------------------------------##

data.annot <- read.table("/home/pyrainer/SVRAW1/prainer/Files/Human/data.annot/Homo_sapiens.GRCh38.92_data.annot.txt")
myseu <- load("0.data/Seurats_objects.rds")

##---------------------------------------------##
##------------------Function-------------------##
##---------------------------------------------##

## Output Summary: 
# p_val : p_val (unadjusted)
# avg_logFC : log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group.
# pct.1 : The percentage of cells where the gene is detected in the first group
# pct.2 : The percentage of cells where the gene is detected in the second group
# p_val_adj : Adjusted p-value, based on BONFERRONI correction using all genes in the dataset.

## test.use:
# roc : Standard AUC classifier
# bimod :  Likelihood-ratio test for single cell gene expression, (McDavid et al., Bioinformatics, 2013)
# wilcox : Wilcoxon rank sum test (default)
# t : Student's t-test
# tobit : Tobit-test for differential gene expression (Trapnell et al., Nature Biotech, 2014)
# poisson, negbinom, MAST, DESeq2 ... 

## Other param:
# min.pct : filter the genes that are not expressed in at least min.pct (min % cells of each pop to compare)
# logfc.threshold : filter the genes that don't have at least logfc.threshold of FC, exemple log(2)
# only.pos : Return only positive markers (default = F)
# thresh.use : Test only genes whose average expression is > thresh.use between cluster (faster for higher values)
# return.thresh = 0.01 : Only return markers that have a p-value (or AUC) > return.tresh

myDEGs_f <- function(d, padj_cutoff = 0.05, log2fc_cutoff = log2(1.2)){
  s <- myseu[[d]]
  s <- SetIdent(s, value = "mySelectedClustering")
  
  #Tests
  bimod <- FindAllMarkers( object = s, test.use = "bimod", only.pos = T, logfc.threshold = log2fc_cutoff)
  wilcox <- FindAllMarkers( object = s, test.use = "wilcox", only.pos = T, logfc.threshold = log2fc_cutoff)
  
  #Significant
  bimod_s <- bimod[which(bimod$p_val_adj < padj_cutoff),]
  wilcox_s <- wilcox[which(wilcox$p_val_adj < padj_cutoff),]
  
  #Common across tests
  print(d)
  common <- list()
  for(i in levels(s$mySelectedClustering)){
    b <- subset(bimod_s, cluster == i)$gene
    w <- subset(wilcox_s, cluster == i)$gene
    common[[i]] <- b[b %in% w]
    print(paste0(i,": ",length(common[[i]])))
  }
  
  #Output
  output <- list()
  for(i in levels(s$mySelectedClustering)){
    output[[i]] <- wilcox_s[wilcox_s$gene %in% common[[i]],]
    output[[i]] <- subset(output[[i]], cluster == i)
    output[[i]]$geneID <- data.annot[output[[i]]$gene, "gene_short_name"]
    output[[i]] <- output[[i]][order(output[[i]]$avg_log2FC, decreasing = T),]
    rownames(output[[i]]) <- output[[i]]$gene
  }
  
  #Save text
  dir.create(paste0(out_dir, d))
  for(i in levels(s$mySelectedClustering)){
    write.table(output[[i]], file = paste0(out_dir, d, "/Markers_", d, "_",i,".txt"))
    write.table(output[[i]]$geneID, file = paste0(out_dir, d, "/Markers_", d, "_",i,"_ID.txt"), 
                quote = F, row.names = F, col.names = F)
  }
  
  #SaveRDS
  saveRDS(output, paste0(out_dir, d, "/Markers_", d,".rds") )
}

myDEGs_f <- function(d, padj_cutoff = 0.05, log2fc_cutoff = log(1.5)){
  s <- myseu[[d]]
  s <- SetIdent(s, value = "mySelectedClustering")
  
  #Tests
  wilcox <- FindAllMarkers( object = s, test.use = "wilcox", only.pos = T, logfc.threshold = log2fc_cutoff )
  
  #Significant
  wilcox_s <- wilcox[wilcox$p_val_adj < padj_cutoff,]

  #Output
  output <- list()
  for(i in levels(s$mySelectedClustering)){
    output[[i]] <- subset(wilcox_s, cluster == i)
    output[[i]]$geneID <- data.annot[output[[i]]$gene, "gene_short_name"]
    output[[i]] <- output[[i]][order(output[[i]]$avg_log2FC, decreasing = T),]
    rownames(output[[i]]) <- output[[i]]$gene
  }
  
  #Save text
  dir.create(paste0(out_dir, d, "_wilcox"))
  for(i in levels(s$mySelectedClustering)){
    write.table(output[[i]], file = paste0(out_dir, d, "_wilcox/Markers_", d, "_",i,".txt"))
    write.table(output[[i]]$geneID, file = paste0(out_dir, d, "_wilcox/Markers_", d, "_",i,"_ID.txt"), 
                quote = F, row.names = F, col.names = F)
  }
  
  #SaveRDS
  saveRDS(output, paste0(out_dir, d, "_wilcox/Markers_", d,".rds") )
}

##---------------------------------------------##
##--------------------DEGs---------------------##
##---------------------------------------------##

myDepots <- c(paste0("EP", c(0,1,7)),
              paste0("SC", c(0,1,7)),
              "MG7", "MK7",
              "GB7",
              paste0("PR", c(3,11,12)))
lapply(myDepots, myDEGs_f)
# [1] "EP0"
# [1] "Meso: 428"
# [1] "ASCs: 461"
# [1] "IGFBP2: 206"
# [1] "PreAs: 335"
# [1] "VSMPs: 310"

# [1] "EP1"
# [1] "PreAs: 478"
# [1] "Meso: 499"
# [1] "ASCs: 433"
# [1] "VSMPs: 504"
# [1] "IGFBP2: 44"

# [1] "EP7"
# [1] "ASCs: 525"
# [1] "Endo: 1367"
# [1] "IGFBP2: 63"
# [1] "Immune: 1008"
# [1] "Meso: 704"
# [1] "PreAs: 506"
# [1] "VSMPs: 275"

# [1] "SC0"
# [1] "PreAs: 198"
# [1] "ASCs: 229"
# [1] "VSMPs: 402"

# [1] "SC1"
# [1] "PreAs: 287"
# [1] "ASCs: 385"
# [1] "VSMPs: 454"

# [1] "SC7"
# [1] "PreAs: 339"
# [1] "ASCs: 511"
# [1] "VSMPs: 860"
# [1] "Endo: 651"
# [1] "Immune: 547"

# [1] "MG7"
# [1] "ASCs: 522"
# [1] "Endo: 1229"
# [1] "Immune: 737"
# [1] "PreAs: 391"
# [1] "VSMPs: 637"

# [1] "MK7"
# [1] "PreAs: 340"
# [1] "ASCs: 393"
# [1] "Endo: 895"
# [1] "VSMPs: 297"
# [1] "Immune: 173"

# [1] "GB7"
# [1] "ASCs: 3"
# [1] "PreAs: 10"

# [1] "PR3"
# [1] "ASCs: 589"
# [1] "PreAs: 348"
# [1] "VSMPs: 192"

# [1] "PR3"
# [1] "ASCs: 589"
# [1] "PreAs: 348"
# [1] "VSMPs: 192"

# [1] "PR11"
# [1] "PreAs: 409"
# [1] "ASCs: 693"
# [1] "VSMPs: 895"
# [1] "Unknown_VSMPs: 176"

# [1] "PR12"
# [1] "ASCs: 371"
# [1] "PreAs: 319"
# [1] "VSMPs: 648"