# Plot look and feel
mashaGgplot2Theme <- list(
  theme_classic(base_size = 14) + 
    theme(text = element_text(size = 14)) +
    theme(axis.line.x = element_line(colour = 'black', size = 0.5,
                                     linetype = 'solid'),
          axis.line.y = element_line(colour = 'black', size=0.5,
                                     linetype ='solid'),
          panel.grid.minor = element_line(colour = "white", size = 0.5,
                                          linetype = 2))
)

# Colors
# myColors <- c("ASCs" = "#33A02C", "PreAs" = "#E31A1C", "IGFBP2" = "#284724", "Meso" = "#810F7C", 
#               "VSMPs" = "#FC8D62", "Endo" = "darkgoldenrod1", "Immune" = "#E889BD", "Unknown_VSMPs" = "gray")
myColors <- c("ASCs" = "#33A02C", "PreAs" = "#E31A1C", "IGFBP2" = "#1f78b4", "Meso" = "#810F7C", 
              "VSMPs" = "#FC8D62", "Endo" = "darkgoldenrod1", "Immune" = "#E889BD", "Unknown_VSMPs" = "gray")

myBatchColors <- c("EP0" = "#cbc9e2", "EP1" = "#9e9ac8", "EP7" = "#6a51a3",
                   "MG7" = "#9ecae1", "MK7" = "#3182bd", 
                   "SC0" = "#F9F167", "SC1" = "gold", "SC7" = "orange1",
                   "PR3" = "#D38A59", "PR11" = "#A05B22", "PR12" = "#673400",
                   "GB7" = "#e78ac3")

myIntegratedColors <- c(myColors[c("ASCs", "PreAs", "IGFBP2")], "CILP" = "#a6d854",  "IFIT" = "gray",
                        "HHIP" = "#8da0cb", "CHI3L1-2" = "#fb9a99", "PR specific" = "sienna4", 
                        myColors[c("Meso", "Endo", "VSMPs", "Immune")],
                        "res.0.2_8" = "black")

myDepotsColors <- c("EP" = "#6F3996", "MG" = "#022AD7", "SC" = "#FEC010", "PR" = "#945200")

#df <- data.frame(x = c(1:length(myBatchColors)), y =rep(1,length(myBatchColors))); df$x <- as.factor(df$x) 
#ggplot(df, aes(x=x, y=y, fill = x)) + geom_bar(stat = "identity") + scale_fill_manual(values = unname(myBatchColors))
#Depots
myDepots <- c(paste0("EP", c(0,1,7)),
              "MG7", "MK7",
              paste0("SC", c(0,1,7)),
              paste0("PR", c(3,11,12)),
              "GB7")

#data.annot
data.annot <- read.table("~/SVRAW1/prainer/Files/Human/data.annot/Homo_sapiens.GRCh38.92_data.annot.txt")

# Functions
#' plot.colored_gene.continuous
#' @description plot reduction and colored by the gene expression of gene of choice
#' @author Vincent & Pernille
#' @param data.tsne tsne coordinates in a matrix i,j i cells, j the two tSNE to use
#' @param data.scale the expression matrix of the genes
#' @param data.annot a vector containing the ensembl names as row names and gene symbole as values
#' @param gene_name the ensembl name of the gene to plot
#' @param gene.text Boolean indicating if should print gene name on the plot
#' @param x.lab label y axis (nothing by default)
#' @param y.lab label y axis (nothing by default)
#' @param axis boolean indicating if should print axis or not 
#' @param col colors to use
#' @param size.point default 0.1, the size of the points
#' @param size.text default 4, size of the title 
#' @param pos.text default -20, 35
#' @param legend position of legend
#' @return ggplot
plot.colored_gene.continuous <- function(data.tsne, data.scale, data.annot, gene_name, 
                                         gene.text = T, x.lab = "", y.lab = "",axis = F,
                                         col = c('#e7e1ef','#d4b9da','#c994c7','#df65b0','#e7298a','#ce1256','#980043','#67001f'), 
                                         size.point = 0.1, size.text = 4, pos.text = c(-20, 35), legend = "right"){
  data.plot = data.frame(data.tsne)
  colnames(data.plot) = c("x","y")
  data.plot$gene = data.scale[gene_name,][rownames(data.plot)]
  data.plot$cell = rownames(data.plot)
  gene_hgnc <- as.character(data.annot[gene_name, "gene_short_name"])
  data.plot = data.plot[with(data.plot, order(gene)), ]
  data.text = data.frame(x=pos.text[1],y=pos.text[2],ens_id=gene_name, gene_id=gene_hgnc)
  p<- ggplot(data=data.plot, aes(x=x, y=y)) + 
    ggrastr::rasterise(geom_point(aes(col=gene), shape=19,  size=size.point), dpi = 650)  +
    scale_colour_gradientn(colours = col, aesthetics = "colour", name = "") +
    theme(legend.position = legend,
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
          axis.text=element_blank(), axis.ticks=element_blank(),  
          plot.background=element_blank(), panel.background=element_blank(),
          axis.line.x = element_line(), axis.line.y = element_line()) +
    xlab(label = x.lab) + ylab(label = y.lab) 
  if(axis == F){
    p <- p + theme( axis.line = element_blank()) + xlab(label = "") + ylab(label = "")
  }
  if(gene.text){
    p <- p + geom_text(data=data.text, aes(x=x, y=y, label=ens_id), inherit.aes = F, size = size.text/3) +
      geom_text(data=data.text, aes(x=x, y=y+5, label=gene_id), inherit.aes = F, size=size.text)
  }
  
  return(p)
}

#' plot.tsne.colored_continuous
#' @description Plot reduc dim colored by a continuous category
#' @author Vincent & Pernille
#' @param data.tsne tsne coordinates in a matrix i,j i cells, j the two tSNE to use
#' @param data.categories category to plot (vector)
#' @param name.category.plot tittle of the plot, default NULL
#' @param label.legend tittle legend
#' @param x.lab label y axis (nothing by default)
#' @param y.lab label y axis (nothing by default)
#' @param axis boolean indicating if should print axis or not 
#' @param col colors to use
#' @param size.point default 0.1, the size of the points
#' @param size.text default 4, size of the title 
#' @param pos.text default 0,0
#' @param legend position of legend
#' @param limits.scale default NULL 
#' @return ggplot
plot.tsne.colored_continuous <- function(data.tsne, data.categories, name.category.plot = "", 
                                         label.legend = NULL, axis = F, col = c('#f7fcf0','#e0f3db','#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#0868ac','#084081'),
                                         size.point = 0.1, size.text = 3, pos.text = c(0,0), legend = "right",x.lab = "", y.lab = "", 
                                         limits.scale = NULL, tittle = ""){
  
  data.plot = data.frame(data.tsne)
  colnames(data.plot) = c("x", "y")
  data.plot$category = data.categories[rownames(data.plot)]
  data.plot = data.plot[with(data.plot, order(category)), ]
  data.plot$cell = rownames(data.plot)
  data.text = data.frame(x=pos.text[1],y=pos.text[2],id=name.category.plot)
  col_to_use = col[1:length(table(data.categories))] #take the number of col depending on the number of categories (max 12)
  if(is.null(label.legend)){ label.legend = as.character(colnames(data.categories))}
  p<-ggplot(data=data.plot, aes(x=x, y=y, label=cell)) + 
    ggrastr::rasterise(geom_point(aes(col=category), shape=19,  size=size.point), dpi = 650) +
    scale_colour_gradientn(colours = col, aesthetics = "colour", name = label.legend, limits = limits.scale)+
    xlab(label = x.lab) + ylab(label = y.lab) +
    theme(axis.line = element_line(colour = "black"),legend.position = legend, 
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),   
          axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
          panel.background=element_blank()) + 
    ggtitle(tittle)
  if(axis == F){ p <- p + theme(axis.line = element_blank())+ xlab(label = "") + ylab(label = "")}
  return(p)
}

#' plot.tsne.colored_continuous_2
#' @description Plot reduc dim colored by a continuous category
#' @author Vincent & Pernille
#' @param data.tsne tsne coordinates in a matrix i,j i cells, j the two tSNE to use
#' @param data.categories category to plot (matrix or data.frame)
#' @param name.category.plot tittle of the plot, default NULL
#' @param label.legend tittle legend
#' @param x.lab label y axis (nothing by default)
#' @param y.lab label y axis (nothing by default)
#' @param axis boolean indicating if should print axis or not 
#' @param col colors to use
#' @param size.point default 0.1, the size of the points
#' @param size.text default 4, size of the title 
#' @param pos.text default 0,0
#' @param legend position of legend
#' @param limits.scale default NULL
#' @return ggplot
plot.tsne.colored_continuous_2 <- function(data.tsne, data.categories, name.category.plot = "", label.legend = NULL, 
                                           axis = F, col = c('#f7fcf0','#e0f3db','#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#0868ac','#084081'),
                                           size.point = 0.1, size.text = 3, pos.text = c(0,0), legend = "right",
                                           x.lab = "", y.lab = "", limits.scale=NULL){
  
  data.plot = data.frame(data.tsne)
  colnames(data.plot) = c("x", "y")
  data.plot$category = data.categories[rownames(data.plot),]
  data.plot = data.plot[with(data.plot, order(category)), ]
  data.plot$cell = rownames(data.plot)
  data.text = data.frame(x=pos.text[1],y=pos.text[2],id=name.category.plot)
  col_to_use = col[1:length(table(data.categories))] #take the number of col depending on the number of categories (max 12)
  if(is.null(label.legend)){ label.legend = as.character(colnames(data.categories))}
  p<-ggplot(data=data.plot, aes(x=x, y=y, label=cell)) + 
    ggrastr::rasterise(geom_point(aes(col=category), shape=19,  size=size.point), dpi = 650) +
    scale_colour_gradientn(colours = col, aesthetics = "colour", name = label.legend, limits = limits.scale)+
    xlab(label = x.lab) + ylab(label = y.lab) +
    theme(axis.line = element_line(colour = "black"),legend.position = legend, 
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),   
          axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
          panel.background=element_blank())
  if(axis == F){ p <- p + theme(axis.line = element_blank())+ xlab(label = "") + ylab(label = "")}
  return(p)
}

#' plot.tsne.colored_discrete
#' @description Plot reduc dim colored by discrete categories
#' @author Vincent & Pernille
#' @param data.tsne tsne coordinates in a matrix i,j i cells, j the two tSNE to use
#' @param data.categories category to plot (vector)
#' @param name.category.plot tittle of the plot, default NULL
#' @param label.legend tittle legend
#' @param x.lab label y axis (nothing by default)
#' @param y.lab label y axis (nothing by default)
#' @param axis boolean indicating if should print axis or not 
#' @param col colors to use
#' @param size.point default 0.1, the size of the points
#' @param size.text default 4, size of the title 
#' @param pos.text default 0,0
#' @param legend position of legend
#' @param limits.scale default NULL
#' @param if.NA indicates if there is some cells attributed to category NA, if yes they are colored in gray, default F
#' @return ggplot
plot.tsne.colored_discrete <- function(data.tsne, data.categories, name.category.plot = "", 
                                       axis = F, x.lab = "", y.lab = "", labels.category = NULL, 
                                       size.point = 0.1, size.text = 3, pos.text = c(-20, 35), 
                                       col = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'),
                                       legend = "right", leg.point.size = 3){
  if(name.category.plot == ""){ pos.text = c(0,0) }
  data.plot = data.frame(data.tsne)
  colnames(data.plot) = c("x", "y")
  data.plot$category = as.integer(data.categories[rownames(data.plot),])
  data.plot$cell = rownames(data.plot)
  data.text = data.frame(x=pos.text[1],y=pos.text[2],id=name.category.plot)
  if( is.null(labels.category)) { labels.category = names(table(data.categories)) }
  col_to_use = col[1:length(table(data.categories))] #take the number of col depending on the number of categories (max 12)
  p<-ggplot(data=data.plot %>% arrange(desc(is.na(data.plot$category))), aes(x=x, y=y, label=cell)) + 
    ggrastr::rasterise(geom_point(aes(col=as.factor(category)), shape=19,  size=size.point), dpi = 650) +
    scale_discrete_manual(values = col_to_use, aesthetics = "colour", label=labels.category, name = as.character(colnames(data.categories)), na.value = "gray") +
    theme(axis.line = element_line(colour = "black"),legend.position = legend, 
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
          axis.text=element_blank(), axis.ticks=element_blank(), 
          plot.background=element_blank(), panel.background=element_blank()) +
    geom_text(data=data.text, aes(x=x, y=y, label=id), inherit.aes = F, size = size.text) +
    xlab(label = x.lab) + ylab(label = y.lab) +
    theme(legend.key=element_blank()) + 
    guides(colour = guide_legend(override.aes = list(size=leg.point.size)))
  if(axis == F){ p <- p + theme(axis.line = element_blank())+ xlab(label = "") + ylab(label = "")}
  
  return(p)
}

#' plot.tsne.colored_discrete_2
#' @description Plot reduc dim colored by discrete categories
#' @author Vincent & Pernille
#' @param data.tsne tsne coordinates in a matrix i,j i cells, j the two tSNE to use
#' @param data.categories category to plot (matrix or data.frame)
#' @param name.category.plot tittle of the plot, default NULL
#' @param label.legend tittle legend
#' @param x.lab label y axis (nothing by default)
#' @param y.lab label y axis (nothing by default)
#' @param axis boolean indicating if should print axis or not 
#' @param col colors to use
#' @param size.point default 0.1, the size of the points
#' @param size.text default 4, size of the title 
#' @param pos.text default 0,0
#' @param legend position of legend
#' @param limits.scale default NULL
#' @return ggplot
plot.tsne.colored_discrete_2 <- function(data.tsne, data.categories, name.category.plot = "", 
                                         axis = F, x.lab = "", y.lab = "", labels.category = NULL, 
                                         size.point = 0.1, size.text = 3, pos.text = c(-20, 35), 
                                         col = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'),
                                         legend = "right", leg.point.size = 3, rasterise_points = T){
  if(name.category.plot == ""){ pos.text = c(0,0) }
  data.plot = data.frame(data.tsne)
  colnames(data.plot) = c("x", "y")
  data.plot$category = as.integer(data.categories[rownames(data.plot)])
  data.plot$cell = rownames(data.plot)
  data.text = data.frame(x=pos.text[1],y=pos.text[2],id=name.category.plot)
  if( is.null(labels.category)){ labels.category = names(table(data.categories)) }
  col_to_use = col[1:length(table(data.categories))]  #take the number of col depending on the number of categories (max 12)
  p<-ggplot(data=data.plot %>% arrange(desc(is.na(data.plot$category))), aes(x=x, y=y, label=cell)) + 
    ggrastr::rasterise(geom_point(aes(col=as.factor(category)), shape=19,  size=size.point), dpi = 650) +
    scale_discrete_manual(values = col_to_use, aesthetics = "colour", label=labels.category, name = as.character(colnames(data.categories)), na.value = "gray") +
    theme(axis.line = element_line(colour = "black"),legend.position = legend, panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), plot.background=element_blank(), panel.background=element_blank()) +
    geom_text(data=data.text, aes(x=x, y=y, label=id), inherit.aes = F, size = size.text) +
    xlab(label = x.lab) + ylab(label = y.lab) +
    theme(legend.key=element_blank()) + 
    guides(colour = guide_legend(override.aes = list(size=leg.point.size)))
  if(axis == F){ p <- p + theme(axis.line = element_blank())+ xlab(label = "") + ylab(label = "")}
  
  return(p)
}


#' convert_geneID_to_data.annot
#' @description convert a vector of gene ID to the equivalent data.annot
#' @param genesID       
#' @param data.annot
#' @return the data.annot of the genes 
convert_geneID_to_data.annot <- function(genesID, data.annot){
  out <- data.annot[data.annot$gene_short_name %in% genesID, ]
  rownames(out) <- out$gene_short_name
  
  missing <- genesID[!genesID %in% as.character(data.annot$gene_short_name)]
  for(m in missing){
    ens <- ""
    while(!ens %in% rownames(data.annot)){
      ens <- readline(prompt=paste("Gene", m, "is not found, please enter ens_id (without quote):"))
    }
    to_add <- data.annot[ens,]
    rownames(to_add) <- m
    out <- rbind(out, to_add)
  }
  #genesID <- genesID[genesID %in% as.character(data.annot$gene_short_name)]
  out <- out[genesID, ]
  rownames(out) <- out$ens_id
  return(out)
}

get_ensid <- function(g){
  d <- data.annot[data.annot$gene_short_name %in% g, ]
  rownames(d) <- d$gene_short_name
  d <- d[g,"ens_id"]
  return(d)
}
