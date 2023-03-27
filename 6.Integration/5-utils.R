#' prepareData
#' @description filter the matrices from list_matrices (keep cells from the corresponding
#' matric in list_original_data) and create seurat object, normalized & set variable genes from the original datasets
#' @param list_matrices list of matrices to prepare 
prepareData <- function(list_matrices, list_original_data){
  for(d in names(list_matrices)){
    print(d)
    list_matrices[[d]] <- list_matrices[[d]][, colnames(list_original_data[[d]])]
    colnames(list_matrices[[d]]) <- paste0(d, "_",colnames(list_matrices[[d]]))
    list_matrices[[d]] <- CreateSeuratObject(list_matrices[[d]])
    list_matrices[[d]] <- NormalizeData(list_matrices[[d]], normalization.method = "LogNormalize", scale.factor = 10000)
    VariableFeatures(list_matrices[[d]]) <- VariableFeatures(list_original_data[[d]])}
  return(list_matrices)
}
