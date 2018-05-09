#' Select the top n cluster markers  
#'
#' This function identifies true label groups between reference groups and clusters.
#' @param clusters A vector of cluster labels.
#' @param markers A data.frame as in the output of the FindAllMarkers Seurat function.
#' @param ntop The number of top markers you want.
#' @return A list of ntop markers per cluster.
#' 
#' 
#' 


cut_markers <- function(clusters,markers,ntop){
  
  #l=length(levels(clusters))-1
  # levels(clusters)=seq(0:l)-1
  lev=levels(clusters)
  gene_cl=lapply(lev,function(x) markers$gene[markers$cluster==x][1:ntop])

  return(gene_cl)
}



