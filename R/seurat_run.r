#' Run Seurat for a simulated data set 
#'
#' This function runs the Seurat package for a data set simulated by using Splatter
#' @param sim A Splatter sim object with simulated cellular groups or paths.
#' @param ntop The absolute number of top ranked markers per cluster. 
#' It can be set according to the number of sim.markers (based on the spec parameter in the markers_by_specificity).
#' @param out_seu The Seurat output of seurat_run. If provided, it recompute the set of group markers, 
#' without rerun Seurat. It is helpful when you need to change the ntop parameter.
#' @param res The value of resolution in the Seurat FindClusters function. Default is set to 0.8.
#' @param dims.use The vector with the dimensions to use in the Seurat FindClusters function. Default is set to 1:2.
#' @param test.de The test used to run the Seurat FindAllMarkers function. 
#' @return A list with clusters, the output of the Seurat FindAllMarkers function, the set of genes per cluster and the running time. 
#' If out_seu is provided, only the set of genes per cluster is returned.
#' @export
#' @examples
#' 

seurat_run = function(sim,ntop,out_seu=NULL,res=NULL,dims.use=NULL,test.de){
 
  fd=rowData(sim)
  genes=fd$Gene
  pd=colData(sim)
  
  if(is.null(res)){
    res=0.5
  }
  
  if(is.null(dims.use)){
    dims.use=c(1:2)
    
  }
  
if(is.null(out_seu)){
  
  
  batch=factor(pd$Batch)
  cd=counts(sim)
  groups=factor(pd$Group)
  rownames(cd)=fd$Gene  
  colnames(cd)=pd$Cell
  k=length(levels(groups)) 
  
  library(Seurat)
  #library(Matrix)
  
  dir.create("out_Seurat/")
  outdir="out_Seurat"  
  cd.sparse=Matrix(data=cd,sparse=T)
  
  start.time <- Sys.time()
  
  data <- CreateSeuratObject(raw.data = cd.sparse, min.cells = 0, min.genes = 0, project = "benchmarking")
  data <- NormalizeData(object = data, normalization.method = "LogNormalize", scale.factor = 10000)
  data <- FindVariableGenes(object = data, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.01, x.high.cutoff = 8, y.cutoff = 1)
  data <- AddMetaData(object = data, metadata = pd[,c(2,3)], col.name = c("Batch","Group"))
  data <- ScaleData(object = data, vars.to.regress = c("Batch"))
  data <- RunPCA(object = data, pc.genes = data@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
  data <- ProjectPCA(object = data, do.print = F)
  #data <- JackStraw(object = data, num.replicate = 100, do.print = FALSE)
  #PCElbowPlot(object = data)
  
  data <- FindClusters(object = data, reduction.type = "pca",resolution = res, dims.use=dims.use,print.output = 0, save.SNN = TRUE,force.recalc=T)
  
  print(table(data@ident))
  data <- RunTSNE(object = data, dims.use = dims.use, do.fast = TRUE,force.recalc=T)
  TSNEPlot(object = data)
  
  markers=FindAllMarkers(data,test.use = test.de,logfc.threshold = 0,only.pos = T)
  
  end.time <- Sys.time()
  time=difftime(end.time,start.time,units="mins")
  
  #save(markers, file = "out_Seurat/markers.RData")
  clusters=data@ident
  l=length(levels(clusters))-1
  gene_cl=lapply(0:l,function(x) which(genes %in% markers$gene[markers$cluster==x])[1:ntop])
  #save(data,file="out_Seurat/data.RData")
  return(list(clusters=clusters,markers_df=markers,gene_cl=gene_cl,min=time))
}else{
  
  markers=out_seu$markers_df
  clusters=out_seu$clusters
  l=length(levels(clusters))-1
  gene_cl=lapply(0:l,function(x) which(genes %in% markers$gene[markers$cluster==x])[1:ntop])
  return(gene_cl)
  
}
  
}

