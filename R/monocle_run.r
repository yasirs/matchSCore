monocle_run = function(sim,ntop,out_monocle){
  
  fd=rowData(sim)
  genes=fd$Gene
  
if(is.null(out_monocle)){
  
  pd=colData(sim)
  
  batch=factor(pd$Batch)
  cd=counts(sim)
  groups=factor(pd$Group)
  rownames(cd)=fd$Gene  
  colnames(cd)=pd$Cell
  
  k=length(levels(groups)) 
  
  #library(monocle)
  require(devtools)
  load_all("~/Downloads/monocle-release-master/")
  dir.create("out_Monocle/")
  outdir="out_Monocle"

  start.time <- Sys.time()

  pheno.data.df <- data.frame(pd)
  pheno <- new('AnnotatedDataFrame', data = pheno.data.df)
  features=data.frame(fd,gene_short_name=rownames(cd))
  rownames(features)=fd$Gene
  features=new('AnnotatedDataFrame', data = features)

  data <- newCellDataSet(cd, featureData = features,phenoData = pheno,expressionFamily = negbinomial())

  data=estimateSizeFactors(data)
  data=estimateDispersions(data)

  disp_table <- dispersionTable(data)

  unsup_clustering_genes <- subset(disp_table, mean_expression>quantile(disp_table$mean_expression,probs=0.05))
  ordering_genes=unsup_clustering_genes$gene_id
  data <- setOrderingFilter(data,ordering_genes)
  #plot_ordering_genes(data)



  data <- reduceDimension(data, max_components = 2,reduction_method = 'tSNE', verbose = T,residualModelFormulaStr="~Batch")
  #data <- clusterCells(data,num_clusters=k,method="DDRTree",verbose = T)

  data <- clusterCells(data,num_clusters=k+1,verbose = T)
  #data <- clusterCells(data,num_clusters=k,verbose = T)
  save(data,file="out_Monocle/data_last.RData")
  pd=pData(data)
  print(table(pd$Cluster))
#plot_cell_clusters(data, 1, 2, color = "Cluster")

# diff_test_res <- differentialGeneTest(data,fullModelFormulaStr = "~Cluster",cores = 4)
# save(diff_test_res,file="out_Monocle/de_test.RData")
# sig_genes <- subset(diff_test_res, qval < 0.1)
# save(sig_genes,file="out_Monocle/sig_genes.RData")
# 
#   library(R.utils)
#   gcDLLs()
  #detach("package:splatter", unload=TRUE)
  #detach("package:scater", unload=TRUE)
  
  library(Seurat)

  seu=exportCDS(data, export_to = c("Seurat"),export_all = T)
  clusters=pd$Cluster
  seu@ident=clusters
  names(seu@ident)=pd$Cell
  table(seu@ident)

  markers=FindAllMarkers(seu,logfc.threshold = 0,latent.vars="Size_Factor",test.use="MAST")
  

  end.time <- Sys.time()
  time=difftime(end.time,start.time,units="mins")


  save(markers, file = "out_Monocle/markers.RData")
  l=length(levels(clusters))
  gene_cl=lapply(1:l,function(x) which(genes %in% markers$gene[markers$cluster==x])[1:ntop])
  save(data,file="out_Monocle/data.RData")
  return(list(clusters=clusters,markers_df=markers,gene_cl=gene_cl,min=time))
  
}else{
  
  markers=out_monocle$markers_df
  clusters=out_monocle$clusters
  l=length(levels(clusters))
  gene_cl=lapply(1:l,function(x) which(genes %in% markers$gene[markers$cluster==x])[1:ntop])
  return(gene_cl)
  
  
  
}

}

