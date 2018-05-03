# load(file="sim_test_de0.5_bcv0.3.RData")
sc3_run = function(sim,ntop,out_sc3){
  

  
  
  fd=rowData(sim)
  
  if(is.null(out_sc3)){
  
    pd=colData(sim)
    batch=factor(pd$Batch)
    cd=counts(sim)
    groups=factor(pd$Group)
    rownames(cd)=fd$Gene  
    colnames(cd)=pd$Cell
    genes=rownames(cd)
    k=length(levels(groups))  
    
    rowData(sim)$feature_symbol=genes
    logcounts(sim)<- log2(cd+1)
    
  library(SC3)
  dir.create("out_SC3")
  out_dir="out_SC3"
  
  start.time <- Sys.time()
  
  sce <- sc3(sim, ks = k, biology = TRUE)
  
  end.time <- Sys.time()
  time=difftime(end.time,start.time,units="mins")
  
  save(sce,file=paste(out_dir,"sce.RData",sep="/"))
  pd_sce=colData(sce)
  j=grep("clusters",names(pd_sce))
  clusters=as.numeric(pd_sce[,j])
  fd_sce=data.frame(rowData(sce))
  m=grep("markers_clusts",names(fd_sce))
  pv=grep("de_padj",names(fd_sce))
  
  ord_fd=order(fd_sce[,m],fd_sce[,pv])
  ord_fd=fd_sce[ord_fd,]
  gene_cl=lapply(1:k,function(x) ord_fd$Gene[which(ord_fd[,m]==x)[1:ntop]])
  gene_cl=lapply(gene_cl,function(x) which(fd$Gene %in% x))

  return(list(clusters=clusters,fd.ord=ord_fd,gene_cl=gene_cl,min=time)) }
  
  else{
    ord_fd=out_sc3$fd.ord
    m=grep("markers_clusts",names(ord_fd))
    pv=grep("de_padj",names(ord_fd))
    k=max(ord_fd[,m],na.rm = T)
    gene_cl=lapply(1:k,function(x) ord_fd$Gene[which(ord_fd[,m]==x)[1:ntop]])
    gene_cl=lapply(gene_cl,function(x) which(fd$Gene %in% x))
    return(gene_cl) }
  
    
}