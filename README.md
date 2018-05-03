## Introduction

matchSCore is a basic R package for the comparison of single cell RNA-seq data across tools and experiments.
The purpose of this vignette is to guide you to the use of the matchSCore with some examples. 
The package and the repository will be regularly updated with other reference data that you can use to annotate your clusters.  

## Installation

This is a development version of matchSCore running on R (>= 3.4.0).
matchSCore can be installed by devtools:

```{r,eval=FALSE}

library("devtools")
install_github('elimereu/matchSCore')


```


## The matchSCore metric

The matchSCore is a Jaccard Index based metric that enables the direct comparison of clusters predicted by a tool for a simulated (1) or real data set (2). It allows to:

1. Track the accuracy trend of a tool in clustering and marker identification compared with the optimal solution provided by a simulated data set. In this case, matchSCore works in combination with the Splatter package - https://github.com/Oshlack/splatter/blob/master/vignettes/splatter.Rmd - . (Benchmarking)

2. Score the matching between your clusters and cell identities from a reference data set. (Cluster annotation) 


##Benchmarking with a simulated data set


```{r,eval=FALSE}

library(matchSCore)
library(splatter)
library(scater)
library(MixSim)

## For example we create a simulated data set with 4 groups of equal proportions. 
sim <- splatSimulate(batchCells=rep(500,2),group.prob=rep(0.25,4),method = "groups")
n_groups <- 4

## Compute the ranking of genes per group
rank_df <- rank_sim(sim)

## Set the proportion of top ranked genes (specificity) we want to output
specificity=0.1
## feature info data frame
fd <- rowData(sim)
##List of markers per group
markers_pos <- markers_by_specificity(rank_df,specificity,n_groups) ## by position in fd$Gene

## List with the top 10% ranked markers per group
sim.markers <- lapply(markers_pos,function(x) fd$Gene[x])
sim.markers


```

For example, we could use Seurat to analyze this data set by running the function seurat_run

```{r,eval=FALSE}

##Run Seurat
out_seu <-seurat_run(sim,ntop=100,out_seu = NULL,res = NULL,dims.use = NULL,test.de = "wilcox")

pd=colData(sim)
lab.sim <- pd$Group
idc <- out_seu$clusters ## Seurat clusters
gene_cl <- out_seu$gene_cl ## Seurat cluster genes

## identify the group labels (defined by the simulated groups) for the clusters 
lab <- compute_labels(lab.sim,idc)


RandIndex(lab.sim,idc) ## we can compute the Rand Index (RI), adjusted Rand Index (ARI) or F score

## Now, we can compute the matchSCore value 
matchSCore(markers,gene_cl,lab)

```

##Clustering Annotation

The matchSCore matrix can be also used to match clusters from two different experiments.
For example, you could use the top100 ranked markers we got (by using Seurat) from the Smart-Seq2 and Chromium bladder sample from the Tabula Muris Atlas.
You can load files directly from the data folder in this repository. 

```{r,eval=FALSE}
## We use Smart-Seq2 data as the reference data
load(file="data/gene_cl.ref_bladder.RData")

## And Chromium data as test data
load(file="data/gene_cl.obs_bladder_droplet.RData")

## The matchSCore2 function computes the clustering comparison and produce the heatmap table with matchSCore values for each group combination
out=matchSCore2(gene_cl.ref,gene_cl.obs,tissue = "Bladder",ylab = "Smart-Seq2",xlab = "Chromium")

## The matchSCore heatmap is stored in the ggplot slot of out.  
out$ggplot

```

