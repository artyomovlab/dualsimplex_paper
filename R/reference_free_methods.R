library(TOAST)
dyn.load('/usr/lib/jvm/java-11-openjdk-amd64/lib/server/libjvm.so') # was not able to make rstudio like this
library(rJava)

## Toast
deconvolve_with_toast <-  function(mat, K, nmarker = dim(mat)[[1]]) {
  data_raw <-  as.matrix(mat)
  refinx <- findRefinx(data_raw, nmarker = nmarker)
  Y <- data_raw[refinx,]
  Y <-  sweep(Y,2,colSums(Y),`/`)
  outT <- myRefFreeCellMix(Y,  mu0=myRefFreeCellMixInitialize(Y, K = K), iters=500, verbose=0)
  solution_H <- t(outT$Omega)
  solution_W <- outT$Mu
  colnames(solution_H) <-  colnames(data_raw)
  rownames(solution_H) <-  paste("toast_cell_type", c(1:K)) 
  colnames(solution_W) <-  paste("toast_cell_type", c(1:K)) 
  #rownames(solution_W) <-  rownames(data_raw)
  return(list(object=outT, proportions=solution_H, basis=solution_W))
}

library(CAM3)
## CAM3
deconvolve_with_CAM3 <-  function(mat, K, thres.low=0.1, thres.high=0.85, radius.thres=0.9,MG.num.thres = 0) {
  data_raw <-  as.matrix(mat)
  
  Y <-  sweep(data_raw,2,colSums(data_raw),`/`)
  
  rCAM3 <- CAM3Run(Y, K = K, dim.rdc = K, thres.low = thres.low , cluster.num = 120, thres.high =thres.high,radius.thres=radius.thres, MG.num.thres = MG.num.thres)
  solution_H <-t(rCAM3@ASestResult[[1]]@Aest)
  solution_W <- rCAM3@ASestResult[[1]]@Sest
  
  colnames(solution_H) <-  colnames(data_raw)
  rownames(solution_H) <-  paste("cam3_cell_type", c(1:K)) 
  
  colnames(solution_W) <-  paste("cam3_cell_type", c(1:K)) 
  #rownames(solution_W) <-  rownames(data_raw)
  return(list(object=rCAM3, proportions=solution_H, basis= solution_W))
}


library(linseed)
## Linseed
deconvolve_with_Linseed <-  function(mat, K, topGenes = 10000) {
  data_raw <-  as.matrix(mat)
  lo <- LinseedObject$new(data_raw,  topGenes=topGenes)
  lo$calculatePairwiseLinearity()
  lo$calculateSpearmanCorrelation()
  lo$calculateSignificanceLevel(100)
  lo$significancePlot(0.01)
  lo$filterDatasetByPval(0.01)
  lo$setCellTypeNumber(K)
  lo$project("filtered")
  lo$smartSearchCorners(dataset="filtered", error="norm")
  solution_H <- lo$proportions
  colnames(solution_H) <-  colnames(data_raw)
  rownames(solution_H) <-  paste("linseed_cell_type", c(1:K)) 
  
  ## infer basis matrix from proportions
  ## Get W for the whole data based on proportions
  res <- DualSimplex::nnls_C__(t(solution_H), t(data_raw))
  solution_W <- t(res)
  colnames(solution_W) <-  paste("cam3_cell_type", c(1:K)) 
  rownames(solution_W) <-  rownames(data_raw)
  
  return(list(object=lo, proportions=solution_H, basis= solution_W))
}


## CellDistinguisher

library("CellDistinguisher")
library(CellMix)
library(NMF)
library(hNMF)
library(matrixStats)
#library(DelayedMatrixStats)

deconvolve_with_CellDistinguisher <-  function(mat, K) {
  data_raw <-  as.matrix(mat)
  inData <-  data_raw
  exprLinear <- inData[ rowSums(inData==0) != ncol(inData), ]
  distinguishers <- gecd_CellDistinguisher(exprLinear, genesymb=NULL, 
                                           numCellClasses=K, 
                                           minDistinguisherAlternatives=20, 
                                           maxDistinguisherAlternatives=100, 
                                           minAlternativesLengthsNormalized=0.5, 
                                           expressionQuantileForFilter=0.995, 
                                           expressionConcentrationRatio=0.333, verbose=0)
  
  
  deconvolutionSSKL <- gecd_DeconvolutionCellMix(exprLinear, distinguishers$bestDistinguishers, method="ssKL", maxIter=5)
  sampleCompositions <- deconvolutionSSKL$sampleCompositions
  
  
  solution_H <- sampleCompositions

  colnames(solution_H) <-  colnames(data_raw)
  rownames(solution_H) <-  paste("cd_cell_type", c(1:K)) 
  
  solution_W <-  deconvolutionSSKL$cellSubclassSignatures
  
  
  return(list(object=deconvolutionSSKL, proportions=solution_H, basis=solution_W))
}


