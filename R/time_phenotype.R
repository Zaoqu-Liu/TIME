# -------------------------------------------------------------------------
#'@title  time_phenotype
#'
#'@description  Get the TIME phenotypes and TIMEscore
#'
#'@details  Use gene expression data to get the corresponding TIME phenotypes and TIMEscore
#'
#'@param data The matrix of gene expression data, rows are samples and columns are genes.
#'
#'@return a list
#'
#'@examples
#'if(!require("GSVA")) BiocManager::install("GSVA",update = F,ask = F)
#'if(!require("clusterRepro")) install.packages("clusterRepro",update = F,ask = F)
#'library(TIME)
#'data('data')
#'result=time_phenotype(data2)
time_phenotype <- function(data=data){
  ssgsea <- t(scale(t((gsva(as.matrix(data2),immune,method = 'ssgsea')))))
  x <- intersect(rownames(ssgsea),rownames(Centroids))
  ssgsea <- ssgsea[x,]
  Centroids <- Centroids[x,]
  Result <- IGP.clusterRepro(ssgsea,Centroids)
  Result2 <- clusterRepro(Centroids = as.matrix(Centroids),New.data = ssgsea,Number.of.permutations = 1000)
  TIMEscore <- gsva(as.matrix(data2),TIME,method = 'ssgsea')
  TIMEscore2 <- data.frame(ID=colnames(TIMEscore),TIMEscore=as.numeric(t(TIMEscore)))
  rr <- list()
  rr[['TIME']] <- data.frame(ID=names(Result$Class),TIME=paste0('TIME-',Result$Class),row.names = NULL)
  rr[['p.value']] <- data.frame(TIME=c('TIME-1','TIME-2','TIME-3'),p.value=Result2$p.value)
  rr[['Number']] <- data.frame(TIME=c('TIME-1','TIME-2','TIME-3'),Number=Result2$Number)
  rr[['Actual.IGP']] <- data.frame(TIME=c('TIME-1','TIME-2','TIME-3'),Actual.IGP=Result2$Actual.IGP)
  rr[['Actual.Size']] <- data.frame(TIME=c('TIME-1','TIME-2','TIME-3'),Actual.Size=Result2$Actual.Size)
  rr[['ssgsea_result']] <- as.data.frame(ssgsea)
  rr[['TIMEscore']] <- TIMEscore2
  return(rr)
}