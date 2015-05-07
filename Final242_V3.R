rm(list=ls())

setwd("~/Documents/courses/stat/stat 242/Final Project/data")
library(elasticnet)

load("small_CTD1.Rdat")
load("small_scaled_micro_array_2.Rdat")
load("gene_to_geneset_matrix.Rdat")
load("tissue_matrix.Rdat")

X <- small_scaled_micro_array_2[-98,]
Y <- small_CTD1[-98,]

ngene <- dim(X)[2]
ndrug <- dim(Y)[2]
nsample <- dim(X)[1]

ntype <- dim(tissue_matrix)[2]
cancertype <- rep(0, 98)
cancertype <- as.numeric(lapply(rownames(X), 
                     function(x) which(tissue_matrix[x,]==1)))

cancernum <- rep(0, ntype)
for (i in 1:98){
  cancernum[cancertype[i]]=cancernum[cancertype[i]]+1
}

bX <- matrix(0, 25, ngene)
bY <- matrix(0, 25, ndrug)
wX <- matrix(0, nsample, ngene)
wY <- matrix(0, nsample, ndrug)

for (i in 1:25){
  if (sum(cancertype==i)>=2){
    bX[i,] <- colMeans(X[cancertype==i,])
    bY[i,] <- colMeans(Y[cancertype==i,])
  }else if(sum(cancertype==i)==1){
    bX[i,] <- X[cancertype==i,]
    bY[i,] <- Y[cancertype==i,]
  }
}

for (i in 1:98){
  wX[i,] <- X[i,]-bX[cancertype[i],]
  wY[i,] <- Y[i,]-bY[cancertype[i],]
}

bzero <- which(rowSums(bY)==0)
bX2 <- bX[-1*bzero,]
bY2 <- bY[-1*bzero,]

wzero <- which(rowSums(wY)==0)
wX2 <- wX[-1*wzero,]
wY2 <- wY[-1*wzero,]

save(bX2, file = "bX2.Rdat")
save(bY2, file = "bY2.Rdat")
save(wX2, file = "wX2.Rdat")
save(wY2, file = "wY2.Rdat")

alpha <- 1.5

bX2.norm <- scale(bX2)
bY2.norm <- scale(bY2)*alpha
bdata <- cbind(bX2.norm, bY2.norm)

wX2.norm <- scale(wX2)
wY2.norm <- scale(wY2)*alpha
wdata <- cbind(wX2.norm, wY2.norm)

k <- 4
para <- rep(1800,k)
sp <- arrayspc(wdata, K=k, para,
               use.corr=FALSE, max.iter=100,trace=FALSE,eps=1e-3)
output <- order(abs(sp$loadings[,1]),decreasing = TRUE)[1:100]

R <- 100

wdrug <- matrix(0, nrow=R, ncol=dim(Y)[2])
wgeneset <- matrix(0, nrow=R, ncol=dim(gene_to_geneset_matrix)[1])

bdrug <- matrix(0, nrow=R, ncol=dim(Y)[2])
bgeneset <- matrix(0, nrow=R, ncol=dim(gene_to_geneset_matrix)[1])

ngene <- dim(X)[2]
nw <- dim(wX2)[1]
nb <- dim(bX2)[1]

## within cancer-type PCA
for (i in 1:R){
  ind <- sample(1:nw, nw, replace=TRUE, prob=NULL)
  sp <- arrayspc(wdata[ind,],K=k,para,
                 use.corr=FALSE, max.iter=100,trace=FALSE,eps=1e-3)
  output <- order(abs(sp$loadings[,1]),decreasing = TRUE)[1:100]
  wdrug[i,output[output>ngene]-ngene] <- 1
  genes <- output[output<=ngene]
  for (j in length(genes)){
    set <- gene_to_geneset_matrix[,genes[j]]==1
    wgeneset[i, set] <- 1
  }
}

save(wdrug, file="wdrug.Rdat")
save(wgeneset, file="wgeneset.Rdat")

## between cancer-type PCA
para2 <- rep(250,k)
for (i in 1:R){
  ind <- sample(1:nb, nb, replace=TRUE, prob=NULL)
  sp <- arrayspc(bdata[ind,],K=k,para2,
                 use.corr=FALSE, max.iter=100,trace=FALSE,eps=1e-3)
  output <- order(abs(sp$loadings[,1]),decreasing = TRUE)[1:100]
  bdrug[i,output[output>ngene]-ngene] <- 1
  genes <- output[output<=ngene]
  for (j in length(genes)){
    set <- gene_to_geneset_matrix[,genes[j]]==1
    bgeneset[i, set] <- 1
  }
}

save(bdrug, file="bdrug.Rdat")
save(bgeneset, file="bgeneset.Rdat")




