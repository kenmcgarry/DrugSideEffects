# reviewers_JBI_complexnetwork.R
# Conduct complex network analysis on the protein ontarget interactions with the candidate drugs, ensure
# constructed network is of good quality by constructing ROC curves based on interactions and associations.

library(igraph)
library(xtable)
library(sand)
library(ROCR)
library(eigenmodel)
library(igraph)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(xtable)

ontargets <- read.csv(file='C://R-files//sider//drugbank-proteins.tsv', header=TRUE, sep="\t")

joint_list <- names_ids(drug_list) # we need DB ids as well as DB names
all_targets <- get_all_drug_targets(joint_list)  
ct <- candidate_targets[,c(2,5)]
ct <- ct %>% drop_na()
targets_ig <- graph.data.frame(ct,directed=FALSE)
gstats <- get_gstatistics(targets_ig)
plot(targets_ig)

gr1 <- sample_gnm(250, 699, directed = FALSE, loops = FALSE) # 250 proteins with 699 connections between them
gr1 <- igraph::simplify(gr1)

nedges<-699
nverts<- 250
avepath <- average.path.length(gr1)
connected<- is.connected(gr1)
avedeg <- mean(degree(gr1))
diam<- diameter(gr1,weights=NA)
modular<- modularity(gr1, membership(cluster_walktrap(gr1)))
transit=transitivity(gr1)
randomg <- data.frame(modular,avepath,nedges,nverts,transit,avedeg,diam,connected)
xtable(randomg) # 


# ----- ROC curves from Kolaczyk book (page 107-108) ----
# ------ multiple models for common practice etc ---------- 
A <- get.adjacency(lazega, sparse=FALSE)
lazega.leig.fit1 <- eigenmodel_mcmc(A, R=2, S=11000,burn=10000)

same.prac.op <- v.attr.lazega$Practice %o% v.attr.lazega$Practice
same.prac <- matrix(as.numeric(same.prac.op  %in% c(1, 4, 9)), 36, 36)
same.prac <- array(same.prac,dim=c(36, 36, 1))
lazega.leig.fit2 <- eigenmodel_mcmc(A, same.prac, R=2,S=11000,burn=10000)
# -------------------------------------

same.off.op <- v.attr.lazega$Office %o% v.attr.lazega$Office
same.off <- matrix(as.numeric(same.off.op %in% c(1, 4, 9)), 36, 36)
same.off <- array(same.off,dim=c(36, 36, 1)) 
lazega.leig.fit3 <- eigenmodel_mcmc(A, same.off,R=2, S=11000, burn=10000)

lat.sp.1 <- eigen(lazega.leig.fit1$ULU_postmean)$vec[, 1:2]
lat.sp.2 <- eigen(lazega.leig.fit2$ULU_postmean)$vec[, 1:2]
lat.sp.3 <- eigen(lazega.leig.fit3$ULU_postmean)$vec[, 1:2]

# CHUNK 1
#A <- get.adjacency(targets_ig,sparse=FALSE)
A <- get.adjacency(lazega, sparse=FALSE)
dim(A)

#v.attrs <- get.data.frame(targets_ig, what="vertices")
v.attrs <- get.data.frame(lazega, what="vertices")
# CHUNK 30
#perm.index <- sample(1:31125) # 250 * 249/2= 31,125 permutations
perm.index <- sample(1:630) #  permutations
nfolds <- 3
#nmiss <- 31125/nfolds
nmiss <- 630/nfolds

Avec <- A[lower.tri(A)]
Avec.pred1 <- numeric(length(Avec))

perf <- list()
  
for (j in 1:3){
  if(j==1) modelversion <- NULL
  if(j==2) modelversion <- same.off
  if(j==3) modelversion <- same.prac
# CHUNK 31 - cross validation.
for(i in seq(1,nfolds)){
  # Index of missing values.
  miss.index <- seq(((i-1) * nmiss + 1),(i*nmiss), 1)
  A.miss.index <- perm.index[miss.index]
  
  # Fill a new Atemp appropriately with NA's.
  Avec.temp <- Avec
  Avec.temp[A.miss.index] <- rep("NA", length(A.miss.index))
  Avec.temp <- as.numeric(Avec.temp)
  Atemp <- matrix(0, 36, 36) # varying these two numbers leads to different accuracies
  Atemp[lower.tri(Atemp)] <- Avec.temp
  Atemp <- Atemp + t(Atemp)
  
  # Now fit model and predict.
  Y <- Atemp
  model1.fit <- eigenmodel_mcmc(Y, modelversion,R=2,S=11000, burn=10000) # change model here
  model1.pred <- model1.fit$Y_postmean
  model1.pred.vec <- model1.pred[lower.tri(model1.pred)]
  Avec.pred1[A.miss.index] <- model1.pred.vec[A.miss.index]
}
  pred1 <- prediction(Avec.pred1, Avec)
  perf[j] <- performance(pred1, "tpr", "fpr")
}
# CHUNK 32
#pred1 <- prediction(Avec.pred1, Avec)
#perf1 <- performance(pred1, "tpr", "fpr")

plot(perf[[1]], col="red", lwd=5)
plot(perf[[2]], add = TRUE, col="blue",lwd=5)
plot(perf[[3]], add = TRUE, col="green",lwd=5)
abline(a=0, b= 1,lty=2 )

# CHUNK 33
perf1.auc <- performance(pred1, "auc")
slot(perf1.auc, "y.values")



