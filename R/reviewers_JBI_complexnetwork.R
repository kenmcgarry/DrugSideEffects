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
library(org.Hs.eg.db)

setwd("C:/R-files/sider")    # point to where my code lives
load("reviewersJBIloadData.RData") # load in the data structures made by reviewers_JBI_LoadData.R
load("reviewers_candidates.RData")
load("12thJune2017.RData")
source("reviewers_JBI_functions.R")  # load in the functions required for finding lists of drugs and side-effects

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
# rebuild the data structure to incorporate parameters for protein-to-drug network
edge.list.candidates <- data.frame(get.edgelist(targets_ig))
data_proteins <- (unique(edge.list.candidates$X1))
data_drugs <- (unique(edge.list.candidates$X2))
data_proteins <- as_tibble(data_proteins)
data_drugs <- as_tibble(data_drugs)
attrib_candidates <- bind_rows(data_drugs,data_proteins)
attrib_candidates[2] <-nrow(attrib_candidates)
attrib_candidates[1:nrow(data_drugs),2] <- 1  # 1 if its a drug
attrib_candidates[nrow(data_drugs)+1:nrow(data_proteins),2] <- 0 # zero if its a protein

attrib_candidates[1:nrow(data_drugs),3] <- substr(candidate_list$atc_codes,start = 1, stop = 3) # 2nd level ATC codes
attrib_candidates[nrow(data_drugs)+1:nrow(data_proteins),3] <- 0   # proteins dont have ATC codes

attrib_candidates[,4] <- sample(1:20, nrow(attrib_candidates), replace=TRUE)
attrib_candidates[,5] <- sample(1:20, nrow(attrib_candidates), replace=TRUE)

names(attrib_candidates)[1] <- "name"
names(attrib_candidates)[2] <- "type"
names(attrib_candidates)[3] <- "atc"
names(attrib_candidates)[4] <- "pathway"
names(attrib_candidates)[5] <- "dose"

# edge list only
A <- get.adjacency(targets_ig, sparse=FALSE)
lazega.leig.fit1 <- eigenmodel_mcmc(A, R=2, S=11000,burn=10000)

# using A
same.prac.op <- attrib_candidates$dose %o% attrib_candidates$dose
same.prac <- matrix(as.numeric(same.prac.op  %in% c(1, 4, 9)), 250, 250)
same.prac <- array(same.prac,dim=c(250, 250, 1))
lazega.leig.fit2 <- eigenmodel_mcmc(A, same.prac, R=2,S=11000,burn=10000)
# -------------------------------------

#using B    
#attrib_candidates$atc <- as.factor(attrib_candidates$atc)
same.off.op <- attrib_candidates$pathway %o% attrib_candidates$pathway
same.off <- matrix(as.numeric(same.off.op %in% c(1, 4, 9)), 250, 250)
same.off <- array(same.off,dim=c(250, 250, 1)) 
lazega.leig.fit3 <- eigenmodel_mcmc(A, same.off,R=2, S=11000, burn=10000)

lat.sp.1 <- eigen(lazega.leig.fit1$ULU_postmean)$vec[, 1:2]
lat.sp.2 <- eigen(lazega.leig.fit2$ULU_postmean)$vec[, 1:2]
lat.sp.3 <- eigen(lazega.leig.fit3$ULU_postmean)$vec[, 1:2]

apply(lazega.leig.fit1$SL_postsamp,2,mean)
apply(lazega.leig.fit2$SL_postsamp,2,mean)
apply(lazega.leig.fit3$SL_postsamp,2,mean)


# CHUNK 1
A <- get.adjacency(targets_ig,sparse=FALSE)
#A <- get.adjacency(lazega, sparse=FALSE)
dim(A)

v.attrs <- get.data.frame(targets_ig, what="vertices")
#v.attrs <- get.data.frame(lazega, what="vertices")
# CHUNK 30
perm.index <- sample(1:31125) # 250 * 249/2= 31,125 permutations
#perm.index <- sample(1:630) #  permutations
nfolds <- 3
nmiss <- 31125/nfolds
#nmiss <- 630/nfolds

Avec <- A[lower.tri(A)]
Avec.pred1 <- numeric(length(Avec))

perf <- list()
auc <- list()
  
for (j in 1:3){
  if(j==1) modelversion <- NULL   # select all variables A + B
  if(j==2) modelversion <- same.off   # select only A
  if(j==3) modelversion <- same.prac  # select only B
# CHUNK 31 - cross validation.
  for(i in seq(1,nfolds)){
    # Index of missing values.
    miss.index <- seq(((i-1) * nmiss + 1),(i*nmiss), 1)
    A.miss.index <- perm.index[miss.index]
  
    # Fill a new Atemp appropriately with NA's.
    Avec.temp <- Avec
    Avec.temp[A.miss.index] <- rep("NA", length(A.miss.index))
    Avec.temp <- as.numeric(Avec.temp)
    Atemp <- matrix(0, 250, 250) # varying these two numbers leads to different accuracies
    Atemp[lower.tri(Atemp)] <- Avec.temp
    Atemp <- Atemp + t(Atemp)
  
    # Now fit model and predict.
    Y <- Atemp
    model1.fit <- eigenmodel_mcmc(Y, modelversion,R=2,S=11000, burn=10000) # change model here
    model1.pred <- model1.fit$Y_postmean
    model1.pred.vec <- model1.pred[lower.tri(model1.pred)]
    Avec.pred1[A.miss.index] <- model1.pred.vec[A.miss.index]
  }
  pred1 <- floor(Avec.pred1)
  Avec[Avec==2] <- 0
  pred1 <- prediction(Avec.pred1, Avec)
  perf[j] <- performance(pred1, "tpr", "fpr")
  auc[j] <- performance(pred1, "auc")
}


plot(perf[[1]], col="red", lwd=5)
plot(perf[[2]], add = TRUE, col="blue",lwd=5)
plot(perf[[3]], add = TRUE, col="green",lwd=5)
abline(a=0, b= 1,lty=2 )

# CHUNK 33
auc[[1]].auc <- performance(pred1, "auc")
slot(auc[[1]], "y.values")



