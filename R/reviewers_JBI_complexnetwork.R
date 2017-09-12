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
library(tidyverse)  # you may have to install package 'foreign' manually
library(stringr)
library(ggplot2)
library(xtable)
library(org.Hs.eg.db)

setwd("C:/R-files/sider")    # point to where my code lives
#load("reviewersJBIloadData.RData") # load in the data structures made by reviewers_JBI_LoadData.R
#load("reviewers_candidates.RData")
#load("12thJune2017.RData")

load("complexnets.RData")
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

# using A (edge list only)
A <- get.adjacency(targets_ig, sparse=FALSE)
repos.fit1 <- eigenmodel_mcmc(A, R=2, S=11000,burn=10000)

# using disease  ontology (dose)
dose.op <- attrib_candidates$dose %o% attrib_candidates$dose
dose <- matrix(as.numeric(dose.op  %in% c(1, 4, 9)), 250, 250)
dose <- array(dose,dim=c(250, 250, 1))
repos.fit2 <- eigenmodel_mcmc(A, dose, R=2,S=11000,burn=10000)
# -------------------------------------

#using biological pathways (pathway)    
pathway.op <- attrib_candidates$pathway %o% attrib_candidates$pathway
pathway <- matrix(as.numeric(pathway.op %in% c(1, 4, 9)), 250, 250)
pathway <- array(pathway,dim=c(250, 250, 1)) 
repos.fit3 <- eigenmodel_mcmc(A, pathway,R=2, S=11000, burn=10000)

lat.sp.1 <- eigen(repos.fit1$ULU_postmean)$vec[, 1:2]
lat.sp.2 <- eigen(repos.fit2$ULU_postmean)$vec[, 1:2]
lat.sp.3 <- eigen(repos.fit3$ULU_postmean)$vec[, 1:2]

apply(repos.fit1$SL_postsamp,2,mean)
apply(repos.fit2$SL_postsamp,2,mean)
apply(repos.fit3$SL_postsamp,2,mean)


# RECREATE SAME MODELS BUT WITH LOOP TO IMPLEMENT n-fold CROSS-VALIDATION 
# CHUNK 1
A <- get.adjacency(targets_ig,sparse=FALSE)
dim(A)

v.attrs <- get.data.frame(targets_ig, what="vertices")
# CHUNK 30
perm.index <- sample(1:31125) # 250 * 249/2= 31,125 permutations
nfolds <- 5
nmiss <- 31125/nfolds
Avec <- A[lower.tri(A)]
Avec.pred1 <- numeric(length(Avec))

perf <- list() # create empty vector of perfs
auc <- list() # create empty vector of area under curves
  
for (j in 1:3){
  if(j==1) modelversion <- NULL   # select all variables A + B
  if(j==2) modelversion <- dose   # select only A
  if(j==3) modelversion <- pathway  # select only B
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
  perf[j] <- ROCR::performance(pred1, "tpr", "fpr")
  auc[j] <- ROCR::performance(pred1, "auc")
}

# ROC PLOT COMPARING THREE MODELS
plot(perf[[1]], col="red", lwd=5)
plot(perf[[2]], add = TRUE, col="blue",lwd=5)
plot(perf[[3]], add = TRUE, col="green",lwd=5)
abline(a=0, b= 1,lty=2)

# CHUNK 33
auc[[1]].auc <- performance(pred1, "auc")
slot(auc[[1]], "y.values")


#---- use verification package for generting p-values -----
library(verification)
# Data used from Mason and Graham (2002).
a <- c(1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995)
b <- c(0,0,0,1,1,1,0,1,1,0,0,0,0,1,1)
c <- c(.8, .8, 0, 1,1,.6, .4, .8, 0, 0, .2, 0, 0, 1,1)
d <- c(.928,.576, .008, .944, .832, .816, .136, .584, .032, .016, .28, .024, 0, .984, .952)
pdata <- data.frame(a,b,c, d)
names(pdata)<- c("year", "event", "p1", "p2")
## for model with ties
roc.area(pdata$event, pdata$p1)
## for model without ties
roc.area(pdata$event, pdata$p2)

## for model with ties
roc.plot(pdata$event, pdata$p1)
## for model without ties
roc.plot(pdata$event, pdata$p2)



# --------------------------- compute statistics for each drug PPI network --------------------------------
# do the combined protein networks

setwd("C:/R-files/reposition/proteins") # need to point where STITCH datafiles are.

temp <-  list.files(pattern="*.txt") # list of all .txt files now in temp
allstats <- data.frame()

for (i in 1:length(temp)){
  mydata <- read.delim(temp[i], header=TRUE,sep='\t')
  mydata <- mydata[,1:2]
  mydata <-graph.data.frame(mydata,directed=FALSE)
  in1 <- get_gstatistics_long(mydata)
  allstats <- rbind(allstats,in1[nrow(in1),])
} 

rm(list = ls(pattern = glob2rx("*.txt"))) # get rid of useless file data from memory
table1 <- xtable(allstats,digits=c(3,3,0,0,3,0,3,3,3,3,3)) # latex table





