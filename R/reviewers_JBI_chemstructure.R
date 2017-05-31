# reviewers_JBI_chemstructure.R

#library(scales)
#library(DT)
#library(caret)
#library(kernlab)
library(dplyr)
library(tidyverse)
#library(stringr)
#library(VennDiagram)
library(ggplot2)
library(xtable)
library(clusterProfiler)
library(igraph)
library(paxtoolsr)
library(rJava)
library(ReactomePA)
library(ChemmineR)
library(ChemmineOB)

setwd("C:/R-files/sider")    # point to where my code lives
load("reviewersJBIloadData.RData") # load in the data structures made by reviewers_JBI_LoadData.R
load("reviewers_candidates.RData")
load("27thMay2017.RData")
source("reviewers_JBI_functions.R")  # load in the functions required for finding lists of drugs and side-effects

# For each candidate drug, get the diseases it treats
inds <-rep("",nrow(candidate_list))
for (i in 1:nrow(candidate_list)){
  disease <- filter(indications, drugbank_name == candidate_list$name[i])
  inds[i] <- disease %>%
    select(meddra_name)
  cat("\ncandidate drug:",i,"is", candidate_list$name[i],"and treats",inds[[i]],"\n")
}
# The [inds] data structure (list) was used manually to create the LaTex table in the article)

# How many indications/diseases on average do our drugs treat, mim and max numbers
x<-base::lapply(inds,length)
x<-as.numeric(unlist(x))
mean(x)
min(x)
max(x)

# Ok, now for chemical structure similarity analysis.
# 

setwd("C:/R-files/bigfiles")
sdfset <- read.SDFset("structures.sdf")
valid <- validSDF(sdfset)
sdfset <- sdfset[valid]

## Assign drugbank IDs from datablock to cid slot
blockmatrix <- datablock2ma(datablocklist=datablock(sdfset))
cid(sdfset) <- as.character(blockmatrix[,"DRUGBANK_ID"])

## Generate APset and FPset (note FPset: has better search performance)
apset <- sdf2ap(sdfset)
fpset <- desc2fp(apset, descnames=2048, type="FPset")

## Subsetting by cid slot using drugbank ids should work now consistently
sdfset["DB00472"]
apset["DB00472"]

# Overwrite drug names with our candidate drugs using drugbank_id
drugs <- candidate_list$drugbank_id
fpdrugs <- fpset[drugs]   # extract our 69 chemical signatures from the many. THINK ABOUT CURRENT DRUGS
params <- genParameters(fpdrugs)  # 

results1 <- fpSim(fpdrugs[[1]], fpdrugs, top=25, parameters=params,method="Tversky") 
results2 <- fpSim(fpdrugs[[2]], fpdrugs, top=25, parameters=params,method="Tanimoto") 
results3 <- fpSim(fpdrugs[[3]], fpdrugs, top=25, parameters=params,method="Tanimoto") 

results1 <- cbind(candidate_list$name[1:25],results1) # USE drug names and NOT drugbank ID's
results2 <- cbind(drugs[rownames(results2)],results2)
results3 <- cbind(drugs[rownames(results3)],results3)

library(xtable)
print.xtable(xtable(results1))
print.xtable(xtable(results2))
print.xtable(xtable(results3))

#print.xtable(results,label='comp',caption = 'Similarity measures drugs based on 2048-length fingerprints')
# Creates the similarity score matrix and clusters them.
simMA <- sapply(cid(fpdrugs), function(x) fpSim(x=fpdrugs[x], fpdrugs, sorted=FALSE)) 
colnames(simMA)<-candidate_list$name
rownames(simMA)<-candidate_list$name

library(ape)
library(sparcl)
library(cluster) # used for kmeans and silhoutte plot

cl <- kmeans(simMA,10,nstart=50)
sk <- silhouette(cl$cl,dist(simMA))
plot(sk)

par(mar=c(3, 3, 3, 3))
hc <- hclust(as.dist(1-simMA), method="complete")
plot(as.phylo(hc), cex = 0.9, label.offset = 0.01)

y <- cutree(hc,6)
ColorDendrogram(hc,y=y,labels=candidate_list$name,branchlength = 0.7,cex = 2)
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=3), horiz=TRUE) #list(col=4, lwd=3)
rect.hclust(hc,k=5,border="red")

























