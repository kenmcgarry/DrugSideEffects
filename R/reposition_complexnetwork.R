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
load("11-09-2017.RData")
source("reviewers_JBI_functions.R")  # load in the functions required for finding lists of drugs and side-effects
ontargets <- read.csv(file='C://R-files//sider//drugbank-proteins.tsv', header=TRUE, sep="\t")

joint_list <- names_ids(drug_list) # we need DB ids as well as DB names
all_targets <- get_all_drug_targets(joint_list)  
ct <- candidate_targets[,c(2,5)]
ct <- ct %>% drop_na()

# these two lines needed further down
ctt <- rbind(ct,all_targets[,c(2,5)])
write.csv(ctt,"drugtoprotein.csv")

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
nfolds <- 2
nmiss <- 31125/nfolds
Avec <- A[lower.tri(A)]
Avec.pred1 <- numeric(length(Avec))

perf1 <- list() # create empty vector of perfs for ROC-Curves
auc1 <- list() # create empty vector of area under curves
pred1 <- list() # create empty vector of predictions for PR-Curves
perf2 <- list() # create empty vector of perfs for ROC-Curves
auc2 <- list() # create empty vector of area under curves
pred2 <- list() #
perf3 <- list() # create empty vector of perfs for ROC-Curves
auc3 <- list() # create empty vector of area under curves
pred3 <- list() #
  
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
  #pred1 <- floor(Avec.pred1)
  Avec[Avec==2] <- 0
  if(j==1){
    pred1 <- ROCR::prediction(Avec.pred1, Avec)
    perf1 <- ROCR::performance(pred1, "tpr", "fpr")
    auc1 <- ROCR::performance(pred1, "auc")}
  if(j==2){
    pred2 <- ROCR::prediction(Avec.pred1, Avec)
    perf2 <- ROCR::performance(pred2, "tpr", "fpr")
    auc2 <- ROCR::performance(pred2, "auc")}
  if(j==3){
    pred3 <- ROCR::prediction(Avec.pred1, Avec)
    perf3 <- ROCR::performance(pred3, "tpr", "fpr")
    auc3 <-  ROCR::performance(pred3, "auc")}
}

# ROC PLOT COMPARING THREE MODELS
#plot(perf1, col="red", lwd=5)
#plot(perf2, add = TRUE, col="blue",lwd=5)
#plot(perf3, add = TRUE, col="green",lwd=5)
 

# PR-CURVE AND ROC PLOT FOR THREE MODELS 
# http://takayasaito.github.io/precrec/articles/introduction.html
library(precrec)
library(ggplot2)

#--------------------- test the package/create four artificial models ------------
#samps2 <- create_sim_samples(1, 100, 100, "all")
# Use a sample dataset created by the create_sim_samples function
#msmdat3 <- mmdata(samps2[["scores"]], samps2[["labels"]], modnames = samps2[["modnames"]])
# Calculate ROC and Precision-Recall curves for multiple models
#mscurves <- evalmod(msmdat3)
# Show ROC and Precision-Recall curves with the ggplot2 package
#autoplot(mscurves)
#----------------------------------------------------------------------

# convert my stuff to precrec format 
score1 <- slot(pred1,"predictions") # must extract slot data from S3 structure
score2 <- slot(pred2,"predictions")
score3 <- slot(pred3,"predictions")
msmdat <- mmdata(scores = c(score1,score2,score3), 
                 labels=c(pred1@labels,pred2@labels,pred3@labels),
                 modnames = c("Full model","chem+ppi","pathway+ontology"))
# Calculate ROC and Precision-Recall curves for multiple models
mscurves <- evalmod(msmdat)
autoplot(mscurves)



#----------------------------------------------------------------------------

# CHUNK 33
auc[[1]].auc <- ROCR::performance(pred1, "auc")
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
  mygraph <-graph.data.frame(mydata,directed=FALSE)
  in1 <- get_gstatistics(mygraph)  #_long
  filen <- temp[i]
  allstats <- rbind(allstats,in1[nrow(in1),],make.row.names = FALSE)
  cat(temp[i],"\n")
} 

rownames(allstats) <- c(sub('\\.txt$', '', temp) ) #get rid of file extension
rm(list = ls(pattern = glob2rx("*.txt"))) # get rid of useless file data from memory
allstats <- allstats[order(allstats$between),]  #sort allstats according to 
table1 <- xtable(allstats,digits=c(3,3,1,1,3,1,3,3,3,3,3,3)) # latex table

# more drug centric oriented analysis
# create the disease to drug to targetprotein network for visual inspection 

setwd("C:/R-files/sider") 
drug25 <- read.csv(file='C://R-files//sider//drugtoprotein25.csv', header=TRUE, sep=",")

graph25 <-graph.data.frame(drug25,directed=FALSE)
graphnames <-V(graph25)$name
ad <- get.adjacency(graph25)

nodecolor=character(ncol(ad))  # create a character for every column in adjaceny matrix
nodelabel<-V(graph25)$name
nodesize <- vector(length=152)
nodeshape <- vector(mode="character",length=152)
nodeshape[1:120]<-"circle"
nodeshape[121:152] <- "triangle"
nodesize[1:120] <- 5
nodesize[121:152] <- 15
nodecolor[121:152] <-"tomato" # figure where the proteins and where the drugs are
nodecolor[1:120] <-"lightgreen" 

graph25 <- igraph::simplify(graph25)#, remove.multiple = T, remove.loops = F, edge.attr.comb=c(weight="sum", type="ignore") )
dg <- igraph::degree(graph25, mode = c("all"),loops = FALSE, normalized = FALSE)

for (i in 1:120){
  if(dg[i] > 1)
    nodecolor[i] <-"gold"
    nodesize[i] <- 9
}

graph25<-as.undirected(graph25); 
plot(graph25, edge.color="darkgray", 
     vertex.color=nodecolor,
     vertex.label=nodelabel,
     vertex.shape=nodeshape,
     vertex.size=nodesize,
     vertex.label.cex=0.6, 
     vertex.label.font=0.5, 
     vertex.frame.color="white",
     #vertex.frame.color="darkgreen",
     vertex.label.color="black", 
     vertex.label.family = "sans",
     layout=layout.kamada.kawai(graph25))

tkplot(graph25,layout = layout.fruchterman.reingold,vertex.label = nodelabel,
       vertex.label.color= "black",
       vertex.size=nodesize, 
       vertex.color=nodecolor,
       vertex.shape=nodeshape,
       vertex.size=nodesize,
       edge.arrow.size=0, edge.curved=FALSE)


tkid <- tkplot(graph25) #tkid is the id of the tkplot that will open, move nodes 
                        # around to your hearts content then...use line below
l <- tkplot.getcoords(tkid) # to grab the coordinates from tkplot


plot(graph25, edge.color="darkgray", 
     vertex.color=nodecolor,
     vertex.label=nodelabel,
     vertex.shape=nodeshape,
     vertex.size=nodesize,
     vertex.label.cex=0.6, 
     vertex.label.font=0.5, 
     vertex.frame.color="white",
     #vertex.frame.color="darkgreen",
     vertex.label.color="black", 
     vertex.label.family = "sans",
     layout=l)


