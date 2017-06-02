# reviewers_JBI_integrate.R
# Integrate all data and knowledge for overall impact in principled way

library(igraph)
library(dplyr)
library(Matrix)
library(tidyverse)
library(stringr)
library(ggplot2)
library(xtable)
library(ReactomePA)
library(clusterProfiler)
library(GOSemSim)
library(DOSE)
library(AnnotationHub)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
#library(biomaRt)

setwd("C:/R-files/sider")    # point to where my code lives
load("reviewersJBIloadData.RData") # load in the data structures made by reviewers_JBI_LoadData.R
load("reviewers_candidates.RData")
load("1stJune2017.RData")
source("reviewers_JBI_functions.R")  # 

# Jaccard similarity coefficient is used to generate an association score to calculate the
# overlap between two sets of features i.e. :
#   1. off-targets & pathway      2. drug chemical similarity & side-effects
#   3. on-targets & pathway       4. disease ontology pairs 
#   5. side-effects & drug ATC
#
# The Jaccard association index then draws together the combined individual scores for each pair and
# providing a value between 0 and 1 which is used to rank each candidate drug.
  
# Calculate on-target promiscuity score for all drugs
on_count <- length(unique(candidate_targets[,2])) # how many ontarget proteins in total?
on_score <- vector()
for (i in 1:nrow(candidate_list)){
  on_temp <- filter(candidate_targets, name == candidate_list[i,2])
  on_score[i] <- length(unique(on_temp$SYMBOL))
# how many ontargets attached to each drug?
}

on_score <- on_score/on_count   # ratio of on-targets to all on-targets for each candidate drug
rm(on_temp,on_count) # tidy up mess

# Calculate side-effect similarity score for all drugs
se_score <- mycandidates[,3]/100  # divide by 100 to remove percentage and get a value between 0 and 1


# Calculate score for number of diseases/symptoms each candidate drug treats
# inds_score
inds_score <- vector()
for (i in 1:length(inds)){
  inds_temp <- unlist(inds[i])
  inds_score[i] <- length(unique(inds_temp))
}
inds_count <- length(unique(unlist(inds)))
inds_score <- inds_score/inds_count


collate_scores <- cbind(inds_score,se_score,on_score)
collate_scores <- jaccard(collate_scores)  # ensure Matrix library is loaded.
collate_scores <- diag(collate_scores)
candidate_list_scores <- cbind(candidate_list,collate_scores)
candidate_list_scores <- arrange(candidate_list_scores,desc(collate_scores))
head(candidate_list_scores)
















