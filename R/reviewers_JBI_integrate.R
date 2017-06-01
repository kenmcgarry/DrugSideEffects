# reviewers_JBI_integrate.R
# Integrate all data and knowledge for overall impact in principled way

library(igraph)
library(dplyr)
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
load("31stMay2017.RData")
source("reviewers_JBI_functions.R")  # 

# Jaccard similarity coefficient is used to generate an association score to calculate the
# overlap between two sets of features i.e. :
#   1. off-targets & pathway      2. drug chemical similarity & side-effects
#   3. on-targets & pathway       4. disease ontology pairs 
#   5. side-effects & drug ATC
#
# The Jaccard association index then draws together the combind individual scores for each pair and
# providing a value between 0 and 1 which is used to rank each candidate drug.
  
# Calculate on-target promiscuity score for all drugs
on_count <- length(unique(candidate_targets[,2])) # how many ontarget proteins in total?
on_score <- vector()
for (i in 1:nrow(candidate_list)){
  on_temp <- filter(candidate_targets, name == candidate_list[i,2])
  on_score[i] <- length(unique(on_temp$SYMBOL))
# how many ontargets attached to each drug?
}
rm(on_temp) # tidy up mess
on_score/on_count   # ratio of on-targets to all on-targets


# Calculate side-effect similarity score for all drugs
se_score <- mycandidates[,3]/100



























