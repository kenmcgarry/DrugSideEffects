# reviewers_JBI_run.R

library(scales)
library(DT)
library(caret)
library(kernlab)
library(dplyr)
library(tidyverse)
library(stringr)
library(VennDiagram)
library(ggplot2)
library(xtable)
library(clusterProfiler)
library(igraph)


setwd("C:/R-files/sider")    # point to where my code lives
load("reviewersJBIloadData.RData") # load in the data structures made by reviewers_JBI_LoadData.R
load("reviewers_candidates.RData")
source("reviewers_JBI_functions.R")  # load in the functions required for finding lists of drugs and side-effects


# Create a structure to hold all known drugs treating your particular disease of interest
# Obviously instantiate these functions before calling them.

drug_list <- get_drugs("C0002395",restrictedlist)  # umls code for Alzheimers but what drugs are used to treat it?
drug_list <- get_drugs("C0020179",restrictedlist)  # umls code for Huntingdon's but what drugs are used to treat it?
drug_list <- get_drugs("C0149925",restrictedlist)  # umls code for small cell lung cancer but what drugs are used to treat it?
drug_list <- get_drugs("C0007131",restrictedlist)  # umls code for non small cell lung cancer but what drugs are used to treat it?
drug_list <- get_drugs("C0024141",restrictedlist)  # umls code for systemic lupus erythematosus
#  C0262584   C0149925   C0278987 C0162296

se_list <- get_sideeffects(drug_list) # what side-effects do these drugs have?

repurposing <- get_repos_sideeffects(drug_list,se_list) # collate common side-effects for repurposing
common_list <- repurposing[[1]]
drugs_used <- repurposing[[2]]


mycandidates <- search_similardrugs(common_list,drug_list) # what drugs could be candidates for repositioning for our disease

drugs_used <- sample(drugs_used)

plot_venn(drugs_used,se_list) # No more than 5 drugs max, or Venn will not work!!!!








