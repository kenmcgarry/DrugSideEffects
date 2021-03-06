# reviewers_JBI_run.R

library(scales)
library(DT)
library(caret)
library(kernlab)
library(dplyr)
#library(tidyverse)
library(stringr)
library(VennDiagram)
library(ggplot2)
library(xtable)
library(clusterProfiler)
library(igraph)
library(ROCR)

setwd("C:/R-files/sider")    # point to where my code lives
load("reviewersJBIloadData.RData") # load in the data structures made by reviewers_JBI_LoadData.R
load("reviewers_candidates.RData")
source("reviewers_JBI_functions.R")  # load in the functions required for finding lists of drugs and side-effects
load("12thJune2017.RData")

# Create a structure to hold all known drugs treating your particular disease of interest
# Obviously instantiate these functions before calling them.

drug_list <- get_drugs("C0002395",restrictedlist)  # umls code for Alzheimers
drug_list <- get_drugs("C0020179",restrictedlist)  # umls code for Huntingdons
drug_list <- get_drugs("C0149925",restrictedlist)  # umls code for small cell lung cancer
drug_list <- get_drugs("C0007131",restrictedlist)  # umls code for non small cell lung cancer
drug_list <- get_drugs("C0024141",restrictedlist)  # umls code for systemic lupus erythematosus
drug_list <- get_drugs("C0029408",restrictedlist)  # umls code for Osteoarthritis
drug_list <- get_drugs("C0003873",restrictedlist)  # umls code for Rheumatoid arthritis

se_list <- get_sideeffects(drug_list) # what side-effects do these drugs have?

repurposing <- get_repos_sideeffects(drug_list,se_list) # collate common side-effects for repurposing
common_list <- repurposing[[1]]
drugs_used <- repurposing[[2]]


mycandidates <- search_similardrugs(common_list,drug_list) # what drugs could be candidates for repositioning for our disease

#drugs_used <- sample(drugs_used)

dindex <- match(drugs_used, drug_list)  # just need side-effects of the drugs we actually will use
plot_venn(drugs_used,se_list,dindex) # No more than 5 drugs max, or Venn will not work!!!!

#OD_RA <- mycandidates
#OD_SLE <- mycandidates
#OD_NSCLC <- mycandidates



