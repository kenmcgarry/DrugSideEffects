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
#drug_list <- c('Donepezil','Galantamine','Rivastigmine')  # this a manual example for Alzheimers

drug_list <- get_drugs("C0002395",restrictedlist)  # using the umls code for Alzheimers what drugs are used to treat it?

se_list <- get_sideeffects(drug_list) # what side-effects do these drugs have?

repos_list <- get_repos_sideeffects(drug_list,se_list) # get collate common side-effects for repurposing

mycandidates <- search_similardrugs(repos_list,drug_list) # what drugs could be candidates for repositioning for our disease

plot_venn(drug_list,se_list) # No more than 5 drugs max, or Venn will not work!!!!







