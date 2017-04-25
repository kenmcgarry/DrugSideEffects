# reviewers_JBI_pathways.R
# analysis of pathway information based on membership proteinS

library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(xtable)
library(ReactomePA)
library(clusterProfiler)
library(bitr)
library(org.Hs.eg.db)
library(biomaRt)

setwd("C:/R-files/sider")    # point to where my code and data files are
load("reviewersJBIloadData.RData") # load in the data structures made by reviewers_JBI_LoadData.R
load("reviewers_candidates.RData") # load in the data structures (drug_list,se_list,repos_list,mycandidates)

source("reviewers_JBI_DrugList.R")  # read in the functions required for finding lists of drugs and side-effects

ontargets <- read.csv(file='C://R-files//sider//drugbank-proteins.tsv', header=TRUE, sep="\t")

#drug_atc <- repair_atc(drugbank_id="DB01175",name="Escitalopram",atc_codes="N06AB10",type="Small Molecule") # add Escitalopram to list
#drug_atc <- repair_atc(drugbank_id="DB07701",name="DB07701",atc_codes="unknown",type="Small Molecule") # add DB07701 to list

joint_list <- names_ids(drug_list) # we need DB ids as well as DB names
all_targets <- get_all_drug_targets(joint_list)  

get_drug_targets <- function(thedrug){
  temp_targets <- filter(ontargets, drugbank_id == thedrug) 
  x <- as.character((temp_targets$uniprot_id))
  y <- select(org.Hs.eg.db, x, "SYMBOL", "UNIPROT")
  drug_targets <- cbind(y,rep(thedrug,nrow(y)))
  names(drug_targets)[3]<-"Drug"
  return(drug_targets)
}

get_all_drug_targets <- function(joint_list){
  ndrugs <- nrow(joint_list)
  drug_targets <- vector()
  for (i in 1:ndrugs){
    details <- get_drug_targets(joint_list$drugbank_id[i])
    drug_targets <- rbind(drug_targets,details)
  }
return(drug_targets)  
}

# names_ids() creates a dataframe that adds drugbank_id to rest of data. 
# collates info for the conventional drugs
names_ids <- function(drug_list){
  temp_atc <- data.frame(drugbank_id=character(0),name=character(0),type=character(0),atc_codes=character(0),stringsAsFactors = FALSE); 
  ndrugs<-length(drug_list)
  for ( i in 1:ndrugs){
    temp_atc[i,] <- filter(drug_atc, name == drug_list[i])
  }
  return(temp_atc)
}


# repair_atc() manually appends an entry to end of the database when you realise a drug is missing,
# stops errors occurring.
repair_atc <- function(drugbank_id,name,atc_codes,type)  {
  drug_atc <- rbind(drug_atc,c(drugbank_id,name,type,atc_codes))
  tail(drug_atc)
  return(drug_atc)
}
  

