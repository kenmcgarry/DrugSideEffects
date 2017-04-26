# reviewers_JBI_pathways.R
# analysis of pathway information based on proteins associated with each drug.
# 1. First get the proteins targetted by the drugs then...
# 2. Get the proteins the drugs inadvertently target (responsible for many side-effects)
# 3. Perform pathway analysis on the two groups

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
library(igraph)

setwd("C:/R-files/sider")    # point to where my code and data files are
load("reviewersJBIloadData.RData") # load in the data structures made by reviewers_JBI_LoadData.R
load("reviewers_candidates.RData") # load in the data structures (drug_list,se_list,repos_list,mycandidates)
ontargets <- read.csv(file='C://R-files//sider//drugbank-proteins.tsv', header=TRUE, sep="\t")  # proteins targetted by our drugs
source("reviewers_JBI_DrugList.R")  # read in the functions required for finding lists of drugs and side-effects

#drug_atc <- repair_atc(drugbank_id="DB01175",name="Escitalopram",atc_codes="N06AB10",type="Small Molecule") # add Escitalopram to list
#drug_atc <- repair_atc(drugbank_id="DB07701",name="DB07701",atc_codes="unknown",type="Small Molecule") # add DB07701 to list

joint_list <- names_ids(drug_list) # we need DB ids as well as DB names
all_targets <- get_all_drug_targets(joint_list)  

# pathway analysis using ReactomePA
entrez_id <- bitr(all_targets$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") #convert to  entrez ids so 
path1 <- enrichPathway(gene=entrez_id$ENTREZID,pvalueCutoff=0.05, readable=TRUE)      # enrichpathway() can use them.
head(as.data.frame(path1))
barplot(path1, showCategory=20)
dotplot(path1, showCategory=20)
enrichMap(path1, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
enrichMap(path1, layout=igraph::layout.fruchterman.reingold,vertex.label.cex = 1)

cnetplot(path1, categorySize="pvalue")

### --------- FUNCTIONS ----------------
# get_drug_targets() receives one drug name and one drugbankid and gets all the proteins it targets
get_drug_targets <- function(drugid,drugname){
  temp_targets <- filter(ontargets, drugbank_id == drugid) 
  x <- as.character((temp_targets$uniprot_id))
  y <- select(org.Hs.eg.db, x, "SYMBOL", "UNIPROT")  # translate uniprots into symbols
  z <- rep(drugname,nrow(y))
  v <- rep(1,nrow(y))  # all 1's because they are all protein targets named by drugbank
  drug_targets <- cbind(y,rep(drugid,nrow(y)),z,v)
  names(drug_targets)[3]<-"drugbank_id"
  names(drug_targets)[4]<-"name"
  names(drug_targets)[5]<-"target"
  return(drug_targets)
}

# get_all_drug_targets() receives drugbank_id to search, finally joins data structure 
get_all_drug_targets <- function(joint_list){
  ndrugs <- nrow(joint_list)
  drug_targets <- vector()
  for (i in 1:ndrugs){
    details <- get_drug_targets(joint_list$drugbank_id[i],joint_list$name[i])
    drug_targets <- rbind(drug_targets,details)
    drug_targets
  }
  #drug_targets <- rbind(drug_targets,target)
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


# repair_atc() manually appends an entry to end of the drug_atc atabase when you realise a drug is missing,
# this will stop errors occurring. Will need four bits of information e.g.
# drug_atc <- repair_atc(drugbank_id="DB07701",name="DB07701",atc_codes="unknown",type="Small Molecule") # add DB07701 to list
repair_atc <- function(drugbank_id,name,atc_codes,type)  {
  drug_atc <- rbind(drug_atc,c(drugbank_id,name,type,atc_codes))
  tail(drug_atc)
  return(drug_atc)
}
  

