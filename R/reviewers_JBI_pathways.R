# reviewers_JBI_pathways.R
# analysis of pathway information based on membership proteinS

library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(xtable)
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
#library(biomaRt)

setwd("C:/R-files/sider")    # point to where my code and data files are
load("reviewersJBIloadData.RData") # load in the data structures made by reviewers_JBI_LoadData.R
load("reviewers_candidates.RData") # load in the data structures (drug_list,se_list,repos_list,mycandidates)

source("reviewers_JBI_functions.R")  # read in the functions required for finding lists of drugs and side-effects

ontargets <- read.csv(file='C://R-files//sider//drugbank-proteins.tsv', header=TRUE, sep="\t")

#commented out as these have already been added to my database.
#drug_atc <- repair_atc(drugbank_id="DB01175",name="Escitalopram",atc_codes="N06AB10",type="Small Molecule") # add Escitalopram to list
#drug_atc <- repair_atc(drugbank_id="DB07701",name="DB07701",atc_codes="unknown",type="Small Molecule") # add DB07701 to list

joint_list <- names_ids(drug_list) # we need DB ids as well as DB names
all_targets <- get_all_drug_targets(joint_list)  
targets_ig <- graph.data.frame(all_targets[,2:3],directed=FALSE)
gstats <- get_gstatistics(targets_ig)
plot(targets_ig)

# Remove entries with no proper drugname i.e. anything like DBXXXX etc
# Also, some drugs do not exist in all of the databases I use, hence causes crashes, so
# we lose 10 drugs.
mycandidates <- mycandidates[- grep("DB0", mycandidates$drugnames),]
mycandidates <- mycandidates[- grep("99mT", mycandidates$drugnames),]
mycandidates <- mycandidates[- grep("SCH", mycandidates$drugnames),]
mycandidates <- mycandidates[- grep("Dextrofloxacine", mycandidates$drugnames),]
mycandidates <- mycandidates[- grep("Iopamidol", mycandidates$drugnames),]
mycandidates <- mycandidates[- grep("Medroxyprogesterone Acetate", mycandidates$drugnames),]
mycandidates <- mycandidates[- grep("Gadobutrol", mycandidates$drugnames),]

candidate_list <- names_ids(mycandidates[,1])
candidate_targets <- get_all_drug_targets(candidate_list) # I need to add drugname and not just DB id
                                                      # matching candidate_target $Drug with candidate_list $name
                                                    # pretty plots of ontarget proteins to drug networks
# http://kateto.net/networks-r-igraph
# --- using Guangchuang Yu reactomepa package for enrichment ------------
# http://www.bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#bitr-biological-id-translator
# http://bioconductor.org/packages/release/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html
data(geneList)
de <- names(geneList)[abs(geneList) > 1.5]
x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
head(as.data.frame(x))
enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1) # interesting plot
barplot(x, showCategory=10)












