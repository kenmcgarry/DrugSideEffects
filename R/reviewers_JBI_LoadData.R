# reviewers_JBI_LoadData.R
# work started: 03/04/2017
# Responses to reviewer comments on repostitioning paper submitted to Journal of Biomedical Informatics
# The code here will load in the necessary files and preprocess them.
#
# source("http://bioconductor.org/biocLite.R")
# biocLite("nameofpackage")

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

# SIDER4.1 database (side-effects)
# WARNING COLUMNS HAVE ALL CHANGED FOR VERSION 4.1, EVEN WORSE - DRUGNAMES ARE NOT USED - MUST LOOKUP STITCH DATABASE FOR CONVERSION
# USING FILES CREATED BY D. HIMMELSTEIN. https://github.com/dhimmel/SIDER4

# SIDER4.1 data downloaded from : http://sideeffects.embl.de/download/
indications  <- file.path('C://R-files//sider', 'indications.tsv') %>% read.delim(na.strings='',header = TRUE,stringsAsFactors=FALSE)
sideeffects <- file.path('C://R-files//sider', 'side-effects.tsv') %>% read.delim(na.strings='',header = TRUE,stringsAsFactors=FALSE)
sideterms <- file.path('C://R-files//sider', 'side-effect-terms.tsv') %>% read.delim(na.strings='',header = TRUE,stringsAsFactors=FALSE)
sidefreq <- file.path('C://R-files//sider', 'meddra_freq.tsv.gz') %>% read.delim(na.strings='',header = FALSE,stringsAsFactors=FALSE)

# It will be useful to compare ATC codes of conventional drugs with any candidate drug.
drug_atc <-file.path('C://R-files//drugbank//','drugbank1.tsv') %>% read.delim(na.strings='',sep='\t',header=TRUE,comment.char="#",stringsAsFactors = FALSE)
drug_atc <- drug_atc[,c(1:2,3,5)]
# we can obviously only use drugs with ATC codes!
drug_atc <- drug_atc %>%
  filter(!is.na(atc_codes))

# SIDER4.1 contains drugs/chemicals with long names such as "(1S,4S)-4-(3,4-dichlorophenyl)-N-methyl-1,2,3,4-tetrahydronaphthalen-1-amine."
# these are causing string searching problems (as well as looking very untidy!) so will replace them with drugbank ID using a regex.
z=0
for (i in 1:nrow(sideeffects)){
  if(str_detect(sideeffects[i,2],"\\([^()]+\\)")){   # regular expression
    sideeffects[i,2] <- sideeffects[i,1] # If drugname is lengthy chemical then replace with drugbank ID
    z=z+1
  }
}
dim(sideeffects)  
paste("Made ",z, " chemical name replacements with DB id's")


# SAme problem with indications so will replace drugs/chemicals with long names such as :
# "1-BENZYL-4-[(5,6-DIMETHOXY-1-INDANON-2-YL)METHYL]PIPERIDINE." with drugbank ID using a regex.
z=0
for (i in 1:nrow(indications)){
  if(str_detect(indications[i,2],"\\([^()]+\\)")){   # regular expression
    indications[i,2] <- indications[i,1] # If drugname is lengthy chemical then replace with drugbank ID
    z=z+1
  }
}
dim(indications)  
paste("Made ",z, " chemical name replacements with DB id's")


# REMOVE COMMON TOP 10% OF SIDE-EFFECTS - IF NOT THEN WE RISK MAKING OUR SEARCH TOO GENERAL
sideTable <- table(sideeffects$side_effect_name)
sideTable <- sort(sideTable,decreasing=TRUE) # 5,734 unique side-effects(SE), with 1,223 unique drugs
tenpercent_drugs <- round(length(unique(sideeffects[,2])) *.10 )# therefore remove top 10% of SE associated with more than 122 drugs
freq_se <- nrow(sideeffects)

drugfreq <- seq(0,0,length.out=freq_se)  # create a vector of empty freqs for each side-effects
sideeffects <- mutate(sideeffects,drugfreq) # add extra column for number of drugs with that side-effect
side_df <- as.data.frame(sideTable, stringAsFactors=FALSE)
names(side_df)[names(side_df)=="Freq"] <- "drugfreq"   # rename Freq and var1 to match the sideffect variables of drugfreq and side_effect_name
names(side_df)[names(side_df)=="Var1"] <- "side_effect_name"

# Using dplyr to join the freq of how many drugs inccur a side-effect
# tempside <- sideeffects
sideeffects <- sideeffects  %>%
  select(drugbank_id, drugbank_name, umls_cui_from_meddra,side_effect_name) %>%
  left_join(side_df, by="side_effect_name")

sideeffects <- filter(sideeffects, drugfreq < tenpercent_drugs) # remove the side-effects that occur in more than 122 drugs (10%)


# This is correct for number of drugs per side-effect
results<-table(sideeffects$side_effect_name)
results<-sort(results,decreasing=TRUE)
plot(results,xlab="Number of side-effects",ylab="Number of drugs",cex.main=1)
abline(v=(seq(0,1000,50)), col="darkgray", lty="dotted")
abline(h=(seq(0,1500,200)), col="darkgray", lty="dotted")

# This is correct for number of side -effects per drug.
results<-table(sideeffects$drugbank_name)
results<-sort(results,decreasing=TRUE)
plot(results,xlab="Number of drugs",ylab="Number of side-effects",cex.main=1)
abline(v=(seq(0,4500,200)), col="darkgray", lty="dotted")
abline(h=(seq(0,2500,200)), col="darkgray", lty="dotted")
















