# reviewers_JBI.R
# work started: 03/04/2017
# Responses to reviwer comments on repostitioning paper submitted to:
# Journal of Biomedical Informatics
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

# SIDER4.1 contains drugs/chemicals with long names such as "(1S,4S)-4-(3,4-dichlorophenyl)-N-methyl-1,2,3,4-tetrahydronaphthalen-1-amine."
# these are causing string searching problems (as well as looking very untidy!) so will replace them with drugbank ID using regex.
z=0
for (i in 1:nrow(sideeffects)){
  if(str_detect(sideeffects[i,2],"\\([^()]+\\)")){
    sideeffects[i,2] <- sideeffects[i,1]
    z=z+1
  }
}
dim(sideeffects)  
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


# Get the four popular drugs used to treat alzheimers: 'donepezil','galantamine','rivastigmine','Memantine'
# and obtain the side effects for each one. NOTE: spellings/case important!!
donepezil <- sideeffects$side_effect_name[grep('Donepezil',sideeffects$drugbank_name)]
galantamine <-sideeffects$side_effect_name[grep('Galantamine',sideeffects$drugbank_name)]
rivastigmine <-sideeffects$side_effect_name[grep('Rivastigmine',sideeffects$drugbank_name)]
#memantine <- sideeffects$side_effect_name[grep('Memantine',sideeffects$drugbank_name)]

# Now get only unique side-effects for each drug as duplications can occur.
donepezil<-unique(donepezil)
galantamine<-unique(galantamine)
rivastigmine<-unique(rivastigmine)
#memantine <- unique(memantine)

#---------------- Now determine combined side effects of the AD drugs ------------
universe <- unique(c(donepezil,galantamine,rivastigmine,memantine))
GroupA <-universe %in% donepezil
GroupB <-universe %in% galantamine
GroupC <-universe %in% rivastigmine
#GroupD <-universe %in% memantine

## All Side effects that are common between all AD drugs (see Venn diagram.)
allSE <- universe[GroupA & GroupB & GroupC]# & GroupD]

## overlaps may be used as I've pruned probably too many side-effects
overlaps1 <- universe[GroupA & GroupB]
overlaps2 <- universe[GroupA & GroupC]
overlaps3 <- universe[GroupB & GroupC]
overlaps4 <- universe[GroupA & GroupB & GroupC]
overlaps <- unique(c(overlaps1,overlaps2,overlaps3,overlaps4))

allSE <- overlaps

#--------- Now determine which drugs have similar side effects to our universe of side-effects -----
# But how many SE's would be deemed 'enough' for a drug to be classed as similar????
# Give a count of SE commonality for each drug in list

NoDrugs <- unique(sideeffects$drugbank_name) # We have 1,223 drugs
DBID <- unique(sideeffects$drugbank_id)

Names <- letters[1:5];
Dates<- 1:length(NoDrugs);
# allDrugs will hold summary information for each drug: Drugname; Nosideffects and %coverage
allDrugs<- data.frame(drugnames=NoDrugs,NoSideEffects=vector(mode="numeric",length=length(Dates)),
                      coverage=vector(mode="numeric",length=length(Dates)),drugbankid=DBID, stringsAsFactors = FALSE); 
allDrugNames <- lapply(NoDrugs, as.character)

# TAKES FIVE MINUTES TO CALCULATE
for (i in 1:length(NoDrugs)){
  allDrugs[i,1]<-as.character(allDrugNames[i])
  TempDrug <- filter(sideeffects,drugbank_name==allDrugNames[i])
  allDrugs[i,2] <- length(unique(TempDrug$side_effect_name)) # side effects found for each drug in database
  TempDrug <- unique(TempDrug$side_effect_name)
  coverage <- TempDrug %in% allSE
  coverage <- coverage[!is.na(coverage)] # get rid of NA when we dont get a matching side-effect
  allDrugs[i,3] <- (length(coverage[coverage==TRUE])/length(allSE))*100 # percentage coverage of joint SE of the four alzheiemers drugs
}

summary(allDrugs[,3])

# Get list of drugs with at least 25% shared side effects with our three Alzheimers drugs. 
# However we do need to remove them as we dont want them to appear in the list
candidates <- subset(allDrugs, allDrugs$coverage > 25.0)
candidates <- filter(candidates, drugnames != "Donepezil")  
candidates <- filter(candidates, drugnames != "Rivastigmine")
candidates <- filter(candidates, drugnames != "Galantamine")

candidates <- candidates[order(-candidates$coverage),]  # sort descending so we get top ten for latex table.
rownames(candidates) <- NULL

# create a LaTex table for my LaTex document using the XTABLE package.
tli.table <- xtable(candidates)
digits(tli.table)[c(2,6)] <- 0
print(tli.table,floating=FALSE)


# Now create the Venn diagram.
plot.new()
venn.plot <- venn.diagram(list(donepezil,galantamine,rivastigmine), 
                          NULL, 
                          fill=c("red", "blue","green"), 
                          alpha=c(0.5,0.5,0.5), 
                          cex = 2, 
                          cat.fontface=2, 
                          margins =c(10,10),
                          cat.cex=2,
                          #main = "Venn Diagram showing shared side effects for donepezil,galantamine,rivastigmine",
                          category.names=c("donepezil", "galantamine","rivastigmine"))
grid.draw(venn.plot)

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


# For each drug identify its on-targets and off-targets

















