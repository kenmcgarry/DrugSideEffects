# reviewers_JBI_DrugList.R
# Work started: 03/04/2017
# Responses to reviewer comments on repostitioning paper submitted to Journal of Biomedical Informatics
# The code here will obtain a list of known drugs that combat the disease of interest.
# The known side-effects of these drugs will be identified and a particular combination
# will be used to search for other drugs that match a certain percentage of similar side-effects.
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

# You can run a previous R file (reviewers_JBI_LoadData.R) to create the data from scratch or simply load
# in the RData session file.
# load("repositionJBI.RData") 

# Create a structure to hold all known drugs treating your particular disease of interest
# Obviously instantiate these functions before calling them.
#drug_list <- c('Donepezil','Galantamine','Rivastigmine')  # this a manual example for Alzheimers

drug_list <- getdrugs("C0002395")  # using the umls code for Alzheimers what drugs are used to treat it?
  
se_list <- getsideeffects(drug_list) # what side-effects do these drugs have?

repos_list <- getreposdrugs(drug_list,se_list) # what drugs could be candidates for repositioning for our disease

plotvenn(drug_list,se_list) # No more than 5 drugs max, or Venn will not work!!!!

## --------------------- FUNCTION DEFINITIONS -----------------------

# getdrugs() assumes that "indications" dataframe is already loaded. You must provide getdrugs() 
# with the "umls_cui_from_meddra" code for your disease. It will return the drugs known to be used...
# e.g. C000239 is the code for Alzheimer's. Using the code is less error prone than typing in disease name.
getdrugs <- function(umls) {
  ilist <- filter(indications, umls_cui_from_meddra == umls)
  if(nrow(ilist) > 0){
    for (j in 1:nrow(ilist)){
      cat("drug",j,"is", ilist$drugbank_name[j],"\n")
    }
  return(ilist$drugbank_name)
  }else{
    cat("\n","Sorry, no drugs found...check umls code is correct for your disease")
    return(NULL)}
}


# getsideeffect() assumes that "sideeffects" and "drug_atc "dataframes are already loaded
getsideeffects <- function(yourdrugs) {
  if(length(yourdrugs) < 1){
    cat("\n","Sorry, you havent any drugs in your list, so no side-effects.")
    return(NULL)}
  ndrugs <- length(yourdrugs)
  se <-rep("",ndrugs)
    for (j in 1:ndrugs){
      tempse <- filter(sideeffects, drugbank_name == yourdrugs[j])
      cat("drug",j,"is", yourdrugs[j],"with", nrow(tempse),"side-effects\n")
      se[j] <- tempse %>%
        select(side_effect_name)
    }
  return(se)  # return a list of size n
}    
   

# getreposdrugs() receives  a list of drugs and a list of side-effects 
getreposdrugs <- function(dlist,selist){
  ndrugs <- length(dlist)
  sidefx <- rep("",ndrugs)  # for every drug in dlist, create a membership entry holding the side-effects
  for (i in 1:ndrugs){
    sidefx[i] <- selist[i] # now populate it from selist
  }
  
  group <- selist # create same sized lists but "group" will contain TRUE or FALSE strings for each side-effect
  universe <- unique(unlist(sidefx))  # universe is all the side-effects for each drug
  
  for (i in 1:ndrugs){
    templist <- unlist(selist[i])   # templist = side-effect names drug by drug
    group[[i]] <- universe %in% templist  # group = TRUE or FALSE values depending on membership in universe
    cat("\n\n",(table(group[i])))
  }
  
  allSE <- universe[Reduce('&',group)] # akrun from stackoverflow.com, 20/04/2017
  
  return(allSE)
}

## The more current treatment drugs are used then less likely to get side-effects common to all of them,
## either use less drugs in inital search or use some combination of pairwise overlaps

overlaps1 <- universe[GroupA & GroupB]
overlaps2 <- universe[GroupA & GroupC]
overlaps3 <- universe[GroupB & GroupC]
overlaps4 <- universe[GroupA & GroupB & GroupC]
overlaps <-  unique(c(overlaps1,overlaps2,overlaps3,overlaps4))

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
  allDrugs[i,3] <- (length(coverage[coverage==TRUE])/length(allSE))*100 # percentage coverage of joint SE of the  drugs
}

summary(allDrugs[,3])

# Get list of drugs with at least 25% shared side effects with our three Alzheimers drugs. 
# However we do need to remove  Alzheimers drugs as we dont want them to appear in the list.
candidates <- subset(allDrugs, allDrugs$coverage > 25.0)
candidates <- filter(candidates, drugnames != "Donepezil")  
candidates <- filter(candidates, drugnames != "Rivastigmine")
candidates <- filter(candidates, drugnames != "Galantamine")

candidates <- candidates[order(-candidates$coverage),]  # sort descending so we get top ten for latex table.
rownames(candidates) <- NULL

# need to rename variables so "key" can work with drugnames - i.e. keep name consistent across data structures
names(drug_atc)[names(drug_atc)=="name"] <- "drugnames"

# we really only need the atc_codes; we keep drugnames for the database key
drug_atc <- drug_atc %>%
  select(drugnames,atc_codes)

# join the atc_codes from the drug_atc data structure to candidates data
newcandidates <- candidates  %>%
  select(drugnames, NoSideEffects, coverage, drugbankid) %>%
  left_join(drug_atc, by="drugnames")

# Create a LaTex table for my LaTex document using the XTABLE package.
# Warning: the file we use for ATC codes does not have them for chemical/compound data. So <NA> will appear
# where DB07701 etc occur, checking DRUGBANK by hand reveals little information.
tli.table <- xtable(newcandidates)
digits(tli.table)[c(2,6)] <- 0
print(tli.table,floating=FALSE)


# Create the Venn diagram but recall that no more than 5 drugs can be plotted - limitation of package.
plotvenn <- function(dlist,sidefx) {
plot.new()
venn.plot <- venn.diagram((sidefx[1:5]), 
                filename=NULL, 
                fill=c("red", "blue","green","pink", "yellow"), 
                alpha=c(0.5,0.5,0.5,0.5,0.5), 
                cex = 2, 
                cat.fontface=2, 
                margins =c(10,10),
                cat.cex=2,
                cat.dist = rep(0.1, 5),
                lty = "blank",
                category.names=(dlist[1:5]))
                grid.draw(venn.plot)}


# ------------ interesting dynamic variable creation -------
#var_names <- paste("v", 1:3, sep="")
#for (v in druglist){ 
#  assign(v, select() %>% filter(sideeffects, drugbank_name == v))
#}

#----- for stackoverflow only -----------
#universe <- c("ted","sara","fred","billy")
#group1 <- as.logical(c("TRUE","TRUE","TRUE","TRUE"))
#group2 <- as.logical(c("FALSE","TRUE","FALSE","TRUE"))
#group3 <- as.logical(c("FALSE","TRUE","TRUE","TRUE"))
#group <- list(group1, group2, group3)





