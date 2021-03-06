# reviewers_JBI_DrugList.R
# Work started: 03/04/2017
# Responses to reviewer comments on repostitioning paper submitted to Journal of Biomedical Informatics
# The code here will obtain a list of known drugs that combat the disease of interest.
# The known side-effects of these drugs will be identified and a particular combination
# will be used to search for other drugs that match a certain percentage of similar side-effects.
#
# source("http://bioconductor.org/biocLite.R")
# biocLite("nameofpackage")

setwd("C:/R-files/sider")


## --------------------- FUNCTION DEFINITIONS -----------------------

# getdrugs() assumes that "indications" dataframe is already loaded. You must provide getdrugs() 
# with the "umls_cui_from_meddra" code for your disease. It will return the drugs known to be used...
# e.g. C000239 is the code for Alzheimer's. Using the code is less error prone than typing in disease name.
# The restricted list of drugs we cant use is passed to this function.
get_drugs <- function(umls,rlist) {
  ilist <- filter(indications, umls_cui_from_meddra == umls)
  ilist <- setdiff(ilist$drugbank_name,rlist$drugbank_name)
  if(length(ilist) > 0){
    for (j in 1:length(ilist)){
      cat("\ndrug",j,"is", ilist[j])
    }
  return(ilist)
  }else{
    cat("\n","Sorry, no drugs found...check umls code is correct for your disease")
    return(NULL)}
}


# getsideeffect() assumes that "sideeffects" and "drug_atc "dataframes are already loaded
get_sideeffects <- function(yourdrugs) {
  if(length(yourdrugs) < 1){
    cat("\n","Sorry, you havent any drugs in your list, so no side-effects.")
    return(NULL)}
  ndrugs <- length(yourdrugs)
  se <-rep("",ndrugs)
    for (j in 1:ndrugs){
      tempse <- filter(sideeffects, drugbank_name == yourdrugs[j])
      cat("\ndrug",j,"is", yourdrugs[j],"with", nrow(tempse),"side-effects")
      se[j] <- tempse %>%
        select(side_effect_name)
    }
  return(se)  # return a list of size n
}    

# getcommonsefx() is the guts   
get_commonsefx <- function(dlist,selist){
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
  cat("\n\n",(table(group[i])))}
  
  commonSE <- universe[Reduce('&',group)] # akrun from stackoverflow.com, 20/04/2017
  return(commonSE)
}

# get_repos_sideeffects() receives  a list of drugs and a list of side-effects and finds SE common to all drugs
get_repos_sideeffects <- function(dlist,selist){
  ndrugs <- length(dlist)
  allSE <- get_commonsefx(dlist,selist)
  # if we have fewer than 3 common side-effects then prune drugs
  if(length(allSE) < 3 & ndrugs > 2){
    cat("\nNo common side effects found for all drugs in your list- pruning search space")
    npdrugs <- round(ndrugs/2)
    cat("\nRandomly selecting",npdrugs,"drugs")
    
    prune_drug_list <- sample(dlist, npdrugs)  # sample the new shorter list of drugs
    cat("\nI have selected",prune_drug_list)
    se_prune_list <- get_sideeffects(prune_drug_list) # get the list of side-effects for our shorter list of drugs
    allSE <- get_commonsefx(prune_drug_list,se_prune_list)
    plotvenn(prune_drug_list,se_prune_list)
  }else{
    cat("\nFound",length(allSE),"common side effects found for all drugs in your list.")
    return(allSE)}
  if(allSE<1){
    cat("\nNo common side effects found for all drugs in your list- Exiting.")
    return(NULL)}
}

get_overlaps <- function(dlist,selist){
# getoverlaps() returns the 'optimum' overlaps of shared side-effects - used only when the
# central common to all is zero or too small.
# NOT WORKING YET!
## The more current treatment drugs are used then less likely to get side-effects common to all of them,
## either use less drugs in inital search or use some combination of pairwise overlaps.
  overlaps1 <- universe[GroupA & GroupB]
  overlaps2 <- universe[GroupA & GroupC]
  overlaps3 <- universe[GroupB & GroupC]
  overlaps4 <- universe[GroupA & GroupB & GroupC]
  overlaps <-  unique(c(overlaps1,overlaps2,overlaps3,overlaps4))
  allSE <- overlaps
  return(overlaps)
}

# search_similardrugs() will determine which drugs have similar side effects to our universe of side-effects 
# But how many SE's would be deemed 'enough' for a drug to be classed as similar?
# Give a count of SE commonality for each drug in list
search_similardrugs <- function(allSE,drug_list){
  NoDrugs <- unique(sideeffects$drugbank_name) # We have 1,223 drugs
  DBID <- unique(sideeffects$drugbank_id)

  Names <- letters[1:5];
  Dates <- 1:length(NoDrugs);
  # allDrugs will hold summary information for each drug: Drugname; Nosideffects and %coverage
  allDrugs<- data.frame(drugnames=NoDrugs,NoSideEffects=vector(mode="numeric",length=length(Dates)),
                      coverage=vector(mode="numeric",length=length(Dates)),drugbankid=DBID, stringsAsFactors = FALSE); 
  allDrugNames <- lapply(NoDrugs, as.character)

# Might take a while to calculate on your computer.
  for (i in 1:length(NoDrugs)){
    allDrugs[i,1]<-as.character(allDrugNames[i])
    TempDrug <- filter(sideeffects,drugbank_name==allDrugNames[i])
    allDrugs[i,2] <- length(unique(TempDrug$side_effect_name)) # side effects found for each drug in database
    TempDrug <- unique(TempDrug$side_effect_name)
    coverage <- TempDrug %in% allSE
    coverage <- coverage[!is.na(coverage)] # get rid of NA when we dont get a matching side-effect
    allDrugs[i,3] <- (length(coverage[coverage==TRUE])/length(allSE))*100 # percentage coverage of joint SE of the  drugs
  }

# Get list of drugs with at least 25% shared side effects with our Alzheimers drugs and remove the 
# drugs used to treat our disease as we dont want them to appear in the list of candidates.
  candidates <- subset(allDrugs, allDrugs$coverage > 25.0)
  for (i in 1:length(drug_list)){
    candidates <- filter(candidates, drugnames != drug_list[i])}
 
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

return(newcandidates)
}



# latex_to_table() will create a LaTex table for my LaTex document using the XTABLE package.
# It expects a dataframe with names to be passed to it.
latex_to_table <- function(thedata){
  tli.table <- xtable(thedata)
  #digits(tli.table)[c(2,6)] <- 0 # uncomment this for those columns that use numbers
  print(tli.table,floating=FALSE)
}

# Create the Venn diagram but recall that no more than 5 drugs can be plotted - limitation of package.
plot_venn <- function(dlist,sidefx) {
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
# stops errors occurring. You need to provide 4 peices of information.
repair_atc <- function(drugbank_id,name,atc_codes,type)  {
  drug_atc <- rbind(drug_atc,c(drugbank_id,name,type,atc_codes))
  tail(drug_atc)
  return(drug_atc)
}

# ------------ interesting dynamic variable names creation -------
#var_names <- paste("v", 1:3, sep="")
#for (v in druglist){ 
#  assign(v, select() %>% filter(sideeffects, drugbank_name == v))
#}

#----- for my stackoverflow plea for help -----------
# good reply from Akrun
# http://stackoverflow.com/questions/43522947/using-r-to-access-a-structure-dynamically-with-pasted-comands
#universe <- c("ted","sara","fred","billy")
#group1 <- as.logical(c("TRUE","TRUE","TRUE","TRUE"))
#group2 <- as.logical(c("FALSE","TRUE","FALSE","TRUE"))
#group3 <- as.logical(c("FALSE","TRUE","TRUE","TRUE"))
#group <- list(group1, group2, group3)





