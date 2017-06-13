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
  sidefx[i] <- selist[i]} # now populate it from selist

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
    cat("\nNo common side effects found for all drugs in your list- pruning search space #1")
    npdrugs <- round(ndrugs/2)
    cat("\nRandomly selecting",npdrugs,"drugs")
    prune_drug_list <- sample(dlist, npdrugs)  # sample the new shorter list of drugs
    cat("\nI have selected",prune_drug_list)
    se_prune_list <- get_sideeffects(prune_drug_list) # get the list of side-effects for our shorter list of drugs
    allSE <- get_commonsefx(prune_drug_list,se_prune_list)}
  
  if(length(allSE) >= 3){
    cat("\nFound",length(allSE),"common side effects for all drugs in your list.")
    repstuff <- list(sideeffects=allSE, drugs=dlist)
    return(repstuff)}
  
  if(length(allSE) < 3 & ndrugs > 2){  ### Try pruning again
    cat("\nAgain, no common side effects found for all drugs in your list- pruning search space again #2")
    npdrugs <- length(prune_drug_list)
    npdrugs <- round(npdrugs/2)
    cat("\nRandomly selecting",npdrugs,"drugs")
    prune_drug_list <- sample(prune_drug_list, npdrugs)  # sample the new shorter list of drugs
    cat("\nI have selected",prune_drug_list)
    se_prune_list <- get_sideeffects(prune_drug_list) # get the list of side-effects for our shorter list of drugs
    allSE <- get_commonsefx(prune_drug_list,se_prune_list)}
  
  if(length(allSE) >= 3){
    cat("\nFound",length(allSE),"common side effects for all drugs in your list.")
    repstuff <- list(sideeffects=allSE, drugs=prune_drug_list)
    return(repstuff)}
  
  if(length(allSE) < 3 & ndrugs > 8){  ### Try pruning again
    cat("\nAgain, no common side effects found for all drugs in your list- pruning search space again #3")
    npdrugs <- length(prune_drug_list)
    npdrugs <- round(npdrugs/2)
    cat("\nRandomly selecting",npdrugs,"drugs")
    prune_drug_list <- sample(prune_drug_list, npdrugs)  # sample the new shorter list of drugs
    cat("\nI have selected",prune_drug_list)
    se_prune_list <- get_sideeffects(prune_drug_list) # get the list of side-effects for our shorter list of drugs
    allSE <- get_commonsefx(prune_drug_list,se_prune_list)}
  
  if(length(allSE) < 3 & ndrugs > 4){  ### Try pruning again
    cat("\nAgain, no common side effects found for all drugs in your list- pruning search space again #4")
    npdrugs <- length(prune_drug_list)
    npdrugs <- round(npdrugs/2)
    cat("\nRandomly selecting",npdrugs,"drugs")
    prune_drug_list <- sample(prune_drug_list, npdrugs)  # sample the new shorter list of drugs
    cat("\nI have selected",prune_drug_list)
    se_prune_list <- get_sideeffects(prune_drug_list) # get the list of side-effects for our shorter list of drugs
    allSE <- get_commonsefx(prune_drug_list,se_prune_list)}
  
  if(length(allSE) >= 3){
    cat("\nFound",length(allSE),"common side effects for all drugs in your list.")
    repstuff <- list(sideeffects=allSE, drugs=prune_drug_list)
    return(repstuff)}
  if(length(allSE) <1){
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

# Create the Venn diagram but recall that no more than 5 drugs can be plotted - limitation of vennpackage.
plot_venn <- function(dlist,sidefx,index) {
  plot.new()
  dlen <- length(dlist)
  if(dlen==2){alpha<-c(0.5,0.5);fill=c("red", "blue"); just=list(c(0.6,1),c(0,0))};
  if(dlen==3){alpha<-c(0.5,0.5,0.5);fill=c("red", "blue","green");just=list(c(0.6,1),c(0,0),c(0,0))};
  if(dlen==4){alpha<-c(0.5,0.5,0.5,0.5);fill=c("red", "blue","green","pink");just=list(c(0.6,1),c(0,0) ,c(0,0),c(1,1))};
  if(dlen==5){alpha<-c(0.5,0.5,0.5,0.5,0.5);fill=c("red", "blue","green","pink","yellow");just=list(c(0.6,1),c(0,0),c(0,0),c(1,1),c(1,0))};
  
  venn.plot <- venn.diagram((sidefx[index]), 
                filename=NULL, 
                fill=fill, 
                alpha=alpha, 
                cex = 2, 
                cat.fontface=2, 
                margins =c(12,12),
                cat.just=just,
                cat.cex=2,
                cat.dist = rep(0.1, dlen),
                lty = "blank",
                category.names=(dlist[1:dlen]))
  grid.newpage()
  grid.draw(venn.plot)}

# get_drug_targets() drug by drug get the known targets from ontargets structure
get_drug_targets <- function(thedrug){
  #str(thedrug) #debug
  temp_targets <- filter(ontargets, drugbank_id == thedrug$drugbank_id) 
  x <- as.character((temp_targets$uniprot_id))
  #str(x)
  y <- select(org.Hs.eg.db, x, columns=c("SYMBOL","ENTREZID"),keytype="UNIPROT")
  #y <- select(org.Hs.eg.db, x, cols=c("SYMBOL", "ENTREZID"))
  #str(y)
  drug_targets <- cbind(y,rep(thedrug$drugbank_id,nrow(y))) # add Symbol & Uniprot ids.
  drug_targets <- cbind(drug_targets,rep(thedrug$name,nrow(y))) # add drug name
  drug_targets <- cbind(drug_targets,rep(thedrug$atc_codes,nrow(y))) # add atc code
  names(drug_targets)[3]<-"entrezid"
  names(drug_targets)[4]<-"drugbank_id"
  names(drug_targets)[5]<-"name"
  names(drug_targets)[6]<-"atc_codes"
  return(drug_targets)
}

# get_all_drug_targets() uses the get_drug_targets() function
get_all_drug_targets <- function(joint_list){
  ndrugs <- nrow(joint_list)
  drug_targets <- vector()
  for (i in 1:ndrugs){
    details <- get_drug_targets(joint_list[i,])
    drug_targets <- rbind(drug_targets,details)
  }
  return(drug_targets)  
}

# names_ids() creates a dataframe that adds drugbank_id to rest of data. 
# collates info for the conventional drugs
names_ids <- function(drug_list){
  temp_atc <- data.frame(drugbank_id=character(0),name=character(0),type=character(0),
                         atc_codes=character(0),stringsAsFactors = FALSE); 
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

# Calculate some statistics about the disease gene network
get_gstatistics <- function(gt) {
  modularity<-modularity(gt, membership(cluster_walktrap(gt)))
  avepath <- average.path.length(gt)
  nedges <- ecount(gt)
  nverts <- vcount(gt)
  transit <- transitivity(gt)
  degree <- igraph::degree(gt)
  diameter <- diameter(gt,weights=NA)
  #connect <- is.connected(gt)
  clos<- closeness(gt)
  between <- betweenness(gt,directed=FALSE)
  dens <-graph.density(gt)
  hubness <-hub_score(gt)$vector
  #author <-authority.score(gt)$vector
  bpow <- mean(igraph::power_centrality(gt))
  gstats <- data.frame(modularity,avepath,nedges,nverts,transit,degree,diameter,
                       clos,between,dens,hubness,bpow )
  return(gstats)
}


# id2name() takes a sequence of drugbank id's and returns the actual drugbank names
# expects to see global variables currentdrugs and drugnames
id2name <- function(dbid) {
  #cat("\n",dbid)
  dlen <- length(dbid)
  drugn <- ""
  for (i in 1:dlen){
    din <- grep(dbid[i],drugs)
    #cat("\ndebug",din)
    drugn[i] <- drugnames[din]
  }
  
  return(drugn) 
}

# obtained from stackoverflow question.
jaccard <- function(m) {
  ## common values:
  A = tcrossprod(m)
  ## indexes for non-zero common values
  im = which(A > 0, arr.ind=TRUE)
  ## counts for each row
  b = rowSums(m)
  
  ## only non-zero values of common
  Aim = A[im]
  
  ## Jacard formula: #common / (#i + #j - #common)
  J = sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
  )
  
  return( J )
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





