# repositionJuly2016.r
# Ken McGarry 27/7/2016
# Use cygwin to download this file as web browsers like firefox are not up to 8GB downloads!
# http://biocluster.ucr.edu/~tbackman/bioassayR/pubchem_protein_only.sqlite
# source("http://bioconductor.org/biocLite.R")
# biocLite("dplyr")

library(DO.db)
library(DOSE)
library(graph)
library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(tidyr)
library(igraph)
library(biomaRt)
library("org.Hs.eg.db") # for ncbi to gene conversion 
library(VennDiagram)
library(paxtoolsr)
library(rJava)
library(ReactomePA)
library(ChemmineR)
library(ChemmineOB)
library(bioassayR)
library(RSQLite)
library(igraph)
library(ROCR)
library(xtable)


database <- connectBioassayDB("C://R-files//bigfiles//pubchem_protein_only.sqlite")
aspirinTargets <- activeTargets(database, 2244) # 2244 = aspirin
targetList <- allTargets(database)

activecompounds <- activeAgainst(database,"166897622")
  
queryCids <- c("2244", "3715", "2662", "3033", "133021","44563999", "44564000", "44564001", "44564002")
queryCids <- c("2244", "3715", "2662", "3033", "133021","44563999", "44564000", "44564001", "44564002")
myAssaySet <- getBioassaySetByCids(database, queryCids)
myFp <- bioactivityFingerprint(bioassaySet=myAssaySet)

# Perform activity profile similarity searching with the FPset object, by comparing the 
# first compound to all compounds.
fpSim(myFp[1], myFp, method="Tanimoto")

# Compute an all-against-all bioactivity fingerprint similarity matrix for these compounds.
simMA <- sapply(cid(myFp), function(x) fpSim(myFp[x], myFp, sorted=FALSE, method="Tanimoto"))

# Convert similarity matrix to a distance matrix and perform hierarcheal clustering.
clusterResults <- hclust(as.dist(1-simMA), method="single")
plot(clusterResults)

# How to conver between the NCBI id and the protein symbol(egSYMBOL);(egENSEMBL);
ids=c("2244", "2662")
genes=c("COX1","COX2")
cox1="5742"

unlist(mget(x=rownames(aspirinTargets),envir=org.Hs.egSYMBOL))


#===================================
# Ming Hao fingerprint example using rcdk package (personal communication 2/6/2016.

# 1. download bioassay file from pubchem, say 540252
# 2. delete first four rows (contains text we dont need)
# 3. read in file and obtain pubchem_cid's
# 4. download the SDF file using the CID number at PubChem Download Service

rawdata <- read.csv("C:\\R-files\\pubchem\\540252.csv",header=TRUE)

head(rawdata[,3])  # pubchem_cid in column 3 for each tested substance (587 in total)
table(rawdata[,4])  # activity_outcome in column 4 for each tested substance (587 in total)

sdfset <- read.SDFset("C:\\R-files\\pubchem\\compounds.sdf") 
fpset<-fingerprintOB(sdfset,fingerprintName="FP2")
view(fpset[1:6]) # view the fingerprints
cid(fpset) <- sdfid(sdfset) # overwrite CMP1 to CMP6 with proper id names e.g. 646043
cid(fpset)

sdf.visualize(sdfset) 
ChemmineR::groups(sdfset, groups="fctgroup", type="countMA") 
propma <- atomcountMA(sdfset, addH=FALSE) 
boxplot(propma, col="blue", main="Atom Frequency") 
data(atomprop)
atomprop[1:8,] 

# download direct from pubchem rather than use browser
compounds <- getIds(as.numeric(rawdata[,3]))
compounds <- getIds(as.numeric(rawdata[1:450,3]))
compounds2 <- getIds(as.numeric(rawdata[450:587,3]))

setwd("C:/R-files/bigfiles")

con <- dbConnect(RSQLite::SQLite(), dbname = "pubchem_protein_only.sqlite")
dbListTables(con)
dbDisconnect(con)

con <- dbConnect(RSQLite::SQLite(), dbname = "chembl_21.db")
dbListTables(con)
dbDisconnect(con)

#-------------------- Drug reposition stuff for Drug Development and Industrial Pharmacy journal.--------------
# DB00472 Fluoxetine
# DB00215 Citalopram
# DB00715 Paroxetine
# DB00193 Tramadol
# DB00230 Pregabalin
# DB01238 Aripiprazole
# DB00334 Olanzapine
# DB01165 Ofloxacin
# DB00268 Ropinirole
# DB00704 Naltrexone

library(ChemmineR)

setwd("C:/R-files/bigfiles")
sdfset <- read.SDFset("structures.sdf")  # Also ensure library(chemmineR) loaded
valid <- validSDF(sdfset); sdfset <- sdfset[valid]   # remove invalid SDFs

drug <- grepSDFset("DB00472", sdfset, field="datablock", mode="subset")

comp <-unique(na.omit(as.numeric(unlist(strsplit(unlist(names(drug)), "[^0-9]+"))))) # gets ridof letters,keeps numbers
plot(sdfset[comp], print=FALSE) 


apset <- sdf2ap(sdfset) # create atom pairs 
apset 
cid(apset) <- sdfid(sdfset) # over stupid CMP names with proper id's
cid(apset) <- makeUnique(cid(apset))
cmp.search(apset,apset["DB00268"], type=3, cutoff = 0.3, quiet=TRUE) # works with the artficial "CMP" number NOT DB id!!
                                                            # you need to find out what CMP number replaces proper id!
cid(apset[1:4]) # Compound IDs 
ap(apset[1:4]) # Atom pair descriptors 
db.explain(apset[1])
cmp.similarity(apset[1], apset[2])
cmp.similarity(apset[1], apset[1]) 
## ---------------- Thomas Girkes code ----------------
# You want to take a look at the short tutorial here: https://goo.gl/iTZbfj
# Especially, the distinction between sdfid and cid is important. The latter is the slot available in most object types 
# in ChemmineR (e.g. SDFset and APset). Drugbank is a specialty case, somehow they decided not to assign their IDs to the 
# SDF ID slot but keep this info nicely hidden in the data block of the SDF. However, with datablock2ma you can get this 
# info and assign it to your SDFset's cid slot. After this you can perform ID based subsetting with the bracket operator 
# (highly preferred) rather than a full text search with grepSDFset. The following demonstrates this for the Drugbank SDF 
# that I just downloaded from their FTP site.
## Import of Drugbank SDF
library(ChemmineR)

# create lookup table to convert between drugbank ids and drug name
drugs <- c("DB00843"= "Donepezil", "DB00674"="Galantamine","DB00989"="Rivastigmine",
           "DB00472"="Fluoxetine","DB00215"="Citalopram","DB00715"="Paroxetine",
           "DB00193"="Tramadol","DB00230"="Pregabalin","DB01238"="Aripiprazole",
           "DB00334"="Olanzapine","DB01165"="Ofloxacin","DB00268"="Ropinirole",
           "DB00704"="Naltrexone")


dnames <- c("Donepezil","Galantamine","Rivastigmine",
            "Fluoxetine","Citalopram","Paroxetine",
            "Tramadol","Pregabalin","Aripiprazole",
            "Olanzapine","Ofloxacin","Ropinirole",
            "Naltrexone")



setwd("C:/R-files/bigfiles")
sdfset <- read.SDFset("structures.sdf")
valid <- validSDF(sdfset)
sdfset <- sdfset[valid]



## Assign drugbank IDs from datablock to cid slot
blockmatrix <- datablock2ma(datablocklist=datablock(sdfset))
cid(sdfset) <- as.character(blockmatrix[,"DRUGBANK_ID"])

## Generate APset and FPset (note FPset: has better search performance)
apset <- sdf2ap(sdfset)
fpset <- desc2fp(apset, descnames=2048, type="FPset")

## Subsetting by cid slot using drugbank ids should work now consistently
sdfset["DB00472"]
apset["DB00472"]

fpdrugs <- fpset[drugs]
params <- genParameters(fpdrugs)  # including Donepezil, Galantamine and Rivastigmine

results1 <- fpSim(fpdrugs[[1]], fpdrugs, top=13, parameters=params,method="Tanimoto") 
results2 <- fpSim(fpdrugs[[2]], fpdrugs, top=13, parameters=params,method="Tanimoto") 
results3 <- fpSim(fpdrugs[[3]], fpdrugs, top=13, parameters=params,method="Tanimoto") 

results1 <- cbind(drugs[rownames(results1)],results1)
results2 <- cbind(drugs[rownames(results2)],results2)
results3 <- cbind(drugs[rownames(results3)],results3)

library(xtable)
print.xtable(xtable(results1))
print.xtable(xtable(results2))
print.xtable(xtable(results3))

print.xtable(results,label='comp',caption = 'Similarity measures for the 10 drugs based on 2048-length, atom-pair fingerprints')

# Creates the similarity score matrix and clusters them.
simMA <- sapply(cid(fpdrugs), function(x) fpSim(x=fpdrugs[x], fpdrugs, sorted=FALSE)) 
colnames(simMA)<-as.character(drugs[cid(fpdrugs)])
rownames(simMA)<-as.character(drugs[cid(fpdrugs)])

library(ape)         
par(mar=c(5, 5, 5, 4))
plot(as.phylo(hc), cex = 0.9, label.offset = 0.01)

hc <- hclust(as.dist(1-simMA), method="single")
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=3), horiz=F) 
rect.hclust(hc,k=4,border="red")

#--------- gene ontology (GO) and disease ontology (DO) integration ------
library(clusterProfiler)
library(DOSE)
library(GOSemSim)
#--------------------- ONTOLOGY INTEGRATION -------------------------------
# Gene ontology (GO) and disease ontology (DO) can be used to enrich the statistical information 
# gained from the graphs by integrating them with biological knowledge.

# alzheimers DOID:10652
# parkinsons DOID:14330
# schizophrenia DOID:5419
# psychosis DOID:8646
# depression DOID:1596
# pain disorder DOID:0060164


# ropinrole -> parkinsons 
# tramadol -> analgesic
# pregabalin -> anticonvulsant
# citalopram -> antidepressants
# paroxetine -> antidepressants
# aripiprazole -> schizophrenia
# olanzapine -> antipsychotics

# Get the semantic similarity between our diseases using doSim() function.
hgncList<-c("APP","PSEN1","PSEN2","APOE");
listDO<-c("DOID:10652","DOID:14330","DOID:5419","DOID:8646","DOID:1596","DOID:0060164","DOID:1826","DOID:750","DOID:14262","DOID:13938")
listDO1<-c("DOID:14330","DOID:5419","DOID:8646","DOID:1596","DOID:0060164")
alzDO="DOID:10652"
s<-doSim(listDO,listDO,measure="Wang")
simplot(s,color.low="white",color.high="blue",labs=TRUE,digits=2,labs.size=5,
        font.size=14, xlab="",ylab="")

# ----- change row and column names from DOID numbers to disease names -----
DoNames<-c("alzheimers","parkinsons","schizophrenia","psychosis","depression","pain_disorder","epilepsy","ulcer","candidiasis","amenorrhea")
colnames(s)<-DoNames
rownames(s)<-DoNames

# -------------- gene ontology --------------------
# ------------ the following are STITCH protein-to-drug interactions for all 13 drugs -----------------
done <- read.delim('C:\\R-files\\reposition\\donepezil2016.txt', header=TRUE,sep='\t') 
gala <- read.delim('C:\\R-files\\reposition\\galantamine2016.txt', header=TRUE,sep='\t') 
riva <- read.delim('C:\\R-files\\reposition\\rivastigmine2016.txt', header=TRUE,sep='\t')

fluo <- read.delim('C:\\R-files\\reposition\\fluoxetine2016.txt', header=TRUE,sep='\t')
cita <- read.delim('C:\\R-files\\reposition\\citalopram2016.txt', header=TRUE,sep='\t')
para <- read.delim('C:\\R-files\\reposition\\paroxetine2016.txt', header=TRUE,sep='\t')
tram <- read.delim('C:\\R-files\\reposition\\tramadol2016.txt', header=TRUE,sep='\t')
preg <- read.delim('C:\\R-files\\reposition\\pregabalin2016.txt', header=TRUE,sep='\t')
arip <- read.delim('C:\\R-files\\reposition\\aripiprazole2016.txt', header=TRUE,sep='\t')
olan <- read.delim('C:\\R-files\\reposition\\olanzapine2016.txt', header=TRUE,sep='\t')
oflo <- read.delim('C:\\R-files\\reposition\\ofloxacin2016.txt', header=TRUE,sep='\t')
ropi <- read.delim('C:\\R-files\\reposition\\ropinirole2016.txt', header=TRUE,sep='\t')
nalt <- read.delim('C:\\R-files\\reposition\\naltrexone2016.txt', header=TRUE,sep='\t')

done <- done[,1:2]  # in all cases only keep protein-to-drug and protein-to-protein interactions
gala <- gala[,1:2]
riva <- riva[,1:2]
fluo <- fluo[,1:2]
cita <- cita[,1:2]
para <- para[,1:2]
tram <- tram[,1:2]
preg <- preg[,1:2]
arip <- arip[,1:2]
olan <- olan[,1:2]
oflo <- oflo[,1:2]
ropi <- ropi[,1:2]
nalt <- nalt[,1:2]

giant<-Reduce(function(x, y) merge(x, y, all=TRUE), list(done,gala,riva,fluo,cita,para,tram,preg,arip,olan,oflo,ropi,nalt))

# first some Gene ID conversion from the clusterprofiler function bitr()
dGene<-as.character((arip[,1]))
gene <- bitr(dGene,fromType = "SYMBOL",toType="ENTREZID",annoDb="org.Hs.eg.db")
ggo <- groupGO(gene=gene$ENTREZID, organism="human", ont="BP",level=3, readable=TRUE) #ont=[BP;MF;CC]
yy <- enrichPathway(gene[,2], pvalueCutoff=0.05)
head(summary(yy))
enrichMap(yy,fixed=FALSE) # "fixed=false" callups tkplot rather than normal plot

# ggo@result[1:5,1:4] # is a S4 object so use '@' instead of '$' 
ggo <- ggo@result[ggo@result$Count > 0,]  # only want results with our proteins!
ggo<-as.data.frame(ggo)           # convert to df
goDone <- ggo[order(-ggo$Count),]   # sort results in descending order i.e. most important first

xtable(goDone[1:10,1:4])


# ------------ Galantamine interacting proteins enrichment
dGene<-gala[,1]
gene <- bitr(dGene,fromType = "SYMBOL",toType="ENTREZID",annoDb="org.Hs.eg.db")
ggo <- groupGO(gene=gene$ENTREZID, organism="human", ont="BP",level=3, readable=TRUE)
#head(summary(ggo))

# ggo@result[1:5,1:4] # is a S4 object so use '@' instead of '$' 
ggo <- ggo@result[ggo@result$Count > 0,]  # only want results with our proteins!
ggo<-as.data.frame(ggo)           # convert to df
goGala <- ggo[order(-ggo$Count),]   # sort results in descending order i.e. most important first

xtable(goGala[1:10,1:4])


# ------------ Rivastigmine interacting proteins enrichment
dGene<-riva[,1]
gene <- bitr(dGene,fromType = "SYMBOL",toType="ENTREZID",annoDb="org.Hs.eg.db")
ggo <- groupGO(gene=gene$ENTREZID, organism="human", ont="BP",level=3, readable=TRUE)
#head(summary(ggo))

# ggo@result[1:5,1:4] # is a S4 object so use '@' instead of '$' 
ggo <- ggo@result[ggo@result$Count > 0,]  # only want results with our proteins!
ggo<-as.data.frame(ggo)           # convert to df
goRiva <- ggo[order(-ggo$Count),]   # sort results in descending order i.e. most important first

xtable(goRiva[1:10,1:4])


# plot a nice combined histogram of the GO enrichments for the three drugs
g1 <- data.frame(len=goGala$Count) # Ensure length is used otherwise ggplot2 wont bloody work.
r1 <- data.frame(len=goRiva$Count)
d1 <- data.frame(len=goDone$Count)

names(g1)[names(g1)=="goGala$Count"] <- "count"  # Replace the stupid names "goXXX$Count" with "count"
names(r1)[names(r1)=="goRiva$Count"] <- "count"
names(d1)[names(d1)=="goDone$Count"] <- "count"

# Give each drug its full name.
g1$drug <- 'Galantamine'
r1$drug <- 'Rivastigmine'
d1$drug <- 'Doneprezil'

#and combine into your new data frame drugLen
drugLen <- rbind(d1,r1,g1)
ggplot(drugLen, aes(len, fill = drug)) + geom_bar(pos="dodge" ) + xlab("GO terms with matching proteins") +  ylab("Frequency")


#----------- DOSE package ----------
alzDO="DOID:13938"
class(alzDO) = "DO"
TERM2NAME(alzDO)

#------- mygene package -------------
dGene<-c("APP","PSEN1","PSEN2","APOE");
t<-getGene(dGene,fields=c("symbol"))
t$symbol


# ------ capture a table of side-effects for tramadol in paper -------
stuff <-db_med_effects[grepl("tramadol", db_med_effects$V4),]
xtable(stuff[40:50,])

# --------------------------- compute statistics for each drug PPI network --------------------------------

setwd("C:/R-files/reposition") # need to point where STITCH datafiles are.

temp <-  list.files(pattern="*.txt") # list of all .txt files now in temp
allstats <- data.frame()

for (i in 1:length(temp)){
  mydata <- read.delim(temp[i], header=TRUE,sep='\t')
  mydata <- mydata[,1:2]
  mydata <-graph.data.frame(mydata,directed=FALSE)
  in1 <- get_gstatistics_long(mydata)
  allstats <- rbind(allstats,in1[nrow(in1),])
} 
  
rm(list = ls(pattern = glob2rx("*.txt"))) # get rid of useless file data from memory
table1 <- xtable(allstats,digits=c(3,3,0,0,3,0,3,3,3,3,3)) # latex table


# just a test to see if aripiprazole has hubness of 1.0
arip <- read.delim('C:\\R-files\\reposition\\aripiprazole2016.txt', header=TRUE,sep='\t')
arip <- arip[,1:2]
arip <-graph.data.frame(arip,directed=FALSE)
in2 <- get_gstatistics_long(arip)

# --- now do combined  network -----
giantnet <-graph.data.frame(giant,directed=FALSE)
g1 <- get_gstatistics_short(giantnet)
# sort structure according to hubness
g2<-g1[157:169,]
g3 <- g2[order(-g2$hubness),] 
table1 <- xtable(g3,digits=c(0,0,3,3))

nodesize=degree(giantnet)
nodelabel<-V(giantnet)$name

ad <- get.adjacency(giantnet)
nodecolor=character(ncol(ad))  # create a character for every column in adjaceny matrix
for (i in 1:length(nodecolor))
  if(nodesize[i] > 19){nodecolor[i] <- "tomato"} else {nodecolor[i]<-"darkgray"}

tkplot(giantnet,layout = layout.fruchterman.reingold,vertex.label = nodelabel,
       vertex.label.color= "black",vertex.size=nodesize, vertex.color=nodecolor,
       edge.arrow.size=0, edge.curved=FALSE)

# -------------- plot Anaums network -------------------
ssri <- read.delim('C:\\R-files\\reposition\\ssri.csv', header=TRUE,sep=',')
disease <- as.matrix(ssri[,c(2,4)])    # V2=drug V4 =targets
g <- graph.edgelist(disease,directed=FALSE)

nodesize=degree(g)*2

ad <- get.adjacency(g)
nodecolor=character(ncol(ad))  # create a character for every column in adjaceny matrix,
x <- 1:ncol(ad)

nodelabel<-V(g)$name
d1<-as.character(disease[,1])
d2<-as.character(disease[,2])

for ( i in 1:ncol(ad)){
  z<-(sapply(nodelabel[i],grep,d2))
  if(is.integer(z)){
    nodecolor[i]<-"pink"}
  y<-(sapply(nodelabel[i],grep,d1)) ## Protein Target
  if(is.integer(y)){   
    nodecolor[i]<-"lightblue"}  ## Drug                            
}

tkplot(g,layout = layout.fruchterman.reingold,vertex.label = nodelabel,
       vertex.label.color= "black",vertex.size=nodesize, vertex.color=nodecolor,
       edge.arrow.size=0, edge.curved=FALSE)


# ------------------------------------------------------------------------------

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
    #bpow<-mean(igraph::power_centrality(gt))
    gstats <- data.frame(modularity,avepath,nedges,nverts,transit,degree,diameter,
                  clos,between,dens,hubness )
  return(gstats)
}


get_gstatistics_long <- function(gt) {
  modularity<-modularity(gt, membership(cluster_walktrap(gt)))
  avepath <- average.path.length(gt)
  nedges <- ecount(gt)
  nverts <- vcount(gt)
  transit <- transitivity(gt)
  degree <- igraph::degree(gt)
  #diameter <- diameter(gt,weights=NA)
  #connect <- is.connected(gt)
  clos<- closeness(gt)
  between <- betweenness(gt,directed=FALSE)
  dens <-graph.density(gt)
  hubness <-hub_score(gt)$vector
  #author <-authority.score(gt)$vector
  #bpow<-mean(igraph::power_centrality(gt))
  gstats <- data.frame(modularity,avepath,nedges,nverts,transit,degree,clos,between,dens,hubness )
  return(gstats)
}

get_gstatistics_short <- function(gt) {
  #modularity<-modularity(gt, membership(cluster_walktrap(gt)))
  #avepath <- average.path.length(gt)
  #nedges <- ecount(gt)
  #nverts <- vcount(gt)
  #transit <- transitivity(gt)
  degree <- igraph::degree(gt)
  #diameter <- diameter(gt,weights=NA)
  #connect <- is.connected(gt)
  #clos<- closeness(gt)
  between <- betweenness(gt,directed=FALSE)
  #dens <-graph.density(gt)
  hubness <-hub_score(gt)$vector
  #author <-authority.score(gt)$vector
  #bpow<-mean(igraph::power_centrality(gt))
  gstats <- data.frame(degree,between,hubness )
  return(gstats)
}

## -------------- Disease Ontology using functions from the defunct DOSim package -------------------

# The diseases that are treated by the top ten drugs as per figure 5 in paper.

# DOID:10652 alzheimers
# DOID:14330 parkinsons
# DOID:5419 schizophrenia
# DOID:8646 psychosis
# DOID:1470 major depressive disorder
# DOID:3223 complex regional pain syndrome
# DOID:1826 epilepsy syndrome
# DOID:10808 gastric ulcer
# DOID:1508 candidiasis
# DOID:13938 amenorrhea

# create  a table or venn diagram to show common terms
diseases <- c("DOID:14330","DOID:10652","DOID:5419","DOID:8646","DOID:1470","DOID:3223","DOID:1826","DOID:10808","DOID:1508","DOID:13938")

terms <- getAncestors(diseases[2])

get(diseases[1], DOANCESTOR) # using library(DO.db)

fullterms <- getDoTerm(terms)

df <- data.frame(matrix(unlist(fullterms),  byrow=T),stringsAsFactors=FALSE)
a<-df[1:24,];b<-df[25:48,]#;c<-df[25:37,]
df2 <- data.frame(a,b)
  
xtable(df2)

g <- getDOGraph("DOID:14330","DOID:10652","DOID:5419","DOID:8646","DOID:1470",
                "DOID:3223","DOID:1826","DOID:10808","DOID:1508","DOID:13938")
plot(g)


venn.plot <- venn.diagram(terms[1:4], 
                          NULL, 
                          fill=c("red", "blue","green","yellow"), 
                          alpha=c(0.5,0.5,0.5,0.5), 
                          cex = 2, 
                          cat.fontface=2, 
                          margins =c(10,10),
                          cat.cex=2,
                          #main = "Venn Diagram showing shared side effects for donepezil,galantamine,rivastigmine",
                          category.names=c("parkinsons", "alzheimers","schizophrenia","psychosis"))
grid.draw(venn.plot)

# --------------- DOSim functions from defunct package ------------------

.First.lib<-function(lib, pkgname){
  #library.dynam(pkgname, pkgname, lib)
  initialize_DOSimEnv()
}


initialize_DOSimEnv <-
  function(){
    #print("initializing DOSim package ...")		
    data("DOSimEnv")
    #load("DOSimEnv.rda");
    #assign("DOSimEnv",DOSimEnv,envir=.GlobalEnv)  	
    #print("finished.")
  }


#-----------------------
getChildren <-
  function(dolist,verbose=TRUE){
    if(is.list(dolist)){
      dolist<-unique(unlist(dolist))
    }else{
      dolist<-unique(dolist)
    }
    if(verbose){
      print("Start to fetch the children")
    }
    if(!exists("DOSimEnv")) initialize_DOSimEnv()	
    children<-get("children",envir=DOSimEnv)
    res<-children[dolist[dolist %in% names(children)]]
    notmatch<-dolist[! dolist %in% names(children)]
    if(length(notmatch)>0){
      for(i in 1:length(notmatch)){
        res[[notmatch[i]]]<-NA
      }
    }
    res
  }
#-------------------
getAncestors <-
  function(dolist,verbose=TRUE){
    if(is.list(dolist)){
      dolist<-unique(unlist(dolist))
    }else{
      dolist<-unique(dolist)
    }
    if(verbose){
      print("Start to fetch the ancestors")
    }
    if(!exists("DOSimEnv")) initialize_DOSimEnv()	
    ancestor<-get("ancestor",envir=DOSimEnv)
    res<-ancestor[dolist[dolist %in% names(ancestor)]]
    notmatch<-dolist[! dolist %in% names(ancestor)]
    if(length(notmatch)>0){
      for(i in 1:length(notmatch)){
        res[[notmatch[i]]]<-NA
      }
    }
    res
  }


getDoTerm <-
  function(dolist){
    if(is.list(dolist)){
      dolist<-unique(unlist(dolist))
    }else{
      dolist<-unique(dolist)
    }
    
    if(!exists("DOSimEnv")) initialize_DOSimEnv()	
    doterm<-get("doterm",envir=DOSimEnv)
    res<-doterm[dolist[dolist %in% names(doterm)]]
    notmatch<-dolist[! dolist %in% names(doterm)]
    if(length(notmatch)>0){
      warning(paste(" ===>", length(notmatch), "of", length(dolist), "DOIDs not mapped to current disease ontology\n"))
    }
    res
  }

getDoAnno <-
  function(dolist){
    if(is.list(dolist)){
      dolist<-unique(unlist(dolist))
    }else{
      dolist<-unique(dolist)
    }
    if(!exists("DOSimEnv")) initialize_DOSimEnv()	
    doanno<-get("doanno",envir=DOSimEnv)
    res<-doanno[dolist[dolist %in% names(doanno)]]
    notmatch<-dolist[! dolist %in% names(doanno)]
    if(length(notmatch)>0){
      warning(paste(" ===>", length(notmatch), "of", length(dolist), "DOIDs not mapped to current disease ontology\n"))
    }
    res
    
  }

DOGraph <-
  function(term, env){
    oldEdges <- vector("list", length = 0)
    oldNodes <- vector("character", length = 0)
    newN <- term
    done <- FALSE
    while (!done) {
      newN <- newN[!(newN %in% oldNodes)]
      if (length(newN) == 0)
        done <- TRUE
      else {
        oldNodes <- c(oldNodes, newN)
        numE <- length(newN)
        #nedges <- AnnotationDbi::mget(newN, env = env, ifnotfound = NA)
        nedges <- mget(newN, env=env,ifnotfound= NA)
        nedges <- nedges[!is.na(nedges)]
        oldEdges <- c(oldEdges, nedges)
        if (length(nedges) > 0)
          newN <- sort(unique(unlist(nedges)))
        else newN <- NULL
      }
    }
    rE <- vector("list", length = length(oldNodes))
    names(rE) <- oldNodes
    rE[names(oldEdges)] <- oldEdges
    rE <- lapply(rE, function(x) match(x, oldNodes))
    names(oldNodes) = oldNodes
    return(new("graphNEL", nodes = oldNodes, edgeL = lapply(rE,function(x) list(edges = x)), edgemode = "directed"))
  }

getDOGraph <-
  function(term, prune=Inf){
    if(!require(graph))
      stop("Package graph is required for function getDOGraph")
    if(!exists("DOSimEnv")) initialize_DOSimEnv()
    ENV_Child2Parent<-get("ENV_Child2Parent",envir=DOSimEnv);
    
    G<-DOGraph(term,ENV_Child2Parent)		
    if(prune != Inf){
      dis = johnson.all.pairs.sp(G)		
      inc = unique(unlist(sapply(term, function(t) names(dis[t,])[dis[t,] < prune])))
      G = subGraph(nodes(G)[inc], G)
    }		
    G
  }








