# pubchem.r
# Ken McGarry 27/5/2016
# Use cygwin to download this file as web browsers like firefox are not up to 8GB downloads!
# http://biocluster.ucr.edu/~tbackman/bioassayR/pubchem_protein_only.sqlite
# source("http://bioconductor.org/biocLite.R")
# biocLite("dplyr")

library(ggplot2)
library(dplyr)
library(tidyr)
library(igraph)
library(biomaRt)
library("org.Hs.eg.db") # for ncbi to gene conversion 

library(paxtoolsr)
library(rJava)

#library(rcdk)
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
fpset <- desc2fp(apset, descnames=1024, type="FPset")

## Subsetting by cid slot using drugbank ids should work now consistently
sdfset["DB00472"]
apset["DB00472"]

fpdrugs <- fpset[names(drugs)]
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

print.xtable(results,label='comp',caption = 'Similarity measures for the 10 drugs based on 1024-length, atom-pair fingerprints')

# Creates the similarity score matrix and clusters them.
simMA <- sapply(cid(fpdrugs), function(x) fpSim(x=fpdrugs[x], fpdrugs, sorted=FALSE)) 
colnames(simMA)<-as.character(drugs[cid(fpdrugs)])
rownames(simMA)<-as.character(drugs[cid(fpdrugs)])

# version 1 of plotting dendro
library(ape)
library("dendextend")
library("dendextendRcpp")  

# using piping to get the dend
par(mar=c(5, 5, 5, 4))
dend <- simMA[,-6] %>% dist %>% hclust %>% as.dendrogram
# plot + color the dend's branches before, based on 3 clusters:
#dend %>% color_branches(k=4) %>% plot(horiz=TRUE,lwd = 3)
dend %>% plot(horiz=TRUE,lwd = 4,col = "blue")
# add horiz rect
dend %>% rect.dendrogram(k=4,horiz=TRUE,border = 2, lty = 3, lwd = 3)
# add horiz (well, vertical) line:
abline(v = heights_per_k.dendrogram(dend)["3"] + .6, lwd = 2, lty = 3, col = "blue")

# version 2 of plotting dendro
par(mar=c(5, 5, 5, 4))
plot(as.phylo(hc), cex = 0.9, label.offset = 0.01)
hc <- hclust(as.dist(1-simMA), method="single")
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=3), horiz=T) 
rect.hclust(hc,k=4,border="red")

# version 3 of plotting dendro
library(cluster)
library(ggplot2)
library(ggdendro)     # for dendro_data(...)
dst   <- daisy(simMA, metric = c("gower"), stand = FALSE)
hca   <- hclust(dst, method = "average")
k     <- 5
clust <- cutree(hca,k=k)  # k clusters

dendr    <- dendro_data(hca, type="rectangle") # convert for ggplot
clust.df <- data.frame(label=rownames(simMA), cluster=factor(clust))
dendr[["labels"]]   <- merge(dendr[["labels"]],clust.df, by="label")
rect <- aggregate(x~cluster,label(dendr),range)
rect <- data.frame(rect$cluster,rect$x)
ymax <- mean(hca$height[length(hca$height)-((k-2):(k-1))])

ggplot() +   geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=cluster), size=5) +
  geom_rect(data=rect, aes(xmin=X1-.3, xmax=X2+.3, ymin=0, ymax=ymax), color="red", fill=NA)+
  coord_flip() + geom_hline(yintercept=0.21, color="blue") + scale_y_reverse(expand=c(0.2, 0))# + theme_dendro()




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





