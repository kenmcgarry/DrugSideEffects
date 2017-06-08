# reviewers_JBI_pathways.R
# 1. analysis of pathway information based on membership proteinS
# 2. also ontological similarity analysis

library(igraph)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(xtable)
library(ReactomePA)
library(clusterProfiler)
library(GOSemSim)
library(DOSE)
library(AnnotationHub)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
#library(biomaRt)

setwd("C:/R-files/sider")    # point to where my code and data files are
load("reviewersJBIloadData.RData") # load in the data structures made by reviewers_JBI_LoadData.R
load("reviewers_candidates.RData") # load in the data structures (drug_list,se_list,repos_list,mycandidates)
source("reviewers_JBI_functions.R")  # read in the functions required for finding lists of drugs and side-effects

# get the ontarget proteins for each drug
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

x <- enrichPathway(gene=candidate_targets$entrezid,pvalueCutoff=0.05, readable=TRUE)
head(as.data.frame(x))
xresult <- as.data.frame(x)
xresult <- as_tibble(x)
print.xtable(xtable(xresult[,c(1:5,7)])) # displays Pathway table for paper.

enrichMap(x, layout=igraph::layout.fruchterman.reingold, fixed=FALSE, vertex.label.cex = 1) # interesting plot
barplot(x, showCategory=15)   # another useful plot

# gene ontology enrichment
hsGO <- godata('org.Hs.eg.db', ont="BP")
genes_to_test <- na.omit(candidate_targets$entrezid) # get rid of NA
genes_to_test <- sort(genes_to_test,decreasing = TRUE)

ggo <- groupGO(gene= genes_to_test, OrgDb = org.Hs.eg.db, ont = "BP", level= 3, readable = TRUE)
goresult <- as.data.frame(ggo)
goresult <- as_tibble(ggo)
print.xtable(xtable(goresult[,1:4])) # displays GO table for paper.

# KEGG enrichment
mkk <- enrichMKEGG(gene = genes_to_test,organism = 'hsa')
head(mkk)

# probably more useful if it works
d <- semData('org.Hs.eg.db', ont="MF")
go1 = c("GO:0004022","GO:0004024","GO:0004174")
go2 = c("GO:0009055","GO:0005515")
mgoSim(go1, go2, semData=hsGO, measure="Wang", combine=NULL)
mgoSim(go1, go2, semData=d, measure="Wang")

# Comparison of diseases using the Disease Ontolgy (DO) http://disease-ontology.org/
# http://bioconductor.org/packages/release/bioc/vignettes/DOSE/inst/doc/DOSE.html
# Alzheimer's disease=DOID:10652 	
# Parkinsons=DOID:14330
# Tachycardia = DOID:0060674 	
# mental disorder= DOID:1561
# depression=DOID:1596
# Fatigue= 	DOID:8544
# bipolar disorder=DOID:3312
# Anxiety = DOID:2030
# Ataxia=DOID:0050951
# Migraine=DOID:6364 	
# Epilepsy=DOID:1826
# Restless legs syndrome=DOID:0050425
# Hypotension=DOID:4723
DOnames <- c("Alzheimers","Parkinsons","Mental disorder","Depression","Fatigue","Bipolar disorder",
             "Anxiety","Migraine","Epilepsy","Restless legs","Hypotension")

listDO1<-c("DOID:10652","DOID:14330","DOID:1561","DOID:1596","DOID:8544",
           "DOID:3312","DOID:2030","DOID:6364","DOID:1826","DOID:0050425","DOID:4723")
s<-doSim(listDO1,listDO1,measure="Wang")
rownames(s) <- DOnames # row and column names are the same
colnames(s) <- DOnames
simplot(s,color.low="white",color.high="lightgreen",labs=TRUE,digits=2,labs.size=5,
        font.size=14, xlab="",ylab="")

# similarity coefficients for gene sets
g1 <- c("84842", "2524", "10590", "3070", "91746")
g2 <- c("84289", "6045", "56999", "9869")
DOSE::geneSim(g1, g2, measure="Wang", combine="BMA")



