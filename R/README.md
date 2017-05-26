### about these files
As the code currently stands it is now generic enough to tackle any disease as long as the users know the [UMLS](https://www.nlm.nih.gov/research/umls/ "Unified Medical Language System") code for their disease of interest e.g. "C0002395" is the code for Alzheimer's disease.

`reviewers_JBI_run.R` - Run this file first. You will need to knmow the UMLS code for your disease. Calls in the functions defined by reviewers_JBI_functions.R, assumes data structures are already loaded in. Then runs the functions and generates the candidate drugs. 

`reviewers_JBI_LoadData.R` - as it suggests loads in the tsv. files and creates the R data structures.

`reviewers_JBI_functions.R` - contains majority (so far) of function defintions for processing to generate the candidates drugs for repurposing.

`reviewers_JBI_pathways.R` - obtains the target proteins for each drug and performs an analysis on their shared pathways.

```
setwd("C:/R-files/sider")    # point to where my code lives
load("reviewersJBIloadData.RData") # load in the data structures made by reviewers_JBI_LoadData.R

source("reviewers_JBI_DrugList.R")  # load in the functions required for finding lists of drugs and side-effects

```

`reviewersJBIloadData.RData` - is the R environment from reviewers_JBI_LoadData.R avoids the need to run it afresh everytime you power up R.

`venn-5.pdf` - example Venn Diagram for Alzheimer's drugs and their common overlapping side-effects. The diagram is limited to five groups (drugs) maximum and so the central value of 10 is in fact 8, because 9 drugs were used to generate the candidate drugs. Just be aware of this limitation when plotting Venn's if you have more than five drugs. Generally, the more drugs you conduct the search with, then the fewer side-effects these will have in common.

### whats next?
I need to write a number functions to handle:

+ **Direct R communications** with STITCH database - so far I have manually created lists of protein/drug interactions using their browser. A more generic and easier option for users would be to interrogate STITCH via an API.

+ **Infer** a given side-effect to the biochemical pathways involved. Requires analysis of the drugs on-targets and especially off-targets.

+ **Integrate** the various the sources of the disparate data in a coherent way. So far the system has been a series of isolated analyses, each important in its own right but lacking an overarching synthesis. The usual way in the literature is to generate a 'fingerprint' that integrates the presence or absence of a parameter using '1s' or '0s', a method that typically generates huge arrays.

+ **Comparison** with existing methods - might be problematic as not all the software is freely available. Maybe at least two will be necessary for an evaluation. So far I have evaluated the usefulness/accuracy of the system by checking if my top ten drugs have appeared in the literature as likely candidates for repurposing. The results so far support our method.

### completed
working functions can now handle:
+ **Pathway analysis** with KEGG and Reactome for groups of gene/proteins

+ **GO analysis** with for each gene/protein

#### `Ken McGarry, Friday, 26th May 2017.`
