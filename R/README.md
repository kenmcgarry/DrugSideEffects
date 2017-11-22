### about these files
As the code currently stands it is now generic enough to tackle any disease as long as the users know the [UMLS](https://www.nlm.nih.gov/research/umls/ "Unified Medical Language System") code for their disease of interest e.g. "C0002395" is the code for Alzheimer's disease.

`reposition_run.R` - Run this file first. You will need to know the UMLS code for your disease. Calls in the functions defined by reviewers_JBI_functions.R, assumes data structures are already loaded in. Then runs the functions and generates the candidate drugs. 

`reposition_LoadData.R` - as it suggests loads in the tsv. files and creates the R data structures.

`reposition_functions.R` - contains majority (so far) of function defintions for processing to generate the candidates drugs for repurposing.

`reposition_pathways.R` - obtains the target proteins for each drug and performs an analysis on their shared pathways.

`reposition_chemstructure.R` - obtains the similarities between the 77 drugs (8 conventional + 69 candidate).

`reposition_integrate.R` - combine all the scores and metrics from pathway analysis, GO and DO, on-targets etc etc.

`reposition_loadData.RData` - is the R environment from reviewers_JBI_LoadData.R avoids the need to run it afresh everytime you power up R.

`venn-5.pdf` - example Venn Diagram for Alzheimer's drugs and their common overlapping side-effects. The diagram is limited to five groups (drugs) maximum and so the central value of 10 is in fact 8, because 9 drugs were used to generate the candidate drugs. Just be aware of this limitation when plotting Venn's if you have more than five drugs. Generally, the more drugs you conduct the search with, then the fewer side-effects these will have in common.

### completed
working functions can now handle:
+ **Pathway analysis** with KEGG and Reactome for groups of gene/proteins

+ **GO analysis** for each gene/protein

+ **DO analysis** for each candidate drugs associated disease, using DOSE we can deduced how similar the diseases are.

+ **Integrate** I used a variation on the Jaccard coefficient to integrate the disparate data in a coherent way. This as stated above needs tinkering with!

+ **Comparison** with existing methods - I compared our system with PREDICT (Gottlieb,2011) and the systems by Zhang (2014) and Wang 2013). Results are not quite the same but a number of drugs have been consistently identified as candidates by the four systems. Previously, I have evaluated the usefulness/accuracy of the system by checking if my top ten drugs have appeared in the literature as likely candidates for repurposing. The results so far support our method.

+ **Other diseases** At the moment I have concentrated on Alzhiemers disease, the system needs to be tested with others. I tackled the diseases investigated by Zhang and Wang.

#### `Ken McGarry, 22nd November 2017.`
