### about these files
`reviewers_JBI_LoadData.R` - as it suggests loads in the tsv. files and creates the R data structures.

`reviewers_JBI_DrugList.R` - performs (so far) the majority of the processing to generate the candidates drugs for repurposing.

`reviewersJBIloadData.RData` - is the R environment from reviewers_JBI_LoadData.R avoids the need to run it everytime you power up R.

`venn-5.pdf` - example Venn Diagram for Alzheimer's drugs and their common overlapping side-effects. The diagram is limited to five groups (drugs) maximum and so the central value of 10 is in fact 8, because 9 drugs were used to generate the candidate drugs. Just be aware of this limitation when plotting Venn's if you have more than five drugs. Generally, the more drugs you conduct the search with, then the fewer side-effects these will have in common.
