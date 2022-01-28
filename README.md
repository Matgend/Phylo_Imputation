# Phylo_Imputation

`main.R` This script runs the frameworks, have to define the number of replicates, the missingness and the variance fraction for phylogenetic informations


Why PhyloMCAR doesn't work:\

Use the function `getTipsNA` which takes as arguments tree (phylogeny)  and minTips (number of species for which NAs should be created). And return the name of the tips selected:\
if(minTips == 2): \
  select a taxa of 2 tips or more (2+-3)\
  etc..\
because we simulate a different tree for each simulated data, sometimes there is no taxa with 2 tips and therefore I can't calculate the mean of the replicates because the taxa selected haven't the same size. 
