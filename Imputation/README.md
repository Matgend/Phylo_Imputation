
## List of scripts

`impute_MissForest.R`
This scripts matches a tree and a trait files and performs imputations using the MissForest package. It uses eigenvalues as summary for the underlying phylogenetic trees using the `PVR` library. 

`impute_mvMORPH.R`
This scripts imputes missing continuous values using the mvMORPH package.

`impute_discrTrait_phytools.R`
This scripts imputes missing states of a single discrete trait using the phytools package.

`imputationApproaches.R`
This script contains these imputation approaches: phytools, MICE, missForest, kNN, phylopars 
