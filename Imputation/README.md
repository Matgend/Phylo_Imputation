
## List of scripts

`impute_MissForest.R`
This scripts matches a tree and a trait files and performs imputations using the MissForest package. It uses eigenvalues as summary for the underlying phylogenetic trees using the `PVR` library. 

`impute_mvMORPH.R`
This scripts imputes missing continuous values using the mvMORPH package.

`impute_discrTrait.R`
This scripts imputes missing states of a single discrete trait using the R package phytools or corHMM.

`imputationApproaches.R`
This script contains these imputation approaches: Rphylopars, corrHMM, MICE, missForest and kNN.

`imputeComparisonV2.R`
This script contain the functions to generate the imputed values and to calculate the imputation accuracy.


