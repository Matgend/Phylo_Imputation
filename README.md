# Phylo_Imputation

`main.R` This script runs the frameworks, have to define the number of replicates, the missingness and the variance fraction for phylogenetic informations

`main.sh` This script runs the frameworks on job arrays(slurm) 

`utils.R` This script contains the error calculation function and function for results representation. 

##########################

## Pipeline

To run the same pipeline that applied for the study, you have to run the following command:
`Rscript --vanilla mainNew.R 0.05 0.4 0.1 100 0.95 NP P 2-step`

 * filename: `mainNew.R`
 * missing rate: `0.05`
 * birth rate : `0.4`
 * death rate: `0.1`
 * number of taxa: `100`
 * amount of variance provide by the eigenvectors: `0.95`
 * strategies: `NP`(no phylogeny), `P` (phylogeny), `2-step`



