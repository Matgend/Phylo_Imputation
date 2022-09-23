# Phylo_Imputation

`main.R` This script runs the frameworks, have to define the number of replicates, the missingness and the variance fraction for phylogenetic informations

`main.sh` This script runs the frameworks on job arrays(slurm) 

`utils.R` This script contains the error calculation function and function for results representation. 

##########################

# Usage

In PDIMP, there are several way to play with. You can replicate the analysis done in the (Manuscript), use the pipeline to gap missing value from an empirical data or just pick some functions for other purposes. Below, the inputs and the command to run to replicate the analysis are presented as well as the command to fill the gaps in an empirical data. 

## Simulation

The inputs to provide to simulate data are the following:

* a CSV file providing all the information for the simulation of the data. The CSV file must contain 10 columns:
  * nbr_traits: number of traits
  * class: non_eq_nominal (nominal), interval, ordinal or continuous
  * model: model of evolution: for the MK model, you have to provide the type of rate matrices (ER, SYM or ARD), for continuous the model is BM1 or OU1 and "Manual" is a way to generate trait that are correlated to a specific trait being not simulated according an evolutionary model.
  * states: the number of states in the trait, for continuous trait you have to write 1. 
  * correlation: integer defining if some trait are correlated with other.
  * uncorr_traits: number of traits which are uncorrelated among the nbr_traits feature
  * fraction_uncorr_traits: percentage of trait which are uncorrelated among the nbr_traits feature
  * lambda: $\lambda$ value
  * kappa: $\kappa$ value
  * highCor: correlation rate between the "manual" traits and the trait of interest.

* parameters to simulate the phylogenetic tree under a birth death model

In the csv folder, you can find some examples of CSV files. 

## Missing data insertion

In case you have an empirical data and you want to insert artificial missing data according to these fourth missing mechanisms (MCAR, MAR, MNAR and PhyloNa), you can use the function **NaNImputationEmp()**. In case you want to generate missing values according a specific missing mechanisms in this case you can use the appropriate functions **myMCAR**, **myMAR**, **myMNAR** or **PhyloNa()**.

## Missing values imputation

Through this package, this imputation approaches are provided:
 * corHMM
 * Rphylopars
 * kNN
 * missForest
 * MICE 
 * GAIN
To apply the phylogenetic imputation methods **corHMM** and **Rphylopars**, a phylogenetic tree must by provided. GAIN is a deep learning method scripted in python which requires a python environment with tensorflow installed. 

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



