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
  * class: non_eq_nominal (nominal), ordinal or continuous
  * model: model of evolution: for the MK model, you have to provide the type of rate matrices (ER, SYM or ARD), for continuous the model is BM1 or OU1 and "Manual" is a way to generate trait that are correlated to a specific trait being not simulated according an evolutionary model.
  * states: the number of states in the trait, for continuous trait you have to write 1. 
  * correlation: integer defining if some trait are correlated with other.
  * uncorr_traits: number of trait which are uncorrelated among the nbr_traits feature
  * fraction_uncorr_traits: percentage of trait which are uncorrelated among the nbr_traits feature
  * lambda: $\lambda$ value
  * kappa: $\kappa$ value
  * highCor: correlation rate between the manual trait and the trait of interest.

In the csv folder, you can find some example of CSV files. 

* parameters to simulate the phylogenetic tree under a birth death model




