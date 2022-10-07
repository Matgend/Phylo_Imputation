# Phylo_Imputation

This repo provides the scripts used to generate all the simulations and the different graphs and tables present in the thesis. 
All the functions used in `pipeline.R` come from the R package `TDIP`

The command line to run the pipeline in order to reproduce the results obtained is the following: 
```{r setup}
Rscript --vanilla pipeline.R 0.05 TRUE 1
```
The arguments:
* `Missing rate`: 0.05
* `GAIN or not`: TRUE or FALSE
* `Replicate`: the index of the replicate

The script automatically installs TDIP, which requires a Python installation with defined packages. Please check the repo of the R package for more information.

To generate the graphs, the user can run this command line 
```{r setup}
Rscript --vanilla tables_plots.R
```
