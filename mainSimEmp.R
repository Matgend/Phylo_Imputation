#install.packages("mvMORPH")
#install.packages("phytools")
#install.packages("Matrix")
#install.packages("castor")
#install.packages("geiger")
#install.packages("tidyverse")
#install.packages("missMethods")
#install.packages("phangorn")
#install.packages("fastDummies")
#install.packages("VIM")
#install.packages("laeken")
#install.packages("mice")
#install.packages("ape")
#install.packages("PVR")
#install.packages("missForest")
#install.packages("tibble")
#install.packages("Rphylopars")
#install.packages("corHMM")

# Attach packages
require(mvMORPH)
require(phytools)
require(Matrix)
require(castor)
require(geiger)
require(tidyverse)
require(missMethods)
require(phangorn)
require(fastDummies)
require(VIM)
require(laeken) #for weigthedMean
require(mice)
require(ape)
require(PVR)
require(missForest)
require(tibble)
require(Rphylopars)
require(corHMM)
require(reticulate)
require(faux)
library(miscset)

setwd("/home/mgendre/Cluster/scripts/")
#setwd("C:/Users/Matthieu/Documents/UNIFR/Master_thesis/Scripts/Phylo_ImputationLocal/")

args = commandArgs(trailingOnly=TRUE)

print(paste(args[1], "start"))

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("One argument must be supplied", call.=FALSE)
} else if (length(args) != 3) {
  stop("3 arguments should be supplied")
}

missingRates <- as.numeric(args[1])
typeOfDataset <- as.character(args[2]) #if simulation step one done, if not start from step 2.
empTree <- as.character(args[3]) #load phyologeny


directoryName <- paste0(gsub("\\.", "", round(missingRates, 2)), "_10T")

# Create directory
dir.create("../Simulation", showWarnings = FALSE)
dir.create(paste0("../Simulation/", directoryName), showWarnings = FALSE)
dir.create(paste0("../Simulation/", directoryName, "/FullData"), showWarnings = FALSE)
dir.create(paste0("../Simulation/", directoryName, "/MissingData"), showWarnings = FALSE)
dir.create(paste0("../Simulation/", directoryName, "/Results"), showWarnings = FALSE)
dir.create(paste0("../Simulation/", directoryName, "/Results/Replicates"), showWarnings = FALSE)
dir.create(paste0("../Simulation/", directoryName, "/Results/Overall"), showWarnings = FALSE)


#load datasets
files <- list.files("../csv/") #No parameter necessary now since you're in the proper directory
datasetList <- list()
for (i in 1:length(files)){
  datasetList[[i]] <- read.csv(paste0("../csv/",files[i]), header = T, sep = ";")
}

#conda environement
#use_condaenv(condaenv = "tfKeras")

#load functions
source(file = "./Simulation/SimData.R")
source(file = "./Simulation/NaNImputation.R")
source(file = "./Imputation/imputeComparisonV3.R")
source(file = "utils.R")

#remove .csv to file name
nameFiles <- gsub("\\.csv", "", files)

#Save session
LoadedPackages = sessionInfo()
saveRDS(LoadedPackages,paste0("../Simulation/",directoryName,"/Results/SessionInfo.rds"))

#tree parameters
tree_arg <- list(Birth = 0.4, Death = 0.1, Ntaxa = 100)

#traits selected for MAR values generation
MARTraits <- 1 #index
MARctrlTraits <- 6 #index

#Run scripts
#############
SlurmID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

for(data in 1:(length(datasetList))){
  set.seed(SlurmID)
  
  nameSimulation <- file.path("../Simulation/", directoryName, "/FullData", sprintf("simulatedData%s_R%d_%s",
                                                                                    nameFiles[data], SlurmID, 
                                                                                    paste(as.character(round(missingRates, 2)),
                                                                                          collapse = "-")))
  if(typeOfDataset == "simulation"){
    #Simulate datasets. Output: SimulatedData
    print("Step 1: Simulate data")
    
    simulatedData <- simData(tree_arg, datasetList[[data]], save = nameSimulation)
    
    #resimulate data in case a trait contains only 1 trait
    if(length(simulatedData) == 1){
      dataOK <- FALSE
      while(!dataOK){
        simulatedData <- simData(tree_arg, datasetList[[data]], save = nameSimulation)
        if(length(simulatedData) != 1){
          dataOK <- TRUE
        }
      }
    }
    
    #Simulate NA in different proportions in full data. Output: NaNImputed
    print("Step 2: Simulate NA in data")
    nameNaNSimulation <- file.path("../Simulation/", directoryName, "/MissingData", 
                                   sprintf("NaNData%s_R%d_%s", 
                                           nameFiles[data], 
                                           SlurmID, paste(as.character(round(missingRates, 2)), collapse = "-")))
    print(nameNaNSimulation)
    partitions <- dataPartition(simulatedData)
    missTraits <- ncol(simulatedData$FinalData)
    NaNData <- NaNImputation(missingRates, partitions, simulatedData, 
                             missTraits, save = nameNaNSimulation)
    #resimulate data in case a in Manual design the discrete trait doesn't contain NA
    if(length(NaNData) == 1){
      NaNdataOK <- FALSE
      while(!dataOK){
        NaNData <- NaNImputation(missingRates, partitions, simulatedData, 
                                 missTraits, save = nameNaNSimulation)
        if(length(NaNData) != 1){
          NaNdataOK <- TRUE
        }
      }
    }
    variance_fractions <- c(0, 0.95, 2)
  }
  
  else{
    print("Step 2: Simulate NA in data")
    nameNaNSimulation <- file.path("../Simulation/", directoryName, "/MissingData", 
                                   sprintf("NaNData%s_R%d_%s", 
                                           nameFiles[data], 
                                           SlurmID, paste(as.character(round(missingRates, 2)), collapse = "-")))
    print(nameNaNSimulation)
    #partitions <- dataPartition(simulatedData)
    
    variance_fractions <- c(0)
      
    #convert csv with empirical data
  
    if(!is.null(empTree)){ #check if possible to do this
      empTree <- ape::read.nexus(empTree) #maybe has to be changed
      variance_fractions <- c(0, 0.95, 2)
    }
    #empTree <- load(empTree)
    dataset <- passInList(datasetList[[data]], empTree = empTree, save = nameSimulation)
    
    missTraits <- ncol(dataset$FinalData)
    
    NaNData <- NaNImputationEmp(missingRate, dataset$FinalData,
                                missingTraits, MARTraits, MARctrlTraits, 
                                tree = dataset$TreeList$`0`, save = nameNaNSimulation)
    
    #resimulate data in case a in Manual design the discrete trait doesn't contain NA
    if(length(NaNData) == 1){
      NaNdataOK <- FALSE
      while(!dataOK){
        NaNData <- NaNImputation(missingRates, partitions, simulatedData, 
                                 missTraits, save = nameNaNSimulation)
        if(length(NaNData) != 1){
          NaNdataOK <- TRUE
        }
      }
    }
  }

  #Impute data where missing values (NA) are. Output: 
  print("Step 3: missing Data imputation")
  nameImputation <- file.path("../Simulation/", directoryName, "/Results/Replicates", sprintf("Results%s_%d_R%d_%s", 
                                                                          nameFiles[data], data, SlurmID, 
                                                                          paste(as.character(round(missingRates, 2)),
                                                                          collapse = "-")))

  ImputationApproachesNames <- c("imputeDiscrete", "imputeContinuous", "imputeMICE", "imputeMissForest", 
                                 "imputeKNN", "gainR")   

  
  generateResults(ImputationApproachesNames, NaNData, simulatedData, 
            		variance_fractions, save = nameImputation, addHint = TRUE)
  
}
print(paste(args, "done"))