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

setwd("/home/mgendre/Cluster/scripts/")

# Create directory
dir.create("../Simulation", showWarnings = FALSE)
dir.create("../Simulation/FullData", showWarnings = FALSE)
dir.create("../Simulation/MissingData", showWarnings = FALSE)
dir.create("../Simulation/Results", showWarnings = FALSE)
dir.create("../Simulation/Results/Replicates", showWarnings = FALSE)
dir.create("../Simulation/Results/Overall", showWarnings = FALSE)


#load datasets
files <- list.files("../csv/") #No parameter necessary now since you're in the proper directory
datasetList <- list()
for (i in 1:length(files)){
  datasetList[[i]] <- read.csv(paste0("../csv/",files[i]), header = T, sep = ";")
}

#remove .csv to file name
nameFiles <- gsub("\\.csv", "", files)

#Save session
LoadedPackages = sessionInfo()
saveRDS(LoadedPackages,"../Simulation/Results/SessionInfo.rds")

#Run scripts
#############
Seed <- seq(1,100,50) #number of replicates

#tree parameters
tree_arg <- list(Birth = 0.4, Death = 0.1, Ntaxa = 100)

for(data in 1:length(datasetList)){

  for(s in Seed){
    set.seed(s)
    SlurmID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
    
    #Simulate datasets. Output: SimulatedData
    print("Step 1: Simulate data")
    source(file = "./Simulation/SimData.R")
    nameSimulation <- file.path("../Simulation/FullData", sprintf("simulatedData%s_R%d",
                                                                  nameFiles[data], s+SlurmID))
    simulatedData <- simData(tree_arg, datasetList[[data]], save = nameSimulation)

    #Simulate NA in different proportions in full data. Output: NaNImputed
    print("Step 2: Simulate NA in data")
    source(file = "./Simulation/NaNImputation.R")
    nameNaNSimulation <- file.path("../Simulation/MissingData", sprintf("NaNData%s_R%d",
                                                                        nameFiles[data], s+SlurmID))
    partitions <- dataPartition(simulatedData)
    missingRates <- 1/3
    missTraits <- ncol(simulatedData$FinalData)
    NaNData <- NaNImputation(missingRates, partitions, simulatedData, 
                             missTraits, save = nameNaNSimulation)
    
    #Impute data where missing values (NA) are. Output: 
    print("Step 3: missing Data imputation")
    source(file = "./Imputation/imputeComparisonV2.R")
    nameImputation <- file.path("../Simulation/Results/Replicates", sprintf("Results%s_%d_R%d", 
                                                                            nameFiles[data], data, s+SlurmID))
    ImputationApproachesNames <- c("imputeDiscrete", "imputeContinuous", "imputeMICE", "imputeMissForest", "imputeKNN")
    variance_fractions <- c(0, 0.95)
    generateResults(ImputationApproachesNames, NaNData, simulatedData, 
                    variance_fractions, save = nameImputation)
  }
  
  #Overall results
  print("Step 4: overall results")
  filesInFolder <- list.files(path = "../Simulation/Results/Replicates", pattern = "*.RData")
  replicatesInFolder <- loadReplicates(data, filesInFolder, "../csv/")
  overallMean(replicatesInFolder, path = "../Simulation/Results/Overall/", 
              pathReplicates = "../Simulation/Results/Replicates/")

}



