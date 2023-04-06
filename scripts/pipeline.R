#load package
if(!require(devtools)){
  install.packages("devtools")
}

if(!require(visdat)){
  install.packages("visdat")
}

if(!require(TDIP)){
  devtools::install_github("Matgend/TDIP")
}

library(devtools)
library(TDIP)

args = commandArgs(trailingOnly=TRUE)

print(paste(args[1], "start"))

args = commandArgs(trailingOnly=TRUE)

print(length(args))
# test if there is at least one argument: if not, return an error
if (length(args) != 3) {
  stop("10 arguments must be supplied", call.=FALSE)
}

missingRates <- as.numeric(args[1])
gain <- as.logical(args[2])
replicateNbr <- as.numeric(args[3])

#tree parameters
tree_arg <- list(Birth = 0.4, Death = 0.1, Ntaxa = 100)
varfrac <- as.numeric(0.95)

#strategies can be "NP", "P" or "2-step"
strategies <- c("NP", "P", "2-step")

ImputationApproachesNames <- c("pi_categorical_traits", "pi_continuous_traits", "mice_phylo", "missForest_phylo",
                                "kNN_phylo")


if(gain){
  ImputationApproachesNames <- c("pi_categorical_traits", "pi_continuous_traits", "mice_phylo", "missForest_phylo",
                                 "kNN_phylo", "gain_phylo")
}

# Create directory
directoryName <- paste0(gsub("\\.", "", round(missingRates, 2)))
dir.create("../Simulations", showWarnings = FALSE)
dir.create(paste0("../Simulations/", directoryName), showWarnings = FALSE)
dir.create(paste0("../Simulations/", directoryName, "/FullData"), showWarnings = FALSE)
dir.create(paste0("../Simulations/", directoryName, "/MissingData"), showWarnings = FALSE)
dir.create(paste0("../Simulations/", directoryName, "/Results"), showWarnings = FALSE)
dir.create(paste0("../Simulations/", directoryName, "/Results/Replicates"), showWarnings = FALSE)
dir.create(paste0("../Simulations/", directoryName, "/Results/Overall"), showWarnings = FALSE)

#load datasets
files <- list.files("../csv/") #No parameter necessary now since you're in the proper directory
datasetList <- list()
for (i in 1:length(files)){
  datasetList[[i]] <- read.csv(paste0("../csv/",files[i]), header = T, sep = ";")
}

#conda environement
#use_condaenv(condaenv = "tfKeras")

#remove .csv to file name
nameFiles <- gsub("\\.csv", "", files)

#Save session
LoadedPackages = sessionInfo()
saveRDS(LoadedPackages,paste0("../Simulations/",directoryName,"/Results/SessionInfo.rds"))

#Run scripts
############
for(data in 1:(length(datasetList))){
  set.seed(replicateNbr)
  #Simulate datasets. Output: SimulatedData
  print("Step 1: Simulate data")
  nameSimulation <- file.path("../Simulations/", directoryName, "/FullData", sprintf("simulatedData%s_R%d_%s",
                                                                nameFiles[data], replicateNbr,
                                                                paste(as.character(round(missingRates, 2)),
                                                                collapse = "-")))

  print(nameFiles[data])
  simulatedData <- TDIP::data_simulator(tree_arg,
                                        datasetList[[data]],
                                        save = nameSimulation)


  #Simulate NA in different proportions in full data. Output: NaNImputed
  print("Step 2: Simulate NA in data")
  nameNaNSimulation <- file.path("../Simulations/", directoryName, "/MissingData", sprintf("NaNData%s_R%d_%s",
                                                                      nameFiles[data], replicateNbr,
                                                                      paste(as.character(round(missingRates, 2)),
                                                                      collapse = "-")))
  print(nameNaNSimulation)


  print(ncol(simulatedData$FinalData))
  print(class(ncol(simulatedData$FinalData)))

  NaNData <- TDIP::na_insertion(missingRates,
                                dataset = simulatedData$FinalData,
                                missingTraits = ncol(simulatedData$FinalData),
                                MARTraits = 1,
                                MARctrlTraits = sample(2:ncol(simulatedData$FinalData) ,1), #in case of correlated traits, MARctrlTraits = NULL provided
                                traitsNoNA = NULL,
                                tree = simulatedData$TreeList$`0`,
                                save = nameNaNSimulation)


  #Impute data where missing values (NA) are.
  print("Step 3: missing Data imputation")
  nameImputation <- file.path("../Simulations/", directoryName, "/Results/Replicates", sprintf("Results%s_%d_R%d_%s",
                                                                          nameFiles[data], data, replicateNbr,
                                                                          paste(as.character(round(missingRates, 2)),
                                                                          collapse = "-")))


  mecList <- vector("list", length(NaNData$DataNaN))
  names(mecList) <- names(NaNData$DataNaN)
  errorList <- vector("list", length(NaNData$DataNaN))
  names(errorList) <- names(NaNData$DataNaN)

  namesMethods <- c("MICE", "MissForest", "KNN") #select the method that we want to ensemble
  for(mC in 1:length(NaNData$DataNaN)){

    for(d in 1:length(NaNData$DataNaN[[mC]])){

      data <- NaNData$DataNaN[[mC]][[d]]
      

      imputedData <- TDIP::missing_data_imputation(ImputationApproachesNames,
                                                   data,
                                                   tree = simulatedData$TreeList$`0`,
                                                   strategies,
                                                   maxit = 5,
                                                   nbrMI = 1,
                                                   k = 3,
                                                   numFun = laeken::weightedMedian,
                                                   catFun = VIM::maxCat,
                                                   mincor = NULL,
                                                   varfrac = varfrac,
                                                   save = NULL)
      

      imputedName <- paste0("imputed_", names(NaNData$DataNaN[[mC]])[d])


      #hard voting
      strategies_HV <- c("NP", "mixed")
      for(s in 1:length(strategies_HV)){

        if(strategies_HV[s] == "mixed"){
          namesToselect <- c("MissForest_NP", "KNN_2-step", "MICE_2-step")
        }
        
        else{
          namesToselect <- paste(namesMethods, strategies_HV[s], sep = "_")
        }
        
        namesImputedData <- names(imputedData)
        datasets <- imputedData[which(namesImputedData %in% namesToselect)]

        print(paste("length datasets:", length(datasets)))
        print(names(datasets))
        HVData <- list(TDIP::hard_voting(datasets))
        names(HVData) <- paste0("HV.ML_", strategies_HV[s])
        imputedData <- c(imputedData, HVData)
      }

      errorData <- vector("list", length(imputedData))
      errorName <- c()

      for(impD in 1:length(imputedData)){

        imputationApproachName <- names(imputedData)[impD]

        error <- TDIP::imputation_error(imputedData[[impD]],
                                        simulatedData$FinalData,
                                        data,
                                        imputationApproachName,
                                        simulatedData$Dataframe)

        errorN <- paste0("error_", imputationApproachName)

        errorData[[impD]] <- error
        errorName <- c(errorName, errorN)

      }

      names(errorData) <- errorName

      mecList[[mC]][[d]] <- imputedData
      names(mecList[[mC]])[d] <- imputedName

      errorList[[mC]][[d]] <- errorData
      names(errorList[[mC]])[d] <- paste0("error_", names(NaNData$DataNaN[[mC]])[d])

    }
  }

  #save imputed object
  saveData <- list(errorData = errorList, imputedData = mecList, missingData = NaNData, simulatedData = simulatedData)

  save(saveData, file = paste0(nameImputation, ".RData"))

}

print(paste(args[1], "done"))
