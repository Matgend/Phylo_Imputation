library(tidyverse)
#library(Metrics)


# load simulated data
setwd("C:/Users/Matthieu/Documents/UNIFR/Master_thesis/Scripts/Phylo_Imputation/")
load("C:/Users/Matthieu/Documents/UNIFR/Master_thesis/Scripts/simulatedData.RData")

# load NaN imputed data
load("DataNaN.RData")

#importe imputation functions
source("Imputation/imputationApproaches.R")

ImputationApproachesNames <- c("imputeDiscrete", "imputeContinuous", "imputeMICE", "imputeMissForest", "imputeKNN")
MixedImputationApproachesNames <- c("imputeMICE", "imputeMissForest", "imputeKNN")

missingData <- NaNImputed$DataNaN$`1`$MCAR$`MCAR/CorrContinuousTraits/1/0.05`
missingData
imputedData <- imputeContinousTraits(missingData)
imputedData
missingData <- NaNImputed$DataNaN$`1`$MCAR$`MCAR/CorrContinuousTraits/4/0.05`
imputedTest <- imputeKNN(missingData, k = 2, weightedMedian, maxCat, variance_fraction = 0.8)

TrueData <- Data$FinalData[, names(missingData), drop = FALSE]
TrueData
imputedData <- imputedTest 
trueData <- TrueData




imputationError <- function(imputedData, trueData, missingData){
  
  #get the ordinal trait reference
  ordinalTraits <- which(Data$dataframe$class == "ordinal")

  errors <- c()
  traitNames <- c() 
  for (c in 1:ncol(missingData)){
    
    #know is NaNs in the columns(trait)
    NaNRowIndex <- which(is.na(missingData[,c]))
    
    if(length(NaNRowIndex != 0)){
    
      traitNames <- c(traitNames, names(trueData[c]))
      missingValues <- missingData[NaNRowIndex, c]
      trueValues <- trueData[NaNRowIndex, c]
      imputedValues <- imputedData[NaNRowIndex, c]
  
      #in case continuous data
      if((grep("F.", names(missingData)[c])) == 1){
        
        #rmse
        error <- sqrt(mean((imputedValues - trueValues)^2))
  
      }
  
      #in case ordinal trait
      else if(grep(paste0("/", ordinalTraits), names(missingData)[i]) == 1){
      
        #imputation error for ordinal traits (absolute error)
        error <- mean(abs((imputedValues - trueValues) / trueValues))
      }
  
      #in case discrete data
      else{
        #imputation error for discrete traits
        error <- mean(length(setdiff(imputedValues, trueValues)))
      }
      errors <- c(errors, error)
    }
  }
  output <- data.frame(trait = traitNames, error = errors)

  #for (i in 1:length(NaNIndex)){
  #   
  #   if(length(NaNIndex[[i]] != 0)){
  #     
  #     for(r in NaNIndex[[i]]){
  #       
  #       #in case continuous data
  #       if(length(grep("F.", names(missingData)[i])) == 1){
  #         imputedContiValues <- c(imputedContiValues, imputedData[r, NaNNames[i]])
  #         trueContiValues <- c(trueContiValues, trueData[r, NaNNames[i]])
  #       }
  #       
  #       #in case ordinal trait
  #       else if(length(grep(paste0("/", ordinalTraits), names(missingData)[i]) == 1)){
  #         imputedOrdinalValues <- c(imputedOrdinalValues, imputedData[r, NaNNames[i]])
  #         trueOrdinalValues <- c(trueOrdinalValues, trueData[r, NaNNames[i]])
  #       }
  #       
  #       #in case discrete data
  #       else if(length(grep("I", names(missingData)[i])) == 1){
  #         imputedDiscValues <- c(imputedContiValues, imputedData[r, NaNNames[i]])
  #         trueDiscValues <- c(trueDiscValues, trueData[r, NaNNames[i]])
  #       }
  #     }
  #   }
  # }
  # 
  # if(length(imputedContiValues) != 0){
  #   #imputation error for continuous traits
  #   rmse <- sqrt(mean((imputedContiValues - trueContiValues)^2))
  # }
  # 
  # if(length(imputedOrdinalValues) != 0){
  #   #imputation error for ordinal traits (absolute error)
  #   absError <- mean(abs((imputedOrdinalValues - trueOrdinalValues) / trueOrdinalValues))
  # }
  # 
  # if(length(imputedDiscValues) != 0){
  #   #imputation error for discrete traits
  #   discError <- mean(length(setdiff(imputedDiscValues, trueDiscValues)))
  # }
  
  
  
  
  #errors <- list(rmse = rmse, absError = absError, discError = discError, NaNNames = NaNNames)
  
  return(output)
}
imputedData
TrueData
missingData
error <- imputationError(imputedData, TrueData, missingData)
error2 <- imputationError(imputedData, TrueData, missingData)
f <- merge(error, error2, by = "trait")
f
# Compare imputation approaches with several datasets
####################################################

datasets <- list(NaNImputed)
start_time <- Sys.time()
FinalOutpouts <- list()
FinalImputed <- list()
for(dataset in 1:length(datasets)){
 
  OutputReplicate <- list()
  for(replicate in 1:length(datasets[[dataset]]$DataNaN)){
    
    OutputDataframes <- list()
    ImputedDataNames <- c()
    ImputedData <- list()
    for(random in 1:(length(datasets[[dataset]]$DataNaN[[replicate]])-1)){
      
      subdatas <- datasets[[dataset]]$DataNaN[[replicate]][[random]]
      subdatasNames <- str_split(names(subdatas), "/", simplify = TRUE)

      #create output dataframe
      dataframe <- data.frame(categories = subdatasNames[13:36 ,2], 
                              nbr_traits = subdatasNames[13:36 ,3], 
                              missingness = subdatasNames[13:36 ,4],
                              rows = vector(mode = "character", length = nrow(subdatasNames) - 12),
                              phytools = rep(NA, nrow(subdatasNames) - 12), 
                              Rphylopars = rep(NA, nrow(subdatasNames)-12)) 
                              #add columns for the different phylo imputations
      
      
      #variables to define
      nbrMI <- 5
      k <- 2
      numFun <- weightedMedian
      catFun <- maxCat
      variance_fractions = c(0.4, 0.6, 0.8)
      errorsMixedVec <- c()
      
      for(subdata in 13:length(subdatas)){ 
        
        #to don't have the error
        if(subdata == 6){
          next
        }
        missingData <- subdatas[[subdata]]
        
        TrueValues <- Data$FinalData[, names(missingData), drop = FALSE]
        
        #rows number
        rows <- paste(sort(unique(as.character(str_extract(colnames(missingData), "(?<=\\/)\\d+")))), collapse = "/")
        dataframe[subdata, 4] <- rows
        
        #only discrete data
        if(length(grep("I.", names(subdatas[[subdata]]))) == ncol(subdatas[[subdata]])){
          imputedData <- imputeDiscreteTraits(subdatas[[subdata]])
          error <- paste(as.character(imputationError(imputedData, TrueValues, subdatas[[subdata]])[c(2,3)]),
                         collapse = "/")
          dataframe[subdata, 5] <- error
          
          #save imputed Data
          ImputedDataNames <- c(ImputedDataNames, paste(names(subdatas)[subdata], "imputeDiscrete", sep = "/"))
          ImputedData <- c(ImputedData, list(imputedData))
        }
        
        #only continuous data
        if(length(grep("F.", names(subdatas[[subdata]]))) == ncol(subdatas[[subdata]])){
          imputedData <- imputeContinousTraits(subdatas[[subdata]])
          error <- imputationError(imputedData, TrueValues, subdatas[[subdata]])$rmse
          dataframe[subdata, 6] <- error
          
          #save imputed Data
          ImputedDataNames <- c(ImputedDataNames, paste(names(subdatas)[subdata], "imputeContinuous", sep = "/"))
          ImputedData <- c(ImputedData, list(imputedData))
        }
        
        #mixed and/or only continuous and/or only discrete data
        approach <- 1
        while(approach < (length(MixedImputationApproachesNames)*2)){
          
          for(variance_fraction in variance_fractions){

            MixedImputationApproaches <- list("imputeMICE", list(missingData, nbrMI, method = "pmm", 
                                                            variance_fraction),
                                         "imputeMissForest", list(missingData, variance_fraction,
                                                                  maxiter = 10, ntree = 100, 
                                                                  mtry = sqrt(ncol(missingData))),
                                         "imputeKNN", list(missingData, k, numFun, catFun, 
                                                           variance_fraction))

            #univariate, missrandomForest and imputeKNN don't work.
            if(MixedImputationApproaches[[approach]] == "imputeMICE" | 
               MixedImputationApproaches[[approach]] == "imputeMissForest" | 
               MixedImputationApproaches[[approach]] == "imputeKNN" & ncol(missingData) == 1){
              next
            }
            

            imputedData <- do.call(MixedImputationApproaches[[approach]], MixedImputationApproaches[[approach + 1]])
            error <- imputationError(imputedData, TrueValues, missingData)$rmse
            errorsMixedVec <- c(errorsMixedVec, error)
            
            #save imputed Data
            ImputedDataNames <- c(ImputedDataNames, paste(names(subdatas)[subdata], 
                                                          MixedImputationApproaches[[approach]], sep = "/"))
            ImputedData <- c(ImputedData, list(imputedData))
            
          }
          approach <- approach + 2
        }
        print(subdata)
      }
      #matrix of mixed imputations
      
      errorsMixed <- as.data.frame(matrix(errorsMixedVec, nrow = nrow(subdatasNames), 
                                          ncol = length(MixedImputationApproachesNames) * length(variance_fractions), 
                                   byrow = TRUE))
      
      errorsMixedNames <- rep(paste(MixedImputationApproaches[c(1,3,5)], variance_fractions), 3)
      
      names(errorsMixed) <- errorsMixedNames
      
      dataframe <- cbind(dataframe, errorsMixed)
      
      OutputDataframes[[i]] <- dataframe
    }
    OutputReplicate[[replicate]] <- OutputDataframes
    ImputedDataRep <- list(ImputedData)
    names(ImputedDataRep) <- as.character(replicate)
  }
  FinalOutpouts[[dataset]] <- OutputReplicate
  FinalImputed <- c(FinalImputed, list(ImputedData))
}
end_time <- Sys.time()

print(end_time - start_time)
errorsMixedVec




# to keep, discrete approach to impute data
############################################


#Discrete traits imputation
###########################

calculateAIC <- function(tree, missingData, model){
  out <- tryCatch(
    {
      message("Use corHMM")
      
      #add the tip names in the dataframe
      missingData <- cbind(species = row.names(missingData), missingData)
      
      #convert missingData as character
      missingData[,2] <- as.character(missingData[,2])
      missingData[,2][which(is.na(missingData[,2]))] <- "?" #replace NA by "?"because corHMM don't like it
      #Define the rate model
      model <- "ER"
      FitCorHMM <- corHMM::corHMM(phy = tree, data = missingData, model = model, 
                                  rate.cat = 1, get.tip.states = TRUE)
      
      #Calculate AIC
      AIC <- FitCorHMM$AIC
      FitCorHMM
      list(AIC = AIC, map = FitCorHMM$tip.states, index = 0)
      
    },
    error = function(cond){
      
      message("Use make.simmap")
      
      missingData[,1] <- NULL
      missingData[,1] <- as.factor(missingData[,1])
      
      #extract number of states
      Nstates <- as.numeric(tail(levels(missingData[,1]), n = 1)) + 1
      
      # We need a matrix of prior probabilities for tip states
      StateMat <- generateDummyVariables(missingData)
      colnames(StateMat) <- 0:(ncol(StateMat)-1)
      NaNrowsIndex <- which(is.na(StateMat[,1]))
      StateMat[NaNrowsIndex, ] <- 1/ncol(StateMat)
      
      #Define the rate model
      SimmapTrees <- make.simmap(tree = tree, x = as.matrix(StateMat),
                                 nsim = 100, model = model, pi = "estimated")
      
      #Calculate AIC
      LogLik <- SimmapTrees[[1]]$logL
      K <- attr(SimmapTrees[[1]]$logL, "df") #get the number of degree of freedom (= the nbr of parameter in the model)
      AIC <- 2 * K - 2 * LogLik
      
      SimmapDescribe <- describe.simmap(SimmapTrees, plot = FALSE)
      
      list(AIC = AIC, map = SimmapDescribe, index = 1)
    }
  )
  return(out)
}


#' @title Imputation of missing data for one discrete trait
#' 
#' @description This function imputes missing data for discrete traits using the R package phytools. The first step is to 
#' run the function make.simmap fits a continuous-time reversible Markov model for the evolution of one trait and then 
#' simulates stochastic character histories using that model and the tip states on the tree. Second step, run make.simmap 
#' for three models(ER, SYM and ARD) and select the one having the smallest AIC. Third step is to run the function 
#' describe.simmap which summerize the reuslts obtained with make.simmap
#'
#' @usage imputeOneDiscreteTrait(trait)
#'
#' @param missingData data.frame of 1 factor column containing NAs

#' @return a data.frame of 1 factor column with the NAs replaced by values. 
#missing <- NaNImputed$DataNaN$`1`$MCAR$`MCAR/IndDiscreteTraits/3/0.35`
#missingData <- missingData[ ,2, drop = F]
#missingData <- NaNImputed$DataNaN$`1`$MCAR$`MCAR/IndDiscreteTraits/3/0.35`
#missingData <- partition[,2, drop = F]
imputeOneDisceteTrait <- function(missingData){
  
  #check if tips in matrix traits are ordered as in the tree
  if(!setequal(Data$TreeList$`0`$tip.label, row.names(missingData))){
    
    #change order of the rows, match the order of the phylogeny
    missingData <- missingData[match(Data$TreeList$`0`$tip.label, row.names(missingData)), drop = FALSE]
  }
  
  colName <- names(missingData)
  #print(missingData)
  
  #if only one state represented in 
  
  
  #get the right tree
  correlationGroup <- as.numeric(str_extract(colnames(missingData), "(?<=\\.)\\d+"))[1]
  tree <- Data$TreeList[[correlationGroup + 1]]
  
  AICndMap <- calculateAIC(tree, missingData, "ER")
  
  #Calculate AIC
  AIC <- AICndMap$AIC
  models <- c("SYM", "ARD")
  model <- "ER"
  for (i in 1:length(models)){
    
    AICndMaPDiffModel <- calculateAIC(tree, missingData, model = models[i])
    #FitCorHMMDiffModel <- corHMM::corHMM(phy = tree, data = missingData, model = models[i], 
    #rate.cat = 1, get.tip.states = TRUE)
    AICndMaPDiffModel
    #Calculate AIC
    AICDiffModel <- AICndMaPDiffModel$AIC
    if(AIC > AICDiffModel){
      AIC <- AICDiffModel
      AICndMap <- AICndMaPDiffModel
      model <- models[i]
    }
  }
  
  # Imputation
  if(AICndMap$index == 0){
    MostLikelyState <- apply(AICndMap$map, 1, which.max)
    MostLikelyState <- MostLikelyState - 1
  }
  #MostLikelyState <- apply(FitCorHMM$tip.states, 1, which.max)
  else{
    MostLikelyState <- apply(AICndMap$map$tips, 1, function(x) which.max(x))
    MostLikelyState <- as.data.frame(as.factor(MostLikelyState - 1))
  }
  
  MostLikelyState <- as.data.frame(MostLikelyState)
  colnames(MostLikelyState) <- colName
  return(MostLikelyState)
}
