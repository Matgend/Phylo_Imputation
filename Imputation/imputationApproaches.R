#install.packages("VIM")
# library(VIM)
# library(laeken) #for weigthedMean
# library(mice)
# library(ape)
# library(PVR)
# library(missForest)
# library(tidyverse)
# library(phytools)
# library(mvMORPH)
# library(geiger)
# library(tibble)
# library(Rphylopars)
# library(corHMM)
#install.packages("snow") #for windows 

#load Gain python class
source_python("./Imputation/GAIN/modelV2R.py")
source("utils.R")

#' @title Imputation of missing data for one discrete trait
#' 
#' @description This function imputes missing data for discrete traits using the R package phytools. The first step is to 
#' run the function make.simmap fits a continuous-time reversible Markov model for the evolution of one trait and then 
#' simulates stochastic character histories using that model and the tip states on the tree. Second step, run make.simmap 
#' for three models(ER, SYM and ARD) and select the one having the smallest AIC. Third step is to run the function 
#' describe.simmap which summarize the results obtained with make.simmap
#'
#' @usage imputeOneDiscreteTrait(trait, Data)
#'
#' @param missingData data frame of 1 factor column containing NAs
#' @param Data simulated Data object
#' @return a data frame of 1 factor column with the NAs replaced by values. 
#' 
imputeOneDiscreteTrait <- function(missingData, Data){

  levelsTrait <- levels(missingData[,1])
  
  #check if tips in matrix traits are ordered as in the tree
  if(!setequal(Data$TreeList$`0`$tip.label, row.names(missingData))){
    
    #change order of the rows, match the order of the phylogeny
    missingData <- missingData[match(Data$TreeList$`0`$tip.label, row.names(missingData)), drop = FALSE]
  }

  colName <- names(missingData)
  #if only one state represented in the trait
  if(sum(!is.na(unique(missingData[, 1]))) == 1){
    #print("one state")
    state <- missingData[which(!is.na(missingData)), ][1]
    missingData[which(is.na(missingData)), ] <- state
    return(list(imputedData = missingData, parameters = list("noModel")))
  }
  
  #if trait is ordinal (use random selection imputation?) (now use continuous method)
  row <- as.numeric(str_extract(colnames(missingData), "(?<=\\/)\\d+"))
  if(Data$dataframe[row,"class"] == "ordinal"){
    print("ordinal trait")
    missingData[, 1] <- as.numeric(missingData[, 1]) - 1
    #print(missingData)
    missingData <- imputeContinousTraits(missingData)
    #print(missingData)
    missingData[, 1] <- round(missingData[, 1], 0)
    return(list(imputedData = missingData, parameters = list("noModel")))
  }
  
  tree <- Data$TreeList$`0`

  #add the tip names in the data frame
  missingData <- cbind(species = row.names(missingData), missingData)
  
  #convert missingData as character
  missingData[,2] <- as.character(missingData[,2])
  missingData[,2][which(is.na(missingData[,2]))] <- "?" #replace NA by "?"because corHMM don't like it
  #Define the rate model
  model <- "ER"
  FitCorHMM <- corHMM::corHMM(phy = tree, data = missingData, model = model, rate.mat = , 
                              rate.cat = 1, get.tip.states = TRUE)
  
  
  #Calculate AIC
  AIC <- FitCorHMM$AIC
  models <- c("SYM", "ARD")
  model <- "ER"
  for (i in 1:length(models)){
    
    FitCorHMMDiffModel <- corHMM::corHMM(phy = tree, data = missingData, model = models[i], 
                                  rate.cat = 1, get.tip.states = TRUE)

    #Calculate AIC
    AICDiffModel <- FitCorHMMDiffModel$AIC
    if(AIC > AICDiffModel){
      AIC <- AICDiffModel
      FitCorHMM <- FitCorHMMDiffModel
      model <- models[i]
    }
  }

  #Imputation
  MostLikelyState <- apply(FitCorHMM$tip.states, 1, which.max)
  
  if(levelsTrait[1] == "0"){
    MostLikelyState <- as.data.frame(as.factor(MostLikelyState - 1)) #-1 because when converted a numeric the value are +1
  }

  if(levelsTrait[1] != "0"){
    MostLikelyState <- as.data.frame(as.factor(MostLikelyState)) #-1 because when converted a numeric the value are +1
  }
  
  colnames(MostLikelyState) <- colName
  
  #Parameters
  parameters <- list(model = model, rate.cat = 1)
  
  #print(MostLikelyState)
  return(list(imputedData = MostLikelyState, parameters = parameters, probabilities = FitCorHMM$tip.states))
}


#' @title Imputation of missing values in discrete traits. 
#' 
#' @description This function applies the function imputeOneDiscretTrait on several columns at one time.
#'
#' @usage imputeDiscreteTrait(trait)
#'
#' @param missingData data frame of 1 or more factor columns containing NAs
#' @param Data simulated Data object
#' @return a data frame of 1 or more factor columns with the NAs replaced by values, parameters used for the imputation, matrix 
#' of probabilities of each states
#' 
imputeDiscreteTraits <- function(missingData, Data){

  #select the columns with missing values
  NaNColIndex <- which(apply(missingData, 2, function(x) any(is.na(x))))
  
  proba <- matrix(NA, nrow = nrow(missingData), ncol = ncol(missingData))
  parameters <- list()
  for(i in 1:length(names(NaNColIndex))){
    imputation <- imputeOneDiscreteTrait(missingData[, names(NaNColIndex)[i], drop = FALSE], Data)
    #print("ok")
    missingData[ ,names(NaNColIndex)[i]] <- imputation$imputedData
    
    if(i == 1){
      proba <- imputation$probabilities
    }
    
    else{
      proba <- cbind(proba, imputation$probabilities)
    }
    parameters <- c(parameters, imputation$parameters)
  }
  proba <- as.data.frame(proba, row.names = rownames(missingData))
  
  return(list(imputedData = missingData, parameters = parameters, probabilities = proba))
}


#' @title Imputation of missing data in continuous traits
#' 
#' @description This function imputes missing data for continuous traits applying the phylopars approach. 
#'
#' @usage imputeContinousTraits(missingData, Data)
#'
#' @param missingData data.frame of 1 or more numeric columns containing NAs
#' @param Data simulated Data object
#' @return a data frame of 1 or more numeric columns with the NAs replaced by values + parameters used for the imputation
#' 
imputeContinousTraits <- function(missingData, Data){
  
  Rphylopars <- tryCatch(
   {
     tree <- Data$TreeList$`0`
      
     if(length(setdiff(rownames(missingData), tree$tip.label)) != 0){
       rownames(missingData) <- tree$tip.label
     }
  
     #add species (tips) in the dataframe
     missingData <- tibble::add_column(missingData, species = row.names(missingData) , .before = 1)

     #impute
     models <- c("BM", "OU")
     AICs <- c()
     imputations <- list()
  
     phylo_correlated <- TRUE 
     pheno_correlated <- FALSE
  
     for(i in 1:length(models)){

        imputeData <- phylopars(trait_data = missingData, tree = tree, model = models[i],
                               phylo_correlated = phylo_correlated, pheno_correlated = pheno_correlated)

        imputations[[i]] <- imputeData
    
     #Calculate AIC
     LogLik <- imputeData$logLik
     K <- imputeData$npars #get the number of degree of freedom (= the nbr of parameter in the model)
     AIC <- 2 * K - 2 * LogLik
     AICs <- c(AICs, AIC)
     }

     #keep only the imputed data
     imputedData <- imputations[[which.min(AICs)]]$anc_recon
  
     #keep only the tip labels and not the node labels
     imputedData <- imputedData[1:nrow(missingData), , drop = FALSE]

     #put again the tips names
     rownames(imputedData) <- tree$tip.label

     parameters <- list(model = imputations[[which.min(AICs)]]$model)
     
     list(imputedData = imputedData, parameters = parameters)     
   },
   error = function(cond){
     message("Error running phylopars")
     list(imputedData = missingData, parameters = "Error in Cholesky decomposition")
   }
 )
 return(Rphylopars)
}
     

# MICE
######
#' @title Imputation of missing values in mixed data by the MICE approach
#' 
#' @description This function imputes missing data for continuous and discrete traits applying the MICE approach using the 
#' predictive mean matching method for continuous traits, the polytomous regression for unordered categorical traits having 
#' >2 states, if binary the algo uses a logistic regression and in case of ordered traits with >2 states, the method uses 
#' proportional odds model.
#'
#' @usage imputeMICE(missingData, nbrMI, variance_fraction = 0, Data, hint = NULL)
#'
#' @param missingData data frame of 1 or more columns containing NAs
#' @param nbrMI integer, mentioning the total number of imputations
#' @param variance_fraction variance_fraction minimum variance (%) explained by the eigenvectors
#' @param Data simulated Data object
#' @param hint data frame already imputed by a comparative methods (Rphylopars for continuous and corHMM for discrete)
#' @return a data frame of 1 or more numeric columns with the NAs replaced by values + parameters used for the imputation
#' 
imputeMICE <- function(missingData, nbrMI, variance_fraction = 0, Data, hint = NULL){
  
  colNames <- names(missingData)
  nativeMissingData <- missingData
  
  if(!is.null(hint) & variance_fraction == 2){
    #change names columns hint
    colnames(hint) <- paste0("H",as.character(1:ncol(hint)))
    
    missingData <- cbind(missingData, hint)
  }
  
  if(variance_fraction != 0 & variance_fraction != 2){
    
    tree <- Data$TreeList$`0`
    
    eigen <- get_eigenvec(tree, variance_fraction)
    missingData <- cbind(missingData[, 1:ncol(missingData), drop = FALSE],
                         eigen[, 1:ncol(eigen), drop = FALSE])
    
  }
  names(missingData) <- paste0("A",as.character(1:ncol(missingData)))
  
  #in case collinearity is too important, and MICE can't impute the data
  ImputedMICE <- tryCatch(
    {
      
      #check the number of levels for factors
      levelsByColumns <- lapply(missingData, levels)
      lengthLevels <- lapply(levelsByColumns, length)
      columnsCatPMM <- as.numeric(lengthLevels > 20) #in literature said that in case nbr state > 20 use pmm instead of polyreg
      columnsLogReg <- as.numeric(lengthLevels == 2)
      
      method <- NULL
      
      if(sum(columnsCatPMM) != 0){
        method <- rep("polyreg", length(columnsCatPMM))
        method[which(columnsCatPMM == 1 | as.numeric(lengthLevels == 0) == 1)] <- "pmm"
        method[which(columnsLogReg == 1)] <- "logreg"
      }
      
      ImputedMICE <- mice(missingData, m = nbrMI, maxit = 5, method = method, printFlag = FALSE)
      
      #choose the first column 
      imputedData <- mice::complete(ImputedMICE, action = 1L)[, 1:length(colNames),  drop = FALSE]
      names(imputedData) <- colNames

      colWithNaN <-  which(colSums(is.na(imputedData)) != 0)
      if(length(colWithNaN) != 0){
        for(v in names(colWithNaN)){
          ry <- !is.na(imputedData[ ,v])
          imputedData[which(!ry) ,v] <- mice.impute.sample(imputedData[ ,v], ry)
        }
      }
      
      parameters <- list(nbrMI = nbrMI)
      list(imputedData = imputedData, parameters = parameters)
    },
    error = function(cond){
      list(imputedData = nativeMissingData, parameters = "Collinearity induces imputation error")
    }
  )
  return(ImputedMICE)
}

#' @title Non-parametric missing values imputation for mixed-type data by missForest(randomForest)
#' 
#' @description This function imputes missing data for continuous and discrete traits applying the missForest approach  
#'
#' @usage imputeMissForest(missingData, variance_fraction = 0, maxiter = 10, ntree = 100,
#' mtry = sqrt(ncol(missingData)), Data, hint = NULL)
#'
#' @param missingData data frame of 1 or more columns containing NAs
#' @param variance_fraction variance_fraction minimum variance (%) explained by the eigenvectors
#' @param maxiter maximum number of iterations to be performed given the stopping criterion is not met beforehand.
#' @param ntree number of trees to grow in each forest.
#' @param mtry number of variables randomly sampled at each split. By default it's the square root of the number of 
#' variables
#' @param hint dataframe already imputed by a comprative methods (Rphylopars for continuous and corHMM for discrete)
#' @return a data frame of 1 or more numeric columns with the NAs replaced by values + parameters used for the imputation
#' 
imputeMissForest <- function(missingData, variance_fraction = 0, maxiter = 10, ntree = 100, 
                             mtry = sqrt(ncol(missingData)), Data, hint = NULL){
  
  Nvariables <- ncol(missingData)
  
  #include imputed data in case 
  if(!is.null(hint) & variance_fraction == 2){
    #change names columns hint
    colnames(hint) <- paste0("H",as.character(1:ncol(hint)))
    
    missingData <- cbind(missingData, hint)
  }

  # want to include phylogeny information
  if(variance_fraction != 0 & variance_fraction != 2){
    
    tree <- Data$TreeList$`0`
    
    eigen <- get_eigenvec(tree, variance_fraction)
    missingData <- cbind(missingData[, 1:ncol(missingData), drop = FALSE], 
                         eigen[row.names(missingData), 1:ncol(eigen), drop = FALSE])
  }

  #run missForest
  missForest_imputation <- missForest(xmis = missingData, maxiter = maxiter, ntree = ntree, mtry = mtry, verbose = TRUE)
  
  #cut eigenvectors columns
  missForest_imputation$ximp <- missForest_imputation$ximp[,1:Nvariables, drop = FALSE]
  
  parameters <- list(maxiter = maxiter, ntree = ntree, mtry = mtry)
  
  return(list(imputedData = missForest_imputation$ximp, parameters = parameters))
}


#' @title non-parametric missing value imputation for mixed-type data by kNN(k-nearest neighbor algorithm)
#' 
#' @description This function imputes missing data for continuous and discrete traits applying the kNN approach  
#'
#' @usage imputeKNN(missingData, k, numFun, catFun, variance_fraction = 0, Data, hint = NULL)
#'
#' @param missingData data frame of 1 or more columns containing NAs
#' @param k integer, number of nearest neighbours used
#' @param numFun numFun: function for aggregating the kNN in the case of a numerical variable
#' @param catFun catFun: function for aggregating the kNN in the case of a categorical variable
#' @param variance_fraction variance_fraction minimum variance (%) explained by the eigenvectors
#' @param Data simulated Data object
#' @param hint data frame already imputed by a comprative methods (Rphylopars for continuous and corHMM for discrete)
#' @return a data frame of 1 or more numeric columns with the NAs replaced by values + parameters used for the imputation
#' 
imputeKNN <- function(missingData, k, numFun, catFun, variance_fraction = 0, Data, hint = NULL){
  
  NbrCol <- ncol(missingData)
  
  #include imputed data in case 
  if(!is.null(hint) & variance_fraction == 2){
    #change names columns hint
    colnames(hint) <- paste0("H",as.character(1:ncol(hint)))
    
    missingData <- cbind(missingData, hint)
  }
    
  if(variance_fraction != 0 & variance_fraction != 2){
    
    tree <- Data$TreeList$`0`
    
    eigen <- get_eigenvec(tree, variance_fraction)
    missingData <- cbind(missingData[, 1:ncol(missingData), drop = FALSE], 
                         eigen[row.names(missingData), 1:ncol(eigen), drop = FALSE])
  }
  
  #get the column name with NaN
  variable <- c()
  for(c in 1:ncol(missingData)){
    if(any(is.na(missingData[,c]))){
      variable <- c(variable, colnames(missingData)[c])
    }
  }
  
  DataImputed <- VIM::kNN(missingData, k = k, variable = variable,
                          numFun = numFun, catFun = catFun, imp_var = FALSE)[,1:NbrCol,  drop = FALSE]
  #keep the tip names
  rownames(DataImputed) <- rownames(missingData)

  parameters <- list(k = k, numFun = numFun, catFun = catFun)
  
  return(list(imputedData = DataImputed, parameters = parameters))
}


#' @title GAIN
#' 
#' @description This function imputes missing data for continuous and discrete traits applying the Generative Adversarial #' 
#' Imputation Networks (GAIN) 
#'
#' @usage gainR(missingData, variance_fraction, Data, batch_size = round(ncol(missingData)*0.2), 
#' hint_rate = 0.9, alpha = 100, epochs = 10000, hint = NULL)
#'
#' @param missingData data.frame of 1 or more columns containing NAs
#' @param variance_fraction total amount of minimum variance to be represented by the eigenvectors which correspond to the 
#' phylogenetic inertia
#' @param Data simulated Data object
#' @param batch_size integer
#' @param hint_rate numerical 
#' @param alpha numerical, hyperparameter
#' @param epochs integer, iterations
#' @param hint dataframe already imputed by a comprative methods (Rphylopars for continuous and corHMM for discrete)
#' @return a list of list containing in the "tab" imputedData the imputed Data and in the "tab" parametersGain, the discriminant 
#' loss values, the generative loss values, the MSE loss values and the iterations correspond to these values (important to 
#' generate a plot).
#' 
gainR <- function(missingData, variance_fraction, Data, batch_size = round(ncol(missingData)*0.2), 
                  hint_rate = 0.9, alpha = 100, epochs = 10000, hint = NULL){
  
  NbrCol <- ncol(missingData)
  
  #get the factors columns
  factorsColumns <- names(Filter(is.factor, missingData))
  
  rNames <- row.names(missingData)
  colNames <- colnames(missingData)
  
  #include imputed data in case
  if(!is.null(hint)){
    #change names columns hint
    colnames(hint) <- paste0("H",as.character(1:ncol(hint)))
    
    missingData <- cbind(missingData, hint)
  }

  # want to include phylogeny information
  if(variance_fraction != 0 & variance_fraction != 2){
    
    tree <- Data$TreeList$`0`
    eigen <- get_eigenvec(tree, variance_fraction)
    missingData <- cbind(missingData[, 1:ncol(missingData), drop = FALSE],
                         eigen[row.names(missingData), 1:ncol(eigen), drop = FALSE])
  }

  #convert factors as dummy
  if(length(factorsColumns) != 0){
    oneHotConvertion <- generateOneHotEncoVariables(missingData)
    missingData <- oneHotConvertion$data
    oneHotColnames <- colnames(missingData)
  }
  
  #convert all the values as numeric (second check)
  missingData <- apply(missingData, 2, as.numeric)


  #run Gain
  obj <- Gain()
  gainOutput <- obj$runGain(missingData, batch_size, hint_rate, alpha, epochs)
  
  #merge list of list
  for (l in 2:length(gainOutput)){
    gainOutput[[l]] <- unlist(gainOutput[[l]], recursive = FALSE)
  }
  
  #convert imputedData in dataframe
  imputedData <- as.data.frame(gainOutput[[1]])
  
  if(length(factorsColumns) != 0){
    names(imputedData) <- oneHotColnames
    imputedData <- convertOneHotEncoIntoFactors(imputedData, colNames)
  }
  
  #remove hint variables or phylogenetic eigenvectors (in case no one hot encoding)
  imputedData <- imputedData[, 1:NbrCol, drop = FALSE]
  names(imputedData) <- colNames
  row.names(imputedData) <- rNames

  #parameters
  lossValues <- gainOutput[c(2:5)]
  names(lossValues) <- c("d_loss", "g_loss", "mse_loss", "epochs")
  
  return(list(imputedData = imputedData, parameters = lossValues))
}

#v <- gainR(missingData, 0.2, simulatedData, batch_size)
#v$imputedData  




