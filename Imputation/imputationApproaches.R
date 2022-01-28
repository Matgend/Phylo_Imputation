#install.packages("VIM")
# #install.packages("fastDummies")
# library(fastDummies)
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

#Categorical in dummy
#####################
#' @title Conversion factors in dummy
#' 
#' @description This function returns a data.frame of dummy corresponding to a data.frame of factors      
#' representing a categorical variable.
#'
#' @usage generateDummyVariables(NaNData)
#'
#' @param NaNData data.frame of one ore several factors columns 
#' @return a data.frame in which each variable are represented as dummies.
generateDummyVariables <- function(NaNData){
  
  Nstates <- as.numeric(tail(levels(NaNData[,1]), n = 1)) + 1
  
  #in case NaNData is a vector
  if(is.null(dim(NaNData))){
    stop("The column(s) shoud be of class data.frame")
  }
  
  #check if all the columns are factors
  if(sum(sapply(NaNData, is.factor)) != ncol(NaNData)){
    stop("The columns should be of class factor")
  }
  
  dummy <- fastDummies::dummy_cols(NaNData)
  row.names(dummy) <- row.names(NaNData)
  columns_to_remove <- c(1:ncol(NaNData), grep("NA", colnames(dummy)))
  dummy <- dummy[, -columns_to_remove]
  
  if(is.integer(dummy) & Nstates != 1){
    stateInData <- as.numeric(levels(NaNData[,1]))+1
    intermediateDummy <- dummy
    intermediateDummy[which(intermediateDummy == 1)] <- 0
    intermediateMatrix <- matrix(rep(intermediateDummy, Nstates), ncol = Nstates)
    intermediateMatrix[, stateInData] <- dummy
    dummy <- intermediateMatrix
    row.names(dummy) <- row.names(NaNData)
  }

  return(dummy)
  
}

# Generating eigenvectors
#########################
#' @title Eigenvectors calculations
#' 
#' @description This function calculates the eigenvectors of a phylo object 
#'
#' @usage get_eigenvec(tree, variance_fraction = 0.8, numEigen = NULL)
#'
#' @param tree phylogenetic tree of class "phylo"
#' @param variance_fraction total amount of minimum variance to be represented by the eigenvectors which correspond to the 
#' phylogenetic inertia
#' @return a data.frame in which each colum represent an eigenvector
get_eigenvec <- function(tree, variance_fraction){
  
  decomp <- PVRdecomp(tree, type = "newick") #produces object of class 'PVR'
  eigvec <- as.data.frame(decomp@Eigen$vectors) ##extract eigenvectors
  
  egval <- decomp@Eigen$values #extract eigenvalues
  eigPerc <- egval/(sum(egval)) #calculate % of variance
  eigPercCum <- t(cumsum(eigPerc)) #cumulated variance
  
  #eigenvectors representing more than X% variance
  numEigen <- sum(eigPercCum < variance_fraction)
  
  eigOK <- eigvec[,1:numEigen, drop = T] 
  
  if(numEigen == 0){
    print(paste("The variance_fraction should at leat be equal to ", 
                eigPercCum[1]))
    eigOK <- eigvec[1]
  }
  # Change 'numEigen' on above line to a number if you want to specify number of eigenvectors
  #Eigenvectors generated in object 'eigenTobind'
  #rename eigenTobind species column so it matches trait dataset species column
  eigOK <- as.data.frame(eigOK)
  names(eigOK)[1] <- "c1"
  row.names(eigOK) <- decomp@phylo$tip.label

  return(eigOK)
}


#' @title Imputation of missing data for one discrete trait
#' 
#' @description This function imputes missing data for discrete traits using the R package phytools. The first step is to 
#' run the function make.simmap fits a continuous-time reversible Markov model for the evolution of one trait and then 
#' simulates stochastic character histories using that model and the tip states on the tree. Second step, run make.simmap 
#' for three models(ER, SYM and ARD) and select the one having the smallest AIC. Third step is to run the function 
#' describe.simmap which summerize the reuslts obtained with make.simmap
#'
#' @usage imputeOneDiscreteTrait(trait, Data)
#'
#' @param missingData data.frame of 1 factor column containing NAs
#' @param Data simulated Data object
#' @return a data.frame of 1 factor column with the NAs replaced by values. 
imputeOneDiscreteTrait <- function(missingData, Data){

  #check if tips in matrix traits are ordered as in the tree
  if(!setequal(Data$TreeList$`0`$tip.label, row.names(missingData))){
    
    #change order of the rows, match the order of the phylogeny
    missingData <- missingData[match(Data$TreeList$`0`$tip.label, row.names(missingData)), drop = FALSE]
  }
  
  colName <- names(missingData)
  #print(missingData)
  
  #if only one state represented in the trait
  if(sum(!is.na(unique(missingData[, 1]))) == 1){
    #print("one state")
    state <- missingData[which(!is.na(missingData)), ][1]
    missingData[which(is.na(missingData)), ] <- state
    return(missingData)
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
    return(missingData)
  }
  
  #get the right tree
  correlationGroup <- as.numeric(str_extract(colnames(missingData), "(?<=\\.)\\d+"))[1]
  tree <- Data$TreeList[[correlationGroup + 1]]

  #add the tip names in the dataframe
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

  # Imputation
  MostLikelyState <- apply(FitCorHMM$tip.states, 1, which.max)
  MostLikelyState <- MostLikelyState - 1
  MostLikelyState <- as.data.frame(MostLikelyState)
  colnames(MostLikelyState) <- colName
  
  return(MostLikelyState)
}


#' @title Imputation of missing data for discrete traits, mandatory phylogenetic information  
#' 
#' @description This function apply the function imputeOneDiscretTrait on several columns at one time.
#'
#' @usage imputeDiscreteTrait(trait)
#'
#' @param missingData data.frame of 1 or more factor columns containing NAs
#' @param Data simulated Data object
#' @return a data.frame of 1 or more factor columns with the NAs replaced by values. 
imputeDiscreteTraits <- function(missingData, Data){


  #select the columns with missing values
  NaNColIndex <- which(apply(missingData, 2, function(x) any(is.na(x))))
  #NaNTraits <- missingData[, NaNColIndex, drop = FALSE]
  
  for(i in NaNColIndex){
    #print(i)
    #print(missingData[,i])
    missingData[, i] <- imputeOneDiscreteTrait(missingData[, i, drop = FALSE], Data)
    #print(missingData[,i])
  }
  
  return(missingData)
}

#impute continuous traits (Rphylopars)
######################################
#' @title Imputation of missing data for continuous traits, mandatory phylogenetic information 
#' 
#' @description This function imputes missing data for continuous traits applying the phylopars approach. 
#'
#' @usage imputeContinous(missingData, tree)
#'
#' @param missingData data.frame of 1 or more numeric columns containing NAs
#' @param Data simulated Data object
#' @return a data.frame of 1 or more numeric columns with the NAs replaced by values.

imputeContinousTraits <- function(missingData, Data){
  
  #get the right tree
  correlationGroup <- as.numeric(str_extract(colnames(missingData), "(?<=\\.)\\d+"))[1]
  tree <- Data$TreeList[[correlationGroup + 1]]
  
  if(length(setdiff(rownames(missingData), tree$tip.label)) != 0){
    rownames(missingData) <- tree$tip.label
  }

  #add species (tips) in the dataframe
  missingData <- tibble::add_column(missingData, species = row.names(missingData) , .before = 1)

  #impute
  models <- c("BM", "OU")
  AICs <- c()
  imputations <- list()
  
  for(i in 1:length(models)){
    
    imputeData <- phylopars(trait_data = missingData, tree = tree, model = models[i])
    imputations[[i]] <- imputeData
    
    #Calculate AIC
    LogLik <- imputeData$logLik
    K <- imputeData$npars #get the number of degree of freedom (= the nbr of parameter in the model)
    AIC <- 2 * K - 2 * LogLik
    AICs <- c(AICs, AIC)
  }

  #keep only the imputed data
  imputeData <- imputations[[which.min(AICs)]]$anc_recon
  
  #keep only the tip labels and not the node labels
  imputeData <- imputeData[1:nrow(missingData), , drop = FALSE]

  return(imputeData)
}

# MICE
######
#' @title Imputation of missing data for continuous and discrete traits by MICE approach
#' 
#' @description This function imputes missing data for continuous and discrete traits applying the MICE approach using the 
#' predictive mean matching method. 
#'
#' @usage imputeMICE(missingData, nbrMI, method, tree)
#'
#' @param missingData data.frame of 1 or more columns containing NAs
#' @param nbrMI integer, mentioning the total number of imputations
#' @param method character, name of the imputation method used, by default is pmm (predictive mean matching)
#' @param variance_fraction total amount of minimum variance to be represented by the eigenvectors which correspond to the 
#' phylogenetic inertia
#' @param Data simulated Data object
#' @return a data.frame of 1 or more numeric columns with the NAs replaced by values.

imputeMICE <- function(missingData, nbrMI, method = "pmm", variance_fraction = 0, Data){
  
  colNames <- names(missingData)

  if(variance_fraction != 0){
    
    #get the right tree
    correlationGroup <- as.numeric(str_extract(colnames(missingData), "(?<=\\.)\\d+"))[1]
    
    tree <- Data$TreeList[[correlationGroup + 1]]
    
    eigen <- get_eigenvec(tree, variance_fraction)
    missingData <- cbind(missingData[, 1:ncol(missingData), drop = FALSE],
                         eigen[row.names(missingData), 1:ncol(eigen), drop = FALSE])
    
  }
  names(missingData) <- paste0("A",as.character(1:ncol(missingData)))

  ImputedMICE <- mice(missingData, m = nbrMI, method = method, maxit = 5, printFlag = FALSE)
  
  #take the columns being closer the to the mean of the variable with missing data.
  NaNColIndex <- which(apply(missingData, 2, function(x) any(is.na(x))))

  #choose the first column 
  imputedData <- mice::complete(ImputedMICE, action = 1)[, 1:length(colNames)]
  names(imputedData) <- colNames
  return(imputedData)
}

# missForest
############
#' @title non-parametric missing value imputation for mixed-type data by missForest(randomForest)
#' 
#' @description This function imputes missing data for continuous and discrete traits applying the missForest approach  
#'
#' @usage imputeMissForest(missingData, variance_fraction = 0.8, maxiter = 10, ntree = 100, 
#' mtry = sqrt(ncol(missingData)), tree = NULL)
#'
#' @param missingData data.frame of 1 or more columns containing NAs
#' @param variance_fraction total amount of minimum variance to be represented by the eigenvectors
#' @param maxiter maximum number of iterations to be performed given the stopping criterion is not met beforehand.
#' @param ntree number of trees to grow in each forest.
#' @param mtry number of variables randomly sampled at each split. By default it's the square root of the number of 
#' variables
#' @param variance_fraction total amount of minimum variance to be represented by the eigenvectors which correspond to the 
#' phylogenetic inertia
#' @return a data.frame of 1 or more numeric columns with the NAs replaced by values.

imputeMissForest <- function(missingData, variance_fraction = 0, maxiter = 10, ntree = 100, 
                             mtry = sqrt(ncol(missingData)), Data){
  
  Nvariables <- ncol(missingData)

  # want to include phylogeny information
  if(variance_fraction != 0){
    
    #get the right tree
    correlationGroup <- as.numeric(str_extract(colnames(missingData), "(?<=\\.)\\d+"))[1]
    tree <- Data$TreeList[[correlationGroup + 1]]
    
    eigen <- get_eigenvec(tree, variance_fraction)
    missingData <- cbind(missingData[, 1:ncol(missingData), drop = FALSE], 
                         eigen[row.names(missingData), 1:ncol(eigen), drop = FALSE])
  }

  #run missForest
  missForest_imputation <- missForest(xmis = missingData, maxiter = maxiter, ntree = ntree, mtry = mtry)
  
  #cut eigenvectors columns
  missForest_imputation$ximp <- missForest_imputation$ximp[,1:Nvariables, drop = FALSE]
  
  return(missForest_imputation$ximp)
}

# KNN
#####
#' @title non-parametric missing value imputation for mixed-type data by kNN(k-nearest neighbor algorithm)
#' 
#' @description This function imputes missing data for continuous and discrete traits applying the kNN approach  
#'
#' @usage imputeKNN(missingData, k, numFun, catFun, variance_fraction = 0.8, tree = NULL)
#'
#' @param missingData data.frame of 1 or morr columns containing NAs
#' @param k integer, number of nearest neighbours used
#' @param numFun numFun: function for aggregating the kNN in the case of a numerical variable
#' @param catFun catFun: function for aggregating the kNN in the case of a categorical variable
#' @param variance_fraction total amount of minimum variance to be represented by the eigenvectors which correspond to the 
#' phylogenetic inertia
#' @param Data simulated Data object
#' @return a data.frame of 1 or more numeric columns with the NAs replaced by values.

imputeKNN <- function(missingData, k, numFun, catFun, variance_fraction = 0, Data){
  
  NbrCol <- ncol(missingData)
    
  if(variance_fraction != 0){
    
    #get the right tree
    correlationGroup <- as.numeric(str_extract(colnames(missingData), "(?<=\\.)\\d+"))[1]
    tree <- Data$TreeList[[correlationGroup + 1]]
    
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
                          numFun = numFun, catFun = catFun, imp_var = FALSE)[,1:NbrCol]
  #keep the tip names
  rownames(DataImputed) <- rownames(missingData)

  return(DataImputed)
}

