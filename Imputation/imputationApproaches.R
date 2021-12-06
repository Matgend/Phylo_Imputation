#install.packages("VIM")
#install.packages("fastDummies")
library(fastDummies)
library(VIM)
library(laeken) #for weigthedMean
library(mice)
library(ape)
library(PVR)
library(missForest)
library(tidyverse)
library(phytools)
library(mvMORPH)
library(geiger)
library(tibble)
library(Rphylopars)
#install.packages("snow") #for windows 


# load simulated data
setwd("C:/Users/Matthieu/Documents/UNIFR/Master_thesis/Scripts/")
load("simulatedData.RData")

# load NaN imputed data
load("DataNaN.RData")
test <- NaNImputed$DataNaN$`1`$MCAR$`MCAR/MixedCorrelatedTraits/3/0.2`
testCat <- NaNImputed$DataNaN$`1`$MAR$`MAR/CorrDiscreteTraits/4/0.35`
testConti <- NaNImputed$DataNaN$`1`$MNAR$`MNAR/CorrContinuousTraits/3/0.2`


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
  #'NaNData: data.frame containing missing NAs
  #'return: the dummy matrix of the categorical variables
  
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
  
  return(dummy)
  
}

onetestCat <- as.data.frame(testCat[,2])
rownames(onetestCat) <- rownames(testCat)
colnames(onetestCat) <- "I2.1"
generateDummyVariables(testCat)
generateDummyVariables(onetestCat)


# Generating eigenvectors
#########################
#' @title Eigenvectors calculations
#' 
#' @description This function calculates the eigenvectors of a phylo object 
#'
#' @usage get_eigenvec(tree, variance_fraction = 0.8, numEigen = NULL)
#'
#' @param tree phylogenetic tree of class "phylo"
#' @param variance_fraction total amount of minimum variance to be represented by the eigenvectors
#' @param numEigen minimal number of eigenvectors
#' @return a data.frame in which each colum represent an eigenvector
get_eigenvec <- function(tree, variance_fraction = 0.8, numEigen = NULL){
  
  decomp <- PVRdecomp(tree, type = "newick") #produces object of class 'PVR'
  
  eigvec <- as.data.frame(decomp@Eigen$vectors) ##extract eigenvectors
  
  if (is.null(numEigen)){
    egval <- decomp@Eigen$values #extract eigenvalues
    eigPerc <- egval/(sum(egval)) #calculate % of variance
    eigPercCum <- t(cumsum(eigPerc)) #cumulated variance
    
    #eigenvectors representing more than X% variance
    numEigen <- sum(eigPercCum < variance_fraction) 
  }
  eigOK <- eigvec[,1:numEigen] 
  # Change 'numEigen' on above line to a number if you want to specify number of eigenvectors
  #Eigenvectors generated in object 'eigenTobind'
  #rename eigenTobind species column so it matches trait dataset species column
  row.names(eigOK) <- decomp@phylo$tip.label
  return(eigOK)
  
}

get_eigenvec(Data$TreeList$`0`)
Data$TreeList$`0`$tip.label

# 
# runMake.SimmapFaster <- function(nbrCore, tree, StateMat, nsim, models){
#   
#   tree.mcs <- list()
#   cl <- makeSOCKcluster(rep("localhost", nbrCore))
#   
#   for(i in 1:length(models)){
#   
#     trees.mc<-clusterApply(cl, x = replicate(nbrCore, as.matrix(StateMat), simplify = FALSE),
#                          fun = make.simmap, tree = tree, nsim = nsim, model = models[i], pi = "estimated")
# 
#     trees.mc <- do.call(c, trees.mc)
#     
#     if(!("multiSimmap" %in% class(trees.mc))){
#       class(trees.mc)<-c("multiSimmap", class(trees.mc))
#     }
#     tree.mcs[[i]] <- trees.mc
#     print("Done")
#   }
#   stopCluster(cl)
#   return(tree.mcs)
# }
# tree <- Data$TreeList$`0`
# nbrCore <- 4
# model <- "SYM"
# system.time(runMake.SimmapFaster(4, Data$TreeList$`0`, StateMat, 100, c("ER","SYM")))


#Discrete traits imputation
###########################
#' @title Imputation of missing data for one discrete trait
#' 
#' @description This function imputes missing data for discrete traits using the R package phytools. The first step is to 
#' run the function make.simmap fits a continuous-time reversible Markov model for the evolution of one trait and then 
#' simulates stochastic character histories using that model and the tip states on the tree. Second step, run make.simmap 
#' for three models(ER, SYM and ARD) and select the one having the smallest AIC. Third step is to run the function 
#' describe.simmap which summerize the reuslts obtained with make.simmap
#'
#' @usage imputeOneDiscreteTrait(trait, nsim)
#'
#' @param trait data.frame of 1 factor column containing NAs
#' @param nsim number of simulation for the Markov Model.

#' @return a data.frame of 1 factor column with the NAs replaced by values. 

imputeOneDisceteTrait <- function(trait, nsim){

  #check if tips in matrix traits are ordered as in the tree
  if(!setequal(Data$TreeList$`0`$tip.label, row.names(trait))){
    
    #change order of the rows, match the order of the phylogeny
    trait <- trait[match(Data$TreeList$`0`$tip.label, row.names(trait)), , drop = FALSE]
  }
  
  #extract number of states
  Nstates <- as.numeric(tail(levels(trait[,1]), n = 1))
  
  # We need a matrix of prior probabilities for tip states
  StateMat <- generateDummyVariables(trait)
  colnames(StateMat) <- 0:(ncol(StateMat)-1)
  NaNrowsIndex <- which(is.na(StateMat[,1]))
  StateMat[NaNrowsIndex, ] <- 1/ncol(StateMat)

  #get the right tree and model
  correlationGroup <- as.numeric(str_extract(colnames(trait), "(?<=\\.)\\d+"), 1)
  tree <- Data$TreeList[[correlationGroup]]

  #Define the rate model
  SimmapTrees <- make.simmap(tree = tree, x = as.matrix(StateMat),
                             nsim = nsim, model = "ER", pi = "estimated")
  #Calculate AIC
  LogLik <- SimmapTrees[[1]]$logL
  K <- attr(SimmapTrees[[1]]$logL, "df") #get the number of degree of freedom (= the nbr of parameter in the model)
  AIC <- 2 * K - 2 * LogLik
  models <- c("SYM", "ARD")
  model <- "ER"
  for (i in 1:length(models)){

    SimmapTreesDiffModel <- make.simmap(tree = tree, x = as.matrix(StateMat), 
                                        nsim = nsim, model = models[i], pi = "estimated")
    #Calculate AIC
    LogLik <- SimmapTreesDiffModel[[1]]$logL
    K <- attr(SimmapTreesDiffModel[[1]]$logL, "df") 
    AICDiffModel <- 2 * K - 2 * LogLik
    
    if(AIC[1] > AICDiffModel[1]){
      AIC <- AICDiffModel
      SimmapTrees <- SimmapTreesDiffModel
      model <- models[i]
    }
  }

  # Imputation
  SimmapDescribe <- describe.simmap(SimmapTrees, plot = FALSE)
  MostLikelyState <- apply(SimmapDescribe$tips, 1, function(x) which.max(x))
  MostLikelyState <- as.data.frame(as.factor(MostLikelyState - 1))
  colnames(MostLikelyState) <- colnames(trait)
  return(MostLikelyState)
}

imputeOneDisceteTrait(onetestCat, nsim = 100)


#' @title Imputation of missing data for discrete traits, mandatory phylogenetic information  
#' 
#' @description This function apply the function imputeOneDiscretTrait on several columns at one time.
#'
#' @usage imputeDiscreteTrait(trait, nsim)
#'
#' @param trait data.frame of 1 or more factor columns containing NAs
#' @param nsim number of simulation for the Markov Model.

#' @return a data.frame of 1 or more factor columns with the NAs replaced by values. 

imputeDiscreteTraits <- function(trait, nsim){


  #select the columns with missing values
  NaNColIndex <- which(apply(trait, 2, function(x) any(is.na(x))))
  NaNTraits <- trait[, NaNColIndex, drop = FALSE]

  
  for(i in NaNColIndex){
    trait[, i] <- imputeOneDisceteTrait(NaNTraits, nsim)
  }
  
  return(trait)
}

imputeDiscreteTraits(trait, 100)


#impute continuous traits (Rphylopars)
######################################
#' @title Imputation of missing data for continuous traits, mandatory phylogenetic information 
#' 
#' @description This function imputes missing data for continuous traits applying the phylopars approach. 
#'
#' @usage imputeContinous(missingData, tree)
#'
#' @param missinData data.frame of 1 or more numeric columns containing NAs
#' @param tree phylogenetic tree of class "phylo"

#' @return a data.frame of 1 or more numeric columns with the NAs replaced by values.

imputeContinous <- function(missingData, tree){
  
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

  return(imputeData)
}

PPE <- imputeContinous(testConti, Data$TreeList$`0`)
PPE

# MICE
######
#' @title Imputation of missing data for continuous and discrete traits by MICE approach
#' 
#' @description This function imputes missing data for continuous and discrete traits applying the MICE approach using the 
#' predictive mean matching method. 
#'
#' @usage imputeMICE(missingData, nbrMI, method, tree)
#'
#' @param missinData data.frame of 1 or morr columns containing NAs
#' @param nbrMI integer, mentioning the total number of imputations
#' @param method character, name of the imputation method used, by default is pmm (predictive mean matching)
#' @param tree, boolean, if true, the phylogenetic inertia is added to the data.frame(missingData) through the computation 
#'of the eigenvectors of the tree

#' @return a data.frame of 1 or more numeric columns with the NAs replaced by values.

imputeMICE <- function(missingData, nbrMI, method = "pmm", tree = NULL){
  
  if(!is.null(tree)){
      eigen <- get_eigenvec(tree, variance_fraction = 0.8)
      missingData <- cbind(missingData[, 1:ncol(missingData), drop = FALSE], 
                           eigen[, 1:ncol(eigen), drop = FALSE])
  }
  
  ImputedMICE <- mice(missingData, m = nbrMI, method = method, maxit = 15, printFlag = FALSE)
  
  #take the columns being closer the to the mean of the variable with missing data.
  NaNColIndex <- which(apply(missingData, 2, function(x) any(is.na(x))))
  meanNaNCol <- mean(missingData[,NaNColIndex], na.rm = T)
  meanMI <- apply(ImputedMICE$imp[[NaNColIndex]], 2 , mean)
  
  imputedData <- mice::complete(ImputedMICE, 
                                action = which(abs(meanMI - meanNaNCol) == min(abs(meanMI - meanNaNCol))))
  
  return(imputedData)
}

imputeMICE(testConti, 10)


# missForest
############
#' @title non-parametric missing value imputation for mixed-type data by missForest(randomForest)
#' 
#' @description This function imputes missing data for continuous and discrete traits applying the missForest approach  
#'
#' @usage imputeMissForest(missingData, variance_fraction = 0.8, maxiter = 10, ntree = 100, 
#' mtry = sqrt(ncol(missingData)), tree = NULL)
#'
#' @param missinData data.frame of 1 or morr columns containing NAs
#' @param variance_fraction total amount of minimum variance to be represented by the eigenvectors
#' @param maxiter maximum number of iterations to be performed given the stopping criterion is not met beforehand.
#' @param ntree number of trees to grow in each forest.
#' @param mtry number of variables randomly sampled at each split. By default it's the square root of the number of 
#' variables
#' @param tree, boolean, if true, the phylogenetic inertia is added to the data.frame(missingData) through the computation 
#'of the eigenvectors of the tree

#' @return a data.frame of 1 or more numeric columns with the NAs replaced by values.

imputeMissForest <- function(missingData, variance_fraction = 0.8, maxiter = 10, ntree = 100, 
                             mtry = sqrt(ncol(missingData)), tree = NULL){
  
  Nvariables <- ncol(missingData)
  # want to include phylogeny information
  if(!is.null(tree)){
    eigen <- get_eigenvec(tree, variance_fraction = 0.8)
    missingData <- cbind(missingData[, 1:ncol(missingData), drop = FALSE], eigen[, 1:ncol(eigen), drop = FALSE])
  }

  #run missForest
  missForest_imputation <- missForest(xmis = missingData, maxiter = maxiter, ntree = ntree, mtry = mtry)
  
  #cut eigenvectors columns
  missForest_imputation$ximp <- missForest_imputation$ximp[,1:Nvariables]
  
  return(missForest_imputation$ximp)
}

imputeMissForest(test, Data$TreeList$`0`)

missForest_imputation <- missForest(xmis = test ,maxiter = 10, ntree = 100)
missForest_imputation$ximp
missForest_imputation$OOBerror


# KNN
#####
#' @title non-parametric missing value imputation for mixed-type data by kNN(k-nearest neighbor algorithm)
#' 
#' @description This function imputes missing data for continuous and discrete traits applying the kNN approach  
#'
#' @usage imputeKNN(missingData, k, numFun, catFun, variance_fraction = 0.8, tree = NULL)
#'
#' @param missinData data.frame of 1 or morr columns containing NAs
#' @param k integer, number of nearest neighbours used
#' @param numFun numFun: function for aggregating the kNN in the case of a numerical variable
#' @param catFun catFun: function for aggregating the kNN in the case of a categorical variable
#' @param variance_fraction total amount of minimum variance to be represented by the eigenvectors
#' @param tree, boolean, if true, the phylogenetic inertia is added to the data.frame(missingData) through the computation 
#'of the eigenvectors of the tree

#' @return a data.frame of 1 or more numeric columns with the NAs replaced by values.

imputeKNN <- function(missingData, k, numFun, catFun, variance_fraction = 0.8, tree = NULL){

  NbrCol <- ncol(missingData)
  
  if(!is.null(tree)){
    eigen <- get_eigenvec(tree, variance_fraction = 0.8)
    missingData <- cbind(missingData[, 1:ncol(missingData), drop = FALSE], eigen[, 1:ncol(eigen), drop = FALSE])
  }
  
  #get the column name with NaN
  variable <- c()
  for(c in 1:ncol(missingData)){
    if(any(is.na(missingData[,c]))){
      variable <- c(variable, colnames(missingData)[c])
    }
  }
  
  DataImputed <- VIM::kNN(missingData, k = k, variable = variable, numFun = numFun, catFun = catFun, imp_var = FALSE)[,1:NbrCol]
  #keep the tip names
  rownames(DataImputed) <- rownames(missingData)

  #keep only the imputed row to save some space (discuss to know if keep this or not)
  #indexKeep <- which(apply(testImputed[,(ncol(test)+1):(ncol(test)+ncol(test))], 1, any))
  #testImputed <- testImputed[indexKeep,1:ncol(test)]
  
  return(DataImputed)
}

imputedTest <- imputeKNN(test, k = 2, weightedMedian, maxCat, tree = Data$TreeList$`0`)
imputedTest
