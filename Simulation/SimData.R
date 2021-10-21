#install.packages("mvMORPH")
#install.packages("phytools")
#install.packages("Matrix")
#install.packages("castor")
#install.packages("geiger")
#install.packages("tidyverse")


library(mvMORPH)
library(phytools)
library(Matrix)
library(castor)
library(geiger)
library(tidyverse)

#' @title Simulate variance-covariance matrix for morphological evolution
#'
#' @description This function generates a variance-covariance matrix for
#' morphological evolution, defining the rate at which traits evolve and their
#' evolutionary correlation.
#'
#' @usage simSigma(Ntraits, Cor = NULL, Sigma2 = NULL, uncovTraits = NULL, FracNocov = NULL)
#'
#' @param Ntraits number of traits
#' @param Cor correlation between traits. Default between -1 and 1.
#' Optional, can be fixed to be equal between all traits by giving one value or
#' Ntraits*(Ntraits-1)/2
#' @param Sigma2 Brownian motion rate. Default between 1e-4 and 0.5.
#' Optional, can be fixed to be equal for all traits or one value per trait
#' @param uncovTraits number of traits having a covariance equal to 0 with the others.
#' @param FracNocov fraction of covariance being zero
#' 
#' @return matrix Ntrait x Ntrait for simulating trait evolution
#' 
simSigma <- function(Ntraits, Cor = NULL, Sigma2 = NULL, uncovTraits = NULL, FracNocov = NULL){
  if(!is.null(Sigma2)){
    if(length(Sigma2) != Ntraits && length(Sigma2) != 1){
      stop("Sigma2 should be of length 1 or Ntraits")
    }else if(length(Sigma2) == 1) {
      Sigma2 <- rep(Sigma2, Ntraits)
    }
  }
  else{
    Sigma2 <- runif(Ntraits, min = 1e-4, max = 0.5)
  }
  
  Sigmas <- matrix(Sigma2, nrow = 1)
  if(Ntraits > 1) {
    Cov <- matrix(1, ncol = Ntraits, nrow = Ntraits)
    Q <- Ntraits*(Ntraits-1)/2 #for me it's +1 and not -1
    if(!is.null(Cor)) {
      if(length(Cor) != Q && length(Cor) != 1) {
        stop("Correlation among traits should be of length 1 or Ntraits*(Ntraits-1)/2")
      }
      SimCov <- Cor
    }
    else{
      SimCov <- runif(Q, min = -1, max = 1) # Trait correlation
    }
    Cov[lower.tri(Cov, diag = FALSE)] <- SimCov
    Cov[upper.tri(Cov, diag = FALSE)] <- SimCov
    Sigmas <- diag(Sigma2)  %*% Cov  %*% diag(Sigma2) # Correlation to covariance 
    
    # Force variance-covariance to be positive definite # can I have some precision?
    Tol <- 1e-6
    Ev <- eigen(Sigmas, symmetric = TRUE)$values
    if(!all( Ev >= -Tol * abs(Ev[1L]))) {
      Sigmas <- as.matrix(nearPD(Sigmas)$mat)
    }
  }
  
  if(!is.null(uncovTraits) && uncovTraits < 0){
    stop("UncovTraits should be equal or greater than 0")
  }
  
  if(!is.null(uncovTraits) && uncovTraits %% 1 != 0){
    stop("UncovTraits should be an integer")
  }
  
  if(length(uncovTraits) == 1 && uncovTraits == 0){
    uncovTraits = NULL
  }
  
  if(!is.null(uncovTraits) && uncovTraits > Ntraits){
    stop("Uncovtraits should be equal or smaller than Ntraits")
  }
  
  if(!is.null(uncovTraits) && Ntraits > 1) {
    columns <- 1:uncovTraits
    for(i in columns){
      valueToKeep <- Sigmas[i, i]
      Sigmas[i,] <- 0
      Sigmas[,i] <- 0
      Sigmas[i, i] <- valueToKeep
    }
  }
  
  if(length(FracNocov) == 1 && FracNocov == 0){
    FracNocov = NULL
  }
  
  if(!is.null(FracNocov) && Ntraits > 1) {
    Q <- Ntraits*(Ntraits-1)/2
    Nzero <- round(Q * FracNocov) #to have the number of variable having no covariance.
    Mask <- matrix(1, nrow = Ntraits, ncol = Ntraits)
    Lt <- which(lower.tri(Mask))
    SetZero <- sample(Lt, Nzero)
    Mask[SetZero] <- 0
    Mask <- t(Mask)
    Mask[SetZero] <- 0
    Sigmas <- Sigmas * Mask
  }
  return(Sigmas)
}
simSigma(3, uncovTraits = 0)
simSigma(5, uncovTraits = 0, FracNocov = 0.6)
simSigma(5, FracNocov = 0.6)
simSigma(5, uncovTraits = 5)

#' @title Simulate alpha matrix for morphological evolution
#'
#' @description This function generates an alpha matrix for
#' morphological evolution, defining the strength of attraction 
#' towards a morphological optimum
#'
#' @usage simSigma(Ntraits, alpha = NULL)
#'
#' @param Ntraits number of traits
#' @param alpha strength of attraction towards morphological optimum. Default between 0.5 and 2.
#' Optional, can be fixed to be equal between all traits by giving one value or
#' Ntraits*(Ntraits-1)/2
#' 
#' @return matrix Ntrait x Ntrait for simulating trait evolution
#'
simAlpha <- function(Ntraits, alpha = NULL) {
  if(is.null(alpha)){
    alpha <- exp(runif(Ntraits, log(0.5), log(2)))
    TreeHeight <- 1 # In case we need to change the age of the phylogenies
    alpha <- alpha * log(2) * 1/TreeHeight 
    AlphaMat <- diag(alpha, nrow = Ntraits, ncol = Ntraits)
  }
  else {
    if(length(alpha) != Ntraits && length(alpha) != 1) {
      stop("alpha should be of length 1 or Ntraits")
    }
    if(any(alpha < 0)) {
      stop("Alpha should be larger or equal to 0")
    }
    if(length(alpha) == 1) {
      AlphaMat <- diag(rep(alpha, Ntraits), nrow = Ntraits, ncol = Ntraits)
    }
    else{
      AlphaMat <- diag(alpha, nrow = Ntraits, ncol = Ntraits)
    }
  }
  return(AlphaMat)
}

#Simulate discrete traits
#########################
#' @title Simulate discrete traits for all species according to the MK model
#'
#' @description This function generates matrix of discrete values (states) of
#' morphological evolution. Define also a matrix of nodes
#'
#' @usage simDiscreteTraits(Ntraits, Nstates, rate_model, max_rate, SimTree)
#'
#' @param Ntraits number of traits 
#' @param Nstates number of states of the defined traits (could be a vector or a value)
#' @param rate_model vector, rate model that the transition matrix must satisfy. Chose between "ER" (=all permitted transitions occur at the same rate), "SYM" (=hat backward & forward transitions occur at the same rate) and "ARD" (=all allowed transitions can occur at different rates). 
#' @param tree stochastic birth-death trees
#' @param equal if TRUE, the transition matrix is equal for all the traits. If FALSE all the traits have a different transition matrix. By default is TRUE. 
#' @param Ordinal simulate ordinal data. Default is FALSE
#' @return list containing morphological evolution for each trait and each species, a matrix Ntrait x Ntrait for simulating trait evolution and sigma matrix Ntraits x Ntraits
#' 

simDiscreteTraits <- function(Ntraits, Nstates, rate_model, max_rate, tree, equal = TRUE, Ordinal = FALSE){
  
  if(rate_model != "ER" & rate_model != "SYM" & rate_model != "ARD"){
    stop("The rate model should be one of the following: ER, SYM or ARD")
  }
  
  if(length(rate_model) != 1 && length(rate_model) != Ntraits){
    stop("rate_model should be of length 1 or Nstates")
  }
  
  if(length(rate_model) == 1 && length(rate_model) != Ntraits) {
    rate_model <- rep(rate_model, Ntraits)
  }
  
  if(length(Nstates) != 1 && length(Nstates) != Ntraits){
    stop("Nstates should be of length 1 or Nstates")
  }
  
  if(length(Nstates) == 1 && length(Nstates) != Ntraits) {
    Nstates <- rep(Nstates, Ntraits)
  }
  
  tip_mat <- matrix(0, nrow = length(tree$tip.label), ncol = Ntraits)
  node_mat <- matrix(0, nrow = tree$Nnode, ncol = Ntraits)
  for(i in 1:Ntraits) {
    # define transition matrix
    if(equal){
      set.seed(1)
      Q <- get_random_mk_transition_matrix(Nstates[i], rate_model = rate_model[i], max_rate = max_rate)
    }
    else{
      Q <- get_random_mk_transition_matrix(Nstates[i], rate_model = rate_model[i], max_rate = max_rate)
    }
    if(Ordinal){ #going to be have more rate of change which are equal to 0.
      for (y in 2:Nstates[i]) {
        Q[row(Q) == (col(Q) - y)] <- 0
        Q[row(Q) == (col(Q) + y)] <- 0
      }
    }
    tip_states <- simulate_mk_model(tree, Q)
    tip_mat[, i] <- tip_states$tip_states
    node_mat[, i] <- tip_states$node_states # Do we need the ancestral states?
  }
  return(list(tip_mat = tip_mat, node_mat = node_mat))
}
simDiscreteTraits(4, c(2,3,3,4), "ER" ,0.9, SimTree, equal = F, Ordinal = F)

# Simulate correlated discrete traits with continuous traits
############################################################

ConvertContinousInDiscreteValues <- function(values, Nstates, subclass){
  #made 2 replacement to don't create an extra matrix
  # values: vector of values
  # Nstates: number of states of the defined traits (could be a vector or a value)
  # subclass:
  # - intervals: states are fairly split
  # - ordinal: ordered states, split is random
  # - nominal_no: split is random, no order (shuffled)
  # - nominal_o: 
  # return a vector of discrete value corresponding of continuous value in a same interval.
  
  if(subclass != "non_eq_nominal" & subclass != "ordinal" & subclass != "interval" & subclass != "eq_nominal"){
    stop("The subclass should be one of the following: nominal, ordinal or interval")
  }
  
  if(subclass == "interval"){
  breaks <- seq(min(values), max(values), length.out = Nstates + 1)
  conversion <- as.character(findInterval(values, breaks[-c(1, length(breaks))]))
  }
  
  if(subclass == "ordinal"){
    breaks <- runif((Nstates - 1), min(values), max(values))
    conversion <- as.character(findInterval(values, sort(breaks)))
  }
  
  if(subclass == "nominal_neq"){
    breaks <- sample(sort(values, Nstates-1))
    nominal_values <- findInterval(values, sort(breaks))
    shuffle <- sample(0:(Nstates-1), Nstates, replace = FALSE)
    conversion <- as.character(1:length(nominal_values)) # vector with shuffling 
    for (i in 1:length(unique(nominal_values))){
      conversion[nominal_values == unique(nominal_values)[i]] <- shuffle[i]
    }
  }
  
  if(subclass == "nominal_eq"){
    breaks <- seq(min(values), max(values), length.out = Nstates + 1)
    nominal_eq_values <- findInterval(values, breaks[-c(1, length(breaks))])
    shuffle <- sample(0:(Nstates-1), Nstates, replace = FALSE)
    conversion <- as.character(1:length(nominal_eq_values)) # vector with shuffling 
    for (i in 1:length(unique(nominal_eq_values))){
      conversion[nominal_eq_values == unique(nominal_eq_values)[i]] <- shuffle[i]
    }
  }
  
  return(conversion)
}

ChangeContinuousTraitInDiscrete <- function(Matrix, columnsIndex, Nstates, subclass){
  # apply the function below to the columns of a matrix which is converted in a dataframe.
  # intervals could be a scalar or a vector of boolean.
  # return a dataframe
  
  if(length(Nstates) >= 1 && length(Nstates) < length(columnsIndex)){
    stop("Nstates length should be equal to 1 or equal to the number of columnsIndex")
  }
  
  if(length(Nstates) == 1){
    Nstates <- rep(Nstates, length(columnsIndex))
  }

  if(length(subclass) == 1){
    subclass <- rep(subclass, length(columnsIndex))
  }
  
  #convert matrix in dataframe
  dataframe <- as.data.frame(Matrix)
  
  for(ci in 1:length(columnsIndex)) {
    dataframe[, columnsIndex[ci]] <- ConvertContinousInDiscreteValues(dataframe[, columnsIndex[ci]], Nstates[ci], subclass[ci])
  }
  return(dataframe)
}

# # Scaled continuous traits
# #########################
# #Scale traits between 0 and 10 
# #(https://github.com/GitTFJ/Handling-missing-values-in-trait-data/blob/main/Script1_SimulateData_V1.0.R)
# rangeScale <- function(x, range){((x-min(x))/(max(x)-min(x))*range)}
# ScaledTraits <- apply(ContinuousData$data, MARGIN = 2, rangeScale, 10)
# ScaledTraits <- as.data.frame(ScaledTraits)
# ScaledTraits


#as first argument, param_tree(Birth, Death, Ntaxa), second arguments is a dataframe, third argument save in rds format
simData <- function(param_tree, dataframe, save = TRUE){
  
  if(ncol(dataframe) != 8){
    stop("The number of columns should be equal to 8.")
  }
  
  #Rename columns of the data frame
  ################################
  newNames <-  c("nbr_traits", "class", "model", "states", "correlation", "uncorr_traits", "fraction_uncorr_traits", "lambda")
  names(dataframe) <- newNames
  #in subclass there are, continuous, ordinal, interval(same quantity), nominal(no order).
  # uncorr_traits >= 0, 0<x<=1 fraction of uncovariant traits among the correlated group of traits, number of uncorrelated should be >=1. 

  
  #Transform the columns (nbrs_traits as integer, uncorr_traits, fraction_uncorr_traits and lambda as numeric, and rest as character)
  #######################################################################
  dataframe[,c(2:3,5)] <- apply(dataframe[,c(2:3,5)], 2, as.character)
  dataframe$nbr_traits <- as.integer(dataframe$nbr_traits) #if the value is a float, converted in integer. 
  dataframe[,c(4,6:8)] <- apply(dataframe[,c(4,6:8)], 2, as.numeric)
  
  
  #Check the data frame
  #####################
  
  #Homogenize the names
  dataframe$class[str_detect(dataframe$class, "^[cC]")] <- "continuous"
  dataframe$class[str_detect(dataframe$class, "^[oO]")] <- "ordinal"
  dataframe$class[str_detect(dataframe$class, "^[iI]")] <- "interval"
  dataframe$class[str_detect(dataframe$class, "^[nN]")] <- "non_eq_nominal"
  dataframe$class[str_detect(dataframe$class, "^[eE]")] <- "eq_nominal"
  dataframe$model[str_detect(dataframe$model, "^[bB]")] <- "BM1"
  dataframe$model[str_detect(dataframe$model, "^[oO]")] <- "OU1"
  dataframe$model[str_detect(dataframe$model, "^[eE]")] <- "ER"
  dataframe$model[str_detect(dataframe$model, "^[sS]")] <- "SYM"
  dataframe$model[str_detect(dataframe$model, "^[aR]")] <- "ARD"
  
  
  
  #Check if all the rows are filled correctly
  
  wrong <- which(dataframe$class == "continuous" &
                   (dataframe$model == "BM1" | dataframe$model == "OU1") & dataframe$states != 1)
  if(length(wrong) != 0){
    stop(paste("This line ", wrong, "is not filled correctly \n"))
  }
  
  wrong <- which(dataframe$class != "continuous" & dataframe$states <= 1)
  if(length(wrong) != 0){
    stop(paste("This line ", wrong, "is not filled correctly \n"))
  }
  
  wrong <- which(dataframe$lambda < 0 & dataframe$lambda > 1)
  if(length(wrong) != 0){
    stop(paste("This line ", wrong, "is not filled correctly \n"))
  }
  
  # Simulate phylogeny
  ####################
  birth = param_tree[[1]]
  death = param_tree[[2]]
  ntaxa = param_tree[[3]]
  
  # Simulating phylogenies fails sometimes. Try until we are successful
  Extant <- FALSE
  while (!Extant) {
    SimTree <- pbtree(b = birth, d = death, n = ntaxa, scale = 1, extant.only = FALSE) #why added extant.only ?
    if(!is.null(SimTree)) {
      Extant <- TRUE
    }
  }
  
  
  #Simulate by correlation
  ########################
  
  #Get the various group of correlated or not traits.
  correlation_values <- unique(dataframe$correlation)
  
  #Final matrix
  FinalData <- as.data.frame(matrix(0, nrow = length(SimTree$tip.label), ncol = 0))
  max_rate <- 0.5
  AlphasList <- list()
  ThetasList <- list()
  SigmasList <- list()
  TreeList <- list()
  for(i in correlation_values){

    subdataTree <- SimTree
    #build a subset of traits being correlated together
    subdata <- subset(dataframe, dataframe$correlation == i)
    
    #rescale phylogeny
    lambdaCheck <- mean(subdata$lambda)
    if(lambdaCheck != 1){
      subdataTree <- rescale(SimTree, "lambda", lambdaCheck)
    }
    
    #check fraction of uncorrelated
    if(sum(subdata$fraction_uncorr_traits) != subdata$fraction_uncorr_traits[1] * nrow(subdata)){
      stop("The fraction of uncorrelated trait should be equal within the same group of correlated traits")
    }
    
    if(nrow(subdata) == 1){ #means no correlated with others group of traits
      Sigmas <- simSigma(subdata$nbr_traits, uncovTraits = subdata$uncorr_traits, FracNocov = subdata$fraction_uncorr_traits)
      Thetas <- runif(subdata$nbr_traits, min = -10, max = 10)
      Alphas <- simAlpha(subdata$nbr_traits)
      if(subdata$class == "continuous"){
        if(subdata$model == "BM1"){
          Alphas <- simAlpha(subdata$nbr_traits, alpha = 1e-5 * log(2))
        }
        #simulate independent continuous traits
        ContinuousData <- mvSIM(tree = subdataTree,
                                model = subdata$model,
                                param = list(ntraits = subdata$nbr_traits,
                                             theta = Thetas, #theta = ancestral states
                                             sigma = Sigmas,
                                             alpha = Alphas))
        #change columns names
        ContinuousData <- as.data.frame(ContinuousData)
        colnames(ContinuousData) <- sprintf("F%s.%s", seq(1:subdata$nbr_traits), i)
        FinalData <- cbind(FinalData, ContinuousData)
      }
      else{
        Ordinal <- FALSE
        if(subdata$class == "ordinal"){
          Ordinal <- TRUE
        }
        
        #simulate independent discrete traits
        DiscreteData <- simDiscreteTraits(subdata$nbr_traits, 
                                          subdata$states, subdata$model, max_rate, subdataTree, equal = FALSE, Ordinal)
        
        #change columns names
        DiscreteData$tip_mat <- as.data.frame(DiscreteData$tip_mat)
        colnames(DiscreteData$tip_mat) <- sprintf("I%s.%s", seq(1:subdata$nbr_traits), i)
        
        FinalData <- cbind(FinalData, DiscreteData$tip_mat)
      }
    }
    if(nrow(subdata) > 1){
      #check if the dataframe is correctly filled
      wrong <- which(subdata$model != "BM1" & subdata$model != "OU1" | subdata$lambda != lambdaCheck)
      if(length(wrong) != 0){
        stop(paste("This line ", wrong, "is not filled correctly \n"))
      }
      
      Ntraits <- sum(subdata$nbr_traits)
      model <- "BM1"
      Thetas <- runif(Ntraits, min = -10, max = 10)
      Alphas <- simAlpha(Ntraits, alpha = 1e-5 * log(2)) #by default, alpha is for a Brownian motion model
      uncorrTraits <- sum(subdata$uncorr_traits)
      #FracNocov <- sum(subdata$nbr_traits * subdata$fraction_uncorr_traits) / Ntraits
      Sigmas <- simSigma(Ntraits, uncovTraits = uncorrTraits, FracNocov = uncorrFractionCheck)
      
      if("OU1" %in% subdata$model && "BM1" %in% subdata$model){
        Alphas <- simAlpha(Ntraits)
        model <- "OU1"
        
        # in case uncovTraits simulated with BM
        indexUncorrTraitsBM <- NULL
        uncorrTraitsBM <- sum(subdata$uncorr_traits[subdata$model == "BM1"])
        uncorrTraitsOU <- uncorrTraits - uncorrTraitsBM #number of uncorrelated traits simulated with OU
        if(uncorrTraitsBM != 0){
          indexUncorrTraitsBM <- 1:uncorrTraitsBM
          indexBM <-c(indexUncorrTraitsBM, sample((1+uncorrTraits):Ntraits, (Ntraits - uncorrTraitsBM - 
                                                    sum(subdata$nbr_trait[subdata$model == "OU1"])), replace = FALSE))
        }
        if(uncorrTraitsBM == 0){
          #in case no uncovTraits in the BM model.
          indexBM <- sample((1 + uncorrTraitsOU):Ntraits, Ntraits - sum(subdata$nbr_trait[subdata$model == "OU1"]), replace = FALSE)
        }
        # set 0 to diagonal for BM1
        diag(Alphas)[indexBM] <- 1e-5
      }
      
      if(sum(subdata$model == "OU1") == nrow(subdata)){
        Alphas <- simAlpha(Ntraits)
        model <- "OU1"
      }

      #Simulate data
      ContinuousData <- mvSIM(tree = subdataTree,
                              model = model,
                              param = list(ntraits = Ntraits,
                                           theta = Thetas, #theta = ancestral states
                                           sigma = Sigmas,
                                           alpha = Alphas))

      
      
      #If correlated discrete traits with continuous traits
      discreteSubdata <- subset(subdata, subdata$class != "continuous")
      if(nrow(discreteSubdata) >= 1){
        Nstates <- rep(discreteSubdata$states, discreteSubdata$nbr_traits)
        class <- rep(discreteSubdata$class, discreteSubdata$nbr_traits)
        
        indexColumConvert <- c()
        for(r in 1:nrow(discreteSubdata)){
          if(discreteSubdata$model[r] == "BM1" & discreteSubdata$uncorr_traits[r] != 0){
            indexColumConvert <- c(indexColumConvert, tail(indexUncorrTraitsBM, n = discreteSubdata$uncorr_traits[r]))
            indexUncorrTraitsBM <- indexUncorrTraitsBM[-tail(indexUncorrTraitsBM, n = discreteSubdata$uncorr_traits[r])] #delete the index selected above
          }
          indexUncorrTraitsOU <- (1 + uncorrTraitsBM):uncorrTraits
          if(discreteSubdata$model[r] == "OU1" & discreteSubdata$uncorr_traits[r] != 0){
            indexColumConvert <- c(indexColumConvert, tail(indexUncorrTraitsOU, n = discreteSubdata$uncorr_traits[r]))
            indexUncorrTraitsOU <- indexUncorrTraitsOU[-tail(indexUncorrTraitsOU, n = discreteSubdata$uncorr_traits[r])] #delete the index selected above
          }
        }
        BMTraitsRemained <- sample(indexBM[(uncorrTraitsBM + 1):length(indexBM)], 
                                   sum(discreteSubdata$nbr_traits[discreteSubdata$model == "BM1"]) - uncorrTraitsBM)
        
        if(length(indexBM[(uncorrTraitsBM + 1):length(indexBM)]) == 1){
          BMTraitsRemained <- indexBM[(uncorrTraitsBM + 1):length(indexBM)]
        }
        
        indexOU <- setdiff(1:Ntraits, indexBM)
        OUTraitsRemained <- sample(indexOU[(uncorrTraitsOU+1):length(indexOU)], 
                                   sum(discreteSubdata$nbr_traits[discreteSubdata$model == "OU1"]) - uncorrTraitsOU)
        
        if(length(indexOU[(uncorrTraitsOU+1):length(indexOU)]) == 1){
          OUTraitsRemained <- indexOU[(uncorrTraitsOU+1):length(indexOU)]
        }
        
        indexColumConvert <- sort(c(indexColumConvert, BMTraitsRemained, OUTraitsRemained))
        DiscreteData <- ChangeContinuousTraitInDiscrete(ContinuousData, indexColumConvert, Nstates, class)
        
        #change columns names
        colnames(DiscreteData) <- sprintf("F%s.%s", seq(1:sum(subdata$nbr_traits)), i)
        colnames(DiscreteData)[indexColumConvert] <- sprintf("I%s.%s", seq(1:sum(discreteSubdata$nbr_traits)), i)
        FinalData <- cbind(FinalData, DiscreteData)
      }

      else{
        #change columns names
        ContinuousData <- as.data.frame(ContinuousData)
        colnames(ContinuousData) <- sprintf("F%s.%s", seq(1:sum(subdata$nbr_traits)), i)
        FinalData <- cbind(FinalData, ContinuousData)
      }
    }
    #Save parameters in an own list
    AlphasList[[i]] <- Alphas
    ThetasList[[i]] <- Thetas
    SigmasList[[i]] <- Sigmas
    TreeList[[i]] <- subdataTree
    
     
  }#close for loop
  
  FinalDiscreteData <- FinalData %>% select(starts_with("I"))
  FinlaContinuousData <- FinalData %>% select(starts_with("F"))
  
  #Define list of object
  ######################
  Data <- list(FinalData = FinalData, ContinuousData = FinlaContinuousData, DiscreteData = FinalDiscreteData, 
               AlphaMatrices = AlphasList, Thetas = ThetasList, SigmaMatrices = SigmasList, TreeList = TreeList, 
               PhyloParam = param_tree, dataframe = dataframe)
  
  #Save data
  ##########
  if(save){
  save(Data, file="simulatedData.RData")
  }
  
  return(Data)
    
} #close function

data <- read.csv("DataTest.csv", header = TRUE, sep = ";")
tree_arg <- list(Birth = 0.4, Death = 0.1, Ntaxa = 30)
new_data <- simData(tree_arg, data, save = FALSE)

