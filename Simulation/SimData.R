#install.packages("mvMORPH")
#install.packages("phytools")
#install.packages("Matrix")
#install.packages("castor")
#install.packages("tidyverse")


library(mvMORPH)
library(phytools)
library(Matrix)
library(castor)
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
  
  if(!is.null(uncovTraits) && Ntraits > 1) {
    columns <- sample(nrow(Sigmas), uncovTraits)
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
    stop("The rate model should one of the following: ER, SYM or ARD")
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

ConvertContinousInDiscreteValues <- function(values, Nstates, intervals = TRUE){
  #made 2 replacement to don't create an extra matrix
  # values: vector of values
  # Nstates: number of states of the defined traits (could be a vector or a value)
  # intervals: boolean arguments, if TRUE want ordered states (0, 1, 2)if FALSE, means unordered (2,0,1).
  #          if TRUE intervals are fairly split, if FALSE intervals are unfairly split. 
  # return a vector of discrete value corresponding of continuous value in a same interval.
  
  #Define the thresholds in case it's ordered
  breaks <- seq(min(values), max(values), length.out = Nstates + 1)
  result <- as.character(findInterval(values, breaks[-c(1, length(breaks))]))

  if(!intervals){
    breaks <- runif((Nstates - 1), min(values), max(values))
    ordered_values <- findInterval(values, sort(breaks))
    shuffle <- sample(0:(Nstates-1), Nstates, replace = FALSE)
    result <- as.character(1:length(ordered_values)) # vector with shuffling 
    for (i in 1:length(unique(ordered_values))){
      result[ordered_values == unique(ordered_values)[i]] <- shuffle[i]
    }
  }
  return(result)
}

ChangeContinuousTraitInDiscrete <- function(Matrix, columnsIndex, Nstates, intervals){
  # apply the function below to the columns of a matrix which is converted in a dataframe.
  # intervals could be a scalar or a vector of boolean.
  # return a dataframe
  
  if(length(Nstates) >= 1 && length(Nstates) < length(columnsIndex)){
    stop("Nstates length should be equal to 1 or equal to the number of columnsIndex")
  }
  
  if(length(Nstates) == 1){
    Nstates <- rep(Nstates, length(columnsIndex))
  }

  if(length(intervals) == 1){
    intervals <- rep(intervals, length(columnsIndex))
  }
  
  #convert matrix in dataframe
  dataframe <- as.data.frame(Matrix)
  
  for(ci in 1:length(columnsIndex)) {
    #print(ConvertContinousInDiscreteValues(matrix[, columnsIndex[ci]], Nstates[ci], intervals[ci]))
    dataframe[, columnsIndex[ci]] <- ConvertContinousInDiscreteValues(dataframe[, columnsIndex[ci]], Nstates[ci], intervals[ci])
  }
  return(dataframe)
}
ContinuousData[,9]<- ConvertContinousInDiscreteValues(ContinuousData[, indexZero[1]], Nstates[1], intervals[1])
ChangeContinuousTraitInDiscrete(ContinuousData, indexZero, Nstates, intervals)

# Scaled continuous traits
#########################
#Scale traits between 0 and 10 
#(https://github.com/GitTFJ/Handling-missing-values-in-trait-data/blob/main/Script1_SimulateData_V1.0.R)
rangeScale <- function(x, range){((x-min(x))/(max(x)-min(x))*range)} #? changer on doit pouvoir d?finir un range pour chaque trait qui est diff?rent!!!!!
ScaledTraits <- apply(ContinuousData$data, MARGIN = 2, rangeScale, 10)
ScaledTraits <- as.data.frame(ScaledTraits)
ScaledTraits


#as first argument, list(Birth, Death, Ntaxa) and second arguments is a dataframe.
simData <- function(list, dataframe){
  
  if(dim(dataframe)[2] != 8){
    stop("The number of columns should be equal to 8.")
  }
  
  #Rename columns of the data frame
  ################################
  newNames <-  c("nbr_traits", "class", "subclass", "model", "states", "correlation", "uncorr_traits", "fraction_uncorr_traits")
  names(dataframe) <- newNames
  #in subclass there are, continuous, ordinal, interval(same quantity), nominative(no order).
  # uncorr_traits >= 0, 0<x<=1 fraction of uncovariant traits, number of uncorrelated should be >=1. 

  
  #Transform the columns (nbrs_traits, uncorr_traits, fraction_uncorr_traits as numeric, and rest as character)
  #######################################################################
  dataframe[,c(2:4,6)] <- apply(dataframe[,c(2:4,6)], 2, as.character)
  dataframe$nbr_traits <- as.integer(dataframe$nbr_traits) #if the value is a float, converted in integer. 
  dataframe[,c(5,7:8)] <- apply(dataframe[,c(5,7:8)], 2, as.numeric)
  
  
  #Check the data frame
  #####################
  
  #Homogenize the names
  
  dataframe$class[str_detect(dataframe$class, "^[cC]")] <- "continuous" #every class have the same name
  dataframe$class[str_detect(dataframe$class, "^[dD]")] <- "discrete"
  dataframe$subclass[str_detect(dataframe$subclass, "^[cC]")] <- "continuous"
  dataframe$subclass[str_detect(dataframe$subclass, "^[oO]")] <- "ordinal"
  dataframe$subclass[str_detect(dataframe$subclass, "^[iI]")] <- "interval"
  dataframe$subclass[str_detect(dataframe$subclass, "^[nN]")] <- "nominative"
  dataframe$model[str_detect(dataframe$model, "^[bB]")] <- "BM1"
  dataframe$model[str_detect(dataframe$model, "^[oO]")] <- "OU1"
  dataframe$model[str_detect(dataframe$model, "^[eE]")] <- "ER"
  dataframe$model[str_detect(dataframe$model, "^[sS]")] <- "SYM"
  dataframe$model[str_detect(dataframe$model, "^[aR]")] <- "ARD"
  
  
  #Check if all the rows are filled correctly
  
  wrong <- which(dataframe$class == "continuous" & dataframe$subclass == "continuous" &
                   (dataframe$model == "BM1" | dataframe$model == "OU1") & dataframe$states != 1)
  
  if(length(wrong) != 0){
    stop(paste("This line ", wrong, "is not filled correctly \n"))
  }
  
  wrong <- which(dataframe$class == "discrete" & (dataframe$subclass == "ordinal" | 
                                                    dataframe$subclass == "interval" | dataframe$subclass == "nominative") &
                   (dataframe$model == "ER" | dataframe$model == "SYM" | dataframe$model == "ARG") & dataframe$states <= 1)
  if(length(wrong) != 0){
    stop(paste("This line ", wrong, "is not filled correctly \n"))
  }
  

  # Simulate phylogeny
  ####################
  birth = list[[1]]
  death = list[[2]]
  ntaxa = list[[3]]
  
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
  max_rate <- 0.6
  for(i in correlation_values){
    #build a subset of traits being correlated together
    subdata <- subset(dataframe, dataframe$correlation == i)
    
    if(nrow(subdata) == 1){ #means no correlated with others group of traits
      Sigmas <- simSigma(subdata$nbr_traits, uncovTraits = subdata$uncorr_traits, FracNocov = subdata$fraction_uncorr_traits)
      Thetas <- runif(subdata$nbr_traits, min = -10, max = 10)
      Alphas <- simAlpha(subdata$nbr_traits)
      if(subdata$class == "continuous"){
        if(subdata$model == "BM1"){
          Alphas <- simAlpha(subdata$nbr_traits, alpha = 1e-5 * log(2))
        }
        #simulate independent continuous traits
        ContinuousData <- mvSIM(tree = SimTree,
                                model = subdata$model,
                                param = list(ntraits = subdata$nbr_traits,
                                             theta = Thetas, #theta = ancestral states
                                             sigma = Sigmas,
                                             alpha = Alphas))
        #change columns names
        ContinuousData <- as.data.frame(ContinuousData)
        colnames(ContinuousData) <- sprintf("T%s.%s", seq(1:subdata$nbr_traits), i)
        FinalData <- cbind(FinalData, ContinuousData)
      }
      else{
        Ordinal <- FALSE
        if(subdata$subclass == "ordinal"){
          Ordinal <- TRUE
        }
        
        #simulate independent discrete traits
        DiscreteData <- simDiscreteTraits(subdata$nbr_traits, 
                                          subdata$states, subdata$model, max_rate, SimTree, equal = FALSE, Ordinal)
        
        #change columns names
        DiscreteData$tip_mat <- as.data.frame(DiscreteData$tip_mat)
        colnames(DiscreteData$tip_mat) <- sprintf("T%s.%s", seq(1:subdata$nbr_traits), i)
        
        FinalData <- cbind(FinalData, DiscreteData$tip_mat)
      }
    }
    if(nrow(subdata) > 1){
      Ntraits <- sum(subdata$nbr_traits)
      model <- "BM1"
      indexZero <- NULL
      Thetas <- runif(Ntraits, min = -10, max = 10)
      Alphas <- simAlpha(Ntraits, alpha = 1e-5 * log(2)) #by default, alpha is for a Brownian motion model
      uncovTraits <- sum(subdata$uncorr_traits)
      FracNocov <- sum(subdata$nbr_traits * subdata$fraction_uncorr_traits) / Ntraits
      Sigmas <- simSigma(Ntraits, uncovTraits = subdata$uncorr_traits, FracNocov = subdata$fraction_uncorr_traits)
      
      if("OU1" %in% subdata$subclass && "BM1" %in% subdata$subclass || "discrete" %in% subdata$class){
        Alphas <- simAlpha(Ntraits)
        model <- "OU1"
        # set 0 to diagonal for BM1 and discrete
        indexZero <- sample(1:Ntraits, Ntraits - sum(subdata$nbr_trait[subdata$model == "OU1"]), replace = FALSE)
        diag(Alphas)[indexZero] <- 1e-5
      }
      
      if(sum(subdata$subclass == "OU1") == nrow(subdata)){
        Alphas <- simAlpha(Ntraits)
        model <- "OU1"
      }
      
      if(sum(subdata$class == "discrete") + sum(subdata$subclass == "BM1") == nrow(subdata)){
        indexZero <- sample(1:Ntraits, sum(subdata$nbr_trait[subdata$class == "discrete"]), replace = FALSE)
      }
      
      #Simulate data
      ContinuousData <- mvSIM(tree = SimTree,
                              model = model,
                              param = list(ntraits = Ntraits,
                                           theta = Thetas, #theta = ancestral states
                                           sigma = Sigmas,
                                           alpha = Alphas))
      
      #rescale data??
     
      if(!is.null(indexZero)){ #mean that there is some discrete correlated data
        Nstates <- rep(subdata$states[subdata$states > 1], subdata$nbr_traits[subdata$states > 1])
        intervals <- c()
        discreteSubdata <- subset(subdata, subdata$class == "discrete")
        for(j in 1:nrow(discreteSubdata)){
          if(discreteSubdata$subclass[j] == "interval"){intervals <- c(intervals,rep(as.logical(TRUE), discreteSubdata$nbr_traits[j]))}
          else{intervals <- c(intervals,rep(as.logical(FALSE),discreteSubdata$nbr_traits[j]))}
        }
        DiscreteData <- ChangeContinuousTraitInDiscrete(ContinuousData, indexZero, Nstates, intervals)
        Nstates
        indexZero
        intervals
        #change columns names
        DiscreteData <- as.data.frame(DiscreteData)
        colnames(DiscreteData) <- sprintf("T%s.%s", seq(1:sum(subdata$nbr_traits)), i)
        FinalData <- cbind(FinalData, DiscreteData)
      }
      else{
        #change columns names
        ContinuousData <- as.data.frame(ContinuousData)
        colnames(ContinuousData) <- sprintf("T%s.%s", seq(1:sum(subdata$nbr_traits)), i)
        FinalData <- cbind(FinalData, ContinuousData)
      }
    }
  }#close for loop
  
  
  #Define list of object
  ######################
  #Data <- list(TotalData = FinalData, )
  
  #Save data
  ##########
  #save(Data, file="simulatedData.RData")
  
  return(FinalData)
    
} #close function

data <- read.csv("DataTest.csv", header = TRUE, sep = ";")
data
tree_arg <- list(Birth = 0.4, Death = 0.1, Ntaxa = 30)
list <- list(Birth = 0.4, Death = 0.1, Ntaxa = 30)
new_data <- simData(tree_arg, data)
new_data
