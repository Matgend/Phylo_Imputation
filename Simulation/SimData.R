library(mvMORPH)
library(phytools)
library(Matrix)
library(castor)

#' @title Simulate variance-covariance matrix for morphological evolution
#'
#' @description This function generates a variance-covariance matrix for
#' morphological evolution, defining the rate at which traits evolve and their
#' evolutionary correlation.
#'
#' @usage simSigma(Ntraits, Cor = NULL, Sigma2 = NULL)
#'
#' @param Ntraits number of traits
#' @param Cor correlation between traits. Default between -1 and 1.
#' Optional, can be fixed to be equal between all traits by giving one value or
#' Ntraits*(Ntraits-1)/2
#' @param Sigma2 Brownian motion rate. Default between 1e-4 and 0.5.
#' Optional, can be fixed to be equal for all traits or one value per trait
#' @param uncovTraits number of traits having a covariance equal to 0 with the others.
#' 
#' @return matrix Ntrait x Ntrait for simulating trait evolution
#' 
simSigma <- function (Ntraits, Cor = NULL, Sigma2 = NULL, uncovTraits = NULL){
  if (!is.null(Sigma2)){
    if (length(Sigma2) != Ntraits && length(Sigma2) != 1){
      stop("Sigma2 should be of length 1 or Ntraits")
    }else if(length(Sigma2) == 1){
      Sigma2 <- rep(Sigma2, Ntraits)
    }
  }
  else {Sigma2 <- runif(Ntraits, min = 1e-4, max = 0.5)} #why max = 0.5?
  
  Sigmas <- matrix(Sigma2, nrow = 1)
  if (Ntraits > 1) {
    Cov <- matrix(1, ncol = Ntraits, nrow = Ntraits)
    Q <- Ntraits*(Ntraits-1)/2 #for me it's +1 and not -1
    if (!is.null(Cor)) {
      if (length(Cor) != Q && length(Cor) != 1) {
        stop("Correlation among traits should be of length 1 or Ntraits*(Ntraits-1)/2")
      }
      SimCov <- Cor
    }
    else {
      SimCov <- runif(Q, min = -1, max = 1) # Trait correlation
    }
    Cov[lower.tri(Cov, diag = FALSE)] <- SimCov
    Cov[upper.tri(Cov, diag = FALSE)] <- SimCov
    Sigmas <- diag(Sigma2)  %*% Cov  %*% diag(Sigma2) # Correlation to covariance 
    
    # Force variance-covariance to be poisitive definite # can I have some precision?
    Tol <- 1e-6
    Ev <- eigen(Sigmas, symmetric = TRUE)$values
    if (!all( Ev >= -Tol * abs(Ev[1L]))) {
      Sigmas <- as.matrix(nearPD(Sigmas)$mat)
    }
  }
  
  if (!is.null(uncovTraits)){
    columns <- floor(runif(uncovTraits, 0, dim(Sigmas)[1]))
    for(i in columns){
      valueToKeep <- Sigmas[i, i]
      Sigmas[i,] <- 0
      Sigmas[,i] <- 0
      Sigmas[i,i] <- valueToKeep
    }
  }
  return(Sigmas)
}
simSigma(3, uncovTraits = 1)
#' @title Simulate variance-covariance matrix for all species according to a chosen model for countinuous traits.
#'
#' @description This function generates a variance-covariance matrix for continuous
#' morphological evolution, defining the rate at which traits evolve, their
#' evolutionary correlation and the model of evolution. Define also a matrix for simulating trait evolution
#'
#' @usage simContinuousData(Ntraits, Cor = NULL, Sigma2 = NULL, Birth = 0, Death = 0, Ntaxa = 2, model = "BM", theta = NULL, alpha = NULL)
#'
#' @param Ntraits number of continuous traits
#' @param Sigmas variance-covariance matrix for morphological evolution, defining the 
#' rate at which traits evolve and their evolutionary correlation.
#' @param Tree stochastic birth-death trees
#' @param Ntaxa number of taxa
#' @param model string precising the evolutionnary model (Brownian, Ornstein Uhlenbeck or MC)
#' @param theta ancestral state (selective optima) ~ long term mean, assume that species traits evolve around this value.
#' @param alpha strength of evolutionnary force that return traits back forward to the long-term mean mu.
#' 
#' @return list containing morphological evolution for each trait and each species, a matrix Ntrait x Ntrait for simulating trait evolution and sigma matrix Ntraits x Ntraits
#' 
# Get matrix of evolutionary rate (diagonal) and covariance (off-diagonals)
###########################################################################
Sigmas <- simSigma(Ntraits, Cor =  Cor, Sigma2 =  Sigma2, uncovTraits =  uncovTraits)
simContinuousData <- function(Ntraits, Sigmas, tree, Ntaxa = 2, model = "BM1", theta = NULL, alpha = NULL){

  # Simulate continuous traits under BM or OU model
  #################################################
  
  if (model == "OU1"){
    
    if (!is.null(theta) && is.null(alpha)){
      if (length(theta) != Ntraits && length(theta) != 1){
        stop("Theta should be of length 1 or Ntraits")}
      else if (length(theta == 1)){
        theta <- rep(theta, Ntraits)}
      alpha <- diag(rexp(Ntraits,rate = 1), nrow = Ntraits, ncol = Ntraits)
    }
    
    else if (is.null(theta) && !is.null(alpha)){
      if (length(alpha) != Ntraits && length(alpha) != 1){
        stop("alpha should be of length 1 or Ntraits")}
      else if (alpha < 0){
        stop("Alpha should be larger or equal to 0")}
      else if (length(alpha) == 1){
        alpha <- diag(rep(alpha, Ntraits), nrow = Ntraits, ncol = Ntraits)}
      theta <- runif(Ntraits, min = -10, max = 10)
    }
    else{
      alpha <- diag(exp(runif(Ntraits,log(0.5),log(2))), nrow = Ntraits, ncol = Ntraits) #should be change with ln(2)/phylogeny half life?
      theta <- runif(Ntraits, min = -10, max = 10)
    }
  }
  
  else if (model == "BM1"){
    theta <- rep(0, Ntraits)
  }
  else{
    stop("The model should be BM1 or OU1")
  }
  
  #Simulate continuous traits
  #################
  SimTraits <- mvSIM(tree = tree, model = model,
                     param = list(ntraits = Ntraits,
                                  theta = theta, #theta = ancestral states
                                  sigma = Sigmas,
                                  alpha = alpha))
  
  
  # #"BMM" for a multi-rate and multi-selective regimes, and "BM1" for a unique rate of evolution per trait.
  # # if want to have a correlation matrix in return
  # if (model == "OU1"){
  # # Fit Ornstein-Uhlenbeck (OU) model
  # ###################################
  # FitModel <- mvOU(tree = tree,
  #               data = SimTraits,
  #               model = "OU1",
  #               method = "rpf",
  #               diagnostic = FALSE,
  #               echo = FALSE)
  # }else{
  # # Fit Brownian motion model
  # ###########################
  # FitModel <- mvBM(tree = tree,
  #               data = SimTraits,
  #               model = "BM1",
  #               method = "rpf",
  #               diagnostic = FALSE,
  #               echo = FALSE)
  # }
  #return(list(data = SimTraits, varMatrix = cov2cor(FitModel$sigma), sigma = Sigmas))
  return(list(data = SimTraits, sigma = Sigmas))
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
#' @param rate_model rate model that the transition matrix must satisfy. Chose between "ER" (=all permitted transitions occur at the same rate), "SYM" (=hat backward & forward transitions occur at the same rate) and "ARD" (=all allowed transitions can occur at different rates)
#' @param tree stochastic birth-death trees
#' @param equal if TRUE, the transition matrix is equal for all the traits. If FALSE all the traits have a different transition matrix. By default is TRUE. 
#' @return list containing morphological evolution for each trait and each species, a matrix Ntrait x Ntrait for simulating trait evolution and sigma matrix Ntraits x Ntraits
#' 

simDiscreteTraits <- function(Ntraits, Nstates, rate_model, max_rate, tree, equal = TRUE){

  if (rate_model != "ER" & rate_model != "SYM" & rate_model != "ARD"){
    stop("The rate model should one of the following: ER, SYM or ARD")
  }
  
  if (length(Nstates) != 1 && length(Nstates) != Ntraits){
    stop("Ntates should be of length 1 or Nstates")
  }
  
  else if (length(Nstates) == 1 && length(Nstates) != Ntraits){
    Nstates <- rep(Nstates, Ntraits)
  }
  
  tip_mat <- matrix(0, nrow = length(tree$tip.label), ncol = Ntraits)
  node_mat <- matrix(0, nrow = tree$Nnode, ncol = Ntraits)
  for (i in 1:Ntraits){
    # define transition matrix
    if(equal){
      set.seed(1)
      Q <- get_random_mk_transition_matrix(Nstates[i], rate_model=rate_model, max_rate = max_rate)}
    
    else{Q <- get_random_mk_transition_matrix(Nstates[i], rate_model=rate_model, max_rate = max_rate)}
    tip_states <- simulate_mk_model(tree, Q)
    tip_mat[, i] <- tip_states$tip_states
    node_mat[, i] <- tip_states$node_states
  }
  return(list(tip_mat = tip_mat, node_mat = node_mat))
}

# Simulate correlated discrete traits with continuous traits
############################################################

ConvertContinousInDiscreteValues <- function(values, Nstates, threshold){
  #made 2 replacement to don't create an extra matrix
  # values: vector of values
  # Nstates: number of states of the defined traits (could be a vector or a value)
  # threshold: vector or scalar of values defining the separation value for the discrete transformation
  # return a vector of discrete value corresponding of continuous value in a same interval.
  
  
  #generate a vector of length Nstates of values larger than the max of the vector
  vectorMax <- max(values)
  transitionVec <- runif(Nstates, vectorMax + 1, vectorMax + 2)
  
  
  if(Nstates != length(threshold) + 1){
    stop("The number of threshold should be equal to Nstates - 1")
  }
  
  Nstates <- 1:Nstates
  
  values[values < threshold[1]] <- transitionVec[1]
  values[values >= tail(threshold, n= 1) & values < transitionVec[1]] <- tail(transitionVec, n= 1)

  # replace values by the Nstates values
  for (i in 2:length(Nstates) - 1){values[values >= threshold[i-1] & values < threshold[i]] <- Nstates[i]}
  
  values[values == transitionVec[1]] <- Nstates[1]
  values[values == tail(transitionVec, n= 1)] <- tail(Nstates, n= 1)
  
  return(values)
}

ChangeContinuousTraitInDiscrete <- function(matrix, columnsIndex, Nstates, threshold){
  # apply the function below to the columns of a matrix.
  # threshold could be a vector or a matrix.
  
  if(length(Nstates) >= 1 && length(Nstates) < columnsIndex){
    stop("Nstates length should be equal to 1 or equal to the number of columnsIndex")
  }
  
  if(length(Nstates) == 1){
    Nstates <- rep(Nstates, length(columnsIndex))
  }
  
  if(is.vector(threshold)){
    threshold <- matrix(rep(threshold, length(Nstates)), nrow = length(Nstates), byrow = TRUE)
  }
  
  if(nrow(threshold) < length(Nstates)){
    stop("The threshold number of rows should be equal to the Nstates length.")
  }
  
  for (c in length(columnsIndex)){
    matrix[,columnsIndex[c]] <- ConvertContinousInDiscreteValues(matrix[,columnsIndex[c]], Nstates[c], threshold[c,])
  }
  return(matrix)
}

# Scaled continuous traits
#########################
#Scale traits between 0 and 10 
#(https://github.com/GitTFJ/Handling-missing-values-in-trait-data/blob/main/Script1_SimulateData_V1.0.R)
RangeScale <- function(x, range){((x-min(x))/(max(x)-min(x))*range)} #à changer on doit pouvoir définir un range pour chaque trait qui est différent!!!!!
ScaledTraits <- apply(ContinuousData$data, MARGIN = 2, RangeScale, 10)
ScaledTraits <- as.data.frame(ScaledTraits)
ScaledTraits

# modifMatrix <- ChangeContinuousTraitInDiscrete(ScaledTraits, c(6,7), 3, c(2,6))
# modifMatrix
# 
# 
# Birth <- 0.4
# Death <- 0.1
# Ntaxa <- 30
# Ntraits <- 7
# 
# 
# #Simulate discrete traits (for independent data)
# #########################
# SimDiscrete<-simDiscreteTraits(2, c(2,3), "ARD", 0.1, SimTree)
# SimTree$tip.label
# class(SimTree)
# SimDiscrete$tip_mat
# 
# # Simulate continuous Data
# ##########################
# ContinuousData <- simContinuousData(Ntraits, Sigma, SimTree, Ntaxa = Ntaxa, model = "BM1", theta = NULL, alpha = NULL, uncovTraits = NULL)
# v<- simContinuousData(Ntraits, Cor = NULL, Sigma2 = NULL, SimTree, Ntaxa = Ntaxa, model = "BM1", theta = NULL, alpha = NULL, uncovTraits = 2)
# ContinuousData$data
# ContinuousData$varMatrix
# 
# # Scaled continuous traits
# #########################
# #Scale traits between 0 and 10 
# #(https://github.com/GitTFJ/Handling-missing-values-in-trait-data/blob/main/Script1_SimulateData_V1.0.R)
# RangeScale <- function(x){((x-min(x))/(max(x)-min(x))*10)}
# ScaledTraits <- apply(ContinuousData$data, MARGIN = 2, RangeScale)
# ScaledTraits <- as.data.frame(ScaledTraits)
# ScaledTraits
# 
# # Compare ground truth and estimated variance-covariance matrix
# ###############################################################
# ContinuousData$data # a convertir 
# #cov2cor(Data$sigma)
# heatmap(ContinuousData$varMatrix) # Evolutionary trait correlation
# ContinuousData$sigma

Birth <- 0.4
Death <- 0.1
Ntaxa <- 30
Ntraits <- 7
NbrConTraitsCor <- 7
NbrConTraitsUnCor <- 1
NbrDisTraitsCor <- 3
NbrDisTraitsUncor <-2
Nstates <- c(2,3,2,2,4) #could be a vector of length equal to NbrNbrDisTraitsCor + NbrDisTraitsUncor or a vector of a lenght of 2. 
rangeConTraits = 10 #a changer ci-dessus
threshold = 

simData <- function(NbrConTraitsCor = NULL, NbrConTraitsUncor = NULL, NbrDisTraitsCor = NULL , NbrDisTraitsUncor = NULL, Birth, Death, Nstates, Cor = NULL, Sigma2 = NULL, Ntaxa = 2, model = "BM1", theta = NULL, alpha = NULL, rangeConTraits = NULL, rate_model = "ER", max_rate = 0.1, threshold, save = TRUE){
  
  Ntraits <- NbrConTraitsCor + NbrConTraitsUncor + NbrDisTraitsCor + NbrDisTraitsUncor
  
  if(length(Nstates) != (NbrNbrDisTraitsCor + NbrDisTraitsUncor) && length(Nstates) != 2){
    stop(sprintf("Nstates should be of length %i or 2", (NbrNbrDisTraitsCor + NbrDisTraitsUncor)))
  }
  
  if(length(Nstates) != 2 && (NbrNbrDisTraitsCor + NbrDisTraitsUncor) > 2){
    Nstates <- c(rep(Nstates[1], NbrNbrDisTraitsCor), rep(Nstates[2], NbrDisTraitsUncor))
  }
  
  # Simulate phylogeny
  ####################
  # Simulating phylogenies fails sometimes. Try until we are successful
  Extant <- FALSE
  while (!Extant) {
    SimTree <- pbtree(b = Birth, d = Death, n = Ntaxa, scale = 1, extant.only = FALSE) #why added extant.only ?
    if (!is.null(SimTree)) {
      Extant <- TRUE
    }
  }
  
  # Simulate sigma matrix
  #####################
  if (!is.null(uncovTraits)){
    Sigmas <- simSigma((Ntraits -  NbrDisTraitsUncor), Cor = Cor, Sigma2 = Sigma2, uncovTraits = NbrConTraitsUncor)
  }
  else{
    Sigmas <- simSigma((Ntraits - NbrDisTraitsUncor), Cor = Cor, Sigma2 = Sigma2)
  }
  
  # Simulate continuous Data
  ########################## 
  ContinuousData <- simContinuousData((Ntraits -  NbrDisTraitsUncor), Sigmas, SimTree, Ntaxa = Ntaxa, model = "BM1", theta = NULL, alpha = NULL)
  
  
  # Column index of uncorrelated continuous traits
  ################################################
  indexUsed <- c()
  for (i in 1:dim(ContinuousData$data)[2]){
    unitVec <- numeric(dim(ContinuousData$data)[1])
    unitVec[i] <- 1
    if (identical(ContinuousData$data[,i], unitVec)){
      indexUsed <- c(indexUsed, i)
    }
  }

  # Scaled the continuous traits in a range
  #########################################
  if (!is.null(rangeConTraits)){
    ContinuousData$data <- apply(ContinuousData$data, MARGIN = 2, RangeScale, rangeConTraits)
  }
  
  # Define correlated discrete traits
  ###################################
  if(!is.null(NbrDisTraitsCor)){
    if(!is.null(uncovTraits)){
      columnIndex <- 1:dim(ContinuousData$data)[2]
      columnIndexFree[!columnIndex %in% indexUsed]
      columnIndexToUse <- sample(columnIndexFree, NbrDisTraitsCor)
      ContinousData$data <- ChangeContinuousTraitInDiscrete(ContinousData$data, columnIndexToUse, Nstates, threshold)
      
    }
    columnIndexToUse <- sample(1:dim(ContinuousData$data[2], NbrDisTraitsCor))
    ContinousData$data <- ChangeContinuousTraitInDiscrete(ContinousData$data, columnIndexToUse, Nstates, threshold)
  }
  
  #Total Data matrix
  ##################
  TotData <- ContinuousData$data
  
  if (!is.null(NbrDisTraitsUncor)){
    DiscreteData <- simDiscreteTraits(NbrDisTraitsUncor, Nstates, rate_model, max_rate, SimTree)
    #Total Data matrix (continuous and independent discrete)
    ########################################################
    TotData <- cbind(ContinuousData$data, DiscreteData)
    if(!is.null(NbrDisTraitsCor)){
      DiscreteData <- cbind(DiscreteData, ContinuousData$data[, columnIndexToUse])
      ContinuousData$data <- ContinuousData$data[, -columnIndexToUse]
    }
  }
  
  
  #Define list of object
  ######################
  Data <- list(TotalData = TotData, corMatrix = cor(TotData), sigmas = Sigmas, DiscreteData = DiscreteData, ContinuousData = ContinuousData$data)
  
  #Save data
  ##########
  if (save){
    save(Data, file="simulatedData.RData")
  }
  
  return(Data)
}

simData()

nrow(threshold)

