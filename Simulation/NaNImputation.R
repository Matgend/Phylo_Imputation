# Simulation of Missing values (NaN) functions
##############################################

#' @title Partition discrete and continuous
#'
#' @description This function split the columns according to its names in a discrete or continuous nested list 
#'
#' @usage splitDiscAndContiColnames(colunmNames)
#'
#' @param columnNames vector of characters (columns names of the dataset)
#' @return 2 nested lists one with discrete and the second with continuous column names
#'
splitDiscAndContiColnames <- function(columnNames){
  
  #discrete traits
  pattern <- "I\\d+\\...."
  DiscreteColnames <- str_extract(columnNames, pattern)
  DiscreteColnames <- DiscreteColnames[!is.na(DiscreteColnames)]

  #continuous traits
  pattern <- "F\\d+\\...."
  ContiColnames <- str_extract(columnNames, pattern)
  ContiColnames <- ContiColnames[!is.na(ContiColnames)]
  
  return(list(Discrete = DiscreteColnames, Continuous = ContiColnames))
}


#Subset generation
##################
#' @title Partition the columns according to a following partition design
#'
#' @description This function split the columns according to its names as following: 
#' # Partition the data:
#   - All traits
#   - Independent traits:
#       - Discrete
#       - Continuous
#
#   - Correlated traits:
#       - Mixed continuous / discrete
#       - Discrete
#       - Continuous 
#'
#' @usage dataPartition(Data)
#' @param Data a nested list with at least the structure of the dataset and the dataset itself
#' @return  a nested list having the column name of the partitionned data as in the description

dataPartition <- function(Data){

  correlation_values <- unique(Data$dataframe$correlation)
  IndDiscreteTraits <- list()
  IndContinousTraits <- list()
  MixedCorrelatedTraits <- list()
  CorrDiscreteTraits <- list()
  CorrContinuousTraits <- list()
  
  for (c in correlation_values){
    
    #extract name of the columns
    pattern <- paste("[:alnum:]*\\.", c,"/[:digit:]", sep = "")
    subColnames <- str_extract(colnames(Data$FinalData), pattern)
    subColnames <- subColnames[!is.na(subColnames)]
    
    #number of independent traits for each correlation group
    uncorrTraitsByCorrGroup <- sum(Data$dataframe$uncorr_traits[Data$dataframe$correlation == c])
    
    #correlated traits
    CorrSubColnames <- subColnames[uncorrTraitsByCorrGroup+1:length(subColnames)]
    CorrSubColnames <- CorrSubColnames[!is.na(CorrSubColnames)]
    ContiTraits <- splitDiscAndContiColnames(CorrSubColnames)
    
    if(length(CorrSubColnames) != 0 & length(grep("F.", CorrSubColnames)) >= 1 & 
       length(grep("I.", CorrSubColnames)) >= 1){
      MixedCorrelatedTraits <- c(MixedCorrelatedTraits, list(CorrSubColnames))
    }
    
    if(length(ContiTraits$Discrete) != 0){
      CorrDiscreteTraits <- c(CorrDiscreteTraits, list(ContiTraits$Discrete))
    }
    
    if(length(ContiTraits$Continuous) != 0){
      CorrContinuousTraits <- c(CorrContinuousTraits, list(ContiTraits$Continuous))
    }
    
    if(uncorrTraitsByCorrGroup != 0){
      
      IndeSubColnames <- subColnames[1:uncorrTraitsByCorrGroup]
      
      IndeTraits <- splitDiscAndContiColnames(IndeSubColnames)
      
      if(length(IndeTraits$Discrete) != 0){
        IndDiscreteTraits <- c(IndDiscreteTraits, list(IndeTraits$Discrete))
      }
      
      if(length(IndeTraits$Continuous) != 0){
        IndContinousTraits <- c(IndContinousTraits, list(IndeTraits$Continuous))
      }
    }
  }
  partitions <- list(AllTraits = list(colnames(Data$FinalData)), IndDiscreteTraits = IndDiscreteTraits, 
                     IndContinousTraits = IndContinousTraits, MixedCorrelatedTraits = MixedCorrelatedTraits,
                     CorrDiscreteTraits = CorrDiscreteTraits, CorrContinuousTraits = CorrContinuousTraits)
  return(partitions)
}


#' @title Identify species for phylogenetic NAs
#'
#' @description This function returns the names of the species that should get NAs for their traits
#'
#' @usage getTipsNA(Tree, MinTips) 
#'
#' @param Tree Phylogeny 
#' @param missingRates numerical vector corresponding to the rate of missing value to introduce in the data
#' @return vector containing the names of the species that should get NAs for their traits

getTipsNA <- function (Tree, missingRates){
  Nodes <- Nnode(Tree)
  Tips <- Ntip(Tree)
  MaxNodeID <- Tips + Nodes
  NodeIDs <- (Tips + 1):MaxNodeID
  MinTips <- round(missingRates * Tips, 1)
  EnoughTips <- FALSE
  while (!EnoughTips) {
    FocalNode <- sample(NodeIDs, 1)
    Desc <- phangorn::Descendants(Tree, FocalNode, type = "tips")[[1]]
    SampledTips <- Tree$tip.label[Desc]
    if(length(SampledTips) >= MinTips){
      EnoughTips <- TRUE
    }
  }
  if(length(SampledTips) > MinTips){
    SampledTips <- sample(SampledTips, MinTips)
  }
  return(SampledTips)
}


#MNAR simulation
################

#'@title Impute NaNs in a simulated data according the MNAR mechanisms
#'
#'@description Simulate NaNs in a partitioned complete data set composed of discrete and continuous traits which are correlated or uncorrelated. The NaNs are imputed according a missing rate and the missing mechanism MNAR.
#'
#'@usage myMNAR(missingRate, partition, col_mis)
#'
#'@param missingRates numerical vector corresponding to the rate of missing value to introduce in the data
#'@param partitions nested list having character vectors corresponding to the data partition of discrete traits 
#'@param col_mis vector of index or columns name where missing value will be created
#'@return a dataset with NA following a pattern of MNAR.

myMNAR <- function(missingRate, partition, col_mis){
  
  if(is.numeric(col_mis)){
    colPartition <- col_mis
  }
  
  else if(is.character(col_mis)){
    colPartition <- match(col_mis, names(partition))
  }
  
  for(col in colPartition){
    
    if(length(grep("F.", names(partition[, col]))) == 1){
      partition[ , col] <- delete_MNAR_censoring(partition[ , col], missingRates, 
                            cols_mis = col, where = "upper") #don't change where arg
    }
    
    else{
      levelsVector <- sort(unique(partition[ ,col]))
      
      nbrValueToRemove <- round(nrow(partition) * missingRate)
      
      maxState <- tail(levelsVector, 1)
      
      #in case not enough value in max state, remove in the other states
      while (nbrValueToRemove > 0){
        
        levelsVector <- levelsVector[-which(levelsVector == maxState)]
        indexLargeState <- which(partition[,col] == maxState)
        
        #keep present 1 value of the state
        indexLargeState <- indexLargeState[-which(indexLargeState == sample(indexLargeState, 1))]
        
        if(length(indexLargeState) >= nbrValueToRemove){
          indexLargeState <- sample(indexLargeState, nbrValueToRemove)
          nbrValueToRemove <- 0
          
        }
        
        else{
          nbrValueToRemove <- nbrValueToRemove - length(indexLargeState)
          indexLargeState <- sample(indexLargeState, length(indexLargeState))
          maxState <- sample(levelsVector, 1)
        }
        
        #apply NA
        partition[indexLargeState, col] <- NA
      }
    }
  }
  return(partition)
}

#NaN imputation
##############

#'@title Impute NaNs in a simulated data according to MCAR, MAR and MNAR
#'
#'@description Simulate NaNs in a partitioned complete data set composed of continuous and discrete traits which are correlated or uncorrelated. The NaNs are imputed according a missing rate and respecting 3 categories of missing data, MCAR, MAR and MNAR.
#'
#'@usage NaNImputation(missingRates, partitions, trees, missTraits, replicates, save = TRUE)
#'
#'@param missingRates numerical vector corresponding to the rate of missing value to introduce in the data
#'@param partitions nested list having character vectors corresponding to the data partition
#'@param data list, containing all the simulations (data, trees and parameters)
#'@param missTraits numerical, number of traits in which there is missing data.
#'@param save character correspond to the name of the saved file in .RData format
#'@return a nested list composed of the partitioned data with the 3 categories of missing data (MCRA, MAR and MNAR) according to a precise missing rate and of list of specific traits with NaNs imputed according to MCAR (phylogeny).

# missingRates <- 0.5
# partitions <- dataPartition(Data)
# data <- Data
# missTraits <- ncol(Data$FinalData)
# data$FinalData
# test <- NaNImputation(missingRates, partitions, data, missTraits)
NaNImputation <- function(missingRates, partitions, data, missTraits, save = NULL){
  DataNaN <- list()
  MCAR <- list()
  MAR <- list()
  MNAR <- list()
  namesMAR <- c()
  namesMCAR <- c()
  namesMNAR <- c()
  mechanism <- c("MCAR", "MAR", "MNAR")
  for(l in seq_along(partitions)){

    for(colList in partitions[[l]]){

      if (length(colList) == 0){
        next
      }

      if (missTraits > length(colList)){
        missTraits <- length(colList)
      }
      
      for(mr in 1:length(missingRates)){

        namesMCAR <- c(namesMCAR, paste("MCAR", names(partitions)[l], length(colList), round(missingRates[mr],2) ,sep = "/"))
        
        namesMNAR <- c(namesMNAR, paste("MNAR", names(partitions)[l], length(colList), round(missingRates[mr],2), sep = "/"))
        
        #univariate
        if(length(colList) == 2){
          
          if(missingRates[1] == missingRates[mr]){
            colMis <- sample(colList, missTraits)
          }

          #MCAR
          missingMCAR <- delete_MCAR(data$FinalData[,colList, drop = FALSE], missingRates[mr], cols_mis = colMis)
          
          #MNAR
          missingMNAR <- myMNAR(missingRates[mr], data$FinalData[,colList, drop = FALSE], colMis)
          
        }
        
        if(length(colList) == 1){
          oneColum <- as.data.frame(data$FinalData[ ,colList, drop = FALSE])
          names(oneColum) <- colList
          missingMCAR <- delete_MCAR(oneColum, missingRates[mr], missTraits)
          missingMNAR <- myMNAR(missingRates[mr], oneColum, missTraits)
        }
        
        if(length(colList) > 2){
          colMis <- sample(colList, missTraits)
          
          #MCAR
          missingMCAR <- delete_MCAR(data$FinalData[,colList], missingRates[mr], 
                                     cols_mis = colMis, p_overall = FALSE)
          
          #MNAR
          missingMNAR <- myMNAR(missingRates[mr], data$FinalData[,colList, drop = FALSE], colMis)
        }
        
        if(l > 3 & (length(colList) != 1)){

          namesMAR <- c(namesMAR, paste("MAR", names(partitions)[l], length(colList), round(missingRates[mr],2), sep = "/"))
          
          if(length(colList) %% 2 == 0){
            
            if(missTraits > length(colList)/2){
              missTraits <- length(colList)/2
            }
            
            if(missingRates[1] == missingRates[mr]){
              cols_misIndex <- sample(length(colList), missTraits)
              cols_ctrlIndex <- c(1:length(colList))[-cols_misIndex]
              cols_ctrlIndex <- sample(cols_ctrlIndex, missTraits)
            }
          }
          if(length(colList) %% 2 != 0 ){
            
            
            if(missTraits > (length(colList)-1)/2){
              missTraits <- (length(colList)-1)/2
            }
            if(missingRates[1] == missingRates[mr]){
              cols_misIndex <- sample(length(colList), missTraits)
              cols_ctrlIndex <- c(1:length(colList))[-cols_misIndex]
              cols_ctrlIndex <- sample(cols_ctrlIndex, missTraits)
            }
          }
          
          if(length(colList) == 2){

            #MAR
            if(missingRates[1] == missingRates[mr]){
              cols_misIndex <- sample(1:2, 1)
              cols_ctrlIndex <- c(1,2)[-cols_misIndex]
            }
            
            missingMAR <- delete_MAR_censoring(data$FinalData[,colList], missingRates[mr], 
                                               cols_mis = cols_misIndex, cols_ctrl = cols_ctrlIndex, 
                                               where = "upper") #don't change where arg

          }
          
          else{

            #MAR
            missingMAR <- delete_MAR_censoring(data$FinalData[,colList], missingRates[mr], 
                                               cols_mis = cols_misIndex, cols_ctrl = cols_ctrlIndex, 
                                               where = "upper") #decide where argument                                                                                                                        we want to.
          }
          MAR <- c(MAR, list(missingMAR))
        }
       MCAR <- c(MCAR, list(missingMCAR))
       MNAR <- c(MNAR, list(missingMNAR))
      }
    }
    #rename the nested lists
    names(MCAR) <- namesMCAR
    names(MAR) <- namesMAR
    names(MNAR) <- namesMNAR
  }
  
  # Drop randomly tips in phylogenetic tree
  PhyloNaN <- list()
  namesTrees <- c()
  
  for(mr in 1:length(missingRates)){
    
    #copy of FinalData
    missingPhylo <- data$FinalData
    
    #in case missingRates larger than 50%
    if(missingRates[mr] > 0.5){
      print("Can't simulate PhyloNaN")
      next
    }
    
    #tips to include NaN
    tips <- getTipsNA(data$TreeList$`0`, missingRates[mr])
    
    namesTrees <- c(namesTrees, paste("PhyloNaN", length(tips), round(missingRates[mr],2), sep = "/"))
  
    missingPhylo[tips, ] <- NA
  
    PhyloNaN <- c(PhyloNaN, list(missingPhylo))
  }
  
  #Rename the nested list
  names(PhyloNaN) <- namesTrees
  
  #Data with NA imputed
  DataNaN <- list(MCAR = MCAR, MAR = MAR, MNAR = MNAR, PhyloNaN = PhyloNaN)
  Parameters <- list(missingRates = missingRates, missTraits = missTraits)
  NaNImputed <- list(DataNaN = DataNaN, Parameters = Parameters)
  
  #Save data
  ##########
  if(!is.null(save)){
    save(NaNImputed, file = paste0(save, ".RData"))
  }
  
  return(NaNImputed)
}



