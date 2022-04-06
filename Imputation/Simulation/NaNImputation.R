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

#MCAR simulation
################

#'@title Impute NaNs in a simulated data according the MCAR mechanisms
#'
#'@description Simulate NaNs in a partitioned complete data set composed of discrete and continuous traits which are correlated or uncorrelated. The NaNs are imputed according a missing rate and the missing mechanism MCAR.
#'
#'@usage myMCAR(missingRate, partition, cols_mis)
#'
#'@param missingRate numerical vector corresponding to the rate of missing value to introduce in the data
#'@param ds dataframe in which NA should be created
#'@param cols_mis vector of index or columns name where missing value will be created
#'@return a dataset with NA following a pattern of MNAR.

# missingRate <- 0.5
# ds <- Data$FinalData
# cols_mis <- names(Data$FinalData)[1]
# cols_ctrl <- "F3.2/3"
# t <- myMCAR(missingRate, ds, cols_mis)

myMCAR <- function(missingRate, ds, cols_mis){
  
  if(is.numeric(cols_mis)){
    cols_mis <- names(ds)[cols_mis]
  }
  
  #call MCAR function
  if(ncol(ds) == 1){
    missingMCAR <- delete_MCAR(ds, missingRate, cols_mis = cols_mis)
  }

  if(ncol(ds) > 1){
    missingMCAR <- delete_MCAR(ds, missingRate, cols_mis = cols_mis, p_overall = FALSE)
  }
  #check if all the states are represented. (discrete traits)
  discIndex <- grep("I.", names(ds[, cols_mis, drop = FALSE]))
  
  if(length(discIndex) != 0){
    
    for(t in 1:length(discIndex)){
      
      nbrStates <- unique(ds[, discIndex[t]])
      nbrStatesMiss <- unique(na.omit(missingMCAR[, discIndex[t]]))
      diff <- setdiff(nbrStates, nbrStatesMiss)
      
      #if > 0 means not all the states are represented. 
      if(length(diff) != 0){
        
        rowIndexSave <- c()
        
        for(s in 1:length(diff)){
          
          rowIndexSate <- which(ds[, discIndex[t]] == diff[s])
          
          indexToKeep <- sample(rowIndexSate, 1)
          
          rowIndexSave <- c(rowIndexSave, indexToKeep)

        }
        
        newMissingRate <- nrow(ds) * missingRate / nrow(ds[-rowIndexSave, ]) 
        
        #call MCAR function
        newMissingMCAR <- delete_MCAR(ds[-rowIndexSave, ], newMissingRate, cols_mis = cols_mis, p_overall = FALSE)
        
        missingMCAR <- ds
        missingMCAR[rownames(newMissingMCAR), ] <- newMissingMCAR
      }
    }
  }
  return(missingMCAR)
}

#MAR simulation
################

#'@title Impute NaNs in a simulated data according the MAR mechanisms
#'
#'@description Simulate NaNs in a partitioned complete data set composed of discrete and continuous traits which are correlated or uncorrelated. The NaNs are imputed according a missing rate and the missing mechanism MAR.
#'
#'@usage myMAR(missingRate, partition, cols_mis, cols_ctrl)
#'
#'@param missingRate numerical vector corresponding to the rate of missing value to introduce in the data
#'@param ds dataframe in which NA should be created
#'@param cols_mis vector of index or columns name where missing value will be created
#'@param cols_ctrl vector of index or columns name which control the creation of missing values (same size of cols_mis)
#'@return a dataset with NA following a pattern of MNAR.

# missingRate <- 0.5
# ds <- Data$FinalData
# cols_mis <- names(Data$FinalData)[1]
# cols_ctrl <- "F3.2/3"
# myMAR(missingRate, ds, cols_mis, cols_ctrl)[, 1]

myMAR <- function(missingRate, ds, cols_mis, cols_ctrl){
  
  if(is.numeric(cols_mis)){
    cols_mis <- names(ds)[cols_mis]
  }

  #call MAR function
  missingMAR <- delete_MAR_censoring(ds, missingRate, 
                                     cols_mis = cols_mis, cols_ctrl = cols_ctrl, 
                                     where = "upper")
  
  #check if all the states are represented. (discrete traits)

  discIndex <- grep("I.", names(ds[, cols_mis, drop = FALSE]))
  
  if(length(discIndex) != 0){
    
    for(t in 1:length(discIndex)){
      nbrStates <- unique(ds[, discIndex[t]])
      
      nbrStatesMiss <- unique(na.omit(missingMAR[, discIndex[t]]))
      
      diff <- setdiff(nbrStates, nbrStatesMiss)
      
      #if > 0 means not all the states are represented. 
      if(length(diff) != 0){
        
        rowIndexSave <- c()
        
        for(s in 1:length(diff)){

          #index missing States
          indexMissStates <- which(ds[, discIndex[t]] == diff[s])
          
          #isolate values of ctrl
          valuesCor <- ds[, cols_ctrl[t]][indexMissStates]

          #if continuous values
          if(is.numeric(valuesCor)){
            
            rowIndexMinVal <- which.min(valuesCor)
            rowIndexSave <- c(rowIndexSave, indexMissStates[rowIndexMinVal])
          }
          
          else{
            rowIndexMinVal <- sample(indexMissStates, 1)
            rowIndexSave <- c(rowIndexSave, rowIndexMinVal)
          }
        }
        
        newMissingRate <- nrow(ds) * missingRate / nrow(ds[-rowIndexSave, ]) 
        #call MAR function
        newMissingMAR <- delete_MAR_censoring(ds[-rowIndexSave, ], newMissingRate, 
                                           cols_mis = cols_mis, cols_ctrl = cols_ctrl, 
                                           where = "upper")
        
        missingMAR <- ds
        missingMAR[rownames(newMissingMAR), ] <- newMissingMAR
      }
    }
  }
  return(missingMAR)
}

#MNAR simulation
################

#'@title Impute NaNs in a simulated data according the MNAR mechanisms
#'
#'@description Simulate NaNs in a partitioned complete data set composed of discrete and continuous traits which are correlated or uncorrelated. The NaNs are imputed according a missing rate and the missing mechanism MNAR.
#'
#'@usage myMNAR(missingRate, partition, cols_mis)
#'
#'@param missingRate numerical vector corresponding to the rate of missing value to introduce in the data
#'@param ds dataframe in which NA should be created
#'@param cols_mis vector of index or columns name where missing value will be created
#'@return a dataset with NA following a pattern of MNAR.
#'

myMNAR <- function(missingRate, ds, cols_mis){
  
  if(is.numeric(cols_mis)){
    colPartition <- cols_mis
  }
  
  else if(is.character(cols_mis)){
    colPartition <- match(cols_mis, names(ds))
  }
  
  for(col in colPartition){
    if(length(grep("F.", names(ds)[col])) == 1){
      ds[ , col] <- delete_MNAR_censoring(ds[ , col, drop = FALSE], missingRate, 
                            cols_mis = cols_mis[col], where = "upper") #don't change where arg
    }
    
    else{
      levelsVector <- sort(unique(ds[ ,col]))
      
      nbrValueToRemove <- round(nrow(ds) * missingRate)
      
      maxState <- tail(levelsVector, 1)
      
      #in case not enough value in max state, remove in the other states
      while (nbrValueToRemove > 0){
        
        levelsVector <- levelsVector[-which(levelsVector == maxState)]
        indexLargeState <- which(ds[,col] == maxState)
        
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
        ds[indexLargeState, col] <- NA
      }
    }
  }
  return(ds)
}

#NaN imputation
##############

#'@title Impute NaNs in a simulated data according to MCAR, MAR and MNAR
#'
#'@description Simulate NaNs in a partitioned complete data set composed of continuous and discrete traits which are correlated or uncorrelated. The NaNs are imputed according a missing rate and respecting 3 categories of missing data, MCAR, MAR and MNAR.
#'
#'@usage NaNImputation(missingRate, partitions, trees, missTraits, replicates, save = TRUE)
#'
#'@param missingRate numerical vector corresponding to the rate of missing value to introduce in the data
#'@param partitions nested list having character vectors corresponding to the data partition
#'@param data list, containing all the simulations (data, trees and parameters)
#'@param missingTraits numerical, number of traits in which there is missing data.
#'@param save character correspond to the name of the saved file in .RData format
#'@return a nested list composed of the partitioned data with the 3 categories of missing data (MCRA, MAR and MNAR) according to a precise missing rate and of list of specific traits with NaNs imputed according to MCAR (phylogeny).

#missingRates <- 0.05
#partitions <- dataPartition(Data)
# partitions
# 
#partitions <- dataPartition(new_data)
# partitions
# # # data <- Data
#missTraits <- ncol(Data$FinalData)
# missTraits <- ncol(new_data$FinalData)
# data$FinalData

#test$DataNaN$MCAR$`MCAR/CorrContinuousTraits/3/0.05`
NaNImputation <- function(missingRates, partitions, data, missingTraits, save = NULL){
  DataNaN <- list()
  MCAR <- list()
  MAR <- list()
  MNAR <- list()
  namesMAR <- c()
  namesMCAR <- c()
  namesMNAR <- c()
  
  for(l in seq_along(partitions)){
    mergedMissgMAR <- list()
    for(colList in partitions[[l]]){
      
      if (length(colList) == 0){
        next
      }
      
      if (missingTraits > length(colList)){
        missTraits <- length(colList)
      }
      
      else if(missingTraits <= length(colList)){
        missTraits <- missingTraits
      }
      
      for(mr in 1:length(missingRates)){

        namesMCAR <- c(namesMCAR, paste("MCAR", names(partitions)[l], length(colList), 
                                        round(missingRates[mr],2) ,sep = "/"))
        
        namesMNAR <- c(namesMNAR, paste("MNAR", names(partitions)[l], length(colList), 
                                        round(missingRates[mr],2), sep = "/"))
        
        #univariate
        if(length(colList) == 2){
          
          if(missingRates[1] == missingRates[mr]){
            colMis <- sample(colList, missTraits)
          }

          #MCAR
          missingMCAR <- myMCAR(missingRates[mr], data$FinalData[,colList, drop = FALSE], cols_mis = colMis)
          
          #MNAR
          missingMNAR <- myMNAR(missingRates[mr], data$FinalData[,colList, drop = FALSE], colMis)
          
        }
        
        if(length(colList) == 1){
          oneColum <- as.data.frame(data$FinalData[ ,colList, drop = FALSE])
          names(oneColum) <- colList
          missingMCAR <- myMCAR(missingRates[mr], oneColum, missTraits)
          missingMNAR <- myMNAR(missingRates[mr], oneColum, missTraits)
        }
        
        
        if(length(colList) > 2){
          
          if(length(colList) == missTraits){
            colMis <- colList
          }
          else{
            colMis <- sample(colList, missTraits)
          }

          #MCAR
          missingMCAR <- myMCAR(missingRates[mr], data$FinalData[,colList], 
                                     cols_mis = colMis)
          
          #MNAR
          missingMNAR <- myMNAR(missingRates[mr], data$FinalData[,colList, drop = FALSE], colMis)
        }
        
        if(l > 3 & (length(colList) != 1)){
          
          #in case manual in experimental design
          corGroup <- unique(str_extract(colList, "\\d+(?=\\/)"))
          if(data$dataframe$model[which(data$dataframe$correlation == corGroup)] == "Manual"){
            
            namesMAR <- c(namesMAR, paste("MAR", names(partitions)[l], length(colList),
                                          round(missingRates[mr],2), sep = "/"))
            
            #sample vector 
            ctrlIndex <- grep(paste0(corGroup, "/"), names(data$FinalData))
            cols_ctrlIndex <- sample(ctrlIndex, 1)
            
            #add discrete trait to data
            missingMAR <- myMAR(missingRates[mr], data$FinalData, 1, cols_ctrlIndex)

           
            if(!any(is.na(missingMAR[,1]))){
              return(0)
            }
             
          }
          
          else{
    
            namesMAR <- c(namesMAR, paste("MAR", names(partitions)[l], length(colList), 
                                          round(missingRates[mr],2), sep = "/"))
            
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
              
              missingMAR <- myMAR(missingRates[mr], data$FinalData[,colList], cols_misIndex, cols_ctrlIndex)

            }
            
            else{
              
              #MAR
              missingMAR <- myMAR(missingRates[mr], data$FinalData[,colList], cols_misIndex, cols_ctrlIndex)
              
            }
            
            if(length(partitions[[l]]) > 1){
              mergedMissgMAR <- c(mergedMissgMAR, list(missingMAR))
            }
            
            #add traits with MAR value to all traits.
            copyData <- data$FinalData
            copyData[, colnames(missingMAR)] <- missingMAR
            MAR <- c(MAR, list(copyData))
            namesMAR <- c(namesMAR, paste("MAR", paste0(names(partitions)[l],"ALL"), ncol(copyData), 
                                          round(missingRates[mr],2), sep = "/"))
          }
          
          MAR <- c(MAR, list(missingMAR))
        }
       MCAR <- c(MCAR, list(missingMCAR))
       MNAR <- c(MNAR, list(missingMNAR))
      }
      
    }
    
    #merge the partition of MAR data
    if(length(mergedMissgMAR) != 0){
      MAR <- c(MAR, list(do.call(cbind, mergedMissgMAR)))
      namesMAR <- c(namesMAR, paste("MAR", paste0(names(partitions)[l],"MERGEDPAR"), ncol(copyData), 
                                    round(missingRates[mr],2), sep = "/"))
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
  DataNaN <- list(MCAR = MCAR[1], MAR = MAR, MNAR = MNAR[1], PhyloNaN = PhyloNaN) 
  #remove if want other partitions for MCAR and MNAR
  Parameters <- list(missingRates = missingRates, missTraits = missTraits)
  NaNImputed <- list(DataNaN = DataNaN, Parameters = Parameters)
  
  #Save data
  ##########
  if(!is.null(save)){
    save(NaNImputed, file = paste0(save, ".RData"))
  }
  
  return(NaNImputed)
}

#partitions <- dataPartition(new_data)
#test <- NaNImputation(missingRates, partitions, Data, missTraits)
# test$DataNaN$MAR$`MAR/CorrContinuousTraits3_2/0.05`
#Data <- get(load("F:\\Master_Thesis\\Simulations\\ARD\\005 _10T\\FullData\\simulatedDataDiscreteCorData10KAP0_R53_V0.05.RData"))

# Data <- get(load("C:/Users/Matthieu/Documents/UNIFR/Master_thesis/Scripts/simulatedDataDiscreteARD1Data_R59_0.05.RData"))
# partitions <- dataPartition(Data)
# test <- NaNImputation(0.5, partitions, Data, ncol(Data$FinalData))
# Data$FinalData[,1]

