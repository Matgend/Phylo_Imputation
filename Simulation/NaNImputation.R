#install.packages("tidyverse")
#install.packages("missMethods")

library(tidyverse)
library(geiger)
library(missMethods)
#https://github.com/torockel/missMethods/blob/master/R/delete_censoring.R


# load simulated data
setwd("C:/Users/Matthieu/Documents/UNIFR/Master_thesis/Scripts/")
load("simulatedData.RData")

# Partition the data:
#   - All traits
#   - Correlated traits:
#       - Mixed continuous / discrete
#       - Discrete
#       - Continuous
#   
#   - Independent traits:
#       - Discrete
#       - Continuous
#
# Missingness rates ....
# Missingness tpyes: MCAR, MAR, MNAR


splitDiscAndContiColnames <- function(columnNames){
  #columnNames: vector of characters (columns names of the dataset)
  #return: 2 nested lists one with discrete and the second with continuous column names 
  
  #discrete traits
  pattern <- "I.\\.."
  DiscreteColnames <- str_extract(columnNames, pattern)
  DiscreteColnames <- DiscreteColnames[!is.na(DiscreteColnames)]

  #continuous traits
  pattern <- "F.\\.."
  ContiColnames <- str_extract(columnNames, pattern)
  ContiColnames <- ContiColnames[!is.na(ContiColnames)]
  
  return(list(Discrete = DiscreteColnames, Continuous = ContiColnames))
}

splitDiscAndContiColnames(colnames(Data$DiscreteData))


#Subset generation
##################
dataPartition <- function(Data){
  # Data: a nested list with at least the structure of the dataset and the dataset itself
  # return: a nested list having the column name of the partitionned data:
    # Partition the data:
    #   - All traits
    #   - Independent traits:
    #       - Discrete
    #       - Continuous
    #
    #   - Correlated traits:
    #       - Mixed continuous / discrete
    #       - Discrete
    #       - Continuous
       
    
  
  correlation_values <- unique(Data$dataframe$correlation)
  IndDiscreteTraits <- list()
  IndContinousTraits <- list()
  MixedCorrelatedTraits <- list()
  CorrDiscreteTraits <- list()
  CorrContinuousTraits <- list()
  
  for (c in correlation_values){
    
    #extract name of the columns
    pattern <- paste("..\\.", c, sep = "")
    subColnames <- str_extract(colnames(Data$FinalData), pattern)
    subColnames <- subColnames[!is.na(subColnames)]
    
    #number of independent traits for each correlation group
    uncorrTraitsByCorrGroup <- sum(Data$dataframe$uncorr_traits[Data$dataframe$correlation == c])
    
    #correlated traits
    CorrSubColnames <- subColnames[uncorrTraitsByCorrGroup+1:length(subColnames)]
    CorrSubColnames <- CorrSubColnames[!is.na(CorrSubColnames)]
    ContiTraits <- splitDiscAndContiColnames(CorrSubColnames)
    
    if(length(CorrSubColnames) != 0){
      #append correlated traits cluster
      #name <- sprintf("CorrDiscrContiTraitsC%s", c)
      #MixedCorrelatedTraits <- c(MixedCorrelatedTraits, list(Data$FinalData[, CorrSubColnames]))
      MixedCorrelatedTraits <- c(MixedCorrelatedTraits, list(CorrSubColnames))
    }
    
    if(length(ContiTraits$Discrete) != 0){
      #name <- sprintf("CorrDiscreteTraitsC%s", c)
      #CorrDiscreteTraits <- c(CorrDiscreteTraits, list(Data$FinalData[,ContiTraits$Discrete]))
      CorrDiscreteTraits <- c(CorrDiscreteTraits, list(ContiTraits$Discrete))
    }
    
    if(length(ContiTraits$Continuous) != 0){
      #name <- sprintf("CorrContinuousTraitsC%s", c)
      #CorrContinuousTraits <- c(CorrContinuousTraits, list(Data$FinalData[,ContiTraits$Continuous]))
      CorrContinuousTraits <- c(CorrContinuousTraits, list(ContiTraits$Continuous))
    }
    
    if(uncorrTraitsByCorrGroup != 0){
      
      IndeSubColnames <- subColnames[1:uncorrTraitsByCorrGroup]
      
      IndeTraits <- splitDiscAndContiColnames(IndeSubColnames)
      
      if(length(IndeTraits$Discrete) != 0){
        #name <- paste0("IndDiscreteTraitsC", c)
        #IndDiscreteTraits <- c(IndDiscreteTraits, list(Data$FinalData[,IndeTraits$Discrete]))
        IndDiscreteTraits <- c(IndDiscreteTraits, list(IndeTraits$Discrete))
      }
      
      if(length(IndeTraits$Continuous) != 0){
        #name <- paste0("IndContinuousTraitsC", c)
        #IndContinousTraits <- c(IndContinousTraits, list(Data$FinalData[,IndeTraits$Continuous]))
        IndContinousTraits <- c(IndContinousTraits, list(IndeTraits$Continuous))
      }
    }
  }
  partitions <- list(AllTraits = list(colnames(Data$FinalData)), IndDiscreteTraits = IndDiscreteTraits, IndContinousTraits = IndContinousTraits, MixedCorrelatedTraits = MixedCorrelatedTraits, CorrDiscreteTraits = CorrDiscreteTraits, CorrContinuousTraits = CorrContinuousTraits)
  return(partitions)
}

#Split data
partitions <- dataPartition(Data)


#NA imputation
##############

NAImputation <- function(missingRates, partitions, trees, save = TRUE){
# missingRates: numerical vector corresponding to the rate of missing value to introduce in the data
# partitions: nested list having character vectors corresponding to the data partition
# trees: list of trees of class "phylo"
# save: booléan: if TRUE, save data in a .RData file
# return: a nested list composed of the partitioned data with the 3 kind of missing data (MCRA, MAR and MNAR) according to a precise missing rate and of phylogenetic trees with missing tips.
  
  MCAR <- list()
  MAR <- list()
  MNAR <- list()
  namesMAR <- c()
  namesMCAR <- c()
  namesMNAR <- c()
  mechanism <- c("MCAR", "MAR", "MNAR")
  for(l in seq_along(partitions)){
    
    for(colList in partitions[[l]]){
  
      for(mr in missingRates){
        
        namesMCAR <- c(namesMCAR, paste("MCAR", names(partitions)[l], length(colList), mr ,sep = "/"))
        
        namesMNAR <- c(namesMNAR, paste("MNAR", names(partitions)[l], length(colList), mr ,sep = "/"))
  
        #univariate
        if(length(colList) == 2){
  
          #MCAR
          missingMCAR <- delete_MCAR(Data$FinalData[,colList], mr, colList)
          
          #MNAR
          OneOrTwoColumns <- sample(1:2, 1) #want MNAR values in both columns or not
          if(OneOrTwoColumns == 2){
            cols_misIndex <- c(1,2)
          }
          cols_misIndex <- sample(2, 1)
          missingMNAR <- delete_MNAR_censoring(Data$FinalData[,colList], mr, cols_mis = cols_misIndex, where = "upper") #don't change where arg
          
        }
        
        if(length(colList) == 1){
          oneColum <- as.data.frame(Data$FinalData[,colList])
          names(oneColum) <- colList
          missingMCAR <- delete_MCAR(oneColum, mr, colList)
          missingMNAR <- delete_MNAR_censoring(oneColum, mr, cols_mis = colList, where = "upper") #don't change where arg
        }
        
        if(length(colList) > 2){
          
          if(length(colList) %% 2 == 0){
            #get randomly the number of columns having missing data 
            randomColumns <- sample(length(colList)/2, 1)
            cols_misIndex <- sample(length(colList), randomColumns)
            cols_ctrlIndex <- c(1:length(colList))[-cols_misIndex]
            cols_ctrlIndex <- sample(cols_ctrlIndex, randomColumns)
          }
          if(length(colList) %% 2 != 0 ){
            #get randomly the number of columns having missing data 
            randomColumns <- sample((length(colList)-1)/2, 1)
            cols_misIndex <- sample(length(colList), randomColumns)
            cols_ctrlIndex <- c(1:length(colList))[-cols_misIndex]
            cols_ctrlIndex <- sample(cols_ctrlIndex, randomColumns)
          }
  
          #MCAR
          missingMCAR <- delete_MCAR(Data$FinalData[,colList], mr, colList, p_overall = TRUE)
          
          #MNAR
          missingMNAR <- delete_MNAR_censoring(Data$FinalData[,colList], mr, cols_mis = cols_misIndex, where = "upper")#decide where argument                                                                                                                       we want to.
          }
        
        if(l > 3 & (length(colList) != 1)){
          namesMAR <- c(namesMAR, paste("MAR", names(partitions)[l], length(colList), mr ,sep = "/"))
          
          if(length(colList) == 2){
            #MAR
            cols_misIndex <- sample(1:2, 1)
            cols_ctrlIndex <- c(1,2)[-cols_misIndex]
            missingMAR <- delete_MAR_censoring(Data$FinalData[,colList], mr, 
                                               cols_mis = cols_misIndex, cols_ctrl = cols_ctrlIndex, where = "upper") #don't change where arg
            
          }
          
          else{
            #MAR
            missingMAR <- delete_MAR_censoring(Data$FinalData[,colList], mr, 
                                               cols_mis = cols_misIndex, cols_ctrl = cols_ctrlIndex, where = "upper") #decide where argument                                                                                                                        we want to.
          }
          
          #Highlight the columns used to impute missing data
          colnames(missingMAR)[cols_misIndex] <- paste0(colnames(missingMAR)[cols_misIndex], "/", 1:length(cols_misIndex))
          colnames(missingMAR)[cols_ctrlIndex] <- paste0(colnames(missingMAR)[cols_ctrlIndex], "/", 1:length(cols_ctrlIndex))
          
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
  
  #Drop randomly tips in phylogenetic tree
  TreeDrop <- list()
  namesTrees <- c()
  for (t in 1:length(trees)){
    Ntips <- length(trees[[t]]$tip.label)
    
    for(mr in missingRates){
    namesTrees <- c(namesTrees, paste("TreeDrop", t-1, mr, sep = "/"))
    
    #drop randomly some tips
    dropTree <- geiger::drop.random(trees[[t]], Ntips * mr)
    TreeDrop <- c(TreeDrop, list(dropTree))
    }
  }
  
  #rename the nested list
  names(TreeDrop) <- namesTrees
  
  #Data with NA imputed
  DataNA <- list(MCAR = MCAR, MAR = MAR, MNAR = MNAR, DropTrees = TreeDrop)
  
  #Save data
  ##########
  if(save){
    save(DataNA, file="DataNA.RData")
  }
  
  return(DataNA)
}

#Split data
partitions <- dataPartition(Data)

missingRates <- seq(0.05, 0.40, 0.15)
missingData <- NAImputation(missingRates, partitions, Data$TreeList, save = FALSE)

