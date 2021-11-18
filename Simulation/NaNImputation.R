#install.packages("tidyverse")
#install.packages("missMethods")

library(tidyverse)
library(geiger)
library(missMethods)
library(phangorn)
#https://github.com/torockel/missMethods/blob/master/R/delete_censoring.R


# load simulated data
setwd("C:/Users/Matthieu/Documents/UNIFR/Master_thesis/Scripts/")
load("simulatedData.RData")


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


#' @title Identify species for phylogenetic NAs
#'
#' @description This function returns the names of the species that should get NAs for their traits
#'
#' @usage getTipsNA(Tree, MinTips) 
#'
#' @param Tree Phylogeny 
#' @param MinTips Minimum number of species for which NAs should be created
#' @return vector containing the names of the species that should get NAs for their traits
#' 
#' 
getTipsNA <- function (Tree, MinTips){
  Nodes <- Nnode(Tree)
  Tips <- Ntip(Tree)
  MaxNodeID <- Tips + Nodes
  NodeIDs <- (Tips + 1):MaxNodeID
  EnoughTips <- FALSE
  while (!EnoughTips) {
    FocalNode <- sample(NodeIDs, 1)
    Desc <- phangorn::Descendants(Tree, FocalNode, type = "tips")[[1]]
    SampledTips <- Tree$tip.label[Desc]
    if (length(SampledTips) >= MinTips & length(SampledTips) <= MinTips + 3){ #give me only a clade of a defined size
      EnoughTips <- TRUE
    }
  }
  return(SampledTips)
}

test <- getTipsNA(Data$TreeList$`0`, 6)
test
plot(Data$TreeList$`0`)

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
#'@param trees phylogenetic tree of class "phylo"
#'@param missTraits numerical, number of traits in which there is missing data.
#'@param save boolean, if TRUE, save data in a .RData file
#'@return a nested list composed of the partitioned data with the 3 categories of missing data (MCRA, MAR and MNAR) according to a precise missing rate and of list of specific traits with NaNs imputed according to MCAR (phylogeny).

NaNImputation <- function(missingRates, partitions, tree, missTraits, replicates, save = TRUE){
  
  DataNaN <- list()
  for (r in 1:replicates){
    MCAR <- list()
    MAR <- list()
    MNAR <- list()
    namesMAR <- c()
    namesMCAR <- c()
    namesMNAR <- c()
    mechanism <- c("MCAR", "MAR", "MNAR")
    for(l in seq_along(partitions)){
      
      for(colList in partitions[[l]]){
    
        if (missTraits > length(colList)){
          missTraits <- length(colList)
        }
        
        for(mr in missingRates){
          
          namesMCAR <- c(namesMCAR, paste("MCAR", names(partitions)[l], length(colList), mr ,sep = "/"))
          
          namesMNAR <- c(namesMNAR, paste("MNAR", names(partitions)[l], length(colList), mr ,sep = "/"))
    
          #univariate
          if(length(colList) == 2){
            
            #MCAR
            missingMCAR <- delete_MCAR(Data$FinalData[,colList], mr, missTraits)
            
            #MNAR
            #OneOrTwoColumns <- sample(1:2, 1) #want MNAR values in both columns or not
            #if(OneOrTwoColumns == 2){
            #  cols_misIndex <- c(1,2)
            #}
            #cols_misIndex <- sample(2, 1)
            missingMNAR <- delete_MNAR_censoring(Data$FinalData[,colList], mr, 
                                                 cols_mis = missTraits, where = "upper") #don't change where arg
            
          }
          
          if(length(colList) == 1){
            oneColum <- as.data.frame(Data$FinalData[,colList])
            names(oneColum) <- colList
            missingMCAR <- delete_MCAR(oneColum, mr, missTraits)
            missingMNAR <- delete_MNAR_censoring(oneColum, mr, cols_mis = missTraits, where = "upper") #don't change where arg
          }
          
          if(length(colList) > 2){
            
            if(length(colList) %% 2 == 0){
              
              if(missTraits > length(colList)/2){
                missTraits <- length(colList/2)
              }
              
              #get randomly the number of columns having missing data 
              #randomColumns <- sample(length(colList)/2, 1)
              
              cols_misIndex <- sample(length(colList), missTraits)
              cols_ctrlIndex <- c(1:length(colList))[-cols_misIndex]
              cols_ctrlIndex <- sample(cols_ctrlIndex, missTraits)
            }
            if(length(colList) %% 2 != 0 ){
              
              if(missTraits > (length(colList)-1)/2){
                missTraits <- (length(colList)-1)/2
              }
              #get randomly the number of columns having missing data 
              #randomColumns <- sample((length(colList)-1)/2, 1)
              cols_misIndex <- sample(length(colList), missTraits)
              cols_ctrlIndex <- c(1:length(colList))[-cols_misIndex]
              cols_ctrlIndex <- sample(cols_ctrlIndex, missTraits)
            }
    
            #MCAR
            missingMCAR <- delete_MCAR(Data$FinalData[,colList], mr, colList, p_overall = TRUE)
            
            #MNAR
            missingMNAR <- delete_MNAR_censoring(Data$FinalData[,colList], mr,
                                                 cols_mis = cols_misIndex, where = "upper")#decide where argument                                                                                                                       we want to.
            }
          
          if(l > 3 & (length(colList) != 1)){
            namesMAR <- c(namesMAR, paste("MAR", names(partitions)[l], length(colList), mr ,sep = "/"))
            
            if(length(colList) == 2){
              #MAR
              cols_misIndex <- sample(1:2, 1)
              cols_ctrlIndex <- c(1,2)[-cols_misIndex]
              missingMAR <- delete_MAR_censoring(Data$FinalData[,colList], mr, 
                                                 cols_mis = cols_misIndex, cols_ctrl = cols_ctrlIndex, 
                                                 where = "upper") #don't change where arg
              
            }
            
            else{
              #MAR
              missingMAR <- delete_MAR_censoring(Data$FinalData[,colList], mr, 
                                                 cols_mis = cols_misIndex, cols_ctrl = cols_ctrlIndex, 
                                                 where = "upper") #decide where argument                                                                                                                        we want to.
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
    
    # Drop randomly tips in phylogenetic tree
    PhyloNaN <- list()
    namesTrees <- c()
    lengthTips <- c()
    for (mt in 2:8){
      tips <- getTipsNA(tree, mt)

      if(length(tips) %in% lengthTips){
        next
      }
      
      lengthTips <- c(lengthTips, length(tips))
      for(mr in missingRates){
      namesTrees <- c(namesTrees, paste("PhyloNaN", mt, mr, sep = "/"))
      
      #MCAR
      #randomColums <- sample(ncol(Data$FinalData), missTraits)
      missingMCAR <- delete_MCAR(Data$FinalData[tips, ], mr, 1:ncol(Data$FinalData), p_overall = TRUE)
      
      PhyloNaN <- c(PhyloNaN, list(missingMCAR))
      }
    }
    
    #rename the nested list
    names(PhyloNaN) <- namesTrees
    
    #Data with NA imputed
    DataNaN[[r]] <- list(MCAR = MCAR, MAR = MAR, MNAR = MNAR, PhyloNaN = PhyloNaN)
    
  }
  names(DataNaN) <- seq(1:replicates)
  Parameters <- list(missingRates = missingRates, missTraits = missTraits, replicates = replicates)
  NaNImputed <- list(DataNaN = DataNaN, Parameters = Parameters)
  #Save data
  ##########
  if(save){
    save(NaNImputed, file="DataNaN.RData")
  }
  
  return(NaNImputed)
}

#Split data
partitions <- dataPartition(Data)

missingRates <- seq(0.05, 0.40, 0.15)
missingData <- NaNImputation(missingRates, partitions, Data$TreeList$`0`, 3, 2,  save = F)
