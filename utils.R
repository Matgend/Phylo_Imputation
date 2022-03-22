#utils
library(tidyverse)
library(dplyr)
library(ggstatsplot)
#install.packages('PMCMRplus')
#library(palmerpenguins) #to get the dataset for example -.-


#' @title Calculate the mean of the error for a partition
#'
#' @description This function calculate the mean of error of a partition for. Mean of all each columns of the matrix containing the error value of each traits obtains with each imputation methods. 
#'
#' @usage meanPartition(ErrorOutputObjects)
#'
#' @param ErrorOutputObjects matrix of imputation error. Row: traits, columns: imputation method/
#' variance fraction/missingness
#' @return return list containing the mean of the error for each random approaches and each partition.
meanPartition <- function(ErrorOutputObject){
  MeanOuput <- list()
  
  for(rdn in 2:length(ErrorOutputObject)){
    partitions <- ErrorOutputObject[[rdn]]$Error
    
    partitionMean <- list()
    for(d in 1:length(partitions)){
      meanCols <- colMeans(partitions[[d]][, -1])
      name <- str_split(names(meanCols), "/", simplify = TRUE)
      rowNames <- unique(name[,1])
      colNames <- unique(paste(name[,2], name[,3], sep = "/"))
      meansVector <- 1:(length(rowNames)*length(colNames))
      
      if(rowNames[1] == "imputeDiscrete" | rowNames[1] == "imputeContinuous"){
        missingValtoKeep <- name[which(name[,1] == rowNames[1]),2]
        indexNA <- seq(1, length(meansVector), by = length(rowNames))
        NAcols <- seq(1, length(meansVector), by = length(rowNames)^2)
        for(i in NAcols){
          indexNA <- indexNA[-which(indexNA == i)]
          indexNA <- c(indexNA, (i+1):(i+length(rowNames) - 1))
        }
        sort(indexNA)
        meansVector[indexNA] <- NA
        meansVector[which(!is.na(meansVector))] <- meanCols
        
        #create a matrix to fill
        meanCols <- matrix(data = meansVector, nrow = length(rowNames), ncol = length(colNames))
      }
      
      else{
        meanCols <- t(meanCols)
        dim(meanCols) <- c(length(rowNames), dim(meanCols)[2]/length(rowNames))
      }
      
      rownames(meanCols) <- rowNames
      colnames(meanCols) <- colNames
      partitionMean[[d]] <- meanCols
      names(partitionMean)[d] <- names(partitions[d])
    }
    MeanOuput[[rdn - 1]] <- partitionMean
  }
  names(MeanOuput) <- names(ErrorOutputObject)[-1]
  
  return(MeanOuput) 
}

#' @title download files
#'
#' @description This function downloads the replicate files of same simulation. 
#'
#' @usage loadReplicates(integer, filesInFolder, pathDataframe)
#'
#' @param integer value corresponding to the position of the csv file that we want to select.
#' @param filesInFolder character vector with name of the replicate files
#' @param phatDataframe path where the .csv file are
#' @return return a character vector with all the replicates of the same type of simulated data

loadReplicates <- function(integer, filesInFolder, pathDataframe){
  simDataFolder <- list.files(path = pathDataframe, pattern = "*.csv")
  simDataname <- unique(gsub("\\.csv", "", simDataFolder))
  pattern <- simDataname[integer]
  replicatesFiles <- grep(pattern, filesInFolder, value = T)
  
  return(replicatesFiles)
}

#' @title Mean of all the replicates
#'
#' @description This function calculates the overall mean of all the replicates for the same type of simulated data
#'
#' @usage overallMean(namesReplicates, path, pathReplicates)
#'
#' @param namesReplicates character string with replicate name files
#' @param path path where want to save the data
#' @param pathReplicates path of the directory containing the replicates
#' @return return a .RData object with a matrix for each random approaches and partition with the overall mean error
overallMean <- function(namesReplicates, path, pathReplicates){
  #use the first replicate as the final output
  replicate1 <- get(load(paste0(pathReplicates, namesReplicates[1])))
  meanRep1 <- meanPartition(replicate1)
  #check if empty element in meanRep1 if yes, element removed
  meanRep1 <- deleteEmptylist(meanRep1)
  
  for (i in 2:length(namesReplicates)){
    replicate <- get(load(paste0(pathReplicates, namesReplicates[i])))
    meanRep <- meanPartition(replicate)
    meanRep <- deleteEmptylist(meanRep)
    for (rdn in 1:(length(meanRep)-1)){
      randomType <- meanRep[[rdn]]
      for(d in 1:length(randomType)){
        meanRep1[[rdn]][[d]] <- meanRep1[[rdn]][[d]] + randomType[[d]]
        
      }
    }
  }
  
  #Divide the sums by the number of replicates
  for(rdn in 1:(length(meanRep1)-1)){
    nbrReplicates <- length(namesReplicates)
    namesPartitions <- names(meanRep1[[rdn]])
    meanRep1[[rdn]] <- lapply(names(meanRep1[[rdn]]), function(x) meanRep1[[rdn]][[x]] / nbrReplicates)
    names(meanRep1[[rdn]]) <- namesPartitions
  }
  digitReplicate <- unique(str_extract(namesReplicates, "(?<=\\_)[:digit:]"))
  nameReplicate <- unique(gsub("\\_..*", "", namesReplicates))
  namefile <- file.path(path, sprintf("overallError%s_%s", nameReplicate, digitReplicate))
  save(meanRep1, file = paste0(namefile, ".RData"))
}

setwd("F:/Master_Thesis/Simulations")
filesInFolder <- list.files(path = "./First/FullData2/Results/", pattern = "*.RData")
replicatesInFolder <- loadReplicates(1, filesInFolder[1:10], "F:/Master_Thesis/Cluster/csv/")
namesReplicates <- replicatesInFolder
path <- "C:/Users/Matthieu/Desktop/"
pathReplicates <- "./First/Results/Replicates/"
replicate1 <- get(load(paste0(pathReplicates, namesReplicates[2])))[-1]


meanSDErrorPerTraitsAndTime <- function(namesReplicates, pathReplicates){
  
  replicate1 <- get(load(paste0(pathReplicates, namesReplicates[2])))
  Approaches <- list()
  #create vector to stock mean and sd per methods and trait type
  meanConti <- c()
  sdConti <- c()
  meanDisc <- c()
  sdDisc <- c()

  #define matrix to stock all the execution times
  timeMatrix <- matrix(0, nrow = nrow(replicate1$TimeDataframe), ncol = length(namesReplicates))
  timeMatrix[ ,1] <- replicate1$TimeDataframe$time
  
  #namesPartitions <- c() #4 partitions (MCAR, MAR, MNAR, PhyloMCAR)
  for(rdn in 2:length(replicate1)){
    partitions <- replicate1[[rdn]]$Error
    matrixPartitionMethods <- list()
    
    for(d in 1:length(partitions)){
      data <- partitions[[d]]
      matricesCollection <- c()
      
      for(col in 2:ncol(data)){ #ncol(data) = 10 
        matrixMethods <- matrix(0, nrow = nrow(data), ncol = length(namesReplicates))
        matrixMethods[,1] <- data[, col]
        
        #add the values of replicates
        for(repli in 2:length(namesReplicates)){
          replicate <- get(load(paste0(pathReplicates, namesReplicates[repli])))
          data <- replicate[[rdn]]$Error[[d]]
          matrixMethods[, repli] <- data[, col]
          timeMatrix[, repli] <- replicate$TimeDataframe$time
        }
        
        #calculate mean and sd for each traits
        matrixMethods <- t(matrixMethods)
        colnames(matrixMethods) <- data$trait
        matrixMethods <- gather(as.data.frame(matrixMethods), factor_key = TRUE)
        names(matrixMethods)[1] <- "traits"
        matrixMethods <- matrixMethods %>% group_by(traits) %>%
          summarise(mean = mean(value), sd = sd(value))
        

        #add overall mean and sd in the last row of dataframe
        #continous traits
        contiTraits <- matrixMethods[which(str_detect(matrixMethods$traits, "F")), ]
        #discrete traits
        discTraits <- matrixMethods[which(str_detect(matrixMethods$traits, "F")), ]
        
        if(length(contiTraits) != 0){
          meanConti <- c(meanConti, mean(contiTraits$mean))
          sdConti <- c(sdConti, sd(contiTraits$mean))
        }
        else if(length(contiTraits) == 0){
          meanConti <- c(meanConti, NA)
          sdConti <- c(sdConti, NA)
        }
        
        else if(length(discTraits) != 0){
          meanDisc <- c(meanDisc, mean(discTraits$mean))
          sdDisc <- c(sdDisc, sd(discTraits$mean))
        }
        
        else{
          meanDisc <- c(meanDisc, NA)
          sdDisc <- c(sdDisc, NA)
        }
        
        matricesCollection <- c(matricesCollection, list(matrixMethods))
        names(matricesCollection)[col-1] <- names(data)[col]
      }
      matrixPartitionMethods[[d]] <- matricesCollection
      names(matrixPartitionMethods)[d] <- names(partitions[d])
    }
    Approaches[[rdn-1]] <- matrixPartitionMethods
    names(Approaches)[rdn-1] <- names(replicate1[rdn])
  }
  
  #mean table
  time <- gather(as.data.frame(t(timeMatrix)), factor_key = TRUE)
  names(time)[1] <- "time"
  time <- time %>% group_by(time) %>%
    summarise(mean = mean(value), sd = sd(value))
  
  FinalTable <- replicate1$TimeDataframe
  FinalTable$time <- NULL
  FinalTable$MeanTime <- time$mean
  FinalTable$SDTime <- time$sd
  FinalTable$MeanConti <- meanConti
  FinalTable$SDConti <- sdConti
  FinalTable$MeanDisc <- meanDisc
  FinalTable$SDDisc <- sdDisc
  
  return(list(meanTable = FinalTable, meanPerTaits = Approaches))
}

v <- meanSDErrorPerTraitsAndTime(namesReplicates, pathReplicates)





#' @title Merge replicates for plotting 
#'
#' @description This function merges the dataset of the replicates and preprocess the merged dataset to allow the 
#' conception of a violin plot
#'
#' @usage processDataViolinPlot(namesReplicates, pathReplicates)
#'
#' @param namesReplicates character string with replicate name files
#' @param pathReplicates path of the directory containing the replicates
#' @return return an object of nested list containing dataframes with 2 columns (methods, and error value)
processDataViolinPlot <- function(namesReplicates, pathReplicates){
  
  replicate1 <- get(load(paste0(pathReplicates, namesReplicates[1])))[-1]
  dataPartitionRandomApproaches <- list()
  #namesPartitions <- c() #4 partitions (MCAR, MAR, MNAR, PhyloMCAR)
  for(rdn in 1:length(replicate1)){
    partitions <- replicate1[[rdn]]$Error
    partitionMerged <- c()
    for(d in 1:length(partitions)){
      
      mergedData <- replicate1[[rdn]]$Error[[d]]
      #add the values of replicates
      for(repli in 2:length(namesReplicates)){
        replicate <- get(load(paste0(pathReplicates, namesReplicates[repli])))[-1]
        dataset <- replicate[[rdn]]$Error[[d]]
        mergedData <- rbind(mergedData, dataset)
      }
      mergedData <- pivot_longer(mergedData, cols = (names(mergedData[-1])))
      names(mergedData)[2] <- "methods"
      partitionMerged <- c(partitionMerged, list(mergedData))
      names(partitionMerged)[d] <- names(replicate1[[rdn]]$Error[d])
    }
    dataPartitionRandomApproaches[[rdn]] <- partitionMerged
    names(dataPartitionRandomApproaches)[rdn] <- names(replicate1[rdn])
  }
  return(dataPartitionRandomApproaches)
}
  

#' @title Plot violin plot 
#'
#' @description This function plots violin plot using the function ggbetweenstats from the ggstatsplot package.
#'
#' @usage violinPlot(datasets, save)
#'
#' @param datasets object containing nested lists composed of datasets
#' @param save character vector defining the directory where to save the plots
#' @param namePDF character vector corresponding to the file name
#' @return return a violin plot are a pdf with all the violin plots in.
violinPlot <- function(datasets, save = NULL, namePDF = NULL){
  
  if(!is.null(save) & is.null(namePDF)){
    stop("Should provide a file name")
  }
  
  if(!is.null(save)){
    pdf(paste0(save, namePDF, ".pdf"))
  }
  
  for(rdn in 1:length(datasets)){
    rdn <- 1
    partition <- datasets[[rdn]]
    for(p in 1:length(partition)){
      title <- names(partition[p])
      data <- partition[[p]]
      plt <- ggbetweenstats(
        data = data,
        x = methods,
        y = value, 
        pairwise.comparisons = FALSE,
        package = "RColorBrewer",
        palette = "Paired", 
        title = title, 
        results.subtitle = FALSE
      )
      plot(plt)
    }
  }
  
  if(!is.null(save)){
    dev.off()
  }
  
  if(is.null(save)){
    plot(plt)
    return(plt)
  }
}

#violinPlot(v)

replicate1 <- get(load(paste0(pathReplicates, namesReplicates[2])))
replicate1$TimeDataframe
