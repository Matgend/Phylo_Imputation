#utils
library(tidyverse)
library(dplyr)
library(ggstatsplot)
library(ggbeeswarm)


#ERROR
######
#' @title Calculate error of imputation
#'
#' @description This function calculates the RMSE for imputed continuous data, the absolute 
#' error for imputed ordinal data and misclassification for the other subcategories of discrete data
#'
#' @usage imputationError(imputedData, trueData, missingData, imputationApproachesName, Data)
#'
#' @param imputedData array of imputed data
#' @param trueData array of true data
#' @param missingData array of data with missing values
#' @param Data simulated Data object
#' @return return a data.frame with in the first column the trait names and in the second the errors
#'
imputationError <- function(imputedData, trueData, missingData, imputationApproachesName, Data){
  
  # #change imputation approach name
  # imputationApproachesName[str_detect(imputationApproachesName, "MICE")] <- "MICE"
  # imputationApproachesName[str_detect(imputationApproachesName, "Miss")] <- "MissForest"
  # imputationApproachesName[str_detect(imputationApproachesName, "KNN")] <- "KNN"
  # imputationApproachesName[str_detect(imputationApproachesName, "gain")] <- "GAIN"
  
  
  #get the ordinal trait reference
  ordinalTraits <- which(Data$dataframe$class == "ordinal") #give the row in dataframe which correspond to /n in data names
  errors <- c()
  traitNames <- c() 
  for (c in 1:ncol(missingData)){
    
    #know is NaNs in the columns(trait)
    NaNRowIndex <- which(is.na(missingData[,c]))
    
    if(length(NaNRowIndex != 0)){
      
      traitNames <- c(traitNames, names(trueData[c]))
      missingValues <- missingData[NaNRowIndex, c]
      trueValues <- trueData[NaNRowIndex, c]
      imputedValues <- imputedData[NaNRowIndex, c]
      
      #in case continuous data
      if(length(grep("F.", names(missingData)[c])) == 1){
        #print("continuous")
        #rmse
        error <- sqrt(mean((as.numeric(imputedValues) - as.numeric(trueValues))^2))
        
      }
      
      #in case ordinal trait
      else if(length(ordinalTraits) != 0 & length(grep(paste0("/", ordinalTraits), names(missingData)[c])) == 1){
        #print("ordinal")
        #imputation error for ordinal traits (absolute error)
        error <- mean(abs((as.numeric(imputedValues) - as.numeric(trueValues)) / as.numeric(trueValues)))
      }
      
      #in case discrete data
      else{
        #print("discrete")
        error <- (sum(as.numeric(imputedValues) != as.numeric(trueValues)) / length(trueValues)) 
        #error <- length(setdiff(imputedValues, trueValues)) / length(trueValues)
      }
      errors <- c(errors, error)
    }
  }
  output <- data.frame(trait = traitNames, c2 = errors)
  names(output)[2] <- imputationApproachesName
  
  return(output)
}

#' @title Delete empty list
#'
#' @description This function delete lists of list that are empty 
#'
#' @usage deleteEmptylist(lists)
#'
#' @param lists list containing lists (lists of list)
#' @return return the list of lists without the empty list

deleteEmptylist <- function(lists){
  saveIndex <- c()
  for (l in 1:length(lists)){
    if(length(lists[[l]]) == 0){
      saveIndex <- c(saveIndex, l)
    }
  }
  if (length(saveIndex) != 0){
    lists[[saveIndex]] <- NULL
  }
  return(lists)
}


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

# setwd("F:/Master_Thesis/Simulations")
# filesInFolder <- list.files(path = "./First/FullData2/Results/", pattern = "*.RData")
# replicatesInFolder <- loadReplicates(1, filesInFolder[1:10], "F:/Master_Thesis/Cluster/csv/")
# namesReplicates <- replicatesInFolder
# path <- "C:/Users/Matthieu/Desktop/"
# pathReplicates <- "./First/Results/Replicates/"
# replicate1 <- get(load(paste0(pathReplicates, namesReplicates[2])))[-1]

#namesReplicates <- replicatesInFolder

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
          #print(namesReplicates[repli])
          data <- replicate[[rdn]]$Error[[d]]
          #print(data)
          matrixMethods[, repli] <- data[, col]
          timeMatrix[, repli] <- as.numeric(replicate$TimeDataframe$time)
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
        discTraits <- matrixMethods[which(str_detect(matrixMethods$traits, "I")), ]
        #print(nrow(discTraits))
        if(nrow(contiTraits) != 0){
          #print("conti")
          meanConti <- c(meanConti, mean(contiTraits$mean))
          sdConti <- c(sdConti, sd(contiTraits$mean))
        }
        if(nrow(contiTraits) == 0){
          #print("noconti")
          meanConti <- c(meanConti, NA)
          sdConti <- c(sdConti, NA)
        }
        
        if(nrow(discTraits) != 0){
          #print("Disc")
          meanDisc <- c(meanDisc, mean(discTraits$mean))
          sdDisc <- c(sdDisc, sd(discTraits$mean))
        }
        
        if(nrow(discTraits) == 0){
          #print("noDisc")
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
  
  #add columns with group of imputation method (class)
  class <- rep(NA, nrow(FinalTable))
  class[grep("^c|^R", FinalTable$impute_Approach)] <- "Phylogenetic imputation"
  class[which(str_detect(FinalTable$impute_Approach, pattern = "^G"))] <- "Deep learning"
  class[which(str_detect(FinalTable$impute_Approach, pattern = "PI"))] <- "2-step"
  class[which(is.na(class))] <- "Machine learning"
  
  FinalTable$Class <- class
  
  return(list(meanTable = FinalTable, meanPerTaits = Approaches))
}

#' @title Change names of methods
#'
#' @description This function changes the names of the methods, adding a "+ P" in case 
#' phylogenetic informations are included in the imputation 
#'
#' @usage changeMethodNames(methodVecotr)
#'
#' @param methodVector string composed of characters with the following structure: "method name/missingness/amount of
#'  phylogenetic information"
#' @return return a string composed of characters with the following structure: if no phylogenetic information, "method
#'  name" if phylogenetic information: "method name+P" 
changeMethodNames <- function(methodVector){
  splitVector <- str_split(methodVector, "/", simplify = TRUE)
  for(n in 1:length(methodVector)){
    if(splitVector[n,3] == 0){
      methodVector[n] <- splitVector[n, 1]
    }
    
    else if(splitVector[n,3] == 2){
      methodVector[n] <- paste0(splitVector[n,1],"+IP")
    }
    else{
      methodVector[n] <- paste0(splitVector[n,1],"+P")
    }
  }
  return(methodVector)
}


# v <- meanSDErrorPerTraitsAndTime(namesReplicates, pathReplicates)
# 
# NaNData <- get(load("005_05/MissingData/NaNDataContinuousCorData10_R32_V0.05-0.5.RData"))
# NaNData$DataNaN$MCAR$`MCAR/AllTraits/10/0.05`
# 
# mymatrix <- NaNData$DataNaN$MAR$`MAR/CorrContinuousTraits/10/0.05`
# colnames(mymatrix)[colSums(is.na(mymatrix)) > 0]
# 
# namesReplicates <- replicatesInFolder[1:2]
# pathReplicates
# p <- processDataViolinPlot(namesReplicates, pathReplicates)


#' @title Remove corrupted objects 
#'
#' @description This function detect simulation that failed and remove the files
#'
#' @usage checkFiles(namesReplicates, pathReplicates)
#'
#' @param namesReplicates character string with replicate name files
#' @param pathReplicates path of the directory containing the replicates
#' @return return an object of nested list containing dataframes with 2 columns (methods, and error value)

checkFiles <- function(namesReplicates, pathReplicates){
  replicate1 <- get(load(paste0(pathReplicates, namesReplicates[1])))[-1]
  namesR <- c()
  
  for(rdn in 1:length(replicate1)){
    partitions <- replicate1[[rdn]]$Error
    for(d in 1:length(partitions)){
      
      print(names(partitions))
      #add the values of replicates
      for(repli in 2:length(namesReplicates)){
        #print(namesReplicates[repli])
        replicate <- get(load(paste0(pathReplicates, namesReplicates[repli])))[-1]
        dataset <- replicate[[rdn]]$Error[[d]]
        
        if(sum(is.na(dataset)) > 0 & rdn != length(replicate)){
          #print(namesReplicates[repli])
          print(names(which(colSums(is.na(dataset)) > 0)))
          if(names(which(colSums(is.na(dataset)) > 0)) != "GAIN/0.5/0"){
            print(namesReplicates[repli])
          }
          
          namesR <- c(namesR, namesReplicates[repli])
        }
        
      }
      # print(names(which(colSums(is.na(mergedData)) > 0)))
      print(length(namesR))
    }
  }
  return(namesR) #list(..., namesR)
}





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
#' 
#namesReplicates <- replicatesInFolder
processDataViolinPlot <- function(namesReplicates, pathReplicates){
  replicate1 <- get(load(paste0(pathReplicates, namesReplicates[1])))[-1]
  dataPartitionRandomApproaches <- list()

  for(rdn in 1:length(replicate1)){

    partitions <- replicate1[[rdn]]$Error
    partitionMerged <- c()
    for(d in 1:length(partitions)){
      
      mergedData <- replicate1[[rdn]]$Error[[d]]
      
      #add the values of replicates
      for(repli in 2:length(namesReplicates)){
        #print(namesReplicates[repli])
        replicate <- get(load(paste0(pathReplicates, namesReplicates[repli])))[-1]
        dataset <- replicate[[rdn]]$Error[[d]]

        #print(namesReplicates[repli])
        mergedData <- rbind(mergedData, dataset)

      }

      mergedData <- pivot_longer(mergedData, cols = (names(mergedData[-1])))
      names(mergedData)[2] <- "methods"
      
      #change name methods
      mergedData$methods <- changeMethodNames(mergedData$methods)
    
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
#' @param fileName character vector corresponding to the filename (name of the csv)
#' @param myPalette vector defining the color of each method
#' @return return a violin plot are a pdf with all the violin plots in.
violinPlot <- function(datasets, save, fileName, myPalette){

  for(rdn in 1:length(datasets)){
  
    partition <- datasets[[rdn]]
    for(p in 1:length(partition)){
  
      title <- names(partition[p])
      dataP <- partition[[p]]
      
      #remove row with NaN in dataP
      indexNA <- which(is.na(dataP$value))
      #print(indexNA)
      if(length(indexNA) != 0){
        dataP <- dataP[-indexNA, ]
      }
    
      #Indentify is mixed data and set the ylab
      continuousTraits <- dataP[grep("F.", dataP$trait), ] 
      discreteTraits <- dataP[grep("I.", dataP$trait), ]
      
      partitions <- list(continuousTraits, discreteTraits)
      errors <- c()
      if(nrow(continuousTraits) != 0){
        errors <- c(errors, "RMSE")
        ggplot.component <- ggplot2::scale_y_continuous(breaks = seq(1, 6, 1), limits = (c(-0.2, 6)))
        
      }
      
      if(nrow(discreteTraits) != 0){
        errors <- c(errors, "misclassification")
        ggplot.component <- ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2), limits = (c(0, 1)))
      }
  
      if(nrow(continuousTraits) == 0){
        partitions[[1]] <- partitions[[2]]
        partitions[[2]] <- NULL
      }
      
      if(nrow(discreteTraits) == 0){
        partitions[[2]] <- NULL
        
      }

      #ask for ordinal traits?
      for(i in 1:length(errors)){
        
        title <- paste(title, errors[i], sep = "_")
        
        # if(length(errors == 2)){
        #   title <- paste(title, errors[i], sep = "_")
        # }
  
        title <- str_replace_all(title, "/", "_")
        namePDF <- paste(paste0(save, fileName), title[i], sep = "_")
        data <- partitions[[i]]
        
        # if(errors[i] == "misclassification"){
        #   limY <- c(0,1)
        # }
        # else{
        #   limY <- c(0, 3)
        # }
        #print(data$value)
        # plt <- (ggplot(data, aes(x = methods, y = value, fill = methods)) +
        #           geom_violin(alpha = 0.5) +
        #           geom_beeswarm(groupOnX = TRUE) +
        #           theme(legend.position = "none") +
        #           scale_fill_manual(values = myPalette, drop = T) +
        #           ggtitle(title)+
        #           ylab(errors[i]))
        
        p <- (ggplot(data, aes(x = methods, y = value)) + 
          geom_violin(aes(fill = methods), trim = FALSE) +
          geom_boxplot(width = 0.3) +
          theme_classic() +
          scale_fill_manual(values = myPalette, drop = T) +
          xlab("Methods") +
          ylab(errors[i]) +
          theme(legend.position = "none") +
          theme(axis.text=element_text(size=25),
                axis.title=element_text(size=25,face="bold")))
          #ylim(limY))
        
        dat <- ggplot_build(p)$data[[2]]
        
        #print(dat)
        
        plt <- p + geom_segment(data=dat, aes(x=xmin, xend=xmax, 
                                              y=middle, yend=middle), colour="red", size=2)
         
        # plt <- (ggplot(data, aes(x = methods, y = value, fill = methods)) +
        #           geom_violin(width=1.4) +
        #           geom_boxplot(width=0.1, color= "grey", alpha=0.2) +
        #           scale_fill_viridis(discrete = TRUE) +
        #           theme_ipsum() +
        #           theme(legend.position="none", plot.title = element_text(size=11)) +
        #           xlab("Methods") +
        #           ylab(errors[i])+
        #           scale_fill_manual(values = myPalette, drop = T))
        
        # plt <- (ggbetweenstats(
        #   data = data,
        #   x = methods,
        #   y = value,
        #   ylab = errors[i],
        #   plot.type = "boxviolin",
        #   centrality.plotting = F,
        #   centrality.type = "nonparametric",
        #   pairwise.comparisons = FALSE,
        #   results.subtitle = FALSE) +
        #   xlab("Methods") +
        #   ylab(errors[i]) +
        #   scale_color_manual(values = myPalette, drop = T) +
        #   theme(axis.text=element_text(size=25),
        #         axis.title=element_text(size=25,face="bold")) +
        #   ggplot.component)

        pdf(paste0(namePDF, ".pdf"), width = 32, height = 18)
        plot(plt)
        dev.off()
        
        #p<-ggplot(data, aes(x=methods, y=log10(value))) +
        #  scale_color_manual(values = myPalette, drop = T) +
        #  geom_violin(trim=FALSE)
        # plot(p)
      }
    }
  }
}


#' @title Process data for barplot
#'
#' @description This function calculate the number of times that a methods imputes a value better than an other
#' 
#' @usage preocessDataBarplot(processedData)
#'
#' @param porcessedData nested list of dataframe processed by the function processDataViolinPlot()
#' @param nameRepliate character defining the name of the saved csv
#' @param pathSaveTable directory path, where to save all the tables
#' @return return a list of nested list containing dataframe composed of 3 columns: methods, score, and accuracy. 
#' Score represents the number of time a methods performs better than the others. Accuracy is the score in frequency. Save the table as a csv.
processDataBarplot <- function(processedData, nameRepliate, pathSaveTable = NULL){
  
  #datasets is processed data
  Approach <- list()
  for(rdn in 1:length(processedData)){
    randomApp <- processedData[[rdn]]
    scorePartition <- list()
    scorePartitionName <- c()
    for(p in 1:length(randomApp)){
      partition <- randomApp[[p]]
      uniqueMethods <- unique(partition$methods)
      dataframeScore <- data.frame(methods = as.factor(uniqueMethods), 
                                   score = 0, accuracy = 0)
      uniqueTraits <- unique(partition$trait)
      for(t in 1:length(uniqueTraits)){
        byTraits <- partition[which(partition$trait == uniqueTraits[t]), ]
        
        #split data
        nRow <- nrow(byTraits)
        step <- length(uniqueMethods)
        for(r in seq(step, nRow, by = step)){
          splittedData <- byTraits[(r-step+1):r, ]
          
          
          #find smallest value
          indexMinVal <- which(splittedData$value == min(splittedData$value, na.rm = TRUE))
          
          #add 1 in result dataframe
          winningMethod <- splittedData$methods[indexMinVal]
          
          add <- 1
          
          #if several "best methods"
          if(length(winningMethod) > 1){
            add <- add/length(winningMethod)
          }
          
          dataframeScore$score[which(dataframeScore$methods %in% winningMethod)] <- 
            dataframeScore$score[which(dataframeScore$methods %in% winningMethod)] + add
          
          #print(add)
          
        }
      }
      dataframeScore$accuracy <- dataframeScore$score / (nrow(partition) / step)
      scorePartition <- c(scorePartition, list(dataframeScore))
      names(scorePartition)[p] <- names(randomApp)[p]
      
      if(!is.null(pathSaveTable)){
        nameTxt <- paste0("Table_Barplot_", nameRepliate,"_", names(scorePartition)[p], ".txt")
        nameTxt <- str_replace_all(nameTxt, "/", "_")
        nameTxt <- paste0(pathSaveTable, nameTxt)
        write.table(dataframeScore, nameTxt, sep = "\t", col.names = TRUE, row.names = FALSE)
      }

    }
    Approach[[rdn]] <- scorePartition
    names(Approach)[rdn] <- names(processedData)[rdn]
  }
  return(Approach)
}

#' @title Plot barplot 
#'
#' @description This function plots barplots, each color correspond to an imputation method. The methods are displayed on 
#' the x axis and the accuracy on the y axis.
#'
#' @usage barplotBestMethod(processData, save, fileName, myPalette)
#'
#' @param porcessData nested list of dataframe processed by the function processDataBarplot()
#' @param save character vector defining the directory where to save the plots
#' @param fileName character vector corresponding to the filename (name of the csv)
#' @param myPalette vector defining the color of each method
#' @return return a violin plot are a pdf with all the violin plots in.
barplotBestMethod <- function(processedData, save, fileName, myPalette){

  methodsNames <- processedData$MCAR[[1]]$methods
  for(rdn in 1:length(processedData)){
    partition <- processedData[[rdn]]
    
    for(p in 1:length(partition)){
      title <- names(partition[p])
      title <- str_replace_all(title, "/", "_")
      namePDF <- paste(paste0(save, fileName), title, sep = "_")
      plt <- (ggplot(partition[[p]], 
                     aes(x = methods, y = accuracy, fill = methods)) + 
                geom_bar(stat = "identity") +
                scale_fill_manual(values = myPalette, drop = T)+
                coord_cartesian(ylim = c(0, 1))+
                #geom_text(aes(label=round(accuracy,2)), vjust=-0.3, size=3.5)+
                xlab("Methods")+
                ylab("Best impuation (n = 10T * 100S)")+
                theme(legend.position="none", axis.text=element_text(size=25),
              axis.title=element_text(size=25,face="bold")))
      pdf(paste0(namePDF, ".pdf"), width = 32, height = 18)
      plot(plt)
      dev.off() 
    }
  }
}


#' @title Process data for barplots, best methods comparison for 1 scenario (one csv file)
#'
#' @description This function saves in a table the best score (RMSE or misclassification) for barplot plotting
#' 
#' @usage tableBestMethodUnique(pathResultTable, nameReplicate)
#'
#' @param csvName character defining the name of the saved csv
#' @param pathResultTable directory path, where to save all the tables and load the tables
#' @return return a table which will be used as input for the bestMethodBarplot function. 
# pathResultTable <- "./Overall_Results/Table/"
# # 
# csvName[6]
# # csvName <- csvName[3]
# # 
# # pathSaveMeanTable
# tableBestMethodUnique(pathResultTable, csvName = csvName[6])

tableBestMethodUnique <- function(pathResultTable, csvName){
  
  tablesToMerge <- list.files(path = pathResultTable, pattern = csvName)
  
  for(f in 1:length(tablesToMerge)){
    #f <- 1
    #merged the tables to have a big table.
    if(f == 1){
      mergedTable <- read.table(paste0(pathResultTable, tablesToMerge[f]), header = TRUE)
    }
    
    else{
      table <- read.table(paste0(pathResultTable, tablesToMerge[f]), header = TRUE)
      mergedTable <- rbind(mergedTable, table)
    }
  }
  
  #loop through random_Type
  randomType <- unique(mergedTable$random_Type)
  for(rdn in 1: length(randomType)){
    #rdn <- 1
    
    #loop through the partition
    partition <- unique(mergedTable$partition[which(mergedTable$random_Type == randomType[rdn])])

    #in case alltraits and the second partition are duplicates.
    if(length(partition) == 2 & (randomType[rdn] != "PhyloNaN")){
      partition <- partition[-1]
    }
    
    else if(randomType[rdn] == "PhyloNaN"){
      partition <- "PhyloNA"
    }
    
    for(p in 1:length(partition)){
      #p <- 1
      
      namesRmatrix <- c("methods", "mc5", "mc33","mc50", "md5", "md33", "md50", "class") 
      Rmatrix <- as.data.frame(matrix(NA, nrow = length(class), ncol = length(namesRmatrix)))
      colnames(Rmatrix) <- namesRmatrix
      #loop through the missing rate
      missingRate <- unique(mergedTable$missing_Degree)
      for(m in 1:length(missingRate)){
        #m <- 1
        #m <- 2
        #loop through class
        class <- unique(mergedTable$Class)
        
        
        
        #namesColFinalData <- c("methods", "partition" , "mc5", "mc30", 
        #                       "mc50", "md5", "md33", "md50", "class", "phy_sig", "miss_mech", "signal_param") 
        
        #add at the end
        
        for (c in 1:length(class)){
          
          #c <- 1
          #create subset
          subdata <- mergedTable[which(mergedTable$random_Type == randomType[rdn] &
                                         mergedTable$partition == partition[p] & 
                                         mergedTable$missing_Degree == missingRate[m] & 
                                         mergedTable$Class == class[c]), ]
          
          if(randomType[rdn] == "PhyloNaN"){
            subdata <- mergedTable[which(mergedTable$random_Type == randomType[rdn] &
                                           mergedTable$missing_Degree == missingRate[m] & 
                                           mergedTable$Class == class[c]), ]
          }
          
          Rmatrix[c, "class"] <- unique(subdata$Class)
          
          if(nrow(subdata) != sum(is.na(subdata$MeanConti))){
            bestConti <- subdata[which(subdata$MeanConti == min(subdata$MeanConti, na.rm = TRUE)), ]
            
            mConti <- c("mc5", "mc33", "mc50")
            
            if(m == 1){
              Rmatrix[c, mConti[m]] <-  bestConti$MeanConti
              Rmatrix[c, "methods"] <- bestConti$impute_Approach
            }
            
            else{

              #print(bestDisc$impute_Approach %in% Rmatrix$methods)
              if(bestConti$impute_Approach %in% Rmatrix$methods){
                indexRow <- which(Rmatrix$methods == bestConti$impute_Approach)
                Rmatrix[indexRow, mConti[m]] <-  bestConti$MeanConti
              }
              
              else{
                Rmatrix[nrow(Rmatrix)+1, "methods"] <- bestConti$impute_Approach
                Rmatrix[nrow(Rmatrix), mConti[m]] <-  bestConti$MeanConti
                Rmatrix[nrow(Rmatrix), "class"] <- bestConti$Class
              }
            }
          }
          
          #else if(nrow(subdata) != sum(is.na(subdata$MeanDisc))){
          else{
            bestDisc <- subdata[which(subdata$MeanDisc == min(subdata$MeanDisc, na.rm = TRUE)), ]
            mDisc <- c("md5", "md33", "md50")
            
            
            if(m == 1){
              Rmatrix[c, mDisc[m]] <-  bestDisc$MeanDisc
              Rmatrix[c, "methods"] <- bestDisc$impute_Approach
            }
            
            else{
              
              #print(bestDisc$impute_Approach %in% Rmatrix$methods)

              if(bestDisc$impute_Approach %in% Rmatrix$methods){
                indexRow <- which(Rmatrix$methods == bestDisc$impute_Approach)
                Rmatrix[indexRow, mDisc[m]] <-  bestDisc$MeanDisc
              }
              
              else{
                Rmatrix[nrow(Rmatrix)+1, "methods"] <- bestDisc$impute_Approach
                Rmatrix[nrow(Rmatrix), mDisc[m]] <-  bestDisc$MeanDisc
                Rmatrix[nrow(Rmatrix), "class"] <- bestDisc$Class
                
              }
            }
          }
        }
      }
      #remove empty rows
      Rmatrix <- Rmatrix[, colSums(is.na(Rmatrix)) != nrow(Rmatrix)]
      Rmatrix$partition <- partition[p]
      
      #add columns
      Rmatrix$miss_mech <- randomType[rdn]
      "phy_sig"
      "signal_param"
      
      if(p == 1 & randomType[rdn] == "MCAR"){
        finalTable <- Rmatrix
      }
      else{
        finalTable <- rbind(finalTable, Rmatrix)
      }
    }
  }
  
  if("L" %in% strsplit(csvName[1], split = "")[[1]]){
    finalTable$phy_sig <- "Weak"
    finalTable$signal_param <- "L"
  }
  else if("K" %in% strsplit(csvName[1], split = "")[[1]]){
    finalTable$phy_sig <- "Weak"
    finalTable$signal_param <- "K"
  }
  else{
    finalTable$phy_sig <- "Strong"
    finalTable$signal_param <- "N"
  }
  
  return(finalTable) 
}



#' @title Process data for barplots, best methods comparison for several scenarios (several csv files)
#'
#' @description This function saves in a table the best score of several scenarios in different txt files (RMSE or misclassification) for barplot plotting

#' @usage tableBestmethodSeveral(pathResultTable, nameReplicate)
#'
#' @param csvName character defining the name of the saved csv
#' @param pathResultTable directory path, where to save all the tables and load the tables
#' @return save a table in .txt format or return which will be used as input for the bestMethodBarplot function. 

#pathResultTable = pathSaveMeanTable

tableBestmethodSeveral <- function(csvName, pathResultTable){
  diffPhyloFiles <- str_extract(csvName, ".+?([0-9])")
  diffPhyloFilesUnique <- unique(diffPhyloFiles)
  #fileCount <- table(diffPhyloFiles)
  
  for(f in 1:length(diffPhyloFilesUnique)){
    
    indexFiles <- grep(diffPhyloFilesUnique[f], csvName)
    
    if(length(indexFiles) > 1){
      
      for(i in 1:length(indexFiles)){
        
        if(i == 1){
          finalTable <- tableBestMethodUnique(pathResultTable, csvName[indexFiles[i]])
        }
        else{
          table <- tableBestMethodUnique(pathResultTable, csvName[indexFiles[i]])
          finalTable <- rbind(finalTable, table)
        }
      }
    }
    
    else{
      finalTable <- tableBestMethodUnique(pathResultTable, csvName[indexFiles])
    }
    
    #Generate TXT with mean + sd of each imputation methods
    filename <- diffPhyloFilesUnique[f]
    fileNameTxt <- unique(paste0(pathResultTable, "TableBestMethod_", filename ,".txt"))
    print(fileNameTxt)
    write.table(finalTable, fileNameTxt, sep = "\t", row.names = FALSE)
  }
}