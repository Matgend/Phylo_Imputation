#' @title Check + conversion columns class. 
#' 
#' @description This function checks if character and integer columns are factors and that the columns composed of float are
#' numeric. If it is not the case the columns are converted.
#' 
#' @usage checkConvert(missingData)
#' 
#' @param missingData array of data with missing values
#' 
#' @return The missing dataset with the columns containing integers or characters as factors and columns containing float as numeric
#' 
checkConvert <- function(missingData){
  
  for(c in 1:ncol(missingData)){
    
    if(class(missingData[ ,c]) == "character" | class(missingData[,c]) == "integer"){
      missingData[, c] <- as.factor(missingData[, c])
    }
    
    else if(class(missingData[,c]) == "numeric"){
      
      #check if composed only of integers 
      nbrIntegers <- (sum(missingData[, c] - floor(missingData[,c]) == 0, na.rm = T))
      
      if(nbrIntegers == length(missingData[,c])){
        missingData[, c] <- as.factor(missingData[, c])
      }
    }
  }
  return(missingData)
}

#' @title Delete empty list
#'
#' @description This function deletes lists of list that are empty 
#'
#' @usage deleteEmptylist(lists)
#'
#' @param lists list containing lists (lists of list)
#' @return A list of lists without the empty list
#' 
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


#' @title Calculation of the average error per scenario (partition).
#'
#' @description This function calculates the average error of a scenario (partition) among all the replicates.
#'
#' @usage meanPartition(ErrorOutputObjects)
#'
#' @param ErrorOutputObjects matrix containing the imputations error values. Row: traits, columns: imputation method
#' /variance fraction/missingness
#' @return A list containing the average error for each random approach and each scenario(partition).
#' 
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
#' @param integer value corresponding to the position of the csv file in the variable "filesInFolder" that we want to select.
#' @param filesInFolder character vector with name of the replicate files
#' @param phatDataframe path where the .csv file are
#' @return A character vector with all the replicates of the same type of simulated data
#' 
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
#' @param trait boolean, if NULL returns the error value for all the traits.
#' @return A list with a matrix for each random approaches and partition with the overall mean error
#' 
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



#' @title Merge methods error by missing mechanisms and replicate 
#'
#' @description This function merges the dataset of the errors
#'
#' @usage mergeErrorData(listOfDatasets)
#'
#' @param listOfDatasets list of datasets to merge
#' @return A dataframe of size number of traits * number of imputation methods
#' 
mergeErrorData <- function(listOfDatasets){
  
  mergedData <- listOfDatasets[[1]]
  
  for(d in 2:length(listOfDatasets)){
    mergedData <- merge(mergedData, listOfDatasets[[d]], by = "trait")
  }
  return(mergedData)
  
}

#' @title Merge replicates for plotting 
#'
#' @description This function merges the dataset of the replicates 
#'
#' @usage processData(namesReplicates, pathReplicates, trait)
#'
#' @param namesReplicates character string with replicate name files
#' @param pathReplicates path of the directory containing the replicates
#' @param trait boolean, if NULL returns the error value for all the traits.
#' @return An object of nested list containing dataframes with 2 columns (methods, and error value)
#' 
processData <- function(namesReplicates, pathReplicates, trait = NULL){
  replicate1 <- get(load(paste0(pathReplicates, namesReplicates[1])))$errorData
  
  dataRandomApproaches <- list()
  
  for(rdn in 1:length(replicate1)){
    
    partitions <- replicate1[[rdn]][[1]]
    
    mergedData <- mergeErrorData(partitions)

    #add the values of replicates
    for(repli in 2:length(namesReplicates)){
      #print(namesReplicates[repli])
      replicate <- get(load(paste0(pathReplicates, namesReplicates[repli])))$errorData
      dataset <- replicate[[rdn]][[1]]
      dataset <- mergeErrorData(dataset)
      
      #print(namesReplicates[repli])
      mergedData <- rbind(mergedData, dataset)
      
    }
      
    mergedData <- pivot_longer(mergedData, cols = (names(mergedData)[-1]))
    names(mergedData)[2] <- "Methods"
    
    
    #change names
    mergedData$Methods <- gsub("_NP", "", mergedData$Methods)
    mergedData$Methods <- gsub("_P", "+P", mergedData$Methods)
    mergedData$Methods <- gsub("_2-step", "+PI", mergedData$Methods)
      
    if(!is.null(trait)){
      mergedData <- mergedData[which(mergedData$trait == trait), ]
    }
    

    dataRandomApproaches[[rdn]] <- mergedData
    names(dataRandomApproaches)[rdn] <- names(replicate1[rdn])
  }
  return(dataRandomApproaches)
}

#' @title Plot boxplot
#'
#' @description This function generates a boxplot using the function ggbetweenstats from the ggstatsplot package.
#'
#' @usage boxPlotByMethods(datasets, fileName, myPalette, returnPlot  = FALSE , save = NULL)
#'
#' @param datasets object containing nested lists composed of datasets
#' @param fileName character vector corresponding to the filename (name of the csv)
#' @param myPalette vector defining the color of each method
#' @param returnPlot boolean specifying if want to return the plot and the file name
#' @param save character vector defining the directory where to save the plots
#' @return A violin plot are a pdf with all the violin plots in.
# 
boxPlotByMethods <- function(datasets, fileName, myPalette, returnPlot  = FALSE , save = NULL){
  
  boxplot.list <- list()
  i <- 1
  for(rdn in 1:length(datasets)){
    
    mmName <- names(datasets)[rdn]
    
    data <- datasets[[rdn]]
    #mmName <- strsplit(names(partition), "_")[[1]][1]
    
    #remove row with NaN in data
    indexNA <- which(is.na(data$value))
    #print(indexNA)
    if(length(indexNA) != 0){
      data <- data[-indexNA, ]
    }
    
    #divide data by methods
    diffMethods <- unique(data$Methods)
    
    noPhyMethods <- diffMethods[which(!grepl("P", diffMethods))]
    
    selectedPI <- grep("^P", data$Methods)
    
    PImedian <- median(data$value[selectedPI])
    PIQuantile <- quantile(data$value[selectedPI], c(0.25, 0.75))
    
    for(met in 1:length(noPhyMethods)){
      
      selectedMet <- grep(noPhyMethods[met], data$Methods)
      
      splitData<- data[selectedMet, ]
      
      #transform data
      splitData$value <- PImedian - splitData$value
      
      namePDF <- paste(paste0(save, fileName),"Boxplot", noPhyMethods[met], mmName, sep = "_")
      
      #to have letters as title
      letters <- c("A)", "B)", "C)", "D)")
      name <- letters[rdn]
      
      p <- (ggplot(splitData, aes(x = Methods, y = value)) +
              geom_boxplot(aes(fill = Methods)) +
              theme_classic() +
              scale_fill_manual(values = myPalette, drop = T) +
              xlab("Methods") +
              ylab("Relative average") +
              ylim(-1, 1)+
              theme(legend.position = "none") +
              theme(axis.text=element_text(size = 30),
                    axis.title=element_text(size = 36, face="bold"),
                    plot.title = element_text(size = 45, face="bold")) +
              labs(title = name) +
              #PI
              geom_hline(yintercept = 0, color ="#663333", linetype = "dashed", size = 2) +
              annotate('ribbon', x = c(-Inf, Inf), ymin = PIQuantile[1] - PImedian, 
                       ymax = PIQuantile[2] - PImedian, alpha = 0.2, fill = "#663333"))
      
      dat <- ggplot_build(p)$data[[1]]
      
      plt <- p + geom_segment(data = dat, aes(x = xmin, xend = xmax,
                                              y = middle, yend = middle), colour = "blue", size = 2)

      if(returnPlot){
        namePlot <- paste(fileName, "Boxplot", noPhyMethods[met], mmName, sep = "_")
        boxplot.list[[i]] <- local(plt)
        names(boxplot.list)[i] <- namePlot 
        i <- i + 1
      }
      
      
      else{
        namePDF <- paste(paste0(save, fileName),"Boxplot", noPhyMethods[met], mmName, sep = "_")
        pdf(paste0(namePDF, ".pdf"), width = 32, height = 18)
        plot(plt)
        dev.off()
      }
    }
  }
  
  if(returnPlot){
    return(boxplot.list)
  }
}

#' @title Create Final Table
#'
#' @description This function generates a dataframe of 3 columns, random mechanisms, missing rate and imputation methods
#'
#' @usage createFinalTable(namesReplicates, pathReplicates)
#'
#' @param namesReplicates character vector containing the names of the replicates
#' @param pathReplicates character vector containing the path where the replicates are stored

#' @return A dataframe of 3 columns
#' 
createFinalTable <- function(namesReplicates, pathReplicates){
  
  replicate1 <- get(load(paste0(pathReplicates, namesReplicates[1])))$errorData
  randomType <- rep(names(replicate1), each = length(replicate1[[1]][[1]]))
  missingDegree <- str_match(namesReplicates[1], "\\_R\\d+\\_(.*?)\\.RData")[,2]
  missingDegree <- rep(missingDegree, times = length(randomType))
  impute_Approach <- rep(gsub("error_", "",names(replicate1[[1]][[1]])), 
                         times = (length(randomType)) / length(names(replicate1[[1]][[1]])))
  
  #change names in impute_Approach
  impute_Approach <- gsub("_NP", "", impute_Approach)
  impute_Approach <- gsub("_P", "+P", impute_Approach)
  impute_Approach <- gsub("_2-step", "+PI", impute_Approach)
  
  FinalTable <- data.frame(random_Type = randomType, missing_Degree = missingDegree, impute_Approach = impute_Approach)
  
  return(FinalTable)
  
}

#' @title Table mean imputation
#'
#' @description This function generates a dataframe of 10 columns, random mechanisms, missing rate, imputation methods, mean error of the continuous imputation same for the discrete imputation as well as the standard deviation and the confidence interval.
#'
#' @usage meanSDErrorPerTraitsAndTime(namesReplicates, pathReplicates, trait = NULL)
#'
#' @param namesReplicates character vector containing the names of the replicates
#' @param pathReplicates character vector containing the path where the replicates are stored
#' @param trait character mentioning the name of the trait that we want to evaluate the imputation. If NULL, the mean of 
#' all the trait is computed.
#' @return A list composed of the dataframe and of a table of mean of traits. 
#' 
meanSDErrorPerTraitsAndTime <- function(namesReplicates, pathReplicates, trait = NULL){
  
  replicate1 <- get(load(paste0(pathReplicates, namesReplicates[1])))$errorData
  Approaches <- list()
  #create vector to stock mean and sd per methods and trait type
  meanConti <- c()
  sdConti <- c()
  meanDisc <- c()
  sdDisc <- c()
  CIdown <- c()
  CIup <- c()
  

  FinalTable <- createFinalTable(namesReplicates[1], pathReplicates)
  
  #namesPartitions <- c() #4 partitions (MCAR, MAR, MNAR, PhyloMCAR)
  for(rdn in 1:length(replicate1)){
    partitions <- replicate1[[rdn]][[1]]
    matrixPartitionMethods <- list()
      
    data <- mergeErrorData(partitions)
      
    matricesCollection <- c()
      
    for(col in 2:ncol(data)){ #ncol(data) = 10 
      matrixMethods <- matrix(0, nrow = nrow(data), ncol = length(namesReplicates))
      matrixMethods[,1] <- data[, col]
      
      #add the values of replicates
      for(repli in 2:length(namesReplicates)){
        replicate <- get(load(paste0(pathReplicates, namesReplicates[repli])))$errorData
        parti <- replicate[[rdn]][[1]]
        data <- mergeErrorData(parti) 
        matrixMethods[, repli] <- data[, col]
      }
      
      #calculate mean and sd for each traits
      matrixMethods <- t(matrixMethods)
      colnames(matrixMethods) <- data$trait
      
      
      if(!is.null(trait)){
        matrixMethods <- matrixMethods[ ,trait, drop = FALSE]
      }
      
      matrixMethods <- gather(as.data.frame(matrixMethods), factor_key = TRUE)
      
      names(matrixMethods)[which(names(matrixMethods) == "key")] <- "traits"
      
      matrixMethods <- matrixMethods %>% group_by(traits) %>%
        summarise(mean = mean(value), sd = sd(value), CIdown = attr(confint(value), "range")[1], 
                  CIup = attr(confint(value), "range")[2])
      
      #add overall mean and sd in the last row of dataframe
      #continous traits
      contiTraits <- matrixMethods[which(str_detect(matrixMethods$traits, "F")), ]
      #discrete traits
      discTraits <- matrixMethods[which(str_detect(matrixMethods$traits, "I")), ]
      
      #save CI borns
      CId <- matrixMethods$CIdown
      CIu <- matrixMethods$CIup
      
      CIdown <- c(CIdown, mean(CId))
      CIup <- c(CIup, mean(CIu))
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
        
        if(nrow(discTraits) == 1){
          meanDisc <- c(meanDisc, discTraits$mean)
          sdDisc <- c(sdDisc, discTraits$sd) 
        }
        
        else{
          meanDisc <- c(meanDisc, mean(discTraits$mean))
          sdDisc <- c(sdDisc, sd(discTraits$mean))
        }
        
      }
      
      if(nrow(discTraits) == 0){
        #print("noDisc")
        meanDisc <- c(meanDisc, NA)
        sdDisc <- c(sdDisc, NA)
      }
      
      matricesCollection <- c(matricesCollection, list(matrixMethods))
      names(matricesCollection)[col-1] <- names(data)[col]
    }

    Approaches[[rdn]] <- matricesCollection
    names(Approaches)[rdn] <- names(replicate1[rdn])
  }
  
  FinalTable$MeanConti <- meanConti
  FinalTable$SDConti <- sdConti
  FinalTable$MeanDisc <- meanDisc
  FinalTable$SDDisc <- sdDisc
  FinalTable$CIdown <- CIdown
  FinalTable$CIup <- CIup
  
  
  #add columns with group of imputation method (class)
  class <- rep(NA, nrow(FinalTable))
  class[grep("^c|^R", FinalTable$impute_Approach)] <- "Phylogenetic imputation"
  class[which(str_detect(FinalTable$impute_Approach, pattern = "^G"))] <- "Deep learning"
  class[which(str_detect(FinalTable$impute_Approach, pattern = "PI"))] <- "2-step"
  class[which(is.na(class))] <- "Machine learning"
  
  FinalTable$Class <- class
  
  return(list(meanTable = FinalTable, meanPerTaits = Approaches))
}