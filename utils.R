#' @title Conversion factors in one hot encoding 
#' 
#' @description This function converts discrete variables in one hot encoding using the function one_hot from the R package 
#' mltools
#'
#' @usage generateOneHotEncoVariables(NaNData)
#'
#' @param NaNData data frame of one or several factors columns 
#' @return a list containing a data frame in which each discrete variable is converted as one hot encoding and vector of
#' characters.
generateOneHotEncoVariables <- function(NaNData){
  
  nativeColNames <- colnames(NaNData)
  nativeRowNames <- rownames(NaNData)
  
  #in case NaNData is a vector
  if(is.null(dim(NaNData))){
    stop("The column(s) shoud be of class data.frame")
  }
  
  #convert automatically the factor or character columns into one hot encoding
  oneHotdata <- newdata <- mltools::one_hot(as.data.table(NaNData))
  oneHotdata <- as.data.frame(oneHotdata)
  row.names(oneHotdata) <- nativeRowNames
  
  return(list(data = oneHotdata, nativeColNames = nativeColNames))
  
}


#' @title Conversion one hot encoding in factors
#' 
#' @description This function converts data frame of one hot encoding    
#' representing a categorical variable.
#'
#' @usage convertOneHotEncoIntoFactors(oneHotData, nativeColNames)
#'
#' @param oneHotData data frame of one one hot encoding
#' @param nativeColNames columns names of the original column names of the one hot encoding(before the conversion in one hot 
#' encoding). 
#' @return a data.frame in which each variable that were one hot encoding are now factors.
#' 
convertOneHotEncoIntoFactors <- function(oneHotData, nativeColNames){
  
  #isolate one hot encoding variables
  colNames <- colnames(oneHotData)
  
  oneHotColNames <- colNames[grep("\\_", colNames)]
  
  #select the corresponding columns in the nativeColNames vector
  factorColumn <- nativeColNames[nativeColNames %in% unique(str_extract(oneHotColNames, "^.*(?=\\_)"))]
  
  factorDataframe <- data.frame(matrix(NA, nrow = nrow(oneHotData),
                                       ncol = length(nativeColNames)))
  names(factorDataframe) <- nativeColNames
  
  for(c in 1:length(nativeColNames)){
    
    if(nativeColNames[c] %in% unique(str_extract(oneHotColNames, "^.*(?=\\_)"))){
      
      #select data by group of one hot encoding
      subgroupOneHotColNames <- oneHotColNames[grep(nativeColNames[c], oneHotColNames)]
      
      #isolate states of the trait
      states <- str_extract(subgroupOneHotColNames, "\\d+$")
      
      discreteTrait <- c()
      #go through the rows
      for(r in 1:nrow(oneHotData)){
        
        valueState <- states[which(unlist(oneHotData[r, subgroupOneHotColNames]) == 1)]
        
        #in case 2 or 3 one in the row
        if(length(valueState) > 1){
          valueState <- sample(valueState, 1)
        }
        
        #in case only 0 in the row
        else if(length(valueState) == 0){
          valueState <- sample(states, 1)
        }
        
        discreteTrait <- c(discreteTrait, valueState)
      }
      
      factorDataframe[ , nativeColNames[c]] <- as.factor(discreteTrait)
      
    }
    else{
      factorDataframe[ , nativeColNames[c]] <- oneHotData[, nativeColNames[c]]
    }
  }
  
  return(factorDataframe)
  
}

# Generating eigenvectors
#########################
#' @title Eigenvectors calculations
#' 
#' @description This function calculates the eigenvectors of a phylo object 
#'
#' @usage get_eigenvec(tree, variance_fraction = 0.8, numEigen = NULL)
#'
#' @param tree phylogenetic tree of class "phylo"
#' @param variance_fraction variance_fraction minimum variance (%) explained by the eigenvectors
#' @return a data frame in which each column represents an eigenvector
get_eigenvec <- function(tree, variance_fraction){
  
  decomp <- PVRdecomp(tree, type = "newick") #produces object of class 'PVR'
  eigvec <- as.data.frame(decomp@Eigen$vectors) ##extract eigenvectors
  
  egval <- decomp@Eigen$values #extract eigenvalues
  eigPerc <- egval/(sum(egval)) #calculate % of variance
  eigPercCum <- t(cumsum(eigPerc)) #cumulated variance
  
  #eigenvectors representing more than X% variance
  numEigen <- sum(eigPercCum < variance_fraction)
  
  eigOK <- eigvec[,1:numEigen, drop = T] 
  
  if(numEigen == 0){
    print(paste("The variance_fraction should at leat be equal to ", 
                eigPercCum[1]))
    eigOK <- eigvec[1]
  }
  # Change 'numEigen' on above line to a number if you want to specify number of eigenvectors
  #Eigenvectors generated in object 'eigenTobind'
  #rename eigenTobind species column so it matches trait dataset species column
  eigOK <- as.data.frame(eigOK)
  names(eigOK)[1] <- "c1"
  row.names(eigOK) <- decomp@phylo$tip.label
  
  return(eigOK)
}


#' @title Check + conversion columns class. 
#' 
#' @description This function checks if character and integer columns are factors and that the columns composed of float are 
#' numeric. If it is not the case the columns are converted.
#' 
#' @usage checkConvert(missingData)
#' 
#' @param missingData array of data with missing values
#' 
#' @return return the missing dataset with the columns containing integers or characters as factors and columns containing 
#' float as numeric
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

#' @title Calculate error of imputation
#'
#' @description This function calculates the RMSE for imputed continuous data, the absolute 
#' error for imputed ordinal data and the Proportion of falsely classified entries (PFC) for the other 
#' subcategories of discrete data
#'
#' @usage imputationError(imputedData, trueData, missingData, imputationApproachesName, Data)
#'
#' @param imputedData array of imputed data
#' @param trueData array of true data
#' @param missingData array of data with missing values
#' @param Data simulated Data object
#' @return return a data frame with in the first column the trait names and in the second the errors
#' 
imputationError <- function(imputedData, trueData, missingData, imputationApproachesName, Data){
  
  #get the ordinal trait reference
  ordinalTraits <- which(Data$dataframe$class == "ordinal") #give the row in dataframe which correspond to /n in data names
  errors <- c()
  traitNames <- c() 
  for (c in 1:ncol(missingData)){

    #know is NaNs in the columns(trait)
    NaNRowIndex <- which(is.na(missingData[,c]))
    
    if(length(NaNRowIndex != 0)){
      
      traitNames <- c(traitNames, names(trueData[c]))
      #missingValues <- missingData[NaNRowIndex, c]
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
#' @description This function deletes lists of list that are empty 
#'
#' @usage deleteEmptylist(lists)
#'
#' @param lists list containing lists (lists of list)
#' @return return the list of lists without the empty list
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
#' @param ErrorOutputObjects matrix containing the imputations error values. Row: traits, columns: imputation method/variance
#' fraction/missingness
#' @return return list containing the average error for each random approach and each scenario(partition).
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
#' @return return a character vector with all the replicates of the same type of simulated data
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
#' @return return a .RData object with a matrix for each random approaches and partition with the overall mean error
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


#' @title Empirical data and phylogeny in R list
#'
#' @description This function saves an empirical dataset and a phylogeny in a list of the same structure than when the data is 
#' simulated
#'
#' @usage passInList(empData, empTree = NULL)
#'
#' @param empData a table of class data.frame. The categorical variables must be factor and the continuous variables, numerical
#' @param empTree a phylogenetic tree, by default empTree = NULL. 
#' @param save path to save the data
#' @return return a list which mimic the structure of the simulated data list. Contain at least the empirical data as data.frame 
#' and the phylogenetic tree as phylo object
#' 

#n <- passInList(Data$FinalData, Data$TreeList$`0`)

passInList <- function(empData, empTree = NULL, save = NULL){
  
  
  
  empData <- as.data.frame(empData)
  
  colNames <- names(empData)
  
  #change names traits
  nominalNames <- names(Filter(is.factor, empData))
  newNominalNames <- paste0("I", 1:length(nominalNames))
  indexNominal <- which(colNames %in% nominalNames)
  
  names(empData)[indexNominal] <- newNominalNames
  
  indexConti <- grep("^I.", names(empData), invert=TRUE)
  newContiNames <- paste0("F", 1:length(indexConti))
  names(empData)[indexConti] <- newContiNames
  
  TreeList <- NULL
  if(!is.null(empTree)){
    TreeList <- list('0' = empTree)
  }
  
  Data <- list(FinalData = empData, nativeColNames = colNames, TreeList = TreeList)

  if(!is.null(save)){
    save(Data, file = paste0(save, ".RData"))
  }
  return(Data)
}




