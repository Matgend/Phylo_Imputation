if(!require(ggpubr)){
  install.packages("ggpubr")
}

if(!require(miscset)){
  install.packages("miscset")
}

if(!require(stringr)){
  install.packages("stringr")
}

if(!require(tidyr)){
  install.packages("tidyverse")
}

library(ggpubr)
library(miscset)
library(stringr)
library(tidyverse)

source("utils.R")

#create directories for plots and tables
dir.create("../Simulations/Overall_Results/", showWarnings = FALSE)
dir.create("../Simulations/Overall_Results/Plots", showWarnings = FALSE)
dir.create("../Simulations/Overall_Results/Table", showWarnings = FALSE)
dir_list <- list.dirs("../Simulations", recursive = F)
dir_list <- dir_list[grep("/0", dir_list)]

#color palette for the plots
############################

myPalette <- c("GAIN" = "#0099FF" , "GAIN+P" = "#0099FF", "GAIN+PI" = "#0099FF", "KNN" = "#66CC00", "KNN+P" = "#66CC00", 
               "KNN+PI" = "#66CC00", "MICE" = "#339900", "MICE+P" = "#339900", "MICE+PI" = "#339900",
               "MissForest" = "#336600", "MissForest+P" = "#336600", "MissForest+PI" = "#336600",
               "PI" = "#993300", "HV.ML" = "darkorchid1", "HV.ML+P" = "darkorchid1", "HV.ML+PI" = "darkorchid1")


#dir_list
#d <- 1

for (d in 1:(length(dir_list))){
  path <- dir_list[d]
  
  #load replicates
  nameDir <- str_split(path, "/", simplify = TRUE)[,3]
  
  pathReplicates <- paste0(path,"/Results/Replicates/")
  #pathSavePlots <- paste0("../Simulations/Overall_Results/Plots", nameDir, "/")
  #pathSaveTable <- paste0("../Simulations/Overall_Results/Table", nameDir, "/")
  pathSaveMeanTable <- "../Simulations/Overall_Results/Table/"

  pathCsv <- "../csv"
  #dir.create(pathSavePlots, showWarnings = FALSE)
  #dir.create(pathSaveTable, showWarnings = FALSE)
  
  replicateFiles <- list.files(path = pathReplicates, pattern = "*.RData")
  csvFiles <- list.files(path = pathCsv, pattern = "*.csv")
  csvName <- unique(gsub("\\.csv", "", csvFiles))
  
  
  boxplot.finalList <- list()
  i <- 1
  for(file in 1:length(csvFiles)){
    replicatesInFolder <- loadReplicates(file, replicateFiles, pathCsv)
    replicatesInFolder <- replicatesInFolder[str_detect(replicatesInFolder, paste0("Results", csvName[file], "_"))] # to add when have 100 replicates!!!! [1:100]
    missingRate <- unique(str_extract(replicatesInFolder, "0\\.\\d+"))

    #Merge all the replicates 
    processedData <- processData(replicatesInFolder, pathReplicates, trait = "I1.0/1")

    #Generate boxplot of comparison
    fileName <- paste0(csvName[file], "_", missingRate)
    boxplots <- boxPlotByMethods(processedData, fileName, myPalette, returnPlot = TRUE)
    boxplot.finalList[[i]] <- boxplots
    names(boxplot.finalList)[i] <- csvName[file]
    i <- i + 1
    
    #Generate TXT with mean + sd of each imputation methods
    fileName <- gsub("\\_..*", "", replicatesInFolder)
    meanData <- meanSDErrorPerTraitsAndTime(replicatesInFolder, pathReplicates,  trait = "I1.0/1")
    missingRate <- gsub("\\.", "", missingRate)

    #create a copy which is saved in another directory for plotting purpose
    fileNameTxt <- unique(paste0(pathSaveMeanTable, fileName,"_Table_", missingRate, ".txt"))
    write.table(meanData$meanTable, fileNameTxt, sep = "\t", row.names = FALSE)
    print("meanData, done")
  }
  
  #Generate Boxplots
  ##################

  pathSaveBoxPlots <- "../Simulations/Overall_Results/Plots/"
  impMethods <- c("KNN", "MICE", "MissForest", "GAIN", "HV.ML")
  
  for(m in 1:length(impMethods)){
    methodsIsolated <- impMethods[m]
    for(l in 1:length(boxplot.finalList)){
      indexPlots <- grep(methodsIsolated, names(boxplot.finalList[[l]]))
      plots <- ggarrange(plotlist = boxplot.finalList[[l]][indexPlots], ncol = 2,
                         nrow = 2, legend = "none")
      
      fileName <- paste0(pathSaveBoxPlots, paste(names(boxplot.finalList[l]), methodsIsolated, "MergedBoxplot",
                                              missingRate, sep = "_"), ".pdf")
      pdf(fileName, width = 32, height = 18)
      plot(plots)
      dev.off()
    }
    graphics.off()
    
  }
}


#Generate Final Tables
######################
#pathSaveMeanTable

txtFiles <- list.files(pathSaveMeanTable)

missingRate <- unique(str_match(txtFiles, "Table(.*?)\\.txt")[,2])

for(mR in 1:length(missingRate)){
  
  list_files <- txtFiles[grep(missingRate[mR], txtFiles)]
  #print(list_files)
  
  randomApproaches <- unique(read.table(paste(pathSaveMeanTable, list_files[1], sep = "/"), 
                                        sep = "\t", header = TRUE)$random_Type)
  
  for(rdn in 1:length(randomApproaches)){
    
    for(f in 1:(length(list_files))){
      
      table <- read.table(paste(pathSaveMeanTable, list_files[f], sep = "/"), sep = "\t", header = TRUE)
      
      subdata <- table[which(table$random_Type == randomApproaches[rdn]), ]
      
      #create final table
      if (f == 1){
        FinalData <- as.data.frame(matrix(NA, nrow = length(unique(table$impute_Approach)), ncol = length(list_files)))
        colNames <- gsub("Results|\\_.*", "", list_files)
        colnames(FinalData) <- colNames
        
        methods <- unique(table$impute_Approach)
      }
        
      if(sum(is.na(subdata$MeanConti)) != nrow(subdata)){
        
        meanP <- subdata$MeanConti
        sdP <- subdata$SDConti
        
      }
      
      else if(sum(is.na(subdata$MeanDisc)) != nrow(subdata)){
        meanP <- subdata$MeanDisc
        sdP <- subdata$SDDisc
        
      }
      
      for(r in 1:nrow(FinalData)){
        
        #accuracy
        meanFill <- round(1 - meanP[which(subdata$impute_Approach == methods[r])], 3)
        sdFill <- round(sdP[which(subdata$impute_Approach == methods[r])], 3)
        numberFill <- paste0(meanFill, "(", sdFill, ")")
        
        if(numberFill == "()" | numberFill == "NA(NA)"){
          numberFill <- NA
        }
        
        else{
          FinalData[r, f] <- numberFill
        }
        
      }
    }
    
    #add name to rows
    FinalData <- cbind(Methods = methods, FinalData)
    
    nameTxt <- paste(paste0(pathSaveMeanTable, "FinalTableAccuracy"), randomApproaches[rdn], 
                     gsub("_", "", missingRate[mR]), sep = "_")
    write.table(FinalData, paste0(nameTxt, ".txt"), sep = "\t", row.names=FALSE)
    
    #nameCSV <- paste(paste0(pathSave, "FinalTableAccuracy"), randomApproaches[rdn], missingRate, sep = "_")
    #print(FinalData)
    #write.table(FinalData, paste0(nameCSV, ".csv"), sep = ";", row.names = FALSE)
  }
}





#Generate Barplots
##################
###barplots with 4x3 (missing mecha x phylo signal)

files <- list.files(path = pathSaveMeanTable, pattern = "Final.*.txt")
missingRate <- unique(str_extract(files, "\\d+[^.txt]+"))

for(mR in 1:length(missingRate)){

  filesMR <- files[grep(paste0("_",missingRate[mR]), files)]
  model <- c("MK", "TM")
  
  #print(list_files)
  
  mergedTable <- c()
  #generate table
  for(mM in 1:length(randomApproaches)){
    
    table <- read.table(paste0(pathSaveMeanTable, "/", filesMR[grep(randomApproaches[mM], filesMR)]), header = T)
    
    #add missing mechanims
    table$m.mechanism <- randomApproaches[mM]
    
    if(mR == 1 & mM == 1){
      mergedTable <- table
      
    }
    else{
      mergedTable <- rbind(mergedTable, table)
    }
  }
  
  mergedTable$m.mechanism[which(mergedTable$m.mechanism == "PhyloNaN")] <- "phyloNa"
  #tun the table
  
  #when data of size 100x13
  #change col names
  names(mergedTable) <- c("Methods", expression(paste("MK, ", kappa, " = ", 1, ", ", lambda, " = ", 1)),
                          expression(paste("MK, ", kappa, " = ", 0, ", ", lambda, " = ", 1)),
                          expression(paste("MK, ", kappa, " = ", 1, ", ", lambda, " ~ ", 0)),
                          expression(paste("TM, ", kappa, " = ", 1, ", ", lambda, " = ", 1)),
                          expression(paste("TM, ", kappa, " = ", 0, ", ", lambda, " = ", 1)),
                          expression(paste("TM, ", kappa, " = ", 1, ", ", lambda, " ~ ", 0)),
                          "m.mechanisms")
  
  colNameChange <- names(mergedTable)[2:(length(mergedTable) - 1)]
  
  mergedTable <- mergedTable %>%
    pivot_longer(cols = all_of(colNameChange), names_to = "scenario", values_to = "acc.sd")
  
  #convert value as numeric and create a sd column
  mean <- str_extract(mergedTable$acc.sd, "[^(]+")
  sd <- str_match(mergedTable$acc.sd, "\\((.*)\\)")[,2]
  
  mergedTable$acc.sd <- as.numeric(mean)
  names(mergedTable)[4] <- "mean"
  mergedTable$sd <- as.numeric(sd)
  
  for(m in 1:length(model)){
    
    for(phy in 1:3){
      
      if(phy == 1){
        rows <- which(!grepl("\\+", mergedTable$Methods) | mergedTable$Methods == "PI")
        name <- "noTree"
      }
      else if(phy == 2){
        rows <- grep("P$", mergedTable$Methods)
        name <- "Tree"
      }
      else{
        rows <- which(grepl("I$", mergedTable$Methods) | mergedTable$Methods == "PI")
        name <- "PI"
      }
      #subset plot
      subdata <- mergedTable[rows, ]
      subdata <- subdata[which(grepl(model[m], subdata$scenario)), ]
      
      subdata$m.mechanisms <- factor(subdata$m.mechanisms, levels=c("MCAR", "MAR", "MNAR", "phyloNa"))
      subdata$scenario <- factor(subdata$scenario, levels = unique(subdata$scenario)[c(1, 2, 3)])
      
      subdata$Methods <- factor(subdata$Methods, levels = unique(subdata$Methods)[c(1:4, 6, 5)])
      
      minBoard <- subdata$mean - subdata$sd
      maxBoard <- subdata$mean + subdata$sd
      
      #check if values smaller than 0
      minBoard[which(minBoard < 0)] <- 0
      maxBoard[which(maxBoard > 1)] <- 1
      
      plotSave <- ggplot(subdata) + 
        geom_bar(aes(x = Methods, y = mean, fill = Methods), stat = "identity") +
        scale_fill_manual(values = myPalette)+
        theme(legend.position = "none")+
        geom_errorbar(aes(x = Methods, ymin =  minBoard, ymax = maxBoard), width = 0.6, colour="black", 
                      alpha = 0.9, size = 0.8) +
        ylab("Accuracy [%]") +
        theme(axis.text.x = element_text(size = 32, angle = 45, hjust = 1),
              axis.text.y = element_text(size = 32, hjust = 1),
              axis.title = element_text(size = 30,face="bold")) +
        ylim(0, 1) +
        geom_hline(yintercept = 1/3, linetype="dashed", color = "red", size = 1.2)+
        facet_grid(vars(m.mechanisms), vars(scenario), labeller = label_parsed)+
        theme(strip.text.x = element_text(size = 40, face="bold"),
              strip.text.y = element_text(size = 36, face="bold"))
      
      pathSavePlots <- "../Simulations/Overall_Results/Plots/"
      fileName <- paste0(pathSavePlots, "BarplotsAllMM_allPS_", missingRate[mR], "_", model[m], "_", name , ".pdf")
      print(fileName)
      pdf(fileName, width = 32, height = 18)
      plot(plotSave)
      dev.off()
    }
    
  }
  
}
graphics.off()

#Lines
######

#Table for NaN plot
dataframelist <- list()
for(mR in 1:length(missingRate)){
  list_files <- list.files(pathSaveMeanTable)
  indexFiles <- grep(paste0("Table_", missingRate[mR]), list_files)
  list_files <- list_files[indexFiles]
  
  for(f in 1:(length(list_files))){ #enlever -2
    table <- read.table(paste(pathSaveMeanTable, list_files[f], sep = "/"), sep = "\t", header = TRUE)
    
    #add variance_fractions column
    table$variance_fractions <- "0"
    
    #change variance_fraction
    table$variance_fractions[grep("+P$", table$impute_Approach)] <- "P"
    table$variance_fractions[grep("+PI", table$impute_Approach)] <- "2-step"
    table$variance_fractions[which(table$variance_fractions == "0")] <- "NP"
    
    #change name of approaches
    table$impute_Approach <- gsub(paste0(c("\\+P", "\\+PI"), collapse = "|"), "", table$impute_Approach)
    
    if(mR == 1){
      #dataframelist[[f]] <- table[,c("random_Type", "partition", "missing_Degree", "variance_fractions", 
      #                               "impute_Approach", "MeanConti", "MeanDisc")]
      dataframelist[[f]] <- table
      nameDataframe <- gsub("_0.*$", "", list_files)
    }
    else{
      # table <- table[,c("random_Type", "partition", "missing_Degree","variance_fractions", 
      #                   "impute_Approach", "MeanConti", "MeanDisc")]
      dataframelist[[f]] <- rbind(dataframelist[[f]], table)
    }
    
    dataframelist[[f]]$impute_Approach[grep(c("^c|R"), dataframelist[[f]]$impute_Approach)] <- "PI"
    
  }
}
names(dataframelist) <- nameDataframe


colorMM <- c("MCAR" = "darkgray", "MAR" = "orange", "MNAR" = "blue", "PhyloNaN" = "purple")
NaNplot.list <- list()
i <- 1
for (d in 1:length(dataframelist)){
  
  data <- dataframelist[[d]]
  methodsName <- unique(data$impute_Approach)
  strategies <- unique(data$variance_fractions)
  
  for(m in 1:length(methodsName)){
    
    for(s in 1:length(strategies)){
      
      subdata <-subset(data, impute_Approach == methodsName[m] & variance_fractions == strategies[s])
      
      figTitle <- paste(methodsName[m], strategies[s], sep = "_")
      save <- pathSavePlots
      figTitle <- str_replace_all(figTitle, "/", "_")
      fileName <- paste(names(dataframelist)[d], figTitle, "lines", sep = "_")
      namePDF <- paste0(save, fileName)

      plt <- (ggplot(subdata, aes(x = missing_Degree, y = MeanDisc, col = random_Type)) +
                geom_point(size = 8)+
                scale_color_manual(values = colorMM, drop = T)+
                geom_line(size = 2)+
                geom_ribbon(aes(y = MeanDisc, ymin = CIdown, ymax = CIup, fill = random_Type), 
                            alpha = 0.1, show.legend = FALSE) +
                xlab("Missing rate [%]") + 
                ylab("Average misclassification [%]") +
                labs(colour = 'Missing mechanisms') +
                ylim(0, 1) +
                theme_minimal()+
                theme(axis.text=element_text(size=25),
                      axis.title=element_text(size=25,face="bold"),
                      legend.key.size = unit(2, "cm"), #change legend key size
                      legend.key.height = unit(2, "cm"), #change legend key height
                      legend.key.width = unit(2, "cm"), #change legend key width
                      legend.title = element_text(size=25), #change legend title font size
                      legend.text = element_text(size=20)))#change legend text font size))

      
      NaNplot.list[[i]] <- local(plt)
      names(NaNplot.list)[i] <- fileName
      i <- i + 1
      
    }
  }
}

#Automatise NaNPlot figure from the list NaNplot.list
#####################################################
pathSavePlots <- "../Simulations/Overall_Results/Plots/"

#names(NaNplot.list)
methodsIsolated <- unique(str_match(names(NaNplot.list), "Table_\\s*(.*?)\\s*_lines")[,2])

models <- c("ARD", "BM")

for(met in 1:length(methodsIsolated)){
  
  indexPlots <- grep(methodsIsolated[met], names(NaNplot.list))
  
  for(m in 1:length(models)){
    
    index.list <- indexPlots[grep(models[m], names(NaNplot.list[indexPlots]))]
    
    plots <- ggarrange(plotlist = NaNplot.list[index.list], ncol = 3,
                       nrow = 1, legend = "bottom", common.legend = TRUE)
    
    fileName <- paste0(pathSavePlots, paste(models[m], methodsIsolated[met], "MergedNaN", 
                                            sep = "_"), ".pdf")
    pdf(fileName, width = 32, height = 18)
    plot(plots)
    dev.off()
  }
  
  graphics.off()
}
