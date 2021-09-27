library(MPSEM)
library(phytools)
library(Rphylopars)

setwd("/Users/dsilvestro/Documents/Projects/Ongoing/Catalina_trait_imputation/")
tree_file = "Phylogeny_NoCon_Henderson_FUGEutd.tre"
trait_file = "PalmTraits_UpToDateTAX_comb2.txt"



tree <- ape::read.nexus(file = tree_file)
trait_data <- read.csv(file = trait_file, header = TRUE, sep = "\t") 
colnames(trait_data)[1] <- "species"

### MATCH TRAITS AND TREE
#Using package MPSEM
#check tip labels match species in dataset - NAs appear if there is no match in data table 
spmatch <- match(tree[[1]]$tip.label, trait_data[,1L])
#drop NAs
if (class(tree) == "phylo"){
        red_tree <- drop.tip(tree, tree$tip.label[is.na(spmatch)])
}else if (class(tree) == "multiPhylo"){
        red_tree <- phytools::drop.tip.multiPhylo(tree, tree[[1]]$tip.label[is.na(spmatch)])
}
# drop data
red_data_tmp <- trait_data[trait_data[,1L] %in% red_tree[[1]]$tip.label, ]

#reorder (if needed) and verify match
spmatch <- match(red_tree[[1]]$tip.label, red_data_tmp[,1L])
red_data <- red_data_tmp[spmatch,] 
all(red_tree[[1]]$tip.label==red_data[,1L])


dtypes = sapply(red_data, class)
drop_categorical = c()
for (i in 2:length(dtypes)){
        if (dtypes[i] == "character"){
                l <- length(unique(red_data[,i]))
                if (l > 10){
                        drop_categorical <- c(drop_categorical, i)
                }else{
                        if (l < 6){
                                red_data[,i] = as.integer(as.factor(red_data[,i]))
                                print(c("Setting to categorical: ", colnames(red_data)[i]))
                        }else{
                                red_data[,i] <- as.numeric(as.factor(red_data[,i]))                        
                        }
                        
                }
        }
        else if (min(red_data[,i], na.rm=T) > 0){
                red_data[,i] = log(red_data[,i]) # log-tranform all continuous traits assuming they are measurements
                print(c("Log transforming: ", colnames(red_data)[i]))
        }
        else if (length(unique(red_data[,i])) < 6){
                red_data[,i] = as.integer(red_data[,i])
                print(c("Setting to categorical: ", colnames(red_data)[i]))
        }      
}
# rm cols with too many categorical states
# also removes taxon column: "TAX_genus_species"
red_data <- subset(red_data, select = -drop_categorical)



### IMPUTE USING MISSFOREST
#loading packages for eigenvectors
library(ape)
library(PVR)
library(missForest)

#generating eigenvectors
get_eigenvec <- function(tree, variance_fraction=0.8, numeigen=NULL){
        decomp <- PVRdecomp(tree, type = "newick") #produces object of class 'PVR'
        label.decomp<-as.data.frame(decomp@phylo$tip.label)
        egvec<-as.data.frame(decomp@Eigen$vectors) ##extract eigenvectors
        if (is.null(numeigen)){
                egval<-decomp@Eigen$values #extract eigenvalues
                eigPerc<-egval/(sum(egval)) #calculate % of variance
                eigPercCum <- t(cumsum(eigPerc)) #cumulated variance
                #eigenvectors representing more than X% variance
                numeigen <- sum(eigPercCum<variance_fraction) 
        }
        egOK <- egvec[,1:numeigen] 
        # Change 'numeigen' on above line to a number if you want to specify number of eigenvectors
        eigenTobind<-cbind(label.decomp,egOK) #add names, these are the eigenvectors to merge with trait database

        #Eigenvectors generated in object 'eigenTobind'
        #rename eigenTobind species column so it matches trait dataset species column
        names(eigenTobind)[1] <- "species"
        return(eigenTobind)
        
}

# turn categorical to factos
ind = which(sapply(red_data, typeof) == "integer")
red_data_cp <- red_data
for ( i in ind){
        red_data_cp[,i] = factor(red_data_cp[,i])
}

res = NULL
for (i in 1:length(red_tree)){
        tree = red_tree[[i]]
        eigen <- get_eigenvec(tree, variance_fraction=0.8)
        #merge trait data and eigenvectors
        traitdat <- merge(red_data_cp, eigen, by= "species")
        #drop first column (with species names) 
        #run missForest
        missForest_imputation <- missForest(xmis = traitdat[- c(1)], maxiter = 1, ntree = 10, mtry = 10)
        head(missForest_imputation$ximp) #missForest printed results
        red.append(missForest_imputation)
}










