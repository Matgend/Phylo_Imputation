library(mvMORPH)
library(phytools)
library(Matrix)

# Simulate phylogeny
####################
Birth <- 0.4
Death <- 0.1
Ntaxa <- 30
Ntraits <- 7
# Simulating phylogenies fails sometimes. Try until we are successful
Extant <- FALSE
while (!Extant) {
  SimTree <- pbtree(b = Birth, d = Death, n = Ntaxa, scale = 1, extant.only = TRUE)
  if (!is.null(SimTree)) {
    Extant <- TRUE
  }
}

# Get matrix of evolutionary rate (diagonal) and covariance (off-diagonals)
###########################################################################
# function simSigma from script Simulation/simContTraits.R
Sigmas <- simSigma(Ntraits)


# Simulate continuous traits
############################
SimTraits <- mvSIM(tree = SimTree, model = "BM1",
                   param = list(ntraits = Ntraits,
                                theta = rep(0, Ntraits),
                                sigma = Sigmas))


# Missing values completely at random
#####################################
MissFrac <- 0.1 # Fraction of missing values
SimTraitsNA <- SimTraits
Obs <- nrow(SimTraits) * ncol(SimTraits)
IdxNA <- sample(Obs, round(Obs * MissFrac))
SimTraitsNA[IdxNA] <- NA 
  

# Fit Brownian motion model
###########################
FitBM <- mvBM(tree = SimTree,
              data = SimTraitsNA,
              model = "BM1",
              method = "rpf",
              diagnostic = FALSE,
              echo = FALSE)

# Impute traits
###############
ImputTraits <- estim(tree = SimTree,
                     data = SimTraitsNA,
                     object = FitBM,
                     asr = FALSE)


# Compare ground truth and imputed values
#########################################
SimTraits[IdxNA]
ImputTraits$estimates[IdxNA]
