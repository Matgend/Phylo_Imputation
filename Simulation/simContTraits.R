library(mvMORPH)
library(phytools)
library(Matrix)


#' @title Simulate variance-covariance matrix for morphological evolution
#'
#' @description This function generates a variance-covariance matrix for
#' morphological evolution, defining the rate at which traits evolve and their
#' evolutionary correlation.
#'
#' @usage simSigma(Ntraits, Cor = NULL, Sigma2 = NULL)
#'
#' @param Ntraits number of traits
#' @param Cor correlation between traits. Default between -1 and 1.
#' Optional, can be fixed to be equal between all traits by giving one value or
#' Ntraits*(Ntraits-1)/2
#' @param Sigma2 Brownian motion rate. Default between 1e-4 and 0.5.
#' Optional, can be fixed to be equal for all traits or one value per trait
#'
#' @return matrix Ntrait x Ntrait for simulating trait evolution
#' 
simSigma <- function (Ntraits, Cor = NULL, Sigma2 = NULL) {
  if (!is.null(Sigma2)) {
    if ( length(Sigma2) != Ntraits && length(Sigma2) != 1) {
      stop("Sigma2 should be of length 1 or Ntraits")
    }
    else {
      if (length(Sigma2) == 1) {
        Sigma2 <- rep(Sigma2, Ntraits)
      }
    }
  }
  else {
    Sigma2 <- runif(Ntraits, min = 1e-4, max = 0.5)
  }
  Sigmas <- matrix(Sigma2, nrow = 1)
  if (Ntraits > 1) {
    Cov <- matrix(1, ncol = Ntraits, nrow = Ntraits)
    Q <- Ntraits*(Ntraits-1)/2
    if (!is.null(Cor)) {
      if (length(Cor) != Q && length(Cor) != 1) {
        stop("Correlation among traits should be of length 1 or Ntraits*(Ntraits-1)/2")
      }
      SimCov <- Cor
    }
    else {
      SimCov <- runif(Q, min = -1, max = 1) # Trait correlation
    }
    Cov[lower.tri(Cov, diag = FALSE)] <- SimCov
    Cov <- t(Cov)
    Cov[lower.tri(Cov, diag = FALSE)] <- SimCov
    Sigmas <- diag(Sigma2)  %*% Cov  %*% diag(Sigma2) # Correlation to covariance
    # Force variance-covariance to be poisitve definite
    Tol <- 1e-6
    Ev <- eigen(Sigmas, symmetric = TRUE)$values
    if ( !all( Ev >= -Tol * abs(Ev[1L]) ) ) {
      Sigmas <- as.matrix(nearPD(Sigmas)$mat)
    }
  }
  return(Sigmas)
}


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
Sigmas <- simSigma(Ntraits)


# Simulate continuous traits
############################
SimTraits <- mvSIM(tree = SimTree, model = "BM1",
                   param = list(ntraits = Ntraits,
                                theta = rep(0, Ntraits),
                                sigma = Sigmas))


# Fit Brownian motion model
###########################
FitBM <- mvBM(tree = SimTree,
              data = SimTraits,
              model = "BM1",
              method = "rpf",
              diagnostic = FALSE,
              echo = FALSE)


# Compare ground truth and estimated variance-covariance matrix
###############################################################
FitBM$sigma
Sigma
cov2cor(FitBM$sigma) # Evolutionary trait correlation
cov2cor(Sigma)