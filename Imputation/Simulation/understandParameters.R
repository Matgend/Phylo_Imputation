evolvTrait <- function (Trait, TimeAdd, Sigma, 
                        Alpha = 0, Theta = 0) {
  TraitAdd <- c(Trait, rep(NA_real_, TimeAdd))
  Lt <- length(Trait) + 1             
  for(i in Lt:length(TraitAdd)){
    TraitAdd[i] <- TraitAdd[i-1] + 
      Alpha * (Theta - TraitAdd[i-1]) + 
      rnorm(1, mean = 0, sd = Sigma)
  }
  return(TraitAdd)
}

resample <- function (x, ...) x[sample.int(length(x), ...)]

simDivDisp <- function (LambdaS1, LambdaS2, ExS1, ExS2,
                        Q, Theta0, Sig2S1, Sig2S2, AlphaS1, AlphaS2, ThetaS1, ThetaS2) {
  Extant <- 1
  Traits <- list()
  Traits[[1]] <- Theta0
  SpEx <- list()
  SpEx[[1]] <- c(1, NA)
  Time <- 1
  State <- 1
  while (Time < 100 & Extant > 0) {
    EverExisted <- length(SpEx)
    Index <- 1:length(SpEx)
    ExtantS1 <- unlist( lapply(Index, function(x) State[x] == 1 & any(is.na(SpEx[[x]]))) )
    ExtantS2 <- unlist( lapply(Index, function(x) State[x] == 2 & any(is.na(SpEx[[x]]))) )
    # Time of events: Speciation, habitat transition, or extinction
    ###############################################################
    NumExtantS1 <- sum(ExtantS1)
    NumExtantS2 <- sum(ExtantS2)
    LambdaS1Tmp <- NumExtantS1 * LambdaS1
    LambdaS2Tmp <- NumExtantS2 * LambdaS2
    ExS1Tmp <- NumExtantS1 * ExS1
    ExS2Tmp <- NumExtantS2 * ExS2
    QS1Tmp <- NumExtantS1 * Q[1]
    QS2Tmp <- NumExtantS2 * Q[2]
    Dt <- round(rexp(1, LambdaS1Tmp + LambdaS2Tmp + QS1Tmp + QS2Tmp + ExS1Tmp + ExS2Tmp), 0)
    #print(Dt)
    if(Dt == 0){
      Dt <- 1
    }
    TimeNew <- Time + Dt
    # What happens
    ##################################################
    # Speciation state 1, 2, transition 1->2, 2->1, Extinction state 1, 2
    Event <- sample(1:6, 1, 
                    prob = c(LambdaS1Tmp, LambdaS2Tmp, QS1Tmp, QS2Tmp, ExS1Tmp, ExS2Tmp))
    if(TimeNew >= 100){
      TimeNew <- 100
      SpEx[ExtantS1] <- lapply(SpEx[ExtantS1], function(x) c(x[1], TimeNew))
      SpEx[ExtantS2] <- lapply(SpEx[ExtantS2], function(x) c(x[1], TimeNew))
    }
    else{
      if(Event %in% c(1, 3, 5)){ # Affects state 1
        W <- resample(which(ExtantS1), 1) # Which lineage?
        SpEx[[W]][2] <- TimeNew
        if(Event == 1){ # Speciation S1
          # Add the new species
          SpEx[[EverExisted + 1]] <- c(TimeNew, NA)
          SpEx[[EverExisted + 2]] <- c(TimeNew, NA)
        }
        if(Event == 3){ # Transition S1 -> S2
          SpEx[[EverExisted + 1]] <- c(TimeNew, NA)
          State <- c(State, 2)
        }
      }
      if(Event %in% c(2, 4, 6)){ # Affects state 2
        W <- resample(which(ExtantS2), 1) # Which lineage?
        SpEx[[W]][2] <- TimeNew
        if(Event == 2){ # Speciation S2
          # Add the new species
          SpEx[[EverExisted + 1]] <- c(TimeNew, NA)
          SpEx[[EverExisted + 2]] <- c(TimeNew, NA)
        }
        if(Event == 4){ # Transition S2 -> S1
          SpEx[[EverExisted + 1]] <- c(TimeNew, NA)
          State <- c(State, 1)
        }
      }
    }
    # Evolve all extant lineages until TimeNew
    ##########################################
    TimeAdd <- TimeNew - Time
    if(TimeAdd > 0){
      Traits[ExtantS1] <- lapply(Traits[ExtantS1], function(x) evolvTrait(Trait = x, 
                                                                          TimeAdd = TimeAdd, 
                                                                          Sigma = Sig2S1,
                                                                          Alpha = AlphaS1, 
                                                                          Theta = ThetaS1))
      Traits[ExtantS2] <- lapply(Traits[ExtantS2], function(x) evolvTrait(Trait = x, 
                                                                          TimeAdd = TimeAdd, 
                                                                          Sigma = Sig2S2,
                                                                          Alpha = AlphaS2, 
                                                                          Theta = ThetaS2))
    }
    # Provide starting value for traits for the descendents
    #######################################################
    if(Event %in% 1:4){
      TraitAncestor <- Traits[[W]]
      TraitAncestor <- TraitAncestor[length(TraitAncestor)]
      Traits[[EverExisted + 1]] <- TraitAncestor
      if(Event %in% 1:2){
        # Add descendents
        Traits[[EverExisted + 2]] <- TraitAncestor
        if(Event == 1){
          State <- c(State, 1, 1)
        }
        else{
          State <- c(State, 2, 2)
        }
      }
    }
    Time <- TimeNew
    Index <- 1:length(SpEx)
    ExtantS1 <- unlist( lapply(Index, function(x) 
      State[x] == 1 & any(is.na(SpEx[[x]]))) )
    ExtantS2 <- unlist( lapply(Index, function(x) 
      State[x] == 2 & any(is.na(SpEx[[x]]))) )
    Extant <- sum(ExtantS1) + sum(ExtantS2)
    Res <- list()
    Res[[1]] <- SpEx
    Res[[2]] <- Traits
    Res[[3]] <- State
  }
  return(Res)
}


# Play around with different values and see what is happening
#############################################################
SimRes <- simDivDisp(LambdaS1 = 0.3, # Speciation rate state 1 
                     LambdaS2 = 0.3, # Speciation rate state 2 
                     ExS1 = 0.1, # Extinction rate state 1 
                     ExS2 = 0.1, # Extinction rate state 2
                     Q = c(0.05, 0.05), # Transition rate state 1->2 2->1
                     Theta0 = 0, # Initial morphological value state 1 
                     # Initial morphological value of state 2 depends on state 1 ancestor!
                     Sig2S1 = 0.1, # Rate of morphological evolution state 1 
                     Sig2S2 = 0.1, # Rate of morphological evolution state 2 
                     AlphaS1 = 0.1, # Strength of stabilizing selection state 1 
                     AlphaS2 = 0.1, # Strength of stabilizing selection state 2 
                     ThetaS1 = 1, # Morphological optima state 1 
                     ThetaS2 = -1) # Morphological optima state 2
# Plot for each lineage the trait values over time
palette(c("dodgerblue", "orange"))
plot(1, 1, type = "n",
     ylab = "Morphology", xlab = "Time",
     xlim = c(1, 100), 
     ylim = range(unlist(SimRes[[2]])))
for(i in 1:length(SimRes[[1]])){
  X <- seq(SimRes[[1]][[i]][1], SimRes[[1]][[i]][2], 1)
  lines(X, SimRes[[2]][[i]], col = SimRes[[3]][i])
}