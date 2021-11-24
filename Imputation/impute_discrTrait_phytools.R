library(phytools)
library(castor)

# Stuff from SimData
####################
birth = 0.5
death = 0.2
ntaxa = 20

# Simulating phylogenies fails sometimes. Try until we are successful
Extant <- FALSE
while (!Extant) {
  SimTree <- pbtree(b = birth, d = death, n = ntaxa, scale = 1, extant.only = FALSE) #why added extant.only ?
  if(!is.null(SimTree)) {
    Extant <- TRUE
  }
}
# Simulate discrete trait
Nstates <- 3 
Q <- get_random_mk_transition_matrix(Nstates = Nstates, rate_model = "ER", max_rate = 0.5)
tip_states <- simulate_mk_model(tree = SimTree, Q)

# ARE THE TRAITS ORDERED ACCORDING TO TO THEIR ORDER IN THE TREE?
Trait <- tip_states$tip_states

# We need a matrix of prior probabilities for tip states
StateMat <- matrix(0, ncol = Nstates, nrow = length(Trait))
colnames(StateMat) <- 1:ncol(StateMat)
rownames(StateMat) <- SimTree$tip.label # Only correct if traits are ordered as in the tree! 
for (i in 1:length(Trait)) {
  StateMat[i, Trait[i]] <- Trait[i]
}
# Some NAs for the example via equal probability for each state
StateMat[c(2, 5, 10), ] <- 1/ncol(StateMat)

# Imputation
SimmapTrees <- make.simmap(tree = SimTree, x = StateMat, nsim = 1000, model = "ER", pi = "estimated")
SimmapDescribe <- describe.simmap(SimmapTrees, plot = FALSE)
MostLikelyState <- apply(SimmapDescribe$tips, 1, function(x) which.max(x))
# Compare ground truth and imputation
Trait[c(2, 5, 10)]
MostLikelyState[c(2, 5, 10)]
