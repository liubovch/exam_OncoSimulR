library(OncoSimulR)

fast_allFitnessEffects <- function(fitnessDF){
  allFitnessEffects(genotFitness = fitnessDF,
                    frequencyDependentFitness = TRUE,
                    frequencyType ='rel')
}

fast_oncoSimulIndiv <- function(fitnessEff){
  oncoSimulIndiv(fitnessEff,
                 model = "McFL",
                 onlyCancer = FALSE,
                 finalTime = 1000,
                 verbosity = 0,
                 mu = 1e-6,
                 initSize = 5000,
                 keepPhylog = TRUE,
                 seed = NULL,
                 detectionProb = NA,
                 detectionSize = NA,
                 errorHitMaxTries = FALSE,
                 errorHitWallTime = FALSE)
}

create_df1 <- function(a, b, c, d, s, t){
  data.frame(Genotype = c("WT", "A", "B", "A, B"),
             Fitness = c(paste(
               as.character(a), " * f_ + ", 
               as.character(b), " * f_1 + ", 
               as.character(s)),
               paste(
                 as.character(c), " * f_ + ",
                 as.character(d), " * f_1 + ",
                 as.character(t)),
               "0", "0"),
             stringsAsFactors = FALSE)
}

# α = a − c, β = b − d, and σ = t − s

#### WT and one type of T cells ####
a <- 1.2   # independent
d <- 1.5   # independent
s <- 1   # idependent

# Simulate stable internal equilibrium: α < σ < β
alpha1 <- -1
beta1 <- 1
sigma1 <- 0
c <- a - alpha1
b <- d + beta1
t <- s + sigma1

r <- create_df1(a, b, c, d, s, t)
afe <- fast_allFitnessEffects(r)
osi <- fast_oncoSimulIndiv(afe)
plot(osi, show = "genotypes", type = "line")   # OK

# Simulate WT dominates T: σ < α, β
alpha1 <- 0
beta1 <- 1
sigma1 <- -1
c <- a - alpha1
b <- d + beta1
t <- s + sigma1

r <- create_df1(a, b, c, d, s, t)
afe <- fast_allFitnessEffects(r)
osi <- fast_oncoSimulIndiv(afe)
plot(osi, show = "genotypes", type = "line")   # OK

# Simulate T dominates WT: σ > α, β
alpha1 <- 0
beta1 <- -1
sigma1 <- 1
c <- a - alpha1
b <- d + beta1
t <- s + sigma1

r <- create_df1(a, b, c, d, s, t)
afe <- fast_allFitnessEffects(r)
osi <- fast_oncoSimulIndiv(afe)
plot(osi, show = "genotypes", type = "line")   # OK

# Simulate unstable internal equilibrium: α > σ > β
alpha1 <- 1
beta1 <- -1
sigma1 <- 0
c <- a - alpha1
b <- d + beta1
t <- s + sigma1

r <- create_df1(a, b, c, d, s, t)
afe <- fast_allFitnessEffects(r)
osi <- fast_oncoSimulIndiv(afe)
plot(osi, show = "genotypes", type = "line")   # The problem is that we cannot 
#  start our simulation from the point where WT and T are equal in quiantity, BUT:
afe2 <- allFitnessEffects(genotFitness = r,
                          frequencyDependentFitness = TRUE,
                          frequencyType ='rel',
                          spPopSizes = c(5000, 5000, 0, 0))
evalAllGenotypes(afe2)   # here we see that fitnesses are equal, so could have had 
#  equilibrium if we had equal amount of WT and T cells at the start




#### Prisoner's Dilemma ??? #### 
#### WT and two types of T ####
# T1 : α1 =−1 β1 = 1 σ1 = 0
# T2 : α2 =−2 β2 = 0 σ2 =−1 (prisoner's dilemma)
# T3 : α3 = 0 β3 = 2 σ3 = 1
# three-player games
alpha1 <- -1
alpha2 <- -2
alpha3 <- 0
beta1 <- 1
beta2 <- 0
beta3 <- 2
sigma1 <- 0
sigma2 <- -1
sigma3 <- 1
a <- 1   # independent
b <- 0   # independent
s <- 0.5   # independent
c1 <- a - alpha1
c2 <- a - alpha2
c3 <- a - alpha3
d1 <- b - beta1
d2 <- b - beta2
d3 <- b - beta3
t1 <- sigma1 + s
t2 <- sigma2 + s
t3 <- sigma3 + s

r <- data.frame(Genotype = c("WT", "A", "B"),
                Fitness = c("f_ + 0.5",
                            "2 * f_ - 0.5",
                            "f_ - 2 * f_2 + 1.5"),
                stringsAsFactors = FALSE)

afe <- allFitnessEffects(genotFitness = r,
                         frequencyDependentFitness = TRUE,
                         frequencyType ='rel')
osi <- oncoSimulIndiv(afe,
                      model = "McFL",
                      onlyCancer = FALSE,
                      finalTime = 5000,
                      verbosity = 0,
                      mu = 1e-6,
                      initSize = 5000,
                      keepPhylog = TRUE,
                      seed = NULL,
                      detectionProb = NA,
                      detectionSize = NA,
                      errorHitMaxTries = FALSE,
                      errorHitWallTime = FALSE)
plot(osi, show = "genotypes", type = "line")

# ω1 = 1 − (b −σ1)(b −σ2), 
# ω2 = 1 − (b −σ1) 
# ω3 = 1 − (b −σ2)

w1 = 1 - (b - sigma1) * (b - sigma2)
w2 = 1 - (b - sigma1)
w3 = 1 - (b - sigma2)
w = w1 + w2 + w3
w1 = w1 / w
w2 = w2 / w
w3 = w3 / w
