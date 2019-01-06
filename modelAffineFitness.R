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
}

create_df1 <- function(a, b, c, d, s, t){
  data.frame(Genotype = c("WT", "A", "B", "A, B"),
             Fitness = c(paste(
               as.character(a), " * f_ +", 
               as.character(b), " * f_1 +", 
               as.character(s)),
               paste(
                 as.character(c), " * f_ +",
                 as.character(d), " * f_1 +",
                 as.character(t)),
               "0", "0"),
             stringsAsFactors = FALSE)
}

# α = a − c, β = b − d, and σ = t − s

#### WT and one type of T cells ####
a <- 1.2   # independent
d <- 1.5   # independent
s <- 1   # idependent

# 1. Simulate stable internal equilibrium: α < σ < β
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

# 2. Simulate WT dominates T: σ < α, β
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

# 3. Simulate T dominates WT: σ > α, β
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

# 4. Simulate unstable internal equilibrium: α > σ > β
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




#### Prisoner's Dilemma ####
# In the presence of a constant fitness advantage, cooperation can arise 
#  if the fitness advantage is higher than the cost to pay for cooperation.
#     C    D
# C  b-c  -c
# D   b    0
#  , where C -- cooperation, D -- defection, b -- benefit, c -- cost.
# So, α = β = −c, and for σ < −c (t − s > c) there is a cooperation.
# If α < σ < β -- stable equilibrium, if β < σ < α -- unstable equilibrium.

#### WT and two types of T (three-player games, no interaction b/w T cells) ####
# T1 : α1 =−1 β1 = 1 σ1 = 0 ---
#   exploitation and attraction of normal cells without additional advantage
# T2 : α2 =−2 β2 = 0 σ2 =−1 ---
#   strong exploitation at the cost of a disadvantage (prisoner's dilemma)
# T3 : α3 = 0 β3 = 2 σ3 = 1 ---
#   strong attraction of normal cells with constant fitness advantage
alpha1 <- -1
alpha2 <- -2
alpha3 <- 0
beta1 <- 1
beta2 <- 0
beta3 <- 2
sigma1 <- 0
sigma2 <- -1
sigma3 <- 1
a <- 0   # independent
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

## 1. WT-T1-T2
#     WT  T1  T2
# WT  a   b   b     s
# T1  c1  d1  0     t1
# T2  c2  0   d2    t2
# fixed point: WT/T1/T2 = 0.5/0.5/0
r <- data.frame(Genotype = c("WT", "A", "B"),
                Fitness = c(paste(as.character(a), "* f_ +",
                                  as.character(b), "* f_1 +",
                                  as.character(b), "* f_2 +",
                                  as.character(s)),
                            paste(as.character(c1), "* f_ +",
                                  as.character(d1), "* f_1 +",
                                  as.character(t1)),
                            paste(as.character(c2), "* f_ +",
                                  as.character(d2), "* f_2 +",
                                  as.character(t2))),
                stringsAsFactors = FALSE)

afe <- fast_allFitnessEffects(r)
osi <- fast_oncoSimulIndiv(afe)
plot(osi, show = "genotypes", type = "line")   # OK

# ω1 = 1 − (b −σ1)(b −σ2), 
# ω2 = 1 − (b −σ1) 
# ω3 = 1 − (b −σ2)

w1 = 1 - (b - sigma1) * (b - sigma2)
w2 = 1 - (b - sigma1)
w3 = 1 - (b - sigma2)
w = w1 + w2 + w3
w1 = w1 / w   # 0.5
w2 = w2 / w   # 0.5
w3 = w3 / w   # 0


## 2. WT-T1-T3
#     WT  T1  T3
# WT  a   b   b     s
# T1  c1  d1  0     t1
# T3  c3  0   d3    t3
# fixed point: WT/T1/T3 = 0.25/0.25/0.5
r <- data.frame(Genotype = c("WT", "A", "B"),
                Fitness = c(paste(as.character(a), "* f_ +",
                                  as.character(b), "* f_1 +",
                                  as.character(b), "* f_2 +",
                                  as.character(s)),
                            paste(as.character(c1), "* f_ +",
                                  as.character(d1), "* f_1 +",
                                  as.character(t1)),
                            paste(as.character(c3), "* f_ +",
                                  as.character(d3), "* f_2 +",
                                  as.character(t3))),
                stringsAsFactors = FALSE)

afe <- fast_allFitnessEffects(r)
osi <- fast_oncoSimulIndiv(afe)
plot(osi, show = "genotypes", type = "line")   # OK

## 3. WT-T2-T3
#     WT  T2  T3
# WT  a   b   b     s
# T2  c2  d2  0     t2
# T3  c3  0   d3    t3
# fixed point: WT/T2/T3 = 0.5/0/0.5
r <- data.frame(Genotype = c("WT", "A", "B"),
                Fitness = c(paste(as.character(a), "* f_ +",
                                  as.character(b), "* f_1 +",
                                  as.character(b), "* f_2 +",
                                  as.character(s)),
                            paste(as.character(c2), "* f_ +",
                                  as.character(d2), "* f_1 +",
                                  as.character(t2)),
                            paste(as.character(c3), "* f_ +",
                                  as.character(d3), "* f_2 +",
                                  as.character(t3))),
                stringsAsFactors = FALSE)

afe <- fast_allFitnessEffects(r)
osi <- fast_oncoSimulIndiv(afe)
plot(osi, show = "genotypes", type = "line")   # OK
