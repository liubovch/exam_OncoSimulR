library(OncoSimulR)
library(ggplot2)

fast_allFitnessEffects <- function(fitnessDF){
  allFitnessEffects(genotFitness = fitnessDF,
                    frequencyDependentFitness = TRUE,
                    frequencyType ='rel')
}

fast_oncoSimulIndiv <- function(fitnessEff){
  oncoSimulIndiv(fitnessEff,
                 model = "McFL",
                 onlyCancer = FALSE,
                 initSize = 5000,
                 finalTime = 1000,
                 detectionProb = NA,
                 detectionSize = NA,
                 errorHitMaxTries = FALSE,
                 errorHitWallTime = FALSE)
}

fast_oncoSimulIndiv2 <- function(fitnessEff){
  oncoSimulIndiv(fitnessEff,
                 model = "McFL",
                 onlyCancer = FALSE,
                 initSize = 5000,
                 finalTime = 5000,
                 sampleEvery = 0.25,
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

create_df2 <- function(a, b, c1, d1, c2, d2, s, t1, t2){
  data.frame(Genotype = c("WT", "A", "B"),
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
}

fast_oncoSimulPop <- function(n, fitnessEff){
  oncoSimulPop(n, fitnessEff, 
               model = "McFL", 
               onlyCancer = FALSE, 
               initSize = 5000,
               finalTime = 5000,
               sampleEvery = 0.5,
               detectionProb = NA, 
               detectionSize = NA, 
               errorHitMaxTries = FALSE, 
               errorHitWallTime = FALSE, 
               mc.cores = 4)   # adapt to your hardware
}

osp_create_df <- function(osp){
  # create an empty dataframe
  df <- data.frame(matrix(ncol = 4, nrow = 0))
  x <- c("", "A", "B", "A, B")
  colnames(df) <- x
  
  # a function to place numbers in the same order as in the df
  get_right_order <- function(model){
    order_by_colname <- match(colnames(df), model$GenotypesLabels) + 1
    last_row <- model$pops.by.time[nrow(model$pops.by.time), ]
    rigth_order <- last_row[order_by_colname]
  }
  
  # interate through models and add results to the df
  for (i in osp) {
    df[nrow(df) + 1, ] <- get_right_order(i)
  }
  return(df)
}

osp_nice_plot <- function(osp) {
  df <- osp_create_df(osp)
  df[is.na(df)] <- 0
  df$n <- 1:nrow(df)
  names(df)[1] <- "WT"
  new_df <- reshape(df[c(1:3, 5)], direction = "long", 
                    varying=names(df[c(1:3)]), v.names=c("values"), 
                    timevar = "genotype", times = c("WT", "A", "B"))
  new_df$genotype <- factor(new_df$genotype, levels=c("WT", "A", "B"))
  # boxplot
  ggplot(new_df, 
         aes(x = genotype, y = values)) +
    geom_point(aes(colour = genotype)) +
    labs(y = "Number of cells") + 
    geom_line(aes(group = n), alpha = 0.5, colour = "grey")
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
plot(osi, show = "genotypes", type = "stacked")   # OK

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
plot(osi, show = "genotypes", type = "stacked")   # OK

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
plot(osi, show = "genotypes", type = "stacked")   # OK

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
plot(osi, show = "genotypes", type = "stacked")   # The problem is that we cannot 
#  start our simulation from the point where WT and T are equal in quiantity, BUT:
afe <- allFitnessEffects(genotFitness = r,
                          frequencyDependentFitness = TRUE,
                          frequencyType ='rel',
                          spPopSizes = c(5000, 5000, 0, 0))
evalAllGenotypes(afe)   # here we see that fitnesses are equal, so could have had 
#  equilibrium if we had equal amount of WT and T cells at the start
# Genotype Fitness
# 1       WT    1.85
# 2        A    1.85
# 3        B    0.00
# 4     A, B    0.00


#### Prisoner's Dilemma ####
# In the presence of a constant fitness advantage, cooperation can arise 
#  if the fitness advantage is higher than the cost to pay for cooperation.
#     C    D
# C  b-c  -c
# D   b    0
#  , where C -- cooperation, D -- defection, b -- benefit, c -- cost.
# So, α = β = −c, and for σ < −c (t − s > c) there is a cooperation.
# If α < σ < β -- stable equilibrium and if β < σ < α -- unstable equilibrium.


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
r1 <- create_df2(a, b, c1, d1, c2, d2, s, t1, t2)
afe1 <- fast_allFitnessEffects(r1)
osi1 <- fast_oncoSimulIndiv2(afe1)
plot(osi1, show = "genotypes", type = "stacked")   # OK

osp1 <- fast_oncoSimulPop(n = 100, fitnessEff = afe1)
osp_nice_plot(osp1)

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
r2 <- create_df2(a, b, c1, d1, c3, d3, s, t1, t3)
afe2 <- fast_allFitnessEffects(r2)
osi2 <- fast_oncoSimulIndiv2(afe2)
plot(osi2, show = "genotypes", type = "stacked")   # OK

osp2 <- fast_oncoSimulPop(n = 100, fitnessEff = afe2)
osp_nice_plot(osp2)

## 3. WT-T2-T3
#     WT  T2  T3
# WT  a   b   b     s
# T2  c2  d2  0     t2
# T3  c3  0   d3    t3
# fixed point: WT/T2/T3 = 0.5/0/0.5
r3 <- create_df2(a, b, c2, d2, c3, d3, s, t2, t3)
afe3 <- fast_allFitnessEffects(r3)
osi3 <- fast_oncoSimulIndiv2(afe3)
plot(osi3, show = "genotypes", type = "stacked")   # OK

osp3 <- fast_oncoSimulPop(n = 100, fitnessEff = afe3)
osp_nice_plot(osp3)

