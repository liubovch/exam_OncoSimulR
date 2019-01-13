library(OncoSimulR)
library(ggplot2)

# used functions
fast_allFitnessEffects <- function(fitnessDF){
  allFitnessEffects(genotFitness = fitnessDF,
                    frequencyDependentFitness = TRUE,
                    frequencyType ='rel')
}

fast_oncoSimulIndiv_t1000 <- function(fitnessEff){
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

fast_oncoSimulIndiv_t5000 <- function(fitnessEff){
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

create_f_coeff1 <- function(a, d, s, alpha, beta, sigma){
  list(a = a, d = d, s = s, c = a - alpha, b = d + beta, t = s + sigma)
}

create_f_df1 <- function(coeff){
  data.frame(Genotype = c("WT", "A", "B", "A, B"),
             Fitness = c(paste(
               as.character(coeff$a), " * f_ +", 
               as.character(coeff$b), " * f_1 +", 
               as.character(coeff$s)),
               paste(
                 as.character(coeff$c), " * f_ +",
                 as.character(coeff$d), " * f_1 +",
                 as.character(coeff$t)),
               "0", "0"),
             stringsAsFactors = FALSE)
}

create_f_coeff2 <- function(a, b, s, strategy1, strategy2){
  list(a = a, b = b, s = s, 
       c1 = a - strategy1$alpha, 
       c2 = a - strategy2$alpha,
       d1 = b - strategy1$beta,
       d2 = b - strategy2$beta,
       t1 = strategy1$sigma + s,
       t2 = strategy2$sigma + s)
}

create_f_df2 <- function(coeff){
  data.frame(Genotype = c("WT", "A", "B"),
             Fitness = c(paste(as.character(coeff$a), "* f_ +",
                               as.character(coeff$b), "* f_1 +",
                               as.character(coeff$b), "* f_2 +",
                               as.character(coeff$s)),
                         paste(as.character(coeff$c1), "* f_ +",
                               as.character(coeff$d1), "* f_1 +",
                               as.character(coeff$t1)),
                         paste(as.character(coeff$c2), "* f_ +",
                               as.character(coeff$d2), "* f_2 +",
                               as.character(coeff$t2))),
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
  
  # reshape df
  new_df <- reshape(df[c(1:3, 5)], direction = "long", 
                    varying=names(df[c(1:3)]), v.names=c("values"), 
                    timevar = "genotype", times = c("WT", "A", "B"))
  new_df$genotype <- factor(new_df$genotype, levels=c("WT", "A", "B"))
  
  # nice plot
  ggplot(new_df, 
         aes(x = genotype, y = values)) +
    geom_point(aes(colour = genotype)) +
    labs(y = "Number of cells") + 
    geom_line(aes(group = n), alpha = 0.5, colour = "grey")
}


# α = a − c, β = b − d, and σ = t − s

#### WT and one type of T cells ####
"Just setting some values here (they are not from the paper)"
a <- 1.2   # independent
d <- 1.5   # independent
s <- 1   # idependent

# 1. Simulate stable internal equilibrium: α < σ < β
"Just setting some values here to satisfy the above condition"
f_coeff <- create_f_coeff1(a, d, s, alpha = -1, beta = 1, sigma = 0)

r <- create_f_df1(f_coeff)
afe <- fast_allFitnessEffects(r)
osi <- fast_oncoSimulIndiv_t1000(afe)
plot(osi, show = "genotypes", type = "stacked")

# 2. Simulate WT dominates T: σ < α, β
f_coeff <- create_f_coeff1(a, d, s, alpha = 0, beta = 1, sigma = -1)

r <- create_f_df1(f_coeff)
afe <- fast_allFitnessEffects(r)
osi <- fast_oncoSimulIndiv_t1000(afe)
plot(osi, show = "genotypes", type = "stacked")

# 3. Simulate T dominates WT: σ > α, β
f_coeff <- create_f_coeff1(a, d, s, alpha = 0, beta = -1, sigma = 1)

r <- create_f_df1(f_coeff)
afe <- fast_allFitnessEffects(r)
osi <- fast_oncoSimulIndiv_t1000(afe)
plot(osi, show = "genotypes", type = "stacked")

# 4. Simulate unstable internal equilibrium: α > σ > β
f_coeff <- create_f_coeff1(a, d, s, alpha = 1, beta = -1, sigma = 0)

r <- create_f_df1(f_coeff)
afe <- fast_allFitnessEffects(r)
osi <- fast_oncoSimulIndiv_t1000(afe)
plot(osi, show = "genotypes", type = "stacked")   # The problem is that we cannot 
#  start our simulation from the point where WT and T are equal in quiantity, BUT:
afe <- allFitnessEffects(genotFitness = r,
                          frequencyDependentFitness = TRUE,
                          frequencyType ='rel',
                          spPopSizes = c(5000, 5000, 0, 0))
evalAllGenotypes(afe)   # here we see that fitnesses are equal, so could have had 
#  the equilibrium
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


#### WT and two types of T (three-player games, no interaction b/w T cells) ####
# T1 : α1 =−1 β1 = 1 σ1 = 0 ---
#   exploitation and attraction of normal cells without additional advantage
# T2 : α2 =−2 β2 = 0 σ2 =−1 ---
#   strong exploitation at the cost of a disadvantage (prisoner's dilemma)
# T3 : α3 = 0 β3 = 2 σ3 = 1 ---
#   strong attraction of normal cells with constant fitness advantage
a <- 0   # independent
b <- 0   # independent
s <- 0.5   # independent
strategy1 = list(alpha = -1, beta = 1, sigma = 0)
strategy2 = list(alpha = -2, beta = 0, sigma = -1)
strategy3 = list(alpha = 0, beta = 2, sigma = 1)

## 1. WT-T1-T2
# fixed point: WT/T1/T2 = 0.5/0.5/0
coeff1 <- create_f_coeff2(a, b, s, strategy1, strategy2)
r1 <- create_f_df2(coeff1)
afe1 <- fast_allFitnessEffects(r1)
set.seed(1)   # fix the behaviour since others are possible
osi1 <- fast_oncoSimulIndiv_t5000(afe1)
plot(osi1, show = "genotypes", type = "stacked")

osp1 <- fast_oncoSimulPop(n = 100, fitnessEff = afe1)
osp_nice_plot(osp1)

## 2. WT-T1-T3
# fixed point: WT/T1/T3 = 0.25/0.25/0.5
coeff2 <- create_f_coeff2(a, b, s, strategy1, strategy3)
r2 <- create_f_df2(coeff2)
afe2 <- fast_allFitnessEffects(r2)
set.seed(1)   # fix the behaviour since others are possible
osi2 <- fast_oncoSimulIndiv_t5000(afe2)
plot(osi2, show = "genotypes", type = "stacked")

osp2 <- fast_oncoSimulPop(n = 100, fitnessEff = afe2)
osp_nice_plot(osp2)

## 3. WT-T2-T3
# fixed point: WT/T2/T3 = 0.5/0/0.5
coeff3 <- create_f_coeff2(a, b, s, strategy2, strategy3)
r3 <- create_f_df2(coeff3)
afe3 <- fast_allFitnessEffects(r3)
set.seed(1)   # fix the behaviour since others are possible
osi3 <- fast_oncoSimulIndiv_t5000(afe3)
plot(osi3, show = "genotypes", type = "stacked")

osp3 <- fast_oncoSimulPop(n = 100, fitnessEff = afe3)
osp_nice_plot(osp3)

