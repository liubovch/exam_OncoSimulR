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
                 mu = 1e-3,
                 initSize = 10000,
                 keepPhylog = TRUE,
                 seed = NULL,
                 detectionProb = NA,
                 detectionSize = NA,
                 errorHitMaxTries = FALSE,
                 errorHitWallTime = FALSE)
}



create_df1 <- function(a, b, c, d, x, y, z){
  data.frame(Genotype = c("WT", "A", "B", "A, B", "C"),
             Fitness = c("0.5", paste("1 +",
               as.character(b*z), " * f_C +", 
               as.character(a*y), " * f_B"),
               paste("1 +",
                 as.character(a*x), " * f_A -",
                 as.character(d*z), " * f_C"),
               "0", 
               paste("1+",
                     as.character(b*x), " * f_A")),
             stringsAsFactors = FALSE)
}

### SYSTEM DEFINITION ###

'''
Cancer cells and stromal cells cooperate by exchanging diffusible factors that 
sustain tumor growth. In the case of multiple myeloma, there are three types
of cells: Malignant plasma cells (MM), osteoblast (OB) and osteoclasts (OC).

Essentially OC and OB cells are in equilibrium in the absence of MM cells, while
MM and OC cells have a stimulatory effect on each other, and MM cells inthibit OB
and OB cells hav little or not effect of MM cells:

Summary of effects of diffusible factors:

     GOC          GOB             GMM
OC  neutral       equilibrium     stimulation
OB  equilibrium   neutral         inhibition
MM  stimulation   neutral         neutral
'''

### Variables definition

# A = OC cells
# B = OB cells
# C = MM cells

# x = Contribution of phenotype A
# y = Contribution of phenotype B
# z = Contribution of phenotype C


# 1. Simulation of y < x < z
'''
In multiple myeloma the contribution of the different cells are different.
In this scenario, the benefit of difusible factos that are secreted by MM cells
is greater than the benefit that OC cells can obtain through the diffusible factos
produced by OB. There exists a polymorphic stable point between MM and OC.
'''

a <- 1
b <- 2.5
d <- 0.3
x <- 1.4
y <- 1.2
z <- 1

r <- create_df1(a, b, c, d, x, y, z)
afe <- fast_allFitnessEffects(r)
osi_1 <- fast_oncoSimulIndiv(afe)
plot(osi_1, show = "genotypes", type = "line", xlim= c(12, 20))



# 2. Simulation of y = x = z

#When a < b an equilibrium OC - OB is created
a <- 1
b <- 0.5
d <- 0.33
x <- 1
y <- 1
z <- 1

r <- create_df1(a, b, c, d, x, y, z)
afe <- fast_allFitnessEffects(r)
osi_2 <- fast_oncoSimulIndiv(afe)
plot(osi_2, show = "genotypes", type = "line", xlim= c(12, 60))

#When b + d > a an equilibrium OC-MM is created

a <- 1
b <- 2
d <- 0
x <- 1
y <- 1
z <- 1

r <- create_df1(a, b, c, d, x, y, z)
afe <- fast_allFitnessEffects(r)
osi_3 <- fast_oncoSimulIndiv(afe)
plot(osi_3, show = "genotypes", type = "line", xlim= c(12, 23))



