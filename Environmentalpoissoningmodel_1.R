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


create_df1 <- function(e, f, g, h){
  data.frame(Genotype = c("WT", "A", "B", "A, B"),
             Fitness = c(paste(as.character(z-f), " * f_A +", 
                               as.character(z)," f_B +",
                               as.character(z)," f_"),
                         paste(as.character(z-e-f+g), " * f_A +",
                               as.character(z-e), " * f_B +",
                               as.character(z-e+g), " * f_"),
                         paste(as.character(z-h), " * f_A +",
                               as.character(z-h), " * f_B +",
                               as.character(z-h), " * f_"),
                         "0"),
             stringsAsFactors = FALSE)
}


### SYSTEM DEFINITION ###

'''
Tumour cells might boost their own replicative potential at the expense of 
other tumour cells by evolving the capability of producing cytotoxic substances.
Strategies(or genotypes)
WT: cells producing neither cytotoxins nor resistance.
A: Cells producing cytotoxic substances against other cells, 
B: cells producing resistance to external cytotoxic substances

z= Baseline fintess
e= Cost of producing cytotoxin
f= Disadvantage of being affected by cytotoxin
g= Advantage conferred after having subjected another cell to the cytotoxin
h= Cost of resistance to cytotoxin

'''
# Model 1: Expected frequencies: WT= 0.604, A=0.396, B=0
z <- 1
e <- 0.1
f <- 0.4
g <- 0.01
h <- 0.25

r <- create_df1(e, f, g, h)
afe <- fast_allFitnessEffects(r)
osi_1 <- fast_oncoSimulIndiv(afe)
plot(osi_1, show = "genotypes", type = "line", xlim= c(12, 20))


# Model 2: Expected frequencies: WT= 1, A=0, B=0
z <- 1
e <- 0.3
f <- 0.4
g <- 0.1
h <- 0.25

r <- create_df1(e, f, g, h)
afe <- fast_allFitnessEffects(r)
osi_2 <- fast_oncoSimulIndiv(afe)
plot(osi_2, show = "genotypes", type = "line", xlim= c(12, 1000))


# Model 3: Expected frequencies: WT= 1, A=0, B=0
z <- 1
e <- 0.12
f <- 0.2
g <- 0.24
h <- 0.05

r <- create_df1(e, f, g, h)
afe <- fast_allFitnessEffects(r)
osi_3 <- fast_oncoSimulIndiv(afe)
plot(osi_3, show = "genotypes", type = "line", xlim= c(12, 20))
