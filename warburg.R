# motility

###############################
# case 1 AG:INV
###############################
library(OncoSimulR)
c <- 0.1 #set cost of moving low
f1 <- paste("1+f_1*",as.character(0.5),"+f_2*",as.character(1),sep="")
f2 <- paste("1+f_1*",as.character(1-c),"+f_2*",as.character(1-c/2),sep="")

r <- data.frame(Genotype = c("WT","AG", "INV"),
                Fitness = c("1",f1,f2),
                stringsAsFactors = FALSE)

afe <- allFitnessEffects(genotFitness = r,
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")
## For reproducibility
set.seed(1)
osi <- oncoSimulIndiv(afe,
                      model = "McFL",
                      onlyCancer = FALSE,
                      finalTime = 5000,
                      verbosity = 0,
                      mu = 1e-5, #any slower and WT will always win
                      initSize = 5000,
                      keepPhylog = TRUE,
                      seed = NULL,
                      detectionProb = NA,
                      detectionSize = NA,
                      errorHitMaxTries = FALSE,
                      errorHitWallTime = FALSE)

osi
plot(osi, show = "genotypes", type = "line")

###########################################
# Case AG:GLY
###########################################

k <- 0.2 #cost of living in acid environment
n <- 0.1 #cost of glycolysis

f1 <- paste("1+f_1*",as.character(0.5),"+f_2*",as.character(0.5-n),sep="")
f2 <- paste("1+f_1*",as.character(0.5+n-k),"+f_2*",as.character(0.5-k),sep="")

r <- data.frame(Genotype = c("WT", "AG", "GLY"),
                Fitness = c("1",f1,  f2),
                stringsAsFactors = FALSE)

afe <- allFitnessEffects(genotFitness = r,
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")

## For reproducibility
set.seed(1)
osi <- oncoSimulIndiv(afe,
                      model = "McFL",
                      onlyCancer = FALSE,
                      finalTime = 5000,
                      verbosity = 0,
                      mu = 1e-5,
                      initSize = 5000,
                      keepPhylog = TRUE,
                      seed = NULL,
                      detectionProb = NA,
                      detectionSize = NA,
                      errorHitMaxTries = FALSE,
                      errorHitWallTime = FALSE)

osi
plot(osi, show = "genotypes", type = "line")

###########################################
# Case AG:INV:GLY
###########################################

k <- 0.2 # the lower, the better for INV cells. Allows GLY cells to compete more with AG cells
c <- 0.6 # Will always win until the cost is greater than 1.95. At 1.95, it matches WT
n <- 0.15 #has no real affect in the INV sells winning as INV will move

f1 <- paste("1+f_1*",as.character(0.5),"+f_2*","+f_3*",as.character(0.5-n),sep="")
f2 <- paste("1+f_1*",as.character(1-c),"+f_2*",as.character(1-c/2),"+f_3*",as.character(1-c),sep="")
f3 <- paste("1+f_1*",as.character(0.5+n-k),"+f_2*",as.character(1-k),"+f_3*",as.character(0.5-k),sep="")

r <- data.frame(Genotype = c("WT",
                             "AG", "INV", "GLY"),
                Fitness = c("1",  f1, f2, f3),
                stringsAsFactors = FALSE)

afe <- allFitnessEffects(genotFitness = r,
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")


## For reproducibility
set.seed(1)
osi <- oncoSimulIndiv(afe,
                      model = "McFL",
                      onlyCancer = FALSE,
                      finalTime = 5000,
                      verbosity = 0,
                      mu = 1e-5,
                      initSize = 5000,
                      keepPhylog = TRUE,
                      seed = NULL,
                      detectionProb = NA,
                      detectionSize = NA,
                      errorHitMaxTries = FALSE,
                      errorHitWallTime = FALSE)

osi
plot(osi, show = "genotypes", type = "line")
