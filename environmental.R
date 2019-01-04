library(OncoSimulR)
library(graph)
library(igraph)

a <- 1
b <- 2.5
-d <- -0.3
c_1 <-1.4
c_2 <- 1.2
c_3 <- 1
N <- 10


f1 <- "f_1*2 + f_2*3 + f_3*4"
f2 <- "f_1*2 + f_2*5 + f_3*1"
f3 <- "f_1*4 +f_2 + f_3"


r <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                Fitness = c("f_",
                            f1,
                            f2,
                            f3),
                stringsAsFactors = FALSE)

r <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                Fitness = c("1",
                            “f_1*(z-e-f+g)+f_2*(z-e)+f_3*(z-e+g)”,
                            “f_1*(z-h)+f_2*(z-h)+f_3*(z-h)”,
                            “f_1*(z-f)+f_2+f_3”),
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
                      mu = 1e-6,
                      initSize = 5000, 
                      keepPhylog = TRUE,
                      seed = NULL,
                      detectionProb = NA,
                      detectionSize = NA,
                      errorHitMaxTries = FALSE, 
                      errorHitWallTime = FALSE)

osi

plot(osi, show = "genotypes", type = "line")

plotClonePhylog(osi, N=0)