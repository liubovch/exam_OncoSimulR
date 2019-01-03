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


f1 <- "((2*f_A_B + 1*f_B) * 1)"
f2 <- "((1*f_A + 3*f_A_B) * 4)"
f3 <- "(2*f_A * 2)/ 2"


r <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                Fitness = c("f_",
                            f1,
                            f2,
                            f3),
                stringsAsFactors = FALSE)

r <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                Fitness = c("1",
                            "((b*c_3*f_A_B + a*c_2*f_B)*(N-1))/(N - c_1)",
                            "((a*c_1*f_A - d*c_3*f_A_B)*(N-1))/(N - c_2)",
                            "((b*c_1*f_A)*(N-1))/(N - c_3)"),
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
