library(OncoSimulR)

"DEFINITION ('ideal' simulation parameters)
WT --> A --> B
r(WT) = 1
r(A) = 1
r(AB) > 1
mu(A) != mu(AB)
if sqrt(mu(AB)) >> 1/N, then tunneling is the leading process
"

"Here the simulations that reproduces the events described in the paper."
r <- data.frame(Genotype = c("WT", "A", "B", "A, B"),
                Fitness = c("f_ / (f_ + 15 * f_A + 200 * f_A_B)",
                            "15 * f_A / (f_ + 15 * f_A + 200 * f_A_B)",
                            "0",
                            "200 * f_A_B / (f_ + 15 * f_A + 200 * f_A_B)"),
                stringsAsFactors = FALSE)

afe <- allFitnessEffects(genotFitness = r,
                         frequencyDependentFitness = TRUE,
                         spPopSizes = c(97, 1, 1, 1), 
                         frequencyType ='rel')
evalAllGenotypes(fitnessEffects = afe)

#### Three types of cells. Genuine two-step process ####
set.seed(5)
osi <- oncoSimulIndiv(afe,
                      model = "McFL",
                      onlyCancer = FALSE,
                      finalTime = 100000,
                      verbosity = 0,
                      mu = c("A" = 1e-4, "B" = 1e-3),   # "A" = 2 * 
                      initSize = 100,
                      keepPhylog = FALSE,
                      seed = NULL,
                      detectionProb = NA,
                      detectionSize = NA,
                      errorHitMaxTries = FALSE,
                      errorHitWallTime = FALSE,
                      sampleEvery=1)
plot(osi, show = "genotypes", type = "line")

#### Three types of cells. Tunneling ####
set.seed(5)
osi <- oncoSimulIndiv(afe,
                      model = "McFL",
                      onlyCancer = FALSE,
                      finalTime = 100000,
                      verbosity = 0,
                      mu = c("A" = 1e-4, "B" = 1e-1),
                      initSize = 100,
                      keepPhylog = FALSE,
                      seed = NULL,
                      detectionProb = NA,
                      detectionSize = NA,
                      errorHitMaxTries = FALSE,
                      errorHitWallTime = FALSE,
                      sampleEvery=1)
plot(osi, show = "genotypes", type = "line")
plot(osi, show = "genotypes", type = "line", xlim = c(48000, 49000))

"Everything can be calculated more precisely for this model, if we know all about McFL."
N <- 100
r <- 15   # This is not the ectual reproductive rate in our model!
rho <- r ^ (N - 1) * (1 - r) / (1 - r ^ N)   # probability that a single A cell 
#  will be fixed in a popilation of N-1 "WT" cells
t <- 100000   # Here also might be something!
mu = 1e-4
rpob_allA_at_t <- 1 - exp(-N * mu * rho * t)

#### A tumor suppressor gene and chromosomal instability ####
# I maybe try to model. 
# Though, I have a doubt that I will be able to repeat their results.

