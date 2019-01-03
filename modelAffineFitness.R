library(OncoSimulR)
# T1 : α1 =−1 β1 = 1 σ1 = 0
# T2 : α2 =−2 β2 = 0 σ2 =−1
# T3 : α3 = 0 β3 = 2 σ3 = 1
# α = a − c, β = b − d, and σ = t − s
a = 1
s = 1
b = 1.5
sigma1 = 0.6
sigma2 = 0.7
alpha1 = -1 + sigma1
alpha2 = -1 + sigma2
c1 = a - alpha1
c2 = a - alpha2
beta1 = 1 + sigma1
beta2 = 1 + sigma2
d1 = b - beta1
d2 = b - beta2
t1 = sigma1 + s
t2 = sigma2 + s

f1 <- paste("f_ *", as.character(a), 
            "+ f_1 *", as.character(b), 
            "+ f_2 *", as.character(b), 
            "+", as.character(s))

f2 <- paste("f_ *", as.character(c1), 
            "+ f_1 *", as.character(d1), 
            "+ f_2 *", as.character(0), 
            "+", as.character(t1))

f3 <- paste("f_ *", as.character(c2), 
            "+ f_1 *", as.character(0), 
            "+ f_2 *", as.character(d2), 
            "+", as.character(t2))

r <- data.frame(Genotype = c("WT", "A", "B"),
                Fitness = c(f1, f2, f3),
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
plotClonePhylog(osi, N = 0)

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


