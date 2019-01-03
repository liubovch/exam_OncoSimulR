#############################################
# Case Study : Evasion of apoptosis
# 2 genes, 3 genotypes
############################################

a <- 1
b <- 2
c <- 3

f1 <- paste("f_1*",as.character(1-a+b),"+f_2*",as.character(1-a),"+f_3*",as.character(1-a),sep="")
f2 <- paste("f_1*",as.character(1+b+c),"+f_2*",as.character(1+c),"+f_3*",as.character(1+c),sep="")
f3 <- paste("f_1*",as.character(1+b),"+f_2+f_3",sep="")



r <- data.frame(Genotype = c("WT",
                             "1", "2", "3"),
                Fitness = c(f_,
                            f1,
                            f2,
                            f3),
                stringsAsFactors = FALSE)
afe <- allFitnessEffects(genotFitness = r,
                         frequencyDependentFitness = TRUE,frequencyType ='rel')

osi <- oncoSimulIndiv(afe,
                      model = "McFL",
                      onlyCancer = FALSE,
                      finalTime = 10000,
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
plot(osi, show = "genotypes", type = "line", xlim = c(-1, 170))
plotClonePhylog(osi, N = 0)