#############################################
# Case Study : Evasion of apoptosis
# 3 genes, 3 genotypes
############################################



apoptosis <- function(a, b, c ){

# Strategies:
#  1. Cells that produce a paracrine growth factor to prevent apoptosis of neighbouring cells. 
#  2. Cells that produce an autocrine growth factor to prevent apoptosis of themselves. 
#  3. Cells susceptible to paracrine growth factors but incapable of production of factors. 
  

f1 <- paste("f_1*",as.character(1-a+b),"+f_2*",as.character(1-a),"+f_3*",as.character(1-a),sep="")
f2 <- paste("f_1*",as.character(1+b+c),"+f_2*",as.character(1+c),"+f_3*",as.character(1+c),sep="")
f3 <- paste("f_1*",as.character(1+b),"+f_2+f_3",sep="")

print(f1)
print(f2)
print(f3)

r <- data.frame(Genotype = c("WT",
                             "1", "2", "3"),
                Fitness = c("f_",
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

plot(osi, show = "genotypes", type = "line", xlim = c(-1, 170))
plotClonePhylog(osi, N = 0)
return(osi)
}

#a is the cost of producing the paracrine factor
#b the benefit of receiving the paracrine factor
#c the benefit of producing the autocrine factor

a <- 5
b <- 1
c <- 0


apoptosis(a,b,c)
osi
