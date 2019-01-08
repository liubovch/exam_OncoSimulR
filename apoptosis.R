library("OncoSimulR")
library("ggplot2")
library("doParallel")

set.seed(1)

####################################
# Functions
####################################

# Creates an object fitnessEffects from a data frame with the genotypes and 
# fitness specification
# If is possible to specify also a vector with some initial population sizes 

wraperFitnessEffects <- function(r,spPS=NULL){
  if (!is.null(spPS)) { 
    afe <- allFitnessEffects(genotFitness = r,
                             frequencyDependentFitness = TRUE,frequencyType ='rel',spPopSizes=spPS) 
  } else{
    afe <- allFitnessEffects(genotFitness = r,
                             frequencyDependentFitness = TRUE,frequencyType ='rel')
  }
  return(afe)
}

# Receives a fitnessEffects object
# Generates a a data frame with two columns, the Genotype and the Fitness
# plots the fitness landscape

genotypes <- function(afe){
  # addwt = TRUE to take into account WT 
  eval <- evalAllGenotypes(afe, addwt = TRUE)
  plotFitnessLandscape(eval)
  eval
}

# Run one oncoSimul trajectory
# receives:
#            afe -  fitnessEffects object 
#            isize - initial size of the population
#            mutrate - mutation rate
#            fintime - maximum number ot time units to run
# Returns an ocosimul object

simulOne <- function(afe,isize,mutrate,fintime){
  
  
  osi <- oncoSimulIndiv(afe,
                        model = "McFL",
                        onlyCancer = FALSE,
                        finalTime = fintime,
                        verbosity = 0,
                        mu = mutrate,
                        initSize = isize,
                        keepPhylog = TRUE,
                        seed = NULL,
                        detectionProb = NA,
                        detectionSize = NA,
                        errorHitMaxTries = FALSE,
                        errorHitWallTime = FALSE)
  
  # Plot genotypes sizes along time
  plot(osi, show = "genotypes", type = "line")
  # Plot a parent-child relationship of the clones
  plotClonePhylog(osi, N = 0)
  
  return(osi)
}

# Run a  number of oncoSImul trajectories defined by the parameter iter
# receives:
#            afe -  fitnessEffects object 
#            isize - initial size of the population
#            mutrate - mutation rate
#            fintime - maximum number ot time units to run
# Returns an ocosimulpop object

simulPop <- function(afe,iter,isize,mutrate,fintime){
  
  
  osp <- oncoSimulPop(iter,afe,
                      model = "McFL",
                      onlyCancer = FALSE,
                      finalTime = fintime,
                      verbosity = 0,
                      mu = mutrate,
                      initSize = isize,
                      keepPhylog = FALSE,
                      seed = NULL,
                      detectionProb = NA,
                      detectionSize = NA,
                      errorHitMaxTries = FALSE,
                      errorHitWallTime = FALSE)
  return(osp)
  
}

# Draw a boxplot with the results of the OncoSimulPop simulation
# Returns a dataframe with the final population for each genotype
# for all simulations

popStats <- function(osp,x){
  
  # a function to place numbers in the same order as in the df
  get_right_order <- function(model,df){
    order_by_colname <- match(colnames(df), model$GenotypesLabels) + 1
    last_row <- model$pops.by.time[nrow(model$pops.by.time), ]
    rigth_order <- last_row[order_by_colname]
  }
  
  # create an empty dataframe
  df <- data.frame(matrix(ncol = length(x), nrow = 0))
  colnames(df) <- x
  # interate through models and add results to the df
  for (i in osp) {
    df[nrow(df) + 1, ] <- get_right_order(i,df)
  }
  
  # boxplot
  ggplot(stack(df), aes(x = ind, y = values)) +
    geom_boxplot()
  
  return(df)
}





############################################
# Case Study : Evasion of apoptosis
############################################

# Creates a dataframe with the genotypes and its
# fitness specification according to the apoptosis
# game model, with the right format to be
# used to create a fitnessEffect object with the 
# allFitnessEffects function

createdf_apoptosis <- function( a, b, c ){

# Strategies:
#  a. Cells that produce a paracrine growth factor to prevent apoptosis of neighbouring cells. 
#  b. Cells that produce an autocrine growth factor to prevent apoptosis of themselves. 
#  c. Cells susceptible to paracrine growth factors but incapable of production of factors. 
  

fa <- paste("1+f_A*(",as.character(1-a+b),
            ")+f_B*(",as.character(1-a),
            ")+f_C*(",as.character(1-a),")",
            sep="")
fb <- paste("1+f_A*(",as.character(1+b+c),
            ")+f_B*(",as.character(1+c),
            ")+f_C*(",as.character(1+c),
            ")",sep="")
fc <- paste("1+f_A*(",as.character(1+b),")+f_B+f_C",sep="")

r <- data.frame(Genotype = c("WT", "A", "B","C"),
                Fitness = c("1", fa, fb, fc),
                stringsAsFactors = FALSE)

return(r)

}

# Creates a dataframe with the genotypes and its
# fitness specification according to the apoptosis
# game model, with the right format to be
# used to create a fitnessEffect object with the 
# allFitnessEffects function

createdf_apoptosis_strategiesAC <- function(a,b){
  
  # Strategies:
  #  a. (1) Cells that produce a paracrine growth factor to prevent apoptosis of neighbouring cells.
  #  c. (2) Cells susceptible to paracrine growth factors but incapable of production of factors. 
  
  
  fa <- paste("1+f_1*(",as.character(1-a+b),
              ")+f_2*(",as.character(1-a),")",
              sep="")
  fc <- paste("1+f_1*(",as.character(1+b),
              ")+f_2",
              sep="")
  
  r <- data.frame(Genotype = c("WT", "1", "2"),
                  Fitness = c("1", fa, fc),
                  stringsAsFactors = FALSE)
  
  return(r)
  
}


# Creates a dataframe with the genotypes and its
# fitness specification according to the apoptosis
# game model, with the right format to be
# used to create a fitnessEffect object with the 
# allFitnessEffects function

createdf_apoptosis_strategiesBC <- function(c){
  
  # Strategies:
  #  b.(1) Cells that produce an autocrine growth factor to prevent apoptosis of themselves. 
  #  c.(2) Cells susceptible to paracrine growth factors but incapable of production of factors. 
  
  
  fb <- paste("1+f_1*(", as.character(1+c),
              ")+f_2*(", as.character(1+c),")",
              sep="")
  fc <- "1+f_1+f_2"
  
  r <- data.frame(Genotype = c("WT", "1", "2"),
                  Fitness = c("1", fb, fc),
                  stringsAsFactors = FALSE)
  
  return(r)
  
}


#a is the cost of producing the paracrine factor
#b the benefit of receiving the paracrine factor
#c the benefit of producing the autocrine factor

a <- 1
b <- 3
c <- 2

isize <- 5000
mutrate <- 1e-5
fintime <- 1000


dfapop_1 <- createdf_apoptosis(a, b, c)
afeapop_1 <- wraperFitnessEffects(dfapop_1)
osi_1 <- simulOne(afeapop_1, isize, mutrate, fintime)
# extended simulations with oncoSimulPop
x <- c("", "A", "B", "C")
osp_1 <- simulPop(afeapop_1, 100, isize, mutrate, fintime)
pstats_1 <- popStats(osp_1, x)

a <- 1
b <- 2
c <- 3

dfapop_2 <- createdf_apoptosis(a, b, c)
afeapop_2 <- wraperFitnessEffects(dfapop_2)
osi_2 <- simulOne(afeapop_2, isize, mutrate, fintime)
# extended simulations with oncoSimulPop
osp_2 <- simulPop(afeapop_2, 100, isize, mutrate, fintime)
pstats_2 <- popStats(osp_2, x)

a <- 0
b <- 2
c <- 3

dfapop_3 <- createdf_apoptosis(a, b, c)
afeapop_3 <- wraperFitnessEffects(dfapop_3)
osi_3 <- simulOne(afeapop_3, isize, mutrate, fintime)
# extended simulations with oncoSimulPop
osp_3 <- simulPop(afeapop_3, 100, isize, mutrate, fintime)
pstats_3 <- popStats(osp_3, x)


a <- 0
b <- 3
c <- 2

dfapop_4 <- createdf_apoptosis(a, b, c)
afeapop_4 <- wraperFitnessEffects(dfapop_4)
osi_4 <- simulOne(afeapop_4, isize, mutrate, fintime)
# extended simulations with oncoSimulPop
osp_4 <- simulPop(afeapop_4, 100, isize, mutrate, fintime)
pstats_4 <- popStats(osp_4, x)


y <- c("", "1", "2")

# Just strategies B(1) an C(2)
c <- 2
dfapop_5 <- createdf_apoptosis_strategiesBC(c)
afeapop_5 <- wraperFitnessEffects(dfapop_5)
osi_5 <- simulOne(afeapop_5, isize, mutrate, fintime)
# extended simulations with oncoSimulPop
osp_5 <- simulPop(afeapop_5, 100, isize, mutrate, fintime)
pstats_5 <- popStats(osp_5, y)


a <- 1
b <- 3
#Just strategies A(1) anc C(2)
dfapop_6 <- createdf_apoptosis_strategiesAC(a, b)
afeapop_6 <- wraperFitnessEffects(dfapop_6)
osi_6 <- simulOne(afeapop_6, isize, mutrate, fintime)
# extended simulations with oncoSimulPop
osp_6 <- simulPop(afeapop_6, 100, isize, mutrate, fintime)
pstats_6 <- popStats(osp_6, y)

# Cost a is bigger that benefit c
a <- 4
b <- 2
c <- 2

dfapop_7 <- createdf_apoptosis(a,b,c)
afeapop_7 <- wraperFitnessEffects(dfapop_7)
osi_7 <- simulOne(afeapop_7, isize, mutrate, fintime)
# extended simulations with oncoSimulPop
osp_7 <- simulPop(afeapop_7, 100, isize, mutrate, fintime)
pstats_7 <- popStats(osp_7, x)
#stripchart(pstats_7 ~ x, vertical = TRUE, pch = 1)

# Just strategies B(1) an C(2)
c <- -1
dfapop_8 <- createdf_apoptosis_strategiesBC(c)
afeapop_8 <- wraperFitnessEffects(dfapop_8)
osi_8 <- simulOne(afeapop_8, isize, mutrate, fintime)
# extended simulations with oncoSimulPop
osp_8 <- simulPop(afeapop_8, 100, isize, mutrate, fintime)
pstats_8 <- popStats(osp_8, y)

# Benefit of paracrine much bigger that benefit of autocrine
a <- 1
b <- 5
c <- 2

dfapop_9 <- createdf_apoptosis(a,b,c)
afeapop_9 <- wraperFitnessEffects(dfapop_9)
osi_9 <- simulOne(afeapop_9, isize, mutrate, fintime)
# extended simulations with oncoSimulPop
osp_9 <- simulPop(afeapop_9, 100, isize, mutrate, fintime)
pstats_9 <- popStats(osp_9, x)

