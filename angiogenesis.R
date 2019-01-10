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
# By default: initial population 5000, mutation rate 1e-5 
# and final time 1000
# Returns an ocosimul object

simulOne <- function(afe,isize=5000,mutrate=1e-5,fintime=1000){
  
  
  osi <- oncoSimulIndiv(afe,
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
  
  # Plot genotypes sizes along time
  plot(osi, show = "genotypes", type = "line")
  
  return(osi)
}


# Run a  number of oncoSimul trajectories defined by the parameter iter
# receives:
#            afe -  fitnessEffects object 
#            isize - initial size of the population
#            mutrate - mutation rate
#            fintime - maximum number ot time units to run
# By default: 100 iterations, initial population 5000, mutation rate 1e-5 
# and final time 1000
# Returns an ocosimulpop object

simulPop <- function(afe,iter=100,isize=5000,mutrate=1e-5,fintime=1000){
  

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


# Receives a dataframe like the one returned by function popStats:
#   Rows are the final samples from each oncoSimulIndiv
#   object inside an OncoSimulPop object
#   Colums correspond to the different genotypes
# Generates a barplot with the number of times each strategy outgrow the others
# considering also polymorphic equilibriums

getPercentages <- function (df_stats,initPopSize=5000){
  
  # remove 5% of initPopsize.
  # we will make the assumption that, if at the end of the 
  # oncoSImulIndiv simulation, population of a genotype 
  # is bigger that roundPop, it has achieved equilibrium
  roundPop <- initPopSize*5/10
  
  # remove NAs from the input data frame
  df_stats_cp <- df_stats[complete.cases(df_stats),] 
  na <- dim(df_stats)[1]-dim(df_stats_cp)[1]
  if ( na > 0 ){
    cat("Removed ",na," colums containing NA\n")
  }
  
  # Cosmetic change to show WT as label in the barplot 
  names(df_stats_cp)[1]<-paste("WT")
  
  #Set colour palette
  pal <- colorRampPalette(colors = c("lightblue", "blue"))(4)
  
  # Create a vector indicating the simulations
  # in which final population was significant for 
  # each genotype
  df_binary <-df_stats_cp > roundPop
  names(df_binary) <- names(df_stats_cp)
  
  winners <- vector('character')
    # Substitute TRUES for the name of the winner genotype 
  for (i in 1:dim(df_binary)[2]){
    df_binary[,i][df_binary[,i] == TRUE] <- names(df_binary)[i]
    df_binary[,i][df_binary[,i] == FALSE] <- ""
  }
  # Create a vector of winners
  for (i in 1:dim(df_binary)[1]){
    winners[i] <- paste(df_binary[i,], collapse = '')
  }
  
  # Plot the winning strategies
  barplot(table(winners), main = "Winning strategies", xlab = "Genotypes",col = pal)
 
  
}



#############################################
# Case Study: Angiogenesis
############################################

# Creates a dataframe with the genotypes and its
# fitness specification according to the angiogenesis
# game model, with the right format to be
# used to create a fitnessEffect object with the 
# allFitnessEffects function

createdf_angiogenesis <- function(i,j){
  
  # Strategies:
  # cells 1 can produce angiogenic factors at a fitness cost i 
  # cells 2  produce no angiogenic factors. 
  # In any case cells will get a benefit j when there is an interaction involving an angiogenic factor producing cell.
  
  
  f1 <- paste("1+f_1*(",as.character(1-i+j),")+f_2*(",as.character(1-i+j),")",sep="")
  f2 <- paste("1+f_1*(",as.character(1+j),")+f_2",sep="")
  
  
  r <- data.frame(Genotype = c("WT",
                               "1", "2"),
                  Fitness = c("1",
                              f1,
                              f2),
                  stringsAsFactors = FALSE)
  
  return(r)
  
}

# i is the cost of producing and angiogenic factor
# and j the benefit 
# i < j polimorphic equilibrium

j <- 3
i <- 1

dfang_1 <- createdf_angiogenesis(i,j)
afeang_11 <- wraperFitnessEffects(dfang_1)
osi_11 <- simulOne(afeang_11)
# extended simulations with oncoSimulPop
x <- c("", "1", "2")
osp_11 <- simulPop(afeang_11)
pstats_1 <- popStats(osp_11, x)
getPercentages(pstats_1)
write.table(pstats_1, file=file.path("/home/maria/Downloads/exam_OncoSimulR/results",'angiogenesis_pstats_1.txt'))

# i >= j 
# Strategy 2 displaces the others

i <- j <- 2
dfang_c2 <- createdf_angiogenesis(i,j)
afeang_c2 <- wraperFitnessEffects(dfang_c2)
osi_c2 <- simulOne(afeang_c2)
# extended simulations with oncoSimulPop
osp_c2 <- simulPop(afeang_c2)
pstats_c2 <- popStats(osp_c2, x)
getPercentages(pstats_c2)
write.table(pstats_c2, file=file.path("/home/maria/Downloads/exam_OncoSimulR/results",'angiogenesis_pstats_2.txt'))

j <- 1
i <- 3
dfang_c3 <- createdf_angiogenesis(i,j)
afeang_c3 <- wraperFitnessEffects(dfang_c3)
osi_c3 <- simulOne(afeang_c3)
# extended simulations with oncoSimulPop
osp_c3 <- simulPop(afeang_c3)
pstats_c3 <- popStats(osp_c3, x)
getPercentages(pstats_c3)
write.table(pstats_c3, file=file.path("/home/maria/Downloads/exam_OncoSimulR/results",'angiogenesis_pstats_3.txt'))



