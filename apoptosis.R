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
# Calculates the number of cases in which polymorfic equilibrium is achieved 
# Generates a barplot with the number of times each strategy outgrow the others

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
  
  # Get % of wins for each strategy
  pal <- colorRampPalette(colors = c("lightblue", "blue"))(4)
  barplot(colSums(df_stats_cp >= apply(df_stats_cp, 1, max)), main = "Winning strategies", 
          xlab = "Genotypes",col = pal)
  
  # get % of polymorfic equilibrium
  polim_eq <-  sum(rowSums(df_stats_cp > roundPop) >=2)
  
  # If we have polimorphic equilibrium, remove this cases to calculate the number of times 
  # each strategy/genotype wins in solitary
  if (polim_eq > 0){
    # creates a results matrix without the polimorphic equilibrium cases
    cat("Number of polymorphic equilibriums achived: "
        , polim_eq,"\n")
    df_stats_nopoleq <- df_stats_cp[rowSums(df_stats_cp > roundPop)==1,]
    # Get % of solely wins for each strategy
    print("Winners in solitary %\n")
    print(colSums(df_stats_nopoleq  >= apply(df_stats_nopoleq , 1, max))) 
  }
  
  
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






############################################
# Case Study : Evasion of apoptosis
############################################
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


createdf_apoptosis_strategiesAC <- function(a,b){
  
  # Strategies:
  #  a. Cells that produce a paracrine growth factor to prevent apoptosis of neighbouring cells. 
  #  c. Cells susceptible to paracrine growth factors but incapable of production of factors. 
  
  
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


createdf_apoptosis_strategiesBC <- function(c){
  
  # Strategies:
  #  b. Cells that produce an autocrine growth factor to prevent apoptosis of themselves. 
  #  c. Cells susceptible to paracrine growth factors but incapable of production of factors. 
  
  
  
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

# equal benefit of paracrine and autocrine (removing cost or producing paracrine)
a <- 1
b <- 3
c <- 2


dfapop_1 <- createdf_apoptosis(a, b, c)
afeapop_1 <- wraperFitnessEffects(dfapop_1)
osi_1 <- simulOne(afeapop_1)
# extended simulations with oncoSimulPop
x <- c("", "A", "B", "C")
osp_1 <- simulPop(afeapop_1)
pstats_1 <- popStats(osp_1, x)
getPercentages(pstats_1)
write.table(pstats_1, file=file.path("/home/maria/Downloads/exam_OncoSimulR/results",'apoptosis_pstats_1.txt'))

# benefit of autocrine factor is much bigger than paracrine
a <- 1
b <- 2
c <- 8

dfapop_2 <- createdf_apoptosis(a, b, c)
afeapop_2 <- wraperFitnessEffects(dfapop_2)
osi_2 <- simulOne(afeapop_2)
# extended simulations with oncoSimulPop
osp_2 <- simulPop(afeapop_2)
pstats_2 <- popStats(osp_2, x)
getPercentages(pstats_2)
write.table(pstats_2, file=file.path("/home/maria/Downloads/exam_OncoSimulR/results",'apoptosis_pstats_2.txt'))


# cost of paracrine is 0, benefit of paracrine bigger than benefit of autocrine 
a <- 0
b <- 2
c <- 3

dfapop_3 <- createdf_apoptosis(a, b, c)
afeapop_3 <- wraperFitnessEffects(dfapop_3)
osi_3 <- simulOne(afeapop_3)
# extended simulations with oncoSimulPop
osp_3 <- simulPop(afeapop_3)
pstats_3 <- popStats(osp_3, x)
getPercentages(pstats_3)
write.table(pstats_3, file=file.path("/home/maria/Downloads/exam_OncoSimulR/results",'apoptosis_pstats_3.txt'))

# benefit of paracrine much bigger than benefit from autocrine
a <- 1
b <- 8
c <- 2

dfapop_4 <- createdf_apoptosis(a, b, c)
afeapop_4 <- wraperFitnessEffects(dfapop_4)
osi_4 <- simulOne(afeapop_4)
# extended simulations with oncoSimulPop
osp_4 <- simulPop(afeapop_4)
pstats_4 <- popStats(osp_4, x)
getPercentages(pstats_4)
write.table(pstats_4, file=file.path("/home/maria/Downloads/exam_OncoSimulR/results",'apoptosis_pstats_4.txt'))

y <- c("", "1", "2")

# Just strategies B(1) an C(2)
c <- 2
dfapop_5 <- createdf_apoptosis_strategiesBC(c)
afeapop_5 <- wraperFitnessEffects(dfapop_5)
osi_5 <- simulOne(afeapop_5)
# extended simulations with oncoSimulPop
osp_5 <- simulPop(afeapop_5)
pstats_5 <- popStats(osp_5, y)
getPercentages(pstats_5)
write.table(pstats_5, file=file.path("/home/maria/Downloads/exam_OncoSimulR/results",'apoptosis_pstats_5.txt'))



a <- 1
b <- 3
#Just strategies A(1) anc C(2)
dfapop_6 <- createdf_apoptosis_strategiesAC(a, b)
afeapop_6 <- wraperFitnessEffects(dfapop_6)
osi_6 <- simulOne(afeapop_6)
# extended simulations with oncoSimulPop
osp_6 <- simulPop(afeapop_6)
pstats_6 <- popStats(osp_6, y)
getPercentages(pstats_6)
write.table(pstats_6, file=file.path("/home/maria/Downloads/exam_OncoSimulR/results",'apoptosis_pstats_6.txt'))

# Cost a is bigger that benefit c
a <- 4
b <- 2
c <- 2

dfapop_7 <- createdf_apoptosis(a,b,c)
afeapop_7 <- wraperFitnessEffects(dfapop_7)
osi_7 <- simulOne(afeapop_7)
# extended simulations with oncoSimulPop
osp_7 <- simulPop(afeapop_7)
pstats_7 <- popStats(osp_7, x)
getPercentages(pstats_7)
write.table(pstats_7, file=file.path("/home/maria/Downloads/exam_OncoSimulR/results",'apoptosis_pstats_7.txt'))

# Just strategies B(1) an C(2)
c <- -1
dfapop_8 <- createdf_apoptosis_strategiesBC(c)
afeapop_8 <- wraperFitnessEffects(dfapop_8)
osi_8 <- simulOne(afeapop_8)
# extended simulations with oncoSimulPop
osp_8 <- simulPop(afeapop_8, 100)
pstats_8 <- popStats(osp_8, y)
getPercentages(pstats_8)
write.table(pstats_8, file=file.path("/home/maria/Downloads/exam_OncoSimulR/results",'apoptosis_pstats_8.txt'))

# Benefit of autocrine negative
a <- 1
b <- 3
c <- -1

dfapop_9 <- createdf_apoptosis(a,b,c)
afeapop_9 <- wraperFitnessEffects(dfapop_9)
osi_9 <- simulOne(afeapop_9)
# extended simulations with oncoSimulPop
osp_9 <- simulPop(afeapop_9)
pstats_9 <- popStats(osp_9,x)
getPercentages(pstats_9)
write.table(pstats_9, file=file.path("/home/maria/Downloads/exam_OncoSimulR/results",'apoptosis_pstats_9.txt'))
