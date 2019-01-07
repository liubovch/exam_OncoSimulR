library("OncoSimulR")
library("ggplot2")
library("doParallel")

set.seed(1)

####################################
# Functions
####################################

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

genotypes <- function(afe,spPS){
  eval <- evalAllGenotypes(afe, addwt = TRUE)
  plotFitnessLandscape(eval)
  eval
}

simulOne <- function(afe,isize,mutrate,fintime){

  ## Run one trayectory simulation
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

  ## Plot genotypes sizes along time
  plot(osi, show = "genotypes", type = "line")
  plotClonePhylog(osi, N = 0)
  return(osi)
}

simulPop <- function(afe,iter,isize,mutrate,fintime){

  ## Run a  number of trayectory simulations defined by the parameter iter
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






#############################################
# Case Study : Evasion of apoptosis
############################################
createdf_apoptosis <- function(a, b, c ){

# Strategies:
#  a. Cells that produce a paracrine growth factor to prevent apoptosis of neighbouring cells. 
#  b. Cells that produce an autocrine growth factor to prevent apoptosis of themselves. 
#  c. Cells susceptible to paracrine growth factors but incapable of production of factors. 
  

fa <- paste("1+f_A*(",as.character(1-a+b),")+f_B*(",as.character(1-a),")+f_C*(",as.character(1-a),")",sep="")
fb <- paste("1+f_A*(",as.character(1+b+c),")+f_B*(",as.character(1+c),")+f_C*(",as.character(1+c),")",sep="")
fc <- paste("1+f_A*(",as.character(1+b),")+f_B+f_C",sep="")

r <- data.frame(Genotype = c("WT",
                             "A", "B","C"),
                Fitness = c("1",
                            fa,
                            fb,
                            fc),
                stringsAsFactors = FALSE)

return(r)

}


createdf_apoptosis_strategiesAC <- function(a,b){
  
  # Strategies:
  #  a. Cells that produce a paracrine growth factor to prevent apoptosis of neighbouring cells. 
  #  b. Cells that produce an autocrine growth factor to prevent apoptosis of themselves. 
  #  c. Cells susceptible to paracrine growth factors but incapable of production of factors. 
  
  
  
  fa <- paste("1+f_1*(",as.character(1-a+b),")+f_2*(",as.character(1-a),")",sep="")
  fc <- paste("1+f_1*(",as.character(1+b),")+f_2",sep="")
  
  r <- data.frame(Genotype = c("WT",
                               "1","2"),
                  Fitness = c("1",
                              fa,
                              fc),
                  stringsAsFactors = FALSE)
  
  return(r)
  
}


createdf_apoptosis_strategiesBC <- function(c){
  
  # Strategies:
  #  a. Cells that produce a paracrine growth factor to prevent apoptosis of neighbouring cells. 
  #  b. Cells that produce an autocrine growth factor to prevent apoptosis of themselves. 
  #  c. Cells susceptible to paracrine growth factors but incapable of production of factors. 
  
  
  
  fb <- paste("1+f_1*(",as.character(1+c),")+f_2*(",as.character(1+c),")",sep="")
  fc <- "1+f_1+f_2"
  
  r <- data.frame(Genotype = c("WT",
                               "1","2"),
                  Fitness = c("1",
                              fb,
                              fc),
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


dfapop_1 <- createdf_apoptosis(a,b,c)
afeapop_1 <- wraperFitnessEffects(dfapop_1)
osi_1 <- simulOne(afeapop_1, isize, mutrate, fintime)
# extended simulations with oncoSimulPop
x <- c("", "A", "B", "C")
osp_1 <- simulPop(afeapop_1, 100, isize, mutrate, fintime)
pstats_1 <- popStats(osp_1, x)

a <- 1
b <- 2
c <- 3

dfapop_2 <- createdf_apoptosis(a,b,c)
afeapop_2 <- wraperFitnessEffects(dfapop_2)
osi_2 <- simulOne(afeapop_2, isize, mutrate, fintime)
# extended simulations with oncoSimulPop
osp_2 <- simulPop(afeapop_2, 100, isize, mutrate, fintime)
pstats_2 <- popStats(osp_2, x)

a <- 0
b <- 2
c <- 3

dfapop_3 <- createdf_apoptosis(a,b,c)
afeapop_3 <- wraperFitnessEffects(dfapop_3)
osi_3 <- simulOne(afeapop_3, isize, mutrate, fintime)
# extended simulations with oncoSimulPop
osp_3 <- simulPop(afeapop_3, 100, isize, mutrate, fintime)
pstats_3 <- popStats(osp_3, x)


a <- 0
b <- 3
c <- 2

dfapop_4 <- createdf_apoptosis(a,b,c)
afeapop_4 <- wraperFitnessEffects(dfapop_4)
osi_4 <- simulOne(afeapop_4,isize,mutrate,fintime)
# extended simulations with oncoSimulPop
osp_4 <- simulPop(afeapop_4, 100, isize, mutrate, fintime)
pstats_4<- popStats(osp_4, x)

# Just strateies B(1) an C(2)
c <- 1
dfapop_5 <- createdf_apoptosis_strategiesBC(c)
afeapop_5 <- wraperFitnessEffects(dfapop_5)
osi_5 <- simulOne(afeapop_5, isize, mutrate, fintime)
# extended simulations with oncoSimulPop
osp_5 <- simulPop(afeapop_5, 100, isize, mutrate, fintime)
pstats_5<- popStats(osp_5, x)


a <- 1
b <- 3
#Just strategies A(1) anc C(2)
dfapop_6 <- createdf_apoptosis_strategiesAC(a, b)
afeapop_6 <- wraperFitnessEffects(dfapop_6)
osi_6 <- simulOne(afeapop_6, isize, mutrate, fintime)
# extended simulations with oncoSimulPop
osp_6 <- simulPop(afeapop_6, 100, isize, mutrate, fintime)
pstats_6<- popStats(osp_6, x)

# Cost a is bigger that benefit c
a <- 4
b <- 2
c <- 2

dfapop_7 <- createdf_apoptosis(a,b,c)
afeapop_7 <- wraperFitnessEffects(dfapop_7)
osi_7 <- simulOne(afeapop_7,isize,mutrate,fintime)
# extended simulations with oncoSimulPop
osp_7 <- simulPop(afeapop_7, 100, isize, mutrate, fintime)
pstats_7<- popStats(osp_7, x)
#stripchart(pstats_7 ~ x, vertical = TRUE, pch = 1)
