library("OncoSimulR")
library("ggplot2")

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


#a is the cost of producing the paracrine factor
#b the benefit of receiving the paracrine factor
#c the benefit of producing the autocrine factor

a <- 2
b <- 6
c <- 3

r <- createdf_apoptosis(a,b,c)

osi<-simulation(r)









#############################################
# Case Study: Angiogenesis
############################################

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
j <- 2
i <- 1
isize <- 5000
mutrate <- 1e-5
fintime <- 1000

dfang_1 <- createdf_angiogenesis(i,j)
afeang_11 <- wraperFitnessEffects(dfang_1)
simulOne(afeang_11,isize,mutrate,fintime)
# extended simulations with oncoSimulPop
x <- c("", "1", "2")
osp_1 <- simulPop(afeang_11, 100, isize=5000, mutrate=1e-5, fintime=1000)
popstats(osp_1, x)

popSizes_1=c(5000,1000,100)
afeagng_12 <- wraperFitnessEffects(dfang_1,popSizes_1)
genotypes(afeang_12)
simulOne(afeang_12,isize=5000,mutrate=1e-5,fintime=1000)
osp_2 <- simulPop(afeang_12, 100, isize=5000, mutrate=1e-5, fintime=1000)
popstats(osp_2, x)

popSizes_2=c(5000,100,1000)
afeang_13 <- wraperFitnessEffects(dfang_1,popSizes_2)
genotypes(afeagng_13)
simulOne(afeang_13,isize=5000,mutrate=1e-5,fintime=1000)
osp_3 <- simulPop(afeang_13, 100, isize=5000, mutrate=1e-5, fintime=1000)
popstats(osp_3, x)


# i > j 
# Stratey 2 displaces the others
j <- 1
i <- 2
df3 <- createdf_angiogenesis(j,i)
simulation(df3)