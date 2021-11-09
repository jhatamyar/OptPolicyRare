#Set to my WD
#setwd("/Users/juliahatamyar/Documents/Research/Indonesia_IMR")
setwd("~/scratch/Indonesia_IMR")

#load required libraries
rm(list=ls(all=TRUE))
vec.pac= c("grf", "causalDML", "BART", "nnet", "policytree",
           "ggplot2", "dplyr", "doParallel")

#install packages if not already installed 
list.of.packages <- vec.pac
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(vec.pac, require, character.only = TRUE) 

# source files
source("simsRO.R")
source("helperfuncs.R")
source("MLfuncs.R")

## set to run in parallel 
set.seed(1211);

## my laptop has 8 cores (n-1 = 7), set to 40 for Viking 
cl   <- makeCluster(40, outfile="")
registerDoParallel(cl)

#------------------------------------------------------------------------
# setting 1: not rare prevalence, moderate overlap, setting 1, n = 800

nsim = 100
#start_time <- Sys.time()

r1 <- foreach(t = 1:nsim, .combine='cbind', .inorder=FALSE, .packages=vec.pac) %dopar% { 
  
  # set seed to t 
  
  ## generate data (change n to fewer obs for testing)
  data <- data_make_het(prevalence = "normal")
  
  X <- data$trtdat
  X <- select(X, -trt_ind)
  Y = data$Yobs 
  W = data$trtdat$trt_ind
  
  # save data every few reps (for examination if needed)
  if (t %% 25 == 0) {
    write.csv(data, paste0("sim1data", t, ".csv"))
  }
  
  ## want to save optimal policy and cates for each sim 
  cates.OR <- data$Y[,2] - data$Y[,1]
  or.ATE <- mean(cates.OR)
  
  # to compare to upper and lower bounds of Athey figure 2
  #upper.OR <- mean(pmin(cates.OR, -cates.OR)) 
  #lower.OR <- min(mean(cates.OR), mean(-cates.OR)) 
  
  # policy where everyone with negative CATE is treated
  # note, this is the same as upper.OR but can be used as input to get_advantage
  policy.OR <- ifelse(cates.OR == -1, 1, 0)
  advantage.OR <- get_advantage(policy.OR, cates.OR)
  
  ############################# RUN METHODS
  cDML <- DML(X, Y, W)
  nDR  <- NDR(X, Y, W) 
  CF   <- forest_OOB(X, Y, W)
  CFTT <- forest_test_train(X, Y, W)
  BART <- BART(X, Y, W, burnin = 200, draws = 200) # change later
  
  
  ## initialize results do I need to do this (freaking ugly code, but need to keep separate for each method)
  cates.DR    <- nDR$drcates
  cates.NDR    <- nDR$ndrcates
  cates.CF    <- CF$cates
  cates.CFTT    <- CFTT$cates
  cates.BART    <- BART$cates
  
  scores.DML <- cDML$scores[,1]
  scores.DR <- nDR$scores
  scores.CF <- CF$scores
  scores.CFTT <- CFTT$scores
  scores.BART <- BART$scores
  
  ################## policytree estimation
  ## get policy from scores
  
  scores <- cbind(scores.DML, scores.DR, scores.CF, scores.CFTT, scores.BART)
  policies <- matrix(NA, length(Y), 5)
  X.safe <- as.matrix(X)
  for (i in 1:5){
    score <- as.matrix(scores[,i])
    policies[,i] <- get_policy(score, X.safe)
  }
  
  policy.DML <- policies[,1]
  policy.DR <- policies[,2]
  policy.CF <- policies[,3]
  policy.CFTT <- policies[,4]
  policy.BART <- policies[,5]
  
  # get advantages
  ## first using estimated cates
  advantages.est <- matrix(NA, 5, 2)
  
  for (i in 1:5) {
    score <- scores[,i]
    policy <- policies[,i]
    advantages.est[i,] <- get_advantage(policy, score)
  }
  
  ## then using ACTUAl cates
  advantages.OR <- matrix(NA, 5, 2)
  
  for (i in 1:5) {
    policy <- policies[,i]
    advantages.OR[i,] <- get_advantage(policy, cates.OR)
  }
  
  
  ## structuring the output: each method will have its own vector (length = 1505), then these combined. 
  ## each column goes: 1. cates, 2. scores 3. ATE 4. policy 5. advantage 
  
  res.OR <- c(cates.OR, cates.OR, or.ATE, policy.OR, advantage.OR, advantage.OR) # note, some repeating to match dims 
  res.DML <- c(scores.DML, scores.DML, mean(scores.DML), policy.DML, advantages.est[1,], advantages.OR[1,])
  res.DR  <- c(cates.DR, scores.DR, mean(scores.DR), policy.DR, advantages.est[2,], advantages.OR[2,])
  res.NDR <- c(cates.NDR, scores.DR, mean(scores.DR), policy.DR, advantages.est[2,], advantages.OR[2,])
  res.CF  <- c(cates.CF, scores.CF, mean(scores.CF), policy.CF, advantages.est[3,], advantages.OR[3,])
  res.CFTT  <- c(cates.CFTT, scores.CFTT, mean(scores.CFTT), policy.CFTT, advantages.est[4,], advantages.OR[4,])
  res.BART <- c(cates.BART, scores.BART, mean(scores.BART), policy.BART, advantages.est[5,], advantages.OR[5,])
  
  res <- c(res.OR, res.DML, res.DR, res.NDR, res.CF, res.CFTT, res.BART)
  
  r1 <- data.frame(res)
  
}

#end_time <- Sys.time()

write.csv(r1, file = "SIM1_est.csv") # because it takes forever 

