#'Julia Hatamyar 26 October 2021
#'
#'This file creates/calls the various ML functions 
#'that will be used in the paper. 
#'
#'

# -----------------------------------------------------------------------------
# Chernozhukov Double-Debiased ML
library(causalDML)

DML = function(X, Y, W,
               cf = 5) {
  
  cDML = DML_aipw(x = X, y = Y, w = W, 
                  cf = cf, quiet = F)
  
  #get scores, and ATE (no cates this method)
  ATE.raw <- c(cDML$ATE$results[,1], cDML$ATE$results[,2])
  
  output = list("scores" = cDML$ATE$delta, "ATE" = ATE.raw)
  
  return(output)
}

# -----------------------------------------------------------------------------
# Dr-learner (Kennedy 2020) and NDR-learner (Knaus 2020)
NDR = function(X, Y, W,
               cf = 5) {
  
  nDR = ndr_learner(x = X, y = Y, w = W, 
                    nfolds = cf, quiet = F)
  
  #get cates
  cates.dr  <- nDR$cates[,1]
  cates.ndr <- nDR$cates[,2]
  
  #get scores, and ATE 
  # take average of the four folds' gammas and ates 
  
  #initialize 
  NDRates       <- matrix(NA,5, 4) # four rows plus avg, plus SE, t, p 
  
  NDRscores      <- matrix(NA, nrow = length(Y), ncol = 5)
  
  for(j in 1:4){
    
    for(i in 1:4){
      
      NDRates[i,j] <- nDR$list[[i]]$ATE$results[j]  
      
    }
  }
  
  
  NDRates[5,] <- colMeans(NDRates[1:4,])
  
  for(i in 1:4) {
    
    NDRscores[!nDR$cf_mat[,i],i] <- nDR$list[[i]]$ATE$delta
    
  }

  NDRscores[,5] <- rowMeans(NDRscores[,1:4], na.rm = TRUE)
  
  
  ATE.raw <- NDRates[5,1:2]
  scores.raw <- NDRscores[,5]
  
  output = list("drcates" = cates.dr,
    "ndrcates" = cates.ndr,
    "scores" = scores.raw, "ATE" = ATE.raw)
  
  return(output)
}

# -----------------------------------------------------------------------------

# Vanilla causal forest using out-of-bag predictions (honest forests)
forest_OOB <- function(X, Y, W, 
                       data = data,
                       tune.parameters = "all", 
                       tune.num.trees = 100,
                       tune.num.reps = 500,
                       num.trees = 2000) {
  
  
  
  Y.forest = regression_forest(X, Y,
                               tune.parameters = tune.parameters, 
                               tune.num.trees = tune.num.trees,
                               tune.num.reps = tune.num.reps,
                               num.trees = num.trees)
  
  Y.hat = predict(Y.forest)$predictions
  
  
  W.forest = regression_forest(X, W, 
                               tune.parameters = tune.parameters, 
                               tune.num.trees = tune.num.trees,
                               tune.num.reps = tune.num.reps,
                               num.trees = num.trees) 
  
  W.hat = predict(W.forest)$predictions
  
  cf <- causal_forest(X, Y, W,
                      Y.hat = Y.hat, 
                      W.hat = W.hat,
                      ci.group.size = 2,
                      tune.parameters = tune.parameters, 
                      tune.num.trees = tune.num.trees,
                      tune.num.reps = tune.num.reps,
                      num.trees = num.trees)
  
  cates.raw <- predict(cf)$predictions
  scores.raw <- get_scores(cf)
  ATE.raw <- average_treatment_effect(cf)
  
  output = list("cates" = cates.raw,
                "scores" = scores.raw, 
                "ATE" = ATE.raw)
  
  return(output)
}

# -----------------------------------------------------------------------------
# Causal forest with test and training sets 
# over 5 folds (default). Forests are estimated with all but i folds, 
# then tau is predicted for the holdout (ith) fold. This is repeated 
# `reps` number of times and the results are averaged across folds.
forest_test_train <- function(X, Y, W, 
                              K = 5, 
                              reps = 4, 
                              num.trees = 2000){
  
  r <- foreach(t = 1:reps, .combine='cbind') %do% { 
    
    print(paste("Repetition", toString(t)))
    
    ## initialize results for this repetition
    res      <- array(NA, dim = c(length(Y),1,5))
    dimnames(res)[[3]] <- c("cates", "errors", "scores", "gamma0", "gamma1")
    
    cates    <- array(NA, dim = length(Y))
    errors   <- array(NA, dim = length(Y))
    scores   <- array(NA, dim = length(Y))
    gamma0   <- array(NA, dim = length(Y))
    gamma1   <- array(NA, dim = length(Y))
    
    ## create fold indicators 
    split             <- runif(nrow(X)) 
    cvgroup           <- as.numeric(cut(split,quantile(split,probs = seq(0, 1, 1/K)),include.lowest = TRUE))  
    
    ## loop over folds
    for(i in 1:K){
      
      ## train using training set !i
      cf.train <- causal_forest(X[cvgroup != i,], Y[cvgroup != i], W[cvgroup != i],
                                num.trees = num.trees, # need more for accurate CIs ?
                                tune.parameters = c("alpha", "honesty.fraction"))
      
      # predict using holdout fold i
      print(paste("CF fold", toString(i)))
      test_predictions <- predict(cf.train, newdata = X[cvgroup == i,], estimate.variance = TRUE)
      
      # We need estimates of Y and W for the score : OK to do this with train/test set?
      Y.forest.train <- regression_forest(X[cvgroup != i,], Y[cvgroup != i]) ## change this later to include other args 
      
      W.forest.train <- regression_forest(X[cvgroup != i,], W[cvgroup != i]) 
      
      ## get nuisances from test set
      Y.hat.test    <- predict(Y.forest.train, X[cvgroup == i,])$predictions 
      W.hat.test    <- predict(W.forest.train, X[cvgroup == i,])$predictions
      
      
      # extract cates and variance
      tau.hat <- test_predictions$predictions
      se.est  <- sqrt(test_predictions$variance.estimates)
      
      ## save for this OOS fold index
      cates[cvgroup == i]  <- tau.hat
      errors[cvgroup == i] <- se.est
      
      ## get the scores (this is why I estimated W.hat and Y.hat again)
      ## NOTE this is adapted from "get_scores" grf source code
      debiasing.weights    <- (W[cvgroup == i] - W.hat.test) / (W.hat.test * (1 - W.hat.test))
      Y.residual           <- Y[cvgroup == i] - (Y.hat.test + tau.hat * (W[cvgroup == i] - W.hat.test))
      scores.raw            <- tau.hat + debiasing.weights * Y.residual
      
      scores[cvgroup == i] <- scores.raw # 
      
      ### now I need to get the gamma scores for policytree
      # conditional means
      Y.hat.0 <- Y.hat.test - W.hat.test * tau.hat
      Y.hat.1 <- Y.hat.test + (1 - W.hat.test) * tau.hat
      
      mu.matrix <- cbind("control" = Y.hat.0, "treated" = Y.hat.1)
      
      ## this is taken from double_robust_scores
      W.hat.matrix <- cbind(1 - W.hat.test, W.hat.test) # [control, treated]
      n.obs <- nrow(W.hat.matrix)
      observed.treatment.idx <- cbind(1:n.obs, W[cvgroup == i] + 1)
      
      YY <- matrix(0, n.obs, 2)
      IPW <- matrix(0, n.obs, 2)
      YY[observed.treatment.idx] <- Y[cvgroup == i]
      IPW[observed.treatment.idx] <- 1 / W.hat.matrix[observed.treatment.idx]
      Gamma.matrix <- (YY - mu.matrix) * IPW + mu.matrix
      
      gamma0[cvgroup == i] <- Gamma.matrix[,1]
      gamma1[cvgroup == i] <- Gamma.matrix[,2]
      
    }
    
    #get ATE and SE for this repetition
    ATE <- mean(scores) 
    SE  <- sqrt(mean((scores - ATE)^2) / length(W))
    
    # save to output arrays
    res[,,1] <- cates
    res[,,2] <- errors
    res[,,3] <- scores
    res[,,4] <- gamma0
    res[,,5] <- gamma1
    
    res <- c(as.vector(res), as.vector(ATE), as.vector(SE))
    
    r <- data.frame(res) # 
  }
  
  
  ### post estimation (ATE, etc)
  
  ## extract from large matrix
  cates_all  <- r[1:nrow(X),]
  errors_all <- r[(nrow(X)+1):(2*nrow(X)),]
  scores_all <- r[((2*nrow(X))+1):(3*nrow(X)),]
  gamma0_all <- r[((3*nrow(X))+1):(4*nrow(X)),]
  gamma1_all <- r[((4*nrow(X))+1):(5*nrow(X)),]
  
  ATE_all    <- r[(length(r[,1])-1),]
  SE_all     <- r[(length(r[,1])),]
  
  ## get column means (average over the repetitions)
  res_cates    <- rowMeans(cates_all)
  res_errors   <- rowMeans(errors_all)
  res_scores   <- rowMeans(scores_all)
  res_gamma0   <- rowMeans(gamma0_all)
  res_gamma1   <- rowMeans(gamma1_all)
  
  res_ATE      <- rowMeans(ATE_all)  # -0.003
  res_SE       <- rowMeans(SE_all)  
  
  res_all      <- list(res_cates, res_errors, res_scores, res_gamma0, res_gamma1, res_ATE, res_SE)
  names(res_all) <- c("cates", "errors", "scores", "gamma0", "gamma1", "ATE", "SE")
  
  return(res_all) 
  
  
}

# -----------------------------------------------------------------------------
# BART. notice that the X matrix here should include Y as a covariate, 
# or will have to be added together later. 
BART <- function(X, Y, W,
                 burnin = 1000,
                 draws = 1000,
                 ntree=50L){
  
  # combine treatment as column in X
  X$trt_ind <- W
  
  # first estimate prop score to add as a covariate
  prop_bart <- pbart(
    x.train = select(X, -trt_ind), 
    y.train = pull(X, trt_ind),  
    nskip = 1000,
    ndpost = 1000
  )
  
  # store propensity score in data
  X$prop_score <-  prop_bart$prob.train.mean
  
  bart_fit = pbart(x.train = X, y.train = Y, 
                   ndpost = draws, nskip = burnin, ntree = ntree)
  
  # get offset for later predictions
  offset <- bart_fit$binaryOffset
  
  # separate draws of potential outcomes for treated 
  Y1_1 = X[X$trt_ind==1,]
  Y1_0 = Y1_1
  
  # switch to counterfactual (treatment = 0)
  Y1_0$trt_ind = ifelse(Y1_0$trt_ind==1, 0, 1)
  
  # get posterior draws (1000 of them)
  bart_pred11 = pwbart(Y1_1, bart_fit$treedraws)
  bart_pred10 = pwbart(Y1_0, bart_fit$treedraws)
  
  # get distribution from posterior draws 
  # note this is adjusted by the offset to avoid 
  # squeezing the probabilities to 0.5
  pred_prop11 = pnorm(bart_pred11+offset)
  pred_prop10 = pnorm(bart_pred10+offset)
  
  
  # separate draws of potential outcomes for control 
  Y0_0 = X[X$trt_ind==0,]
  Y0_1 = Y0_0
  
  # switch to counterfactual (treatment = 0)
  Y0_1$trt_ind = ifelse(Y0_0$trt_ind==0, 1, 0)
  
  # get posterior draws (1000 of them)
  bart_pred00 = pwbart(Y0_0, bart_fit$treedraws)
  bart_pred01 = pwbart(Y0_1, bart_fit$treedraws)
  
  # get distribution from posterior draws 
  pred_prop00 = pnorm(bart_pred00+offset)
  pred_prop01 = pnorm(bart_pred01+offset)
  
  n1 = sum(X$trt_ind==1) # treatment 
  n2 = sum(X$trt_ind==0) # control
  
  # initialize values
  RD.est = NULL
  RR.est = NULL
  tau.est = NULL
  y1 = NULL
  y0 = NULL
  y1.draws <- matrix(NA,length(Y),draws)
  y0.draws <- matrix(NA,length(Y),draws)
  
  for (m in 1:draws) {
    
    # Estimate E(Y1), E(Y0),
    y1[X$trt_ind==1] = rbinom(n1, 1, pred_prop11[m,])
    y1[X$trt_ind!=1] = rbinom(n2, 1, pred_prop01[m,])
    y0[X$trt_ind==0] = rbinom(n2, 1, pred_prop00[m,])
    y0[X$trt_ind!=0] = rbinom(n1, 1, pred_prop10[m,])
    
    y1.draws[,m] = (y1)
    y0.draws[,m] = (y0)
    
    y1.pred = mean(y1)
    y0.pred = mean(y0)
    
    # Calculate risk difference (RD): this is the ATE per draw?
    RD.est[m] = y1.pred - y0.pred
    
  }
  
  ate.raw <- mean(RD.est)
  se.raw  <- sd(RD.est)
  
  y1.avg <- rowMeans(y1.draws)
  y0.avg <- rowMeans(y0.draws)
  cates.raw <- y1.avg - y0.avg
  
  # get scores
  scores.raw <- get_scores_BART(W = W, W.hat = X$prop_score,
                            Y = Y, Y.hat.0 = y0.avg, Y.hat.1 = y1.avg)
  
  
  output = list("cates" = cates.raw,
                "scores" = scores.raw, "ATE" = c(ate.raw, se.raw))
  
  return(output)
  
}


# testing 
#library(dplyr)
#tdata <- data_make_het(n = 500)
#x.train <- tdata$trtdat
#x.train <- select(x.train, -trt_ind)
#x.bart <- data$trtdat

#bart <- BART(X = x.train, Y = data$Yobs, W = data$trtdat$trt_ind)

#bart_fit = pbart(x.train, y.train = data$Yobs, 
#                  k = 3,  
#                  ndpost = 1000, nskip = 1000)

#nDR = ndr_learner(x = x.train, y = data$Yobs, w = data$trtdat$trt_ind)

#nDR = NDR(X = x.train, Y = data$Yobs, W = data$trtdat$trt_ind)

#cDML = DML(X = x.train, Y = tdata$Yobs, W = tdata$trtdat$trt_ind)

#cDML = DML_aipw(x = x.train, y = data$Yobs, w = data$trtdat$trt_ind)

#hist(cDML$ATE$delta, breaks = 100)
#max(cDML$ATE$delta)
#summary(cDML$ATE$results)
#printCoefmat(cDML$ATE$results,has.Pvalue = TRUE)

#setwd("~/Research/IMR")

