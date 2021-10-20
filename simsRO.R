#'
#'Simulation for rare binary outcome data 
#'with varying levels of overlap and heterogeneity
#'Julia Hatamayar Oct, 2021
#'
#'========================================
#'
#' @param n The number of observations
#' @param overlap The degree of covariate overlap in the treatment assignment mechanism. 
#' This will impact the distribution of the true propensity score. Options are "strong", "moderate", "weak". 
#' @param prevalence The degree of outcome prevalence. "low" is a 2-3% prevalence consistent with 
#' the IFLS data. "normal" is a 40-50% outcome prevalence. 
#' @param setting The treatment heterogeneity. "1" is nonlinear covariate effect in treated group only.
#' "2" is nonlinear covariate effect in both treatment and control groups. "3" is as 2, but with covariate het effect
#' in the treatment group generated from some different covariates than that in the control group. 
#' 
#' Adapted from: https://github.com/liangyuanhu/CIMTx/tree/master/R 
#' 
#' @return A list containing the simulated data, the observed outcome, tru propensity score, and number of obs

expit = function(x) {exp(x)/(1+exp(x))}

data_make_het = function(n = 10000,  
                         overlap = "moderate", prevalence = "low",
                         setting = 1) {
  # default of 10 confounders but can increase
  p = 10
  
  #continuous covariates
  Xcon = matrix(rnorm(p*n), nrow=n, ncol=p/2)
  x1 = Xcon[,1]; x2 = Xcon[,2]; x3 = Xcon[,3]; x4 = Xcon[,4]; x5 = Xcon[,5]
  
  #categorical covariates
  Xcat = matrix(sample(0:1, n*p, replace = T, prob = c(.5,.5)), nrow=n, ncol=p/2)
  x6 = Xcat[,1]; x7 = Xcat[,2]; x8 = Xcat[,3]; x9 = Xcat[,4]; x10 = Xcat[,5]
  
  #generate treatment label W depending on overlap argument 

  if (overlap == "strong") {
    p1 <- expit(
      0.8 -0.1 * x1 - .9 * x2 - .3 * x3  - .1 * x5 - 0.1 * x6 -0.2*x7 - 0.4 * x9 + .5 * x10
    )
  }
  
  
  if (overlap == "moderate") {
    p1 <- expit(
      1.6-0.1 * 2 * x1 - .9 * 2 * x3 - .3 * 2 * x3  - .1 * 2 * x5 - 0.1 * 2 * x6 -0.3 * 2* x7 - 0.4 * 2 * x9 + .5 * 2 * x10
    )
  }
  
  if (overlap == "weak") {
    p1 <- expit(
      3.1 - 0.1 * 5 * x1 - .9 * 5 * x3 - .3 * 5 * x3  - .1 * 5 * x5 - 0.1 * 5 * x6 -0.3 * 5* x7 - 0.4 * 5 * x9 + .5 * 5 * x10
    )
  }
  
  
  p2 <- 1 - p1
  
  W = NULL
  for (i in 1:n) {
    W[i] <- sample(c(0, 1),
                   size = 1,
                   replace = TRUE,
                   prob = c(p1[i], p2[i]))
  }
  #table(W)
  
  #parallel response surface model: note, outcomes under treatment/control specified differently to create het
  
  if (prevalence == "normal") {
    alpha1 = alpha2 = 0
  } 
  
  if (prevalence == "low") {
    if (setting == 1){
      alpha1 = -5
      alpha2 = -5.6
    }
    
    if (setting == 2 | 3) {
      alpha1 = -4
      alpha2 = -4.2
    }
  }
  
  # nonlinear covariate effect in treated group only
  if (setting ==1) {
    Yp1 = expit(alpha1 - 0.5 * x1 - 0.8 * x3 - 1.8 * x5 - 0.9 * x6 - 0.1 * x7)
    Yp2 = expit(alpha2 + 0.1 * expit(x1) + 0.1 * sin(x3) - 0.1 * x5 ^ 2   - 0.3 * x6 - 0.2 * x7)
    
  }
  
  # nonlinear covariate effect in both groups 
  if (setting ==2) {
    Yp1 = expit(alpha1 + 0.1 * x1 ^ 2 - 0.2 * sin(x3) + 0.2 * expit(x5) + 0.2 * x6 -0.3 * x7)
    Yp2 = expit(alpha2 + 0.1 * expit(x1) + 0.1 * sin(x3) - 0.1 * x5 ^ 2   - 0.3 * x6 - 0.2 * x7)
    
  }
  
  # nonlinear covariate effect in both groups from different confounders 
  if (setting ==3) {
    Yp1 = expit(alpha1 + 0.1 * x2 ^ 2 - 0.2 * sin(x4) + 0.2 * expit(x5) + 0.2 * x6 -0.3 * x8)
    Yp2 = expit(alpha2 + 0.1 * expit(x1) + 0.1 * sin(x3) - 0.1 * x5 ^ 2   - 0.3 * x6 - 0.2 * x7)
    
  }
  
  
  
  #potential outcomes
  Y1 = Y2 = NULL
  for (i in 1:n) {
    Y1[i] = rbinom(1, 1, Yp1[i])
    Y2[i] = rbinom(1, 1, Yp2[i])
  }
  #table(Y1)[2]/sum(table(Y1)) 
  #table(Y2)[2]/sum(table(Y2)) 
  
  Y = cbind(Y1, Y2)
  YW = cbind(Y, W)
  
  Yobs <- apply(YW, 1, function(x) x[1:2][x[3]+1])
  #table(Yobs,W)
  
  trtdat = data.frame(W, Xcon, Xcat)
  colnames(trtdat) = c("trt_ind", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10")
  
  return(list(trtdat=trtdat, Y=Y, Yobs=Yobs, pscore = p2, n=n))
  
}

#uncomment to test 
#test <- data_make_het(overlap = "strong", prevalence = "normal", setting = 3)
#table(test$trtdat$trt_ind)
#table(test$Y[,1])
#table(test$Y[,2])
#hist(test$pscore)