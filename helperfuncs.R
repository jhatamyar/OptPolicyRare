#'
#'Helper functions for Optimal Policy with Rare Outcomes
#'Julia Hatamayar Oct, 2021
#'
#'========================================
#'


#' This function calculates the double-robust score as in 
#' Athey and Wager 2021. 
#' Note: this is implicitly using the following
# estimates for the regression surfaces E[Y|X, W=0/1]:
# Y.hat.0 <- Y.hat - W.hat * tau.hat
# Y.hat.1 <- Y.hat + (1 - W.hat) * tau.hat
#' 
#' @param W original treatment vector
#' @param Y observed outcome vector
#' @param Y.hat expected conditional outcome (MARGINALIZED OVER TX)
#' What does this mean? The Gamma0 amd Gamma1 in double_robust_scores subtract to this value
#' @param W.hat estimated propensity score 
#' @param tau.hat Estimated CATE
#'
dr_scores <- function(W, Y, Y.hat, W.hat, tau.hat){
  debiasing.weights    <- (W - W.hat) / (W.hat * (1 - W.hat))
  Y.residual           <- Y - (Y.hat + tau.hat * (W - W.hat))
  scores.raw            <- tau.hat + debiasing.weights * Y.residual
  
  return(scores = scores.raw)
}

#'========================================
#'
#'This function calculates the DR-scores from a BART model.  
#' We cannot marginalize over treatment to get Y.hat so use the separate estimates for the 
#' conditional outcomes then subtract to get overall score. 
#' 
#' Code adapted from 'double_robust_scores' in the policytree package. 
#' 

get_scores_BART <- function(W, W.hat, Y, Y.hat.0, Y.hat.1) {
  mu.matrix <- cbind("control" = Y.hat.0, "treated" = Y.hat.1)
  
  W.hat.matrix <- cbind(1 - W.hat, W.hat) # [control, treated]
  n.obs <- nrow(W.hat.matrix)
  observed.treatment.idx <- cbind(1:n.obs, W + 1)
  
  ## initialize values 
  YY <- matrix(0, n.obs, 2)
  IPW <- matrix(0, n.obs, 2)
  
  YY[observed.treatment.idx] <- Y
  IPW[observed.treatment.idx] <- 1 / W.hat.matrix[observed.treatment.idx]
  Gamma.matrix <- (YY - mu.matrix) * IPW + mu.matrix
  
  scores.raw <- Gamma.matrix[,2] - Gamma.matrix[,1]
  return(scores = scores.raw)
}



#'========================================
#'
#'This function calculates the advantage of a policy
#'compared to the mysterious "random baseline" 
#'
#'@param policy a vector of treatment assignment from a policy, i.e. 
#'1 = treat, 0 = don't treat
#'@param scores a vector of double-robust scores 

get_advantage = function(policy, scores) {
  mu = mean((2 * policy - 1) * (scores), na.rm = T)
  se = sqrt(var((2 * policy - 1) * (scores), na.rm = T) / length(scores))
  c(point.estimate=mu, std.err=se)
}


