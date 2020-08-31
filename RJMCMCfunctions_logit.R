#---------------------------------------------------------------------------------------------------------------------------------------
# RJMCMC 
# functions 
#---------------------------------------------------------------------------------------------------------------------------------------
# MH step
# this performs the within models update using a Metropolis-Hastings update
updateparam_logit <- function(beta_vec, sigma, llikhood, prior_beta = "norm", prior_beta_par, prop_beta = "norm", prop_beta_par)
{

  # update each beta coefficient that is not equal to 0 at the current iteration 
  for(j in 1:length(beta_vec))
  {
    beta_st <- beta_vec
    
    if(beta_st[j] != 0)
    {
      # random walk update N(current value, chosen sigma)
      beta_st[j] <- get(paste("r", prop_beta, sep=""))(1, beta_vec[j], prop_beta_par)
    
      # log-likelihood evaluated at proposed beta values
      newllikhood <- sum(dbinom(y, 1, plogis(X%*%beta_st), log = T))
    
      # log of numerator of acceptance probability 
      num <- newllikhood + get(paste("d", prior_beta, sep=""))(beta_st[j], prior_beta_par[1], prior_beta_par[2], log=TRUE) + get(paste("d", prop_beta, sep=""))(beta_vec[j], beta_st[j], prop_beta_par, log=TRUE) 
          
      den <- llikhood + get(paste("d", prior_beta, sep=""))(beta_vec[j], prior_beta_par[1], prior_beta_par[2], log=TRUE) + get(paste("d", prop_beta, sep=""))(beta_st[j], beta_vec[j], prop_beta_par, log=TRUE)
          
      u <- runif(1)
    
      if(is.na(exp(num-den)))  browser()
        
      # accept/reject step
      # if move is accepted then set beta to the proposed beta and update the log-likelihood value
      if (u < exp(num-den)) 
      {
        #print("accepted")
        beta_vec <- beta_st
        llikhood <- newllikhood
      }
    }
  } 
  
  sigma <- 3 # placeholder
  
  # the current beta, sigma and log-likelihood values
  list("beta_vec" = beta_vec, "sigma" = sigma, "llikhood"= llikhood)
  
}

#---------------------------------------------------------------------------------------------------------------------------------------
# RJ step
# this performs the between models update
updatemodel_logit <- function(beta_vec, sigma, llikhood, prior_beta = "norm", prior_beta_par, prop_beta = "norm", prop_beta_par)
{
  
  beta_st <- beta_vec
  
  # which coefficient to update (exclude the intercept)
  r <- sample((1:length(beta_vec))[-1], 1)
  
  # if at the current iteration the rth coefficient is equal to 0 propose a new value for it
  if(beta_vec[r] == 0)
  {
    beta_st[r] <- get(paste("r", prop_beta, sep=""))(1, prop_beta_par[1], prop_beta_par[2]) 
    
    #calculate acceptance probability - this is just the ratio of the log-likelihoods here
    # because we have chosen proposal distributions equal to the prior distributions
    # all prior model probabilities are equal so they cancel out
    
    # value of the log-likelihood for the proposed model
    newllikhood <- sum(dbinom(y,  1, plogis(X%*%beta_st), log = T))
    
    num <- newllikhood + log(get(paste("d", prior_beta, sep=""))(beta_st[r], prior_beta_par[1], prior_beta_par[2]))
    
    den <- llikhood + log(get(paste("d", prop_beta, sep=""))(beta_st[r], prop_beta_par[1], prop_beta_par[2])) 
    
  }
  # if at the current iteration the rth coefficient is not equal to 0 propose to set it equal to 0
  else
  {
    beta_st[r] <- 0  
    
    #calculate acceptance probability
    
    # value of the log-likelihood for the proposed model
    newllikhood <- sum(dbinom(y, 1,  plogis(X%*%beta_st), log = T))
    
    num <- newllikhood + log(get(paste("d", prop_beta, sep=""))(beta_vec[r], prop_beta_par[1], prop_beta_par[2]))
    
    den <- llikhood + log(get(paste("d", prior_beta, sep=""))(beta_vec[r], prior_beta_par[1], prior_beta_par[2]))
  }
  
  A <- min(1, exp(num-den))
  
  u <- runif(1)
  
  # accept/reject step
  # if move is accepted then beta parameter to the proposed ones
  # and update the log-likelihood value
  
  if(u <= A)
  {
    llikhood <- newllikhood
    beta_vec <- beta_st
  } 
  
  list("beta_vec" = beta_vec, "llikhood"= llikhood)
} 
  