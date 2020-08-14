### Converse/Gardner Lab Meeting
### 17 August 2020
### Model Averaging Exercise (Dormann et al. 2018)

#  load libraries
library(MuMIn)
library(BayesianTools)
library(BayesVarSel)

# Normal data
normdat <- read.csv("normalSimulatedData1April2020.csv")[,-c(1, 7:10)] 

# Fit the global model
fm1 <- glm(y~., data = normdat, na.action = "na.fail", family = "gaussian")

# Use dredge to fit all possible models
dfm1 <- dredge(fm1, rank="AIC") 
# one model with most of the weight (mod 16)

# Get all of the model fit results
fm1mods <- get.models(dfm1, subset = TRUE)

# example using demo(BayesVarSel.Hald)
# demo(BayesVarSel.Hald)
# data(Hald)

#Bvs() looks at model probabilities and variable inclusion
mod <- Bvs(y ~ V1+V2+V3+V4, data = normdat)
summary(mod)
plot(mod, option = 'conditional')

#Btest() computes bayes factors
#amelia - maybe figure out how to get a bigger list of models - all the combinations, and maybe expand to all covariates if we want...?
mods <- list("y ~ V1+V2+V3+V4", "y ~ V1+V2+V3", 'y ~ V1+V2')
Btest(
  models = mods,
  data = normdat,
  prior.betas = "Robust",
  prior.models = "Constant",
  priorprobs = NULL,
  null.model = NULL
)

# example from Dormann code - didn't fully execute
# Bayes Factor-based weights
# first set up an indicator matrix for all models:
# varIndMatrix <- ifelse(is.na(dfm1[, 1:5]), 0, 1)
# varIndMatrix <- varIndMatrix[order(as.numeric(rownames(varIndMatrix))),] 
# M <- length(fm1mods)
# N <- nrow(normdat)
# X <- cbind(1, normdat[1:N, -1]) # only the predictors, plus intercept as first column
# Y <- normdat[1:N, 1]
# 
# likelihood <- function(x, option = 1){
#   # a function to switch between the four models and compute the respective likelihood
#   # x is a vector with betas, plus a parameter for sd, i.e. NCOL(X)+1
#   res = as.matrix(X) %*% (as.numeric(x[-(NCOL(X)+1)]) * varIndMatrix[option, ]) 
#   # sets all parameters to 0 that are not in the model
#   ll = sum(dnorm(res - Y), sd = x[NCOL(X)+1], log = T)
#   return(ll)
# }
# 
# prior <- function(x){
#   # double-exponential/Laplace prior, with lambda set to 10
#   # this is a lasso-shrinkage prior 
#   ll = sum(dexp(abs(10*x)), log = T)
#   return(ll)
# }
# 
# # For illustration only: setting up two models to be fitted (for
# # details see package BayesianTools):
# #setup1 = createBayesianSetup(likelihood = function(x) likelihood(x, option  = 1), 
# #  prior=createPrior(density = prior, lower = c(rep(-5, NCOL(X)), 0.0001), upper = 
# #  c(rep(5, NCOL(X)),5)))
# #setup2 = createBayesianSetup(likelihood = function(x) likelihood(x, option  = 2), 
# #  prior=createPrior(density = prior, lower = c(rep(-5, NCOL(X)), 0.0001), upper = 
# #  c(rep(5, NCOL(X)),5)))
# # For illustration only: fitting the two models:
# #res1 <- runMCMC(setup1, sampler = "Metropolis", settings = list(iterations = 30000, 
# #                                          startValue = c(rep(0, NCOL(X)), 0.01)))
# #res2 <- runMCMC(setup2, sampler = "Metropolis", settings = list(iterations = 30000, 
# #                                          startValue = c(rep(0, NCOL(X)),0.01)))
# 
# # Instead, we loop the setup, and lapply the run:
# setups <- list()
# for (m in 1:M){
#   setups[[m]] <- createBayesianSetup(likelihood = function(x) likelihood(x, option  = m), 
#                                      prior=createPrior(density = prior, lower = 
#                                                          c(rep(-5, NCOL(X)), 0.0001), upper = c(rep(5, NCOL(X)),5)))
# }
# 
# system.time({
#   resBayesFits <- lapply(setups, runMCMC, sampler = "Metropolis", settings = list(iterations = 
#                                                                                     40000, startValue = c(rep(0, NCOL(X)), 0.01)))
# })
# 
# # extract the marginal likelihoods:
# ML <- unlist(sapply(resBayesFits, function(x) marginalLikelihood(x)[1]))
# names(ML) <- paste0("m", 1:16)
# 
# # compute Bayes Factor weights:
# (weightsBF <- exp(ML) / sum(exp(ML)))
# 
# # Although we may in a typical setup use the MCMC-samples to derive an averaged model prediction, 
# # we here use Bayes factors solely to derive Bayesian model weights, and compute predictions from them for 
# # comparison with the other methods. The code above is a good starting place for computing any other posterior parameter of interest.
# 
# 
# #weightsBF
# weightedPredsBF <- preds %*% weightsBF
# (RMSEBF <- sqrt(mean((weightedPredsBF - truth)^2)))


# Note that rmMCMC and Bayes factor approximate the thing: the probability of model $i$. 
# That these two approaches do not converge here to the same values is a consequence of different priors 
# and different ways to compute model weights. It is beyond the scope of this study to provide a detailed 
# discussion of how to implement different Bayesian fitting procedures and prior-choices.
