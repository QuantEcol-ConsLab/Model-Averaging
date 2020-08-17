#### Converse/Gardner Lab Meeting
## 17 August 2020
## Model Averaging Exercise (Dormann et al. 2018)
## Amelia DuVall and Amanda Warlick

#  load libraries
library(MuMIn)
library(BayesianTools)
library(BayesVarSel)

#### Normal Data ####
# Load normal data
normdat <- read.csv("normalSimulatedData1April2020.csv")[,-c(1, 7:10)] ## select 4 predictors

# Split into train and test data
train.norm <- normdat[1:50,]
test.norm <- normdat[51:100,]

# Fit the global model
fm1.norm <- glm(y~., data = normdat, na.action = "na.fail", family = "gaussian")

# Use dredge to fit all possible models
dfm1.norm <- dredge(fm1.norm, rank="AIC") 
## one model with most of the weight

# Get all of the model fit results
model.list <- get.models(dfm1.norm, subset = TRUE) ## do all of these models have an intercept?

# Use BayesVarSel package to get Bayes Factors

## Manual method for fewer variables: 
# Bvs() looks at model probabilities and variable inclusion
mod.norm <- Bvs(y ~ V1 + V2 + V3 + V4, data = normdat)
summary(mod.norm)
plot(mod.norm, option = 'conditional') 

# automatic testing of all possible models
# Create grid of all model combinations
predictors <- c("V1", "V2", "V3", "V4")
regMat <- expand.grid(c(TRUE,FALSE), c(TRUE,FALSE),
                      c(TRUE,FALSE), c(TRUE,FALSE)) 
names(regMat) <- paste("V", 1:4, sep="")
allModelsList <- apply(regMat, 1, function(x) as.formula(
  paste(c("y ~ 1", predictors[x]),
        collapse=" + ")) )

# Use Btest() to compute Bayes Factors
weights <- Btest(
  models = allModelsList,
  data = normdat,
  prior.betas = "Robust", # other options, how does this influence results?
  prior.models = "Constant", # other options, how does this influence results?
  priorprobs = NULL,
  null.model = NULL)
## Need to be mindful of priors in model
## Test other priors to see if it changes results 

# Extract Bayes Factors and model weights 
BFs.norm <- weights[[1]]
BFweights.norm <- weights[[2]]

# Calculate weighted predictions
preds <- sapply(model.list, predict, newdata = test.norm)
weightedPredsBF.norm <- preds %*% BFweights.norm

# Calculate RMSE
(RMSEBF <- sqrt(mean((weightedPredsBF.norm - test.norm$y)^2)))

#### END Normal Data ####

#### Logit Data ####
# Load logit data
logitdat <- read.csv("logitSimulatedData1April2020.csv")[,-c(1, 7:10)] ## select 4 predictors

# Split into train and test data
train.logit <- logitdat[1:50,] ## we don't use this train dataset?
test.logit <- logitdat[51:100,]

# Fit the global model
fm1.logit <- glm(y~., data = logitdat, na.action = "na.fail", family = "binomial")

# Use dredge to fit all possible models
dfm1.logit <- dredge(fm1.logit, rank="AIC") 
## one model with most of the weight

# Get all of the model fit results
model.list <- get.models(dfm1.logit, subset = TRUE) ## do all of these models have an intercept?

# Use BayesVarSel package to get Bayes Factors
## Bvs() looks at model probabilities and variable inclusion
mod.logit <- Bvs(y ~ V1 + V2 + V3 + V4, data = logitdat)
summary(mod.logit)
plot(mod.logit)

# Create grid of all model combinations
predictors <- c("V1", "V2", "V3", "V4")
regMat <- expand.grid(c(TRUE,FALSE), c(TRUE,FALSE),
                      c(TRUE,FALSE), c(TRUE,FALSE))
names(regMat) <- paste("V", 1:4, sep="")
allModelsList <- apply(regMat, 1, function(x) as.formula(
  paste(c("y ~ 1", predictors[x]),
        collapse=" + ")) )

# Use Btest() to compute Bayes Factors
weights <- Btest(
  models = allModelsList,
  data = logitdat,
  prior.betas = "Robust", # other options, how does this influence results?
  prior.models = "Constant", # other options, how does this influence results?
  priorprobs = NULL,
  null.model = NULL)
## Need to be mindful of priors in model
## Test other priors to see if it changes results 

# Extract Bayes Factors and model weights 
BFs.logit <- weights[[1]]
BFweights.logit <- weights[[2]]

# Calculate weighted predictions
preds <- sapply(model.list, predict, newdata = test.logit)
weightedPredsBF.logit <- preds %*% BFweights.logit

# Calculate RMSE
(RMSEBF <- sqrt(mean((weightedPredsBF.logit - test.logit$y)^2)))

#### END Logit Data ####

#### Dormann et al. code ####
## Did not fully execute -- couldn't get the code to work.

# Bayes Factor-based weights
# first set up an indicator matrix for all models:
# varIndMatrix <- ifelse(is.na(mytab[, 1:5]), 0, 1)
# varIndMatrix <- varIndMatrix[order(as.numeric(rownames(varIndMatrix))),] 
# # sort sequence, same as for model.list!
# 
# X <- cbind(1, train[1:N, -1]) # only the predictors, plus intercept as first column
# Y <- train[1:N, 1]
# library(BayesianTools)
# likelihood <- function(x, option = 1){ ## what is option 1? x?
#   # a function to switch between the four models and compute the respective likelihood
#   # x is a vector with betas, plus a parameter for sd, i.e. NCOL(X)+1
#   res = as.matrix(X) %*% (x[-(NCOL(X)+1)] * varIndMatrix[option, ]) 
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
# system.time({
#   resBayesFits <- lapply(setups, runMCMC, sampler = "Metropolis", settings = list(iterations = 
#                                                                                     40000, startValue = c(rep(0, NCOL(X)), 0.01)))
# })
# # extract the marginal likelihoods:
# ML <- unlist(sapply(resBayesFits, function(x) marginalLikelihood(x)[1]))
# names(ML) <- paste0("m", 1:16)
# 
# # compute Bayes Factor weights
# (weightsBF <- exp(ML) / sum(exp(ML)))
# 
# # Although we may in a typical setup use the MCMC-samples to derive an averaged model prediction, 
# # we here use Bayes factors solely to derive Bayesian model weights, and compute predictions from 
# # them for comparison with the other methods. The code above is a good starting place for computing 
# # any other posterior parameter of interest.
# 
# #weightsBF
# weightedPredsBF <- preds %*% weightsBF
# (RMSEBF <- sqrt(mean((weightedPredsBF - truth)^2)))
# 
# # Note that rmMCMC and Bayes factor approximate the thing: the probability of model $i$. That these 
# # two approaches do not converge here to the same values is a consequence of different priors and 
# # different ways to compute model weights. It is beyond the scope of this study to provide a detailed 
# # discussion of how to implement different Bayesian fitting procedures and prior-choices.

#### END Dormann et al. code ####