#### Normal first then logit
library(MuMIn)
setwd("C:/Users/Robbie/Desktop/Model-Averaging")
dat.norm <-  read.csv("normalSimulatedData1April2020.csv")

#  Fit global model
fm1.norm <- glm(y~., data = dat.norm, na.action = "na.fail", family = "gaussian")

mytab <- dredge(fm1.norm, rank=AIC)
model.list <- get.models(mytab, subset=1:20)
# Reorder model list?

# Maybe get only top 20 models vel sim.

## Implement stacking function from Dormann et al. 2018 Supplementary Information
stacking <- function(test.preds, test.obs){
  # this function computes the optimal weight for a single train/test split;
  # from the models fitted to the training data it uses the predictions to the test;
  # then it optimises the weight vector across the models for combining these 
  # predictions to the observed data in the test;
  # trick 1: each weight is between 0 and 1: w <- exp(-w)
  # trick 2: weights sum to 1: w <- w/sum(w)
  #
  # weights are weights for each model, between -infty and +infty!
  # preds are predictions from each of the models
  
  if (NCOL(test.preds) >= length(test.obs)) stop("Increase the test set! More models 
	   than test points.")
  
  # now do an internal splitting into "folds" data sets:
  weightsopt <- function(ww){ 
    # function to compute RMSE on test data
    w <- c(1, exp(ww)); w <- w/sum(w) ## w all in (0,1) SIMON; 
    # set weight1 always to 1, other weights are scaled accordingly 
    # (this leads to a tiny dependence of optimal weights on whether model1 is any 
    # good or utter rubbish; 
    # see by moving the 1 to the end instead -> 3rd digit changes)
    pred <- as.vector(test.preds %*% w)
    return(sqrt(mean((pred - test.obs)^2)))
  }
  
  # This implements the equation at top right of p. 495 of Dormann et al. 2018
  ops <- optim(par=runif(NCOL(test.preds)-1), weightsopt, method="BFGS")
  if (ops$convergence != 0) stop("Optimisation not converged!")
  round(c(1, exp(ops$par))/sum(c(1, exp(ops$par))), 4)
}

Nstack <- 500 # This should be a high number like 1000
N <- nrow(dat.norm)
Ntrain <- N/2
datsplit <- sample(rep(c(T,F),floor(N/2)))
train <- dat.norm[datsplit,]
testdat <- dat.norm[!datsplit,]
M <- 20 # This should be < Ntrain/2
weightsStack <- matrix(NA, ncol=M, nrow=Nstack)
colnames(weightsStack) <- paste0("m", 1:M)
i = 0
a <- proc.time()[3]
while (i < Nstack){
  trainsplit <- sample(rep(c(T, F), floor(Ntrain/2)))
  trainsub1 <- train[trainsplit, ]
  trainsub2 <- train[!trainsplit, ]
  trainsub1fits <- lapply(model.list, update, .~. , data=trainsub1) # re-fit models on train1
  stackpreds <- sapply(trainsub1fits, predict, newdata=trainsub2) # predict them to train2
  optres <- try(stacking(test.preds=stackpreds, test.obs=trainsub2$y))
  if (inherits(optres, "try-error"))  next;
  i = i + 1
  weightsStack[i,] <- optres
  rm(trainsplit, trainsub1, trainsub2, trainsub1fits, stackpreds)
  #print(i)
}
b <- proc.time()[3]
b-a

## Calculate weights
(weightsStacking <- colSums(weightsStack)/sum(weightsStack))
preds <- sapply(model.list, predict, newdata=testdat)
truth <- testdat$y
weightedPredsStack <- preds %*% weightsStacking # Add back in preds
(RMSEstack <- sqrt(mean((weightedPredsStack - truth)^2))) # Add back in truth

## Compare to AIC weights in terms of RMSE and predictions
avgm.norm <- model.avg(model.list)
weightedPredsAIC <- predict(avgm.norm, testdat, se.fit=TRUE, type="response")
(RMSE_AIC <- sqrt(mean((weightedPredsAIC$fit-truth)^2)))

#### Logit data ####
dat.logit <-  read.csv("logitSimulatedData1April2020.csv")

#  Fit global model
fm1.logit <- glm(y~., data = dat.logit, na.action = "na.fail", family = "binomial")

mytab.logit <- dredge(fm1.logit, rank=AIC)
model.list.logit <- get.models(mytab.logit, subset=1:20)
# Reorder model list?

Nstack <- 500
N <- nrow(dat.logit)
Ntrain <- N/2
datsplit <- sample(rep(c(T,F),floor(N/2)))
train <- dat.logit[datsplit,]
testdat <- dat.logit[!datsplit,]
M <- 20 # This should be < Ntrain/2
weightsStack <- matrix(NA, ncol=M, nrow=Nstack)
colnames(weightsStack) <- paste0("m", 1:M)
i = 0
while (i < Nstack){
  trainsplit <- sample(rep(c(T, F), floor(Ntrain/2)))
  trainsub1 <- train[trainsplit, ]
  trainsub2 <- train[!trainsplit, ]
  trainsub1fits <- lapply(model.list.logit, update, .~. , data=trainsub1) # re-fit models on train1
  stackpreds <- sapply(trainsub1fits, predict, newdata=trainsub2, type="response") # predict them to train2
  optres <- try(stacking(test.preds=stackpreds, test.obs=trainsub2$y))
  if (inherits(optres, "try-error"))  next;
  i = i + 1
  weightsStack[i,] <- optres
  rm(trainsplit, trainsub1, trainsub2, trainsub1fits, stackpreds)
  #print(i)
}

## Calculate weights
(weightsStacking <- colSums(weightsStack)/sum(weightsStack))
preds <- sapply(model.list.logit, predict, newdata=testdat,type="response")
truth <- testdat$y
weightedPredsStack <- preds %*% weightsStacking # Add back in preds
(RMSEstack <- sqrt(mean((weightedPredsStack - truth)^2))) # Add back in truth

avgm.logit <- model.avg(model.list.logit)
weightedPredsAIC_logit <- predict(avgm.logit, testdat, se.fit=TRUE, type="response")
(RMSE_AIC_logit <- sqrt(mean((weightedPredsAIC_logit$fit-truth)^2)))