### Converse/Gardner Lab Meeting
### 18 May 2020
### Model Averaging Exercise (Dormann et al. 2018)
#This code does jackknife model averageing for two simulated data sets, plots fits vs truth, and calculates RMSE.

#  load libraries
library(MuMIn)
library(here)


# Function to conduct jackknife model averaging. Return RMSE, predictions, and weights
jackknife_MA_func<-function(data,fam){
# Fit the global model
fm1 <- glm(y~., data = data, na.action = "na.fail", family = fam)

# Use dredge to fit all possible models
dfm1.jack <- dredge(fm1, rank="AIC") 

# Get all of the model fit results
fm1mods <- get.models(dfm1.jack, subset = TRUE)

# Jackknife 
# The jackknife model averaging optimises the fit of the prediction onto an omitted data point. 
# It is in a way similar to stacking, but requires only N steps (N = number of data points).
# Code adapted from Dormann et al. 2018

# Fit the candidate models, omitting one data point at a time.
N <- nrow(data)
M <- length(fm1mods)
train <- data # we will use the entire dataset here, rather than a subset
J <- matrix(NA, N, M) # matrix with jackknifed predictions
for (i in 1:N){
  # re-fit models on train with one less data point:
  jfits <- lapply(fm1mods, update, .~. , data=train[-i,],family=fam) 
  # predict them to omitted data point:
  J[i,] <- sapply(jfits, predict, newdata=train[i, , drop=F],family=fam, type = "response") 
  rm(jfits)
}

# Compute RMSE for a value of w, given J.
weightsopt <- function(ww, J){ # ww is a vector of model weights on transformed scale, J is our matrix of predicted points
  # function to compute RMSE on test data
  # at some point to also use likelihood instead of RMSE, but primarily for 0/1 data
  w <- c(1, exp(ww)); w <- w/sum(w) # sums to weight
  Jpred <- J %*% w # weighting predictions from all 256 models so we get 1 prediction for each point (100)
  return(sqrt(mean((Jpred - train$y)^2))) # RMSE calculation
}
ops <- optim(par=runif(NCOL(J)-1), weightsopt, method="BFGS", control=list(maxit=5000), J=J) # finds vector of weights which makes best prediction, closet to actual data
if (ops$convergence != 0) stop("Not converged!") 
round(weightsJMA <- c(1, exp(ops$par))/sum(c(1, exp(ops$par))),3) # using "par" from ops, optimizing vector of weights (transformed)
sum(weightsJMA) # check to see if sums to 1
weightedPredsJMA <- J %*% weightsJMA # optimized predictions (same as line 46 in function)
(RMSEjma <- sqrt(mean((weightedPredsJMA - train$y)^2))) # RMSE with best vector weights
names(RMSEjma)<-"RMSE_jma" # naming object
print(RMSEjma) # root mean square error, measuring how close predictions to the actual data

return(list(RMSEjma=RMSEjma, 
            weightedPredsJMA=weightedPredsJMA,
            weightsJMA=weightsJMA,
            LOO_preds=J))
}

#function to plot jackknife predictions vs truth
plot_jackknifeM<-function(truth,preds){
  # Plot predictions versus true data
  plot(x = preds, y = truth, xlab = "Predicted", ylab = "Observations")
  abline(0,1)
}


# Read-in data
normdat<-read.csv(here("normalSimulatedData1April2020.csv"))[,2:10]
logitdat<-read.csv(here("logitSimulatedData1April2020.csv"))[,2:10]

#fit normal
norm_jma<-jackknife_MA_func(data=normdat,fam="gaussian")
plot_jackknifeM(truth=normdat$y,preds=norm_jma$weightedPredsJMA)#plot fits 
hist(norm_jma$weightsJMA) #plot weights
cor_preds<-cov((norm_jma$LOO_preds)) # is this how we 
hist(cor_preds[upper.tri(cor_preds)]) #covariance between leave one out predictions

#fit logit
logit_jma<-jackknife_MA_func(data=logitdat,fam="binomial")
hist(logit_jma$weightsJMA) #plot weights
