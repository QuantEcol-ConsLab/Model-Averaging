### Converse/Gardner Lab Meeting
### 18 May 2020
### Model Averaging Exercise (Dormann et al. 2018)
#This code does jackknife model averageing for two simulated data sets, plots fits vs truth, and calculates RMSE.

#  load libraries
library(MuMIn)
library(here)


#function to conduct jackknife model averaging. Return RMSE, predictions, and weights
jackknife_MA_func<-function(data){

## Normal data ##
# Fit the global model
fm1.norm <- glm(y~., data = normdat, na.action = "na.fail", family = "gaussian")

# Use dredge to fit all possible models
dfm1.norm.jack <- dredge(fm1.norm, rank="AIC") # do these need to be ranked?

# Get all of the model fit results
fm1mods.norm <- get.models(dfm1.norm.jack, subset = TRUE)

# Jackknife 
# The jackknife model averaging optimises the fit of the prediction onto an omitted data point. 
# It is in a way similar to stacking, but requires only N steps (N = number of data points).
# Code adapted from Dormann et al. 2018

# Fit the candidate models, omitting one data point at a time.
N <- nrow(normdat)
M <- length(fm1mods.norm)
train <- normdat # we will use the entire dataset here, rather than a subset
J <- matrix(NA, N, M) # matrix with jackknifed predictions
for (i in 1:N){
  # re-fit models on train with one less data point:
  jfits <- lapply(fm1mods.norm, update, .~. , data=train[-i,]) 
  # predict them to omitted data point:
  J[i,] <- sapply(jfits, predict, newdata=train[i, , drop=F]) 
  rm(jfits)
}

# Compute RMSE for a value of w, given J.
weightsopt <- function(ww, J){ 
  # function to compute RMSE on test data
  # at some point to also use likelihood instead of RMSE, but primarily for 0/1 data
  w <- c(1, exp(ww)); w <- w/sum(w) 
  Jpred <- J %*% w
  return(sqrt(mean((Jpred - train$y)^2)))
}
ops <- optim(par=runif(NCOL(J)-1), weightsopt, method="BFGS", control=list(maxit=5000), J=J)
if (ops$convergence != 0) stop("Not converged!")
round(weightsJMA <- c(1, exp(ops$par))/sum(c(1, exp(ops$par))),3)
sum(weightsJMA) # check to see if sums to 1
weightedPredsJMA <- J %*% weightsJMA
(RMSEjma <- sqrt(mean((weightedPredsJMA - train$y)^2)))
names(RMSEjma)<-"RMSE_jma"
print(RMSEjma)

return(list(RMSEjma=RMSEjma,
            weightedPredsJMA=weightedPredsJMA,
            weightsJMA=weightsJMA))
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
norm_jma<-jackknife_MA_func(normdat)
plot_jackknifeM(truth=normdat$y,preds=norm_jma$weightedPredsJMA)#plot fits 
hist(norm_jma$weightsJMA) #plot weights

#fit logit
logit_jma<-jackknife_MA_func(logitdat)
hist(logit_jma$weightsJMA) #plot weights
