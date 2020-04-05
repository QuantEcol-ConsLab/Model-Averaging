## Converse Gardner lab meeting - April 6 2020
##
## Script generates data (linear regression)
##
## compares models using AIC and BIC
##
##


#rm(list=ls())
set.seed(2020)

### Define stuff

## samples
N <- 100	

## Covariates and effects
beta0 <- 5	# intercept
Beta <- c(0.6, -0.5, -0.4, 0.3, -0.2, 0.1, -0.05, 0.0)	# covariate effects
sig <- 0.5
nCovs <- length(Beta)	# number of covariates

## covariates 
X <- matrix( rnorm(N*nCovs, 0, 1), nc=nCovs, nr=N, byrow=T)
X <- scale(X)[1:N, 1:nCovs]
colnames(X) <- paste("V",1:nCovs, sep="")

## realizations
mu <- ( as.vector( beta0 + Beta%*%t(X) ) ) 	# expected
y <- rnorm(N, mu, sd=sig) 			# realized

hist(y)
pairs(cbind(y, X))

## The steps above generated the dataset for size-varying nudibranchs


###########################################################################################
##
## Option 1. AIC
## Compare all possible models using AIC.
## Generate predicted probablities for observed data
##
###########################################################################################

#### Easy way using packages ####

library(MuMIn)
x<-data.frame(cbind(y, X))

## fit global model
fm1 <- glm(y~., data=x, na.action = "na.fail", family="gaussian")

## fun all possible models
dd <- dredge(fm1, rank = "AIC")

## keep all or subset models to those with most support
ddMods <- get.models(dd, subset = TRUE) 
ddModsTop <- get.models(dd, subset = delta < 2)


#### Predict back to original dataset ####

## Predicted Pr(success) -  predict responses then average
avgm <- model.avg(ddMods )
meanPredY <- (predict(avgm , as.data.frame(X), se.fit = TRUE, type="response"))

#### Alternative to dredge function ####
## Pull sigma from glm model and use to calculate loglikelihood to feed into AIC calculation
## sigma(model) needs to be corrected for bias

## bias correction from sigma to sigma_hat is sqrt( (N-nparameters)/N )
sigma.ML <- sigma(fm1) * sqrt((N-length(fm1$coefficients))/N)
LL <- sum(dnorm(y, mean=fitted(fm1), sd=sigma.ML, log=TRUE)) 
## Can check using logLik(fm1)

## -2loglikelihood + 2*nparameters
AICfm1 <- -2*LL + 2*(length(fm1$coefficients)+1)
## Can check using AIC(fm1)



###########################################################################################
##
## Option 2. BIC (same as above, but with BIC)
## Compare all possible models using BIC.
## Generate predicted probablities for observed data
##
###########################################################################################

#### Easy way using packages ####

library(MuMIn)
x<-data.frame(cbind(y, X))

## fit global model
fm1 <- glm(y~., data=x, na.action = "na.fail", family="gaussian")

## fun all possible models
bd <- dredge(fm1, rank = "BIC")

## keep all or subset models to those with most support
bdMods <- get.models(bd, subset = TRUE) 
bdModsTop <- get.models(bd, subset = delta < 2)


#### Predict back to original dataset ####

## Predicted Pr(success) -  predict responses then average
avgmb <- model.avg(bdMods )
meanPredYb <- (predict(avgmb , as.data.frame(X), se.fit = TRUE, type="response"))

#### Alternative to dredge function ####
## Pull sigma from glm model and use to calculate loglikelihood to feed into BIC calculation
## sigma(model) needs to be corrected for bias

## bias correction from sigma to sigma_hat is sqrt( (N-nparameters)/N )
sigmaMLE <- sigma(fm1) * sqrt((N-length(fm1$coefficients))/N)
LL <- sum(dnorm(y, mean=fitted(fm1), sd=sigmaMLE, log=TRUE))
## Can check using logLik(fm1)

## -2loglikelihood + nparameters*log(number observations in model)
BICfm1 <- -2*LL + (length(fm1$coefficients)+1)*log(nrow(x))
## Can check using BIC(fm1)


##
##
## AIC vs BIC predictions
##
##

plot(meanPredY$fit, y, pch=16, xlab="Predicted")
abline(0,1)
points(meanPredYb$fit, y, col="red")

head(dd)
head(bd)

## END