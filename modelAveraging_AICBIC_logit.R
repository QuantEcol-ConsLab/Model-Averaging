## Converse Gardner lab meeting - April 6 2020
##
## Script generates data (logistic regression)
##
## compares models using AIC and BIC
##
##


#rm(list=ls())
set.seed(2020)

### Define stuff (very technical term) to create dataset

## samples
N <- 100	

## Covariates and effects
beta0 <- qlogis(.6)	# intercept
Beta <- c(0.6, -0.5, -0.4, 0.3, -0.2, 0.1, -0.05, 0.0)	# covariate effects
nCovs <- length(Beta)	# number of covariates


## covariates 
X <- matrix( rnorm(N*nCovs, 0, 1), nc = nCovs, nr=N, byrow=T)
X <- scale(X)[1:N, 1:nCovs]
colnames(X) <- paste("V",1:nCovs, sep="")

## realizations
p <- plogis( as.vector( beta0 + Beta%*%t(X) ) ) 	# expected
y <- rbinom(N, 1, p) 			# realized

## visualizing
hist(p, breaks=seq(0,1,.05))
hist(qlogis(p))
table(y)
pairs(cbind(y, X))

## The steps above generated the dataset for maze-solving nudibranchs


####################################################################
##
## Objective 1. AIC
## Compare all possible models using AIC.
## Generate predicted probablities for observed data
##
####################################################################

#### Easy way using packages ####

library(MuMIn)
x<-data.frame(cbind(y, X))

## fit global model
fm1 <- glm(y~., data=x, na.action = "na.fail", family="binomial")

## find all possible models
dd <- dredge(fm1, rank = "AIC")

## keep all or subset models to those with most support
ddMods <- get.models(dd, subset = TRUE) 
ddModsTop <- get.models(dd, subset = delta < 2)


#### Predict back to original dataset ####

## Predicted Pr(success) -  predict responses then average
avgm <- model.avg(ddMods)
meanPredY <- (predict(avgm , as.data.frame(X), se.fit = TRUE, type="response"))

 ## Can manually verify 
 # predY.avg<- (sapply(ddMods, predict, newdata = as.data.frame(X), family="binomial", type="response")) 
 # meanPredY.avg <- (apply(predY.avg, 1, function(x) weighted.mean(x, Weights(dd)[1:length(ddMods)])))
 ## Or (as in Cade 2015)
 # predY.avg%*%Weights(dd)[1:length(ddMods)]

## OR Predicted Pr(success) -  average logit scale parameters then predict responses 
lpredAvgm <- (predict(avgm , as.data.frame(X), se.fit = TRUE, type="link"))
predAvgm <- plogis(lpredAvgm$fit)

 ## can manually verify (careful of indexing)
 # pred.avg <- plogis( as.vector( avgm $coeff[1,1] + avgm $coeff[1,2:(K+1)][order(names(avgm $coeff[1,2:(K+1)])) ]%*%t(X) ) )
 # range(predAvgm-pred.avg)


#### Plot predictions relative to truth (they're not quite the same) ####
plot(meanPredY$fit, p, pch=16, xlab="Predicted")
abline(0,1)
points(predAvgm  , p, col="red")

## differences between prediction approaches
hist(meanPredY$fit-predAvgm )
range(meanPredY$fit-predAvgm)


#### Alternative to dredge function ####
## Use to calculate loglikelihood to feed into AIC calculation

pMLE <- fm1$ fitted.values # fitted values on probability scale
LL<-sum(y*log(pMLE) + (1-y)*log(1-pMLE))
## note: previous line just writes out the binomial pmf... sum(dbinom(y, 1, pMLE, log=TRUE))
## Can check using logLik(fm1)

## -2loglikelihood + 2*nparameters
AICfm1 <- -2*LL + 2*length(fm1$coefficients)
## Can check using AIC(fm1)

## More of a glimpse under the hood
get_logLik <- function(s_model, family = binomial(logit)) {
        n <- length(s_model$y)
        wt <- rep(1, n)
        deviance <- sum(family$dev.resids(s_model$y, s_model$fitted.values, wt))
        mod_rank <- sum(!is.na(s_model$coefficients))
        
        aic <- family$aic(s_model$y, rep(1, n), s_model$fitted.values, wt, deviance) + 2 * mod_rank
        log_lik <- mod_rank - aic/2
        return(log_lik)
}

fm1_model <- fm1[c("theta", "coefficients", "linear.predictors", "fitted.values", "y", "P", "residuals", "cov", "converged")]
LL <- get_logLik(fm1_model)
## -2loglikelihood + 2*K
altAICfm1 <- -2*LL + 2*K
## Can check using AIC(fm1)


####################################################################
##
## Objective 2. BIC (same as above, but with BIC)
## Compare all possible models using BIC.
## Generate predicted probablities for observed data
##
####################################################################

#### Easy way using packages ####

library(MuMIn)
x<-data.frame(cbind(y, X))

## fit global model
fm1 <- glm(y~., data=x, na.action = "na.fail", family="binomial")

## fun all possible models
bd <- dredge(fm1, rank = "BIC")

## keep all or subset models to those with most support
bdMods <- get.models(bd, subset = TRUE) 
bdModsTop <- get.models(bd, subset = delta < 2)


#### Predict back to original dataset ####

## Predicted Pr(success) -  predict responses then average
avgmb <- model.avg(bdMods)
meanPredYb <- (predict(avgmb , as.data.frame(X), se.fit = TRUE, type="response"))

## OR Predicted Pr(success) -  average logit scale parameters then predict responses 
lpredAvgmb <- predict(avgmb , as.data.frame(X), se.fit = TRUE, type="link")
predAvgmb <- plogis(lpredAvgmb$fit)


#### Plot predictions relative to truth (they're not quite the same) ####
plot(meanPredYb$fit, p, pch=16, xlab="Predicted")
abline(0,1)
points(predAvgmb , p, col="red")

## differences between prediction approaches
hist(meanPredYb$fit-predAvgmb )
range(meanPredYb$fit-predAvgmb )

#### Alternative to dredge function ####
## Use to calculate loglikelihood to feed into BIC calculation

pMLE <- fm1$ fitted.values # fitted values on probability scale
LL<-sum(y*log(pMLE) + (1-y)*log(1-pMLE)) # log(bernoulli pmf)
## Can check using logLik(fm1)

## -2loglikelihood + nparameters*log(number observations in model)
BICfm1 <- -2*LL + length(fm1$coefficients)*log(nrow(x))
## Can check using BIC(fm1)



##
##
## AIC vs BIC predictions
##
##

plot(meanPredY$fit, p, pch=16, xlab="Predicted")
abline(0,1)
points(meanPredYb$fit, p, col="red")

head(dd)
head(bd)

## END