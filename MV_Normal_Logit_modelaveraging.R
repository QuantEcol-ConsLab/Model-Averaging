setwd("/Users/sipeh/Desktop/Lab Meetings/Model averaging")

#load in the data that Nathan and Staci created
normdat<-read.csv("normalSimulatedData1April2020.csv")[,2:10]

#####################
#'*Normal Dataset using MuMIn package*
#Minimal Variance/bates-granger weights using the MuMIn package
fm1.norm<-glm(y~., data=normdat, na.action="na.fail", family="gaussian")
#use dredge
dfm1.norm<-dredge(fm1.norm)
#all models
fm1mods.norm<-get.models(dfm1.norm, subset = TRUE)
#Bates-granger/minimal variance weights
BGweights<-(BGWeights(fm1mods.norm, data=normdat))
#predict for all models
preds_ <-sapply(fm1mods.norm, predict, newdata = normdat)
#weigh the models by the BG/MV weights
wtpreds_BG<-preds_%*%BGweights
#RMSE 
rmse_norm.1<-sqrt(mean((wtpreds_BG-normdat$y)^2))
rmse_norm.1

#plot with the observations
plot(wtpreds_BG, normdat$y)
abline(0,1)

#########normal data MV/BG weights by hand#########
#'*Normal Dataset WITHOUT using MuMIn package*

set.seed(21)
#split the data into a testing and training set
tsplit<-sample(rep(c(T,F), floor((nobs(fm1.norm))/2)))
t1<-normdat[tsplit,]
t2<-normdat[!tsplit,]
#fit models using the training data
t1fits<-lapply(fm1mods.norm, update, .~., data=t1)
#predict using the testing data
t2preds<-sapply(t1fits, predict, newdata=t2)
#get the residuals from the predictions
tresid<-t2$y-t2preds

#get the covariance matrix of the residuals
sigma<-cov(tresid)
#create a vector of ones
one<-rep(1, ncol(sigma))

library(MASS)
#Equation 14 in Dormann et al., 2018
#calculates the weights for each model based on the covariance matrix
bgweights<-ginv(t(one)%*%ginv(sigma)%*%one)%*%one%*%ginv(sigma)
#predict responses for all models
preds_ <-sapply(fm1mods.norm, predict, newdata = normdat)
#model averaging
bgwtpred<-preds_%*%t(bgweights)

plot(bgwtpred, normdat$y)
rmse_norm.2<-sqrt(mean((bgwtpred-normdat$y)^2))
rmse_norm.1 #using the MuMIn package
rmse_norm.2 #without the package
#similar results, but depends on how the data is split
# would multiple test/train sets averaged be better?


###############################################
###############################################
#'*Logit Dataset using MuMIn package*
logitdat<-read.csv("logitSimulatedData1April2020.csv")[,2:10]
#Minimal Variance/bates-granger weights using the MuMIn package
fm1.logit<-glm(y~., data=logitdat, na.action="na.fail", family="binomial")
#use dredge to fit all possible models
dfm1.logit<-dredge(fm1.logit)
#get all of the model fit results
fm1mods.logit<-get.models(dfm1.logit, subset = TRUE)
#use MuMIn to get BG/MV weights 
BGweights_logit<-Weights(BGWeights(fm1mods.logit, data=logitdat))
#predict on the response scale
preds_logit <-sapply(fm1mods.logit, predict, newdata = logitdat, type="response") 
#model averaged results
wtpreds_BG_logit<-((preds_logit)%*%BGweights_logit)
#RMSE
rmse_logit.1<-sqrt(mean((wtpreds_BG_logit-logitdat$y)^2))
rmse_logit.1
#predict on the link scale
preds_logit2<-sapply(fm1mods.logit, predict, newdata = logitdat, type="link")
#model averaged weights
wtpreds_BG_logit2<-(plogis(preds_logit2)%*%BGweights_logit)
#RMSE
rmse_logit.11<-sqrt(mean((wtpreds_BG_logit2-logitdat$y)^2))
rmse_logit.11


###MV by hand
set.seed(5)
#split the data into training and testing
tsplit<-sample(rep(c(T,F), floor(100/2)))
t1<-logitdat[tsplit,]
t2<-logitdat[!tsplit,]
#get all models for the training data
t1fits<-lapply(fm1mods.logit, update, .~., data=t1)
#predict to the test data
t2preds<-sapply(t1fits, predict, newdata=t2)
#get the residuals
tresid<-t2$y-t2preds
#get the covariance matrix for the residuals
sigma<-cov(tresid)
#column of 1s
one<-rep(1, ncol(sigma))

library(MASS)
#Equation 14 Dormann et al., 2018
#model weights
bgweights<-ginv(t(one)%*%ginv(sigma)%*%one)%*%one%*%ginv(sigma)
#model averaging
#as predicted on the response scale
bgwtpred<-(preds_logit)%*%t(bgweights)
#on the link scale
bgwtpred2<-plogis(preds_logit2)%*%t(bgweights)
#RMSE for both
rmse_bg.1<-sqrt(mean((bgwtpred-logitdat$y)^2))
rmse_bg.1
rmse_bg.2<-sqrt(mean((bgwtpred2-logitdat$y)^2))
rmse_bg.2



