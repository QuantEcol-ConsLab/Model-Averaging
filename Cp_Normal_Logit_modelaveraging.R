setwd("/Users/sipeh/Desktop/Lab Meetings/Model averaging")

#load in the data that Nathan and Staci created
normdat<-read.csv("normalSimulatedData1April2020.csv")[,2:10]

#'*First do Mallow's Cp using the MuMIn package for both the normal and logit data*
library(MuMIn)
#############
#normal data Cp model averaging using MuMIn package#
##############
#####can do it a couple of ways
#'*1st Way*#####

#Fit the global model
fm1.norm<-glm(y~., data=normdat, na.action="na.fail", family="gaussian")
#use dredge to fit all possible models
dfm1.norm.cp<-dredge(fm1.norm, rank="Cp")
#get all of the model fit results
fm1mods.norm<-get.models(dfm1.norm.cp, subset = TRUE)
#grab the model averaging weights based on Mallow's Cp
avgfm1.norm<-model.avg(fm1mods.norm)
###then predict back to data, using those weights
pred.norm.1<-predict(avgfm1.norm, normdat, se.fit=TRUE, type="response")
#Use RMSE to compare
rmse.norm_1.1<-sqrt(mean((pred.norm.1$fit-normdat$y)^2))

####plot the predictions from this approach vs. truth
plot(pred.norm.1$fit,normdat$y, xlab="predicted", ylab="observations")
abline(0,1)

#'*2nd Way*#####
#get all possible models
dfm1.norm<-dredge(fm1.norm)
#get all of the model output
fm1.mods.all.norm<-get.models(dfm1.norm, subset=TRUE)
#Get the Mallow's CP for all models
Cps.norm<-sapply(fm1.mods.all.norm, Cp)
# Calculate the weights for each model based on Cp
cp.weights.norm<-exp(-0.5*(Cps.norm-min(Cps.norm)))/sum(exp(-0.5*(Cps.norm-min(Cps.norm))))
# then predict 
pred.norm.2<-sapply(fm1.mods.all.norm, predict, newdata=normdat, type="response")
#Cp averaged predictions
wtpreds.cp.norm<-pred.norm.2%*%cp.weights.norm
#RMSE this way and compare with the other way
rmse.norm_1.2<-sqrt(mean((wtpreds.cp.norm-normdat$y)^2))
rmse.norm_1.1 #first way
rmse.norm_1.2 #this way
#they are the same

#########################
#'*Mallow's Cp by hand*
########################
#'*Normal data*

####write a function to calculate the cp for each dredge model
# then do the weights as above and similar 
sigma.ML.norm<-sigma(fm1.norm)*sqrt((100-length(fm1.norm$coefficients))/100)
LL.norm<-sum(dnorm(normdat$y, mean=fitted(fm1.norm), sd=(sigma.ML.norm), log=TRUE))
#check it
LL.norm
logLik(fm1.norm)
#what the equation is in Dormann et al., 2018
cpfm1.norm<- -2*LL.norm -100 + 2*(length(fm1.norm$coefficients)+1)
Cp(fm1.norm) #how MuMIn calculates it
cpfm1.norm
#doesnt match?


###how MuMIn Cp calculates it
Cpmumin<-deviance(fm1.norm)+2*(deviance(fm1.norm)/df.residual(fm1.norm))*(nobs(fm1.norm)-df.residual(fm1.norm))
Cpmumin
Cp(fm1.norm)

#then can calculate for all models
dfm1.norm<-dredge(fm1.norm)
#get all of the model output
fm1.mods.all.norm<-get.models(dfm1.norm, subset=TRUE)
deviance.norm<-sapply(fm1.mods.all.norm, deviance)
dfr.norm<-sapply(fm1.mods.all.norm, df.residual)  
Cps.norm.again<-deviance.norm+2*(deviance.norm/dfr.norm)*(100-dfr.norm)
#then calculate the weights based on cp
cp.weights.norm.again<-exp(-0.5*(Cps.norm.again-min(Cps.norm.again)))/sum(exp(-0.5*(Cps.norm.again-min(Cps.norm.again))))
# then predict 
norm.preds.cp.again<-sapply(fm1.mods.all.norm, predict, newdata=normdat, type="response")
#Cp averaged predictions
wtpreds.cp.norm.again<-norm.preds.cp.again%*%cp.weights.norm.again
#RMSE this way and compare with the other way
rmse.norm_1.3<-sqrt(mean((wtpreds.cp.norm.again-normdat$y)^2))
rmse.norm_1.1 #first way
rmse.norm_1.2 #second way
rmse.norm_1.3 #by hand


#############
#'*logit data Cp model averaging using MuMIn package*
##############
#load in the data that Nathan and Staci created
logitdat<-read.csv("logitSimulatedData1April2020.csv")[,2:10]

#'*1st Way*#####
#Fit the global model
fm1.logit<-glm(y~., data=logitdat, na.action="na.fail", family="binomial")
#use dredge to fit all possible models
dfm1.logit<-dredge(fm1.logit, rank="Cp")
#get all of the model fit results
fm1mods.logit<-get.models(dfm1.logit, subset = TRUE)
#grab the model averaging weights based on Mallow's Cp
avgfm1.logit<-model.avg(fm1mods.logit)
###then predict back to data, using those weights
#predicted responses then average
meanPredY.logit<-predict(avgfm1.logit, logitdat, se.fit=TRUE, type="response")
#Use RMSE to compare
rmse.logit_1.1<-sqrt(mean((meanPredY.logit$fit-logitdat$y)^2))

#average logit scale then predict 
lPredY.logit<-predict(avgfm1.logit, logitdat, se.fit=TRUE, type="link")
lp.logit<-plogis(lPredY.logit$fit)
rmse.logit_2.1<-sqrt(mean((lp.logit-logitdat$y)^2))

p<-read.csv("knownp.csv")[,3]
####plot the predictions from this approach vs. truth
plot(meanPredY.logit$fit,p, xlab="predicted", ylab="observations")
points(lp.logit, p, col="red", pch=16)
abline(0,1)

#'*2nd Way*#####
#get all possible models
dfm1.logit<-dredge(fm1.logit)
#get all of the model output
fm1.mods.all.logit<-get.models(dfm1.logit, subset=TRUE)
#Get the Mallow's CP for all models
Cps.logit<-sapply(fm1.mods.all.logit, Cp)
# Calculate the weights for each model based on Cp
cp.weights.logit<-exp(-0.5*(Cps.logit-min(Cps.logit)))/sum(exp(-0.5*(Cps.logit-min(Cps.logit))))
# then predict 
#predict responses then average
logit.preds.cp<-sapply(fm1.mods.all.logit, predict, newdata=logitdat, type="response")
#Cp averaged predictions
wtpreds.cp.logit<-logit.preds.cp%*%cp.weights.logit
#RMSE this way and compare with the other way
rmse.logit_1.2<-sqrt(mean((wtpreds.cp.logit-logitdat$y)^2))
rmse.logit_1.1 #first way
rmse.logit_1.2 #this way

#  average logit scale parameters then predict responses
llogit.preds.cp<-(sapply(fm1.mods.all.logit, predict, newdata=logitdat, type="link"))
#Cp averaged predictions
lwtpreds.cp.logit<-plogis(llogit.preds.cp%*%cp.weights.logit)
rmse.logit_2.2<-sqrt(mean((lwtpreds.cp.logit-logitdat$y)^2))
rmse.logit_2.1 #first way
rmse.logit_2.2 #this way
#both the same

#########
#'*logit data Cp by hand*
#'########
Cpmumin<-deviance(fm1.logit)+2*(nobs(fm1.logit)-df.residual(fm1.logit))
Cpmumin
Cp(fm1.logit)

#then can calculate for all models
dfm1.logit<-dredge(fm1.logit)
#get all of the model output
fm1.mods.all.logit<-get.models(dfm1.logit, subset=TRUE)
deviance.logit<-sapply(fm1.mods.all.logit, deviance)
dfr.logit<-sapply(fm1.mods.all.logit, df.residual)  
Cps.logit.again<-deviance.logit+2*(100-dfr.logit)
#then calculate the weights based on cp
cp.weights.logit.again<-exp(-0.5*(Cps.logit.again-min(Cps.logit.again)))/sum(exp(-0.5*(Cps.logit.again-min(Cps.logit.again))))
# then predict 
logit.preds.cp.again<-sapply(fm1.mods.all.logit, predict, newdata=logitdat, type="response")
#Cp averaged predictions
wtpreds.cp.logit.again<-logit.preds.cp.again%*%cp.weights.logit.again
#RMSE this way and compare with the other way
rmse.logit_1.3<-sqrt(mean((wtpreds.cp.logit.again-logitdat$y)^2))
rmse.logit_1.1 #first way
rmse.logit_1.2 #second way
rmse.logit_1.3 #by hand

##  average logit scale parameters then predict responses
ll.preds.cp<-(sapply(fm1.mods.all.logit, predict, newdata=logitdat, type="link"))
#Cp averaged predictions
llwtpreds.cp.logit<-plogis(ll.preds.cp%*%cp.weights.logit.again)
rmse.logit_2.3<-sqrt(mean((llwtpreds.cp.logit-logitdat$y)^2))
rmse.logit_2.1
rmse.logit_2.2
rmse.logit_2.3




