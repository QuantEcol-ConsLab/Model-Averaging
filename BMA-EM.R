# This script estimates Bayesian Model Averaging - Expectation Maximazation model weights, and is build off the code provided in appendix S9 of Dormall et al. 2017.

#packages
library(EBMAforecast)
library(MuMIn)
library(tidyverse)

#load data
dat.norm <-  read.csv("normalSimulatedData1April2020.csv") %>% #normal
  select(2:6)                                                  #trim columns to response and first 4 predictors
dat.logit <- read.csv("logitSimulatedData1April2020.csv")   %>%#logit
  select(2:6)                                                  #trim columns to response and first 4 predictors

#begin function for BMA-EM
BMA_EM_func<-function(dat                          # model data (response and predictors)
                      ,model){                     # model family ("normal" or "logit")

  N<-nrow(dat)                                     # number of data points
fam<-ifelse(model=="normal","gaussian","binomial") # model family for linear model 
model.list<-get.models(dredge(glm(y~., data = dat,
                                  na.action = "na.fail", family = fam))
                       ,subset=TRUE)               # get all possible models
preds <- sapply(model.list, predict,type="response")# predictions for all data points with all models

M<-length(model.list)                              # number of models

set.seed(1)                                        # set seed for consistant results
trainsplit <- sample(rep(c(T, F), N/2))            # random assignment of data rows to test and train
trainsub1 <- dat[trainsplit, ]                     # make trainsub1
trainsub2 <- dat[!trainsplit, ]                    # make trainsub2
trainsub1fits <- lapply(model.list, update, .~. ,
                        data=trainsub1,family=fam) # re-fit models on trainsub1
bmafits <- sapply(trainsub1fits, predict,
                  newdata=trainsub2, 
                  type="response")                 # predict them to trainsub2
bmaY <- trainsub2$y                                #to train the EMA-algorithm
EBMAdata <- makeForecastData(.predCalibration=bmafits,
                             .outcomeCalibration=bmaY,
                            .modelNames=paste0("m", 1:M)) # make dataset
EBMAfit <- calibrateEnsemble(EBMAdata, model=model) # compute weights

print(" ")                                          # empty space for printing
 
print(EBMAfit@modelWeights)                         # weights

weightedPredsEBMA <- preds %*% EBMAfit@modelWeights # weighted predictions
RMSEebma <- sqrt(mean((weightedPredsEBMA - dat$y)^2))#RMSE
print(paste0("RMSE= ",round(RMSEebma,3)))           # print RMSE
if(model=="normal"){                                # plot predicted vs observed for normal data
plot(dat.norm$y,weightedPredsEBMA)
abline(0,1,col="firebrick3")}
}                                                   #end of function

#call function for normal data
BMA_EM_func(dat=dat.norm,model="normal")
#call function for logit data
BMA_EM_func(dat=dat.logit,model="logit")
