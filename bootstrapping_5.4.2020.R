#### Converse-Gardner Joint Lab Group Meeting ####
### May 4, 2020 ###
#### May the Fourth Be With You ####

#### BOOTSTRAPPING APPROACH TO COMPUTING MODEL WEIGHTS ####

#load libraries
library(MuMIn)

#set working directory
setwd("/Users/owner/Desktop")

################

## Let's start with the normal distribution. ##

#read in the example data set
normdat <- read.csv("normalSimulatedData1April2020.csv")[ , 2:10]
head(normdat)

#fit the global model
m.all <- glm(y~., data = normdat, na.action = "na.fail", family = "gaussian")

#find all possible models
mytab <- dredge(m.all, rank = AIC)
mytab

#sort the model list from from simplest to full model:
model.list <- get.models(mytab, subset = NA)
model.list <- model.list[order(as.numeric(rownames(mytab)))]
M <- length(model.list) #note that the intercept may or may not be included in each of these models

#set up the bootstrap

N <- 50
train <- normdat[1:50,]
test <- normdat[(N+1):(2*N),]

Nboots <- 100 #ideally this should be at least 10,000
bestCounter <- rep(0, M)

for (i in 1:Nboots){
  #re-fit models for each bootstrap using the training data
  bsfits <- lapply(model.list, update, .~. , data = train[sample(nrow(train), nrow(train), 
                                                                 replace = T), ])
  bsAICs <- sapply(bsfits, AIC)
  bestCounter[which.min(bsAICs)] <- bestCounter[which.min(bsAICs)] + 1
}
weightsboot <- bestCounter/sum(bestCounter)

#compute predictions
preds <- sapply(model.list, predict, newdata = test)
dim(preds)

#compute weighted predction
weightedPredsBoot <- preds %*% weightsboot
RMSEboot <- sqrt(mean((weightedPredsBoot - test$y)^2))
RMSEboot #not bad, but check out what the authors wrote....

################

## Now let's turn our attention to the binomial data. ##

#read in the example data set
logdat <- read.csv("logitSimulatedData1April2020.csv")[ , 2:10]
head(logdat)

#fit the global model
m.all.log <- glm(y~., data = logdat, na.action = "na.fail", family = "binomial")

#find all possible models
mytab.log <- dredge(m.all.log, rank = AIC)
mytab.log

#sort the model list from from simplest to full model:
model.list.log <- get.models(mytab.log, subset = NA)
model.list.log <- model.list.log[order(as.numeric(rownames(mytab.log)))]
M.log <- length(model.list.log) #note that the intercept may or may not be included in each of these models

#set up the bootstrap

N.log <- 50
train.log <- logdat[1:N.log,]
test.log <- logdat[(N.log+1):(2*N.log),]

Nboots.log <- 100 #ideally this should be at least 10,000
bestCounter.log <- rep(0, M.log)

for (i in Nboots.log){
  #re-fit models for each bootstrap using the training data
  bsfits <- lapply(model.list.log, update, .~. , data = train.log[sample(nrow(train.log), nrow(train.log), 
                                                                 replace = T), ])
  bsAICs <- sapply(bsfits, AIC)
  bestCounter.log[which.min(bsAICs)] <- bestCounter.log[which.min(bsAICs)] + 1
}
weightsboot.log <- bestCounter.log/sum(bestCounter.log)

#compute predictions
preds.log <- sapply(model.list.log, predict, newdata = test.log)
dim(preds.log)

#compute weighted predction
weightedPredsBoot.log <- preds.log %*% weightsboot.log
RMSEboot.log <- sqrt(mean((weightedPredsBoot.log - test.log$y)^2))
RMSEboot.log

