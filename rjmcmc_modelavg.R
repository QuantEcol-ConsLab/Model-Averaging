# Read in normal data first
setwd("C:/Users/Robbie/Desktop/Model-Averaging")
source("RJMCMCfunctions.R")
dat.norm <- read.csv("normalSimulatedData1April2020.csv")

dat.norm <- dat.norm[,-1]
set.seed(18)
trainsplit <- sample(1:nrow(dat.norm), nrow(dat.norm)/2, replace=FALSE)
train <- dat.norm[trainsplit,]
testdat <- dat.norm[-trainsplit,]

X <- as.matrix(cbind(1, train[,-1])) # include 1 for the intercept
y <- train$y

# starting values for coefficients
beta_vec <- runif(ncol(X), -1, 1)
# starting value for residual noise estimate sigma
sigma <- rgamma(1, 1, 1)

llikhood <- sum(dnorm(y, X%*%beta_vec, sigma, log = T))

niter <- 20000 # originally 100000 
burnin <- niter/2
thin <- 20

(niter-burnin)/thin
# store values
beta_mat <- matrix(0, nrow = (niter-burnin)/thin, ncol = length(beta_vec))
sigma_vec <- c()

#choose prior parameters for coefficients
prior_sigma_beta <- 3
prior_pars <- c(0, prior_sigma_beta)

#choose parameters of proposal distributions for updates
prop_sigma_beta <- 3
prop_pars <- c(0, prop_sigma_beta)
# now run rjMCMC:
system.time({
  for (iter in 1:niter){
    temp1 <- updateparam(beta_vec, sigma, llikhood, prior_beta = "norm", prior_beta_par = 
                           prior_pars, prop_beta = "norm", prop_beta_par = prop_pars[2])
    beta_vec <- temp1$beta_vec
    sigma <- temp1$sigma
    llikhood <- temp1$llikhood
    # between models update
    temp2 <- updatemodel(beta_vec, sigma, llikhood,  prior_beta = "norm", prior_beta_par = 
                           prior_pars, prop_beta = "norm", prop_beta_par = prop_pars)
    # read current values for beta, current model and log-likelihood value
    beta_vec <- temp2$beta_vec
    llikhood <- temp2$llikhood
    # Store output 
    if (iter>burnin & (iter-burnin) %% thin == 0){
      ind_st  <-  ceiling((iter-burnin)/thin)
      beta_mat[ind_st, ] <- beta_vec
      sigma_vec[ind_st] <- sigma
    } 
    #there are 16 possible models, count the number of times each model appears
    model_mat <- matrix(0, nrow = nrow(beta_mat), ncol = 1)
    for(sim in 1:nrow(beta_mat)){
      model_mat[sim,] <- paste(c((1:length(beta_vec))[beta_mat[sim,1:length(beta_vec)]!=0]), 
                               sep="", collapse = "+")
    }
  }  
})

model_weights <- sort(round(table(model_mat)/nrow(model_mat), 5), TRUE)
allmodelnames <- c("1", paste0("1+", unlist(sapply(1:8, function(n) 
  apply(combn(2:9, n), 2, paste0, collapse="+")))))
weightsrjMCMC <- rep(0, M) # M is what? 256?
names(weightsrjMCMC) <- allmodelnames
# put values into vector of M models:
for (i in 1:length(model_weights)){# loop through weights of rjMCMC
  ind <- which(allmodelnames == names(model_weights[i]) )
  weightsrjMCMC[ind] <- model_weights[i]
}

m.all <- lm(y~., data=train)
library(MuMIn)
options(na.action = "na.fail")
mytab <- dredge(m.all, rank=AIC)
mytab
model.list <- get.models(mytab, subset=NA)
model.list <- model.list[order(as.numeric(rownames(mytab)))]
M <- length(model.list)

# Get model list sequence to order rjMCMC weights
library(stringr)
model.list.modnames <- rep("", M)
for(i in 1:M){
  varnumbers <- str_extract(names(model.list[[i]]$coefficients), "[1-8]")
  varnumbers <- varnumbers[which(!is.na(varnumbers))]
  varnumbers <- sort(as.numeric(varnumbers)+1)
  model.list.modnames[i] <- paste("1",paste(as.character(varnumbers),collapse="+"),sep="+")
  
}

model.list.modnames[1] <- "1"
rightSequence <- rep(0, M)
for(i in 1:M){
  rightSequence[i] <- which(allmodelnames==model.list.modnames[i])
}

truth <- testdat[,1] # we use test as stand-in
preds <- sapply(model.list, predict, newdata=testdat)

(weightsrjMCMC <- weightsrjMCMC[rightSequence])

weightedPredsrjMCMC <- preds %*% weightsrjMCMC
(RMSErjMCMC <- sqrt(mean((weightedPredsrjMCMC - truth)^2)))


#### Logit ####
source("RJMCMCfunctions_logit.R")
dat.logit <- read.csv("logitSimulatedData1April2020.csv")

dat.logit <- dat.logit[,-1]
set.seed(18)
trainsplit <- sample(1:nrow(dat.logit), nrow(dat.logit)/2, replace=FALSE)
train <- dat.logit[trainsplit,]
testdat <- dat.logit[-trainsplit,]

X <- as.matrix(cbind(1, train[,-1])) # include 1 for the intercept
y <- train$y

# starting values for coefficients
beta_vec <- runif(ncol(X), -1, 1)
# starting value for residual noise estimate sigma - not necessary for binom
sigma <- rgamma(1, 1, 1)

llikhood <- sum(dbinom(y, 1, plogis(X%*%beta_vec), log = T))

niter <- 20000 # originally 100000 
burnin <- niter/2
thin <- 20

(niter-burnin)/thin
# store values
beta_mat <- matrix(0, nrow = (niter-burnin)/thin, ncol = length(beta_vec))
sigma_vec <- c()

#choose prior parameters for coefficients
prior_sigma_beta <- 3
prior_pars <- c(0, prior_sigma_beta)

#choose parameters of proposal distributions for updates
prop_sigma_beta <- 3
prop_pars <- c(0, prop_sigma_beta)
# now run rjMCMC:
system.time({
  for (iter in 1:niter){
    temp1 <- updateparam_logit(beta_vec, sigma, llikhood, prior_beta = "norm", prior_beta_par = 
                           prior_pars, prop_beta = "norm", prop_beta_par = prop_pars[2])
    beta_vec <- temp1$beta_vec
    sigma <- temp1$sigma
    llikhood <- temp1$llikhood
    # between models update
    temp2 <- updatemodel_logit(beta_vec, sigma, llikhood,  prior_beta = "norm", prior_beta_par = 
                           prior_pars, prop_beta = "norm", prop_beta_par = prop_pars)
    # read current values for beta, current model and log-likelihood value
    beta_vec <- temp2$beta_vec
    llikhood <- temp2$llikhood
    # Store output 
    if (iter>burnin & (iter-burnin) %% thin == 0){
      ind_st  <-  ceiling((iter-burnin)/thin)
      beta_mat[ind_st, ] <- beta_vec
      sigma_vec[ind_st] <- sigma
    } 
    #there are 16 possible models, count the number of times each model appears
    model_mat <- matrix(0, nrow = nrow(beta_mat), ncol = 1)
    for(sim in 1:nrow(beta_mat)){
      model_mat[sim,] <- paste(c((1:length(beta_vec))[beta_mat[sim,1:length(beta_vec)]!=0]), 
                               sep="", collapse = "+")
    }
  }  
})

model_weights <- sort(round(table(model_mat)/nrow(model_mat), 5), TRUE)
allmodelnames <- c("1", paste0("1+", unlist(sapply(1:8, function(n) 
  apply(combn(2:9, n), 2, paste0, collapse="+")))))
weightsrjMCMC <- rep(0, M) # M is what? 256?
names(weightsrjMCMC) <- allmodelnames
# put values into vector of M models:
for (i in 1:length(model_weights)){# loop through weights of rjMCMC
  ind <- which(allmodelnames == names(model_weights[i]) )
  weightsrjMCMC[ind] <- model_weights[i]
}

m.all <- glm(y~., data=train, family="binomial")
library(MuMIn)
options(na.action = "na.fail")
mytab <- dredge(m.all, rank=AIC)
mytab
model.list <- get.models(mytab, subset=NA)
model.list <- model.list[order(as.numeric(rownames(mytab)))]
M <- length(model.list)

# Get model list sequence to order rjMCMC weights
library(stringr)
model.list.modnames <- rep("", M)
for(i in 1:M){
  varnumbers <- str_extract(names(model.list[[i]]$coefficients), "[1-8]")
  varnumbers <- varnumbers[which(!is.na(varnumbers))]
  varnumbers <- sort(as.numeric(varnumbers)+1)
  model.list.modnames[i] <- paste("1",paste(as.character(varnumbers),collapse="+"),sep="+")
  
}

model.list.modnames[1] <- "1"
rightSequence <- rep(0, M)
for(i in 1:M){
  rightSequence[i] <- which(allmodelnames==model.list.modnames[i])
}

truth <- testdat[,1] # we use test as stand-in
preds <- sapply(model.list, predict, newdata=testdat)

(weightsrjMCMC <- weightsrjMCMC[rightSequence])

weightedPredsrjMCMC <- plogis(preds %*% weightsrjMCMC)
(RMSErjMCMC <- sqrt(mean((weightedPredsrjMCMC - truth)^2)))
