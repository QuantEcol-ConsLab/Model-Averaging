# Demoing tools for RJMCMC in NIMBLE
# AEB
# August 27, 2020
# Converse Gardner Lab Meeting
# Model Averaging exercises from Dormann et al 2018

# Fake data story (logit): 
# We tested N = 100 nudibranchs to see if they failed (0) or succeeded (1) in solving a maze. 
# There are 8 covariates associated to varying degrees with nudibranch maze-solving ability

# Fake data story (normal): 
# We measured N = 100 nudibranchs to determine how much 8 covariates could help explain their size (in cm). 

# Reference
# https://r-nimble.org/variable-selection-in-nimble-using-reversible-jump-mcmc

##### LIBRARIES ##### 
library(nimble)
library(coda)

##### DATA ##### 
logitdat <- read.csv("logitSimulatedData1April2020.csv")[,-1]
normaldat <- read.csv("normalSimulatedData1April2020.csv")[,-1]

##### RJMCMC IN NIMBLE ##### 

# there are two approaches
# (1) do you need to estimate inclusion nodes? use indicator variables
  # covered here
# (2) do you know a prior probability for variable inclusion?
  # covered at reference above

##### WITH INDICATOR VARIABLES ##### 

# NORMAL MODEL
normalindicator <- nimbleCode({
  sigma ~ dunif(0, 20)  ## uniform prior per Gelman (2006)
  psi ~ dunif(0,1)    ## prior on inclusion probability
  
  for(i in 1:numVars) {
    z[i] ~ dbern(psi) ## indicator variable for each coefficient
    beta[i] ~ dnorm(0, sd = 100)
    zbeta[i] <- z[i] * beta[i]  ## indicator * beta
  }
  for(i in 1:N) {
    pred.y[i] <- inprod(X[i, 1:numVars], zbeta[1:numVars])
    y[i] ~ dnorm(pred.y[i], sd = sigma)
  }
})

normalindicatorConstants <- list(N = 100, numVars = 8)
normalindicatorInits <- list(sigma = 1, psi = 0.5,
                         beta = rnorm(normalindicatorConstants$numVars),
                         z = sample(0:1, normalindicatorConstants$numVars, 0.5))

normalindicatorData  <- list(y = normaldat[, 1], X = normaldat[, 2:9])
normalindicatorModel <- nimbleModel(code = normalindicator, constants = normalindicatorConstants,
                                inits = normalindicatorInits, data = normalindicatorData)

normalindicatorConf <- configureMCMC(normalindicatorModel)
normalindicatorConf$addMonitors('z')

configureRJ(normalindicatorConf,
            targetNodes = 'beta',
            indicatorNodes = 'z',
            control = list(mean = 0, scale = .2)) # mean and scale for proposal distrib (should coeff go back in model)

normalindicatorConf$printSamplers(c("z[1]", "beta[1]"))

mcmcnormalindicatorRJ <- buildMCMC(normalindicatorConf)

cnormalindicatorModel <- compileNimble(normalindicatorModel)
CMCMCnormalindicatorRJ <- compileNimble(mcmcnormalindicatorRJ, project = normalindicatorModel)

system.time(samplesnormalindicator <- runMCMC(CMCMCnormalindicatorRJ, niter = 1000000, nburnin = 100000))

summary(samplesnormalindicator)
colSums(samplesnormalindicator)
plot(as.mcmc(samplesnormalindicator))

# TODO
# check with robbie to make sure these results make sense
# TRUTH
# beta0 <- qlogis(.6)	# intercept
# Beta <- c(0.6, -0.5, -0.4, 0.3, -0.2, 0.1, -0.05, 0.0)	# covariate effects

# LOGIT MODEL
# TODO
logitindicator <- nimbleCode({
  psi ~ dunif(0,1)    ## prior on inclusion probability
  
  for(i in 1:numVars) {
    z[i] ~ dbern(psi) ## indicator variable for each coefficient
    beta[i] ~ dnorm(0, sd = 100)
    zbeta[i] <- z[i] * beta[i]  ## indicator * beta
  }
  for(i in 1:N) {
    logit(p[i]) <- inprod(X[i, 1:numVars], zbeta[1:numVars])
    y[i] ~ dbern(p[i])
  }
})

logitindicatorConstants <- list(N = 100, numVars = 8)
logitindicatorInits <- list(psi = 0.5,
                             beta = rnorm(logitindicatorConstants$numVars),
                             z = sample(0:1, logitindicatorConstants$numVars, 0.5))

logitindicatorData  <- list(y = logitdat[, 1], X = logitdat[, 2:9])
logitindicatorModel <- nimbleModel(code = logitindicator, constants = logitindicatorConstants,
                                    inits = logitindicatorInits, data = logitindicatorData)

logitindicatorConf <- configureMCMC(logitindicatorModel)
logitindicatorConf$addMonitors('z')

configureRJ(logitindicatorConf,
            targetNodes = 'beta',
            indicatorNodes = 'z',
            control = list(mean = 0, scale = .2)) # mean and scale for proposal distrib (should coeff go back in model)

logitindicatorConf$printSamplers(c("z[1]", "beta[1]"))

mcmclogitindicatorRJ <- buildMCMC(logitindicatorConf)

clogitindicatorModel <- compileNimble(logitindicatorModel)
CMCMClogitindicatorRJ <- compileNimble(mcmclogitindicatorRJ, project = logitindicatorModel)

system.time(sampleslogitindicator <- runMCMC(CMCMClogitindicatorRJ, niter = 1000000, nburnin = 100000))

summary(sampleslogitindicator)
colSums(sampleslogitindicator)
plot(as.mcmc(sampleslogitindicator))

# TODO
# check with robbie to make sure these results make sense
# they do not...what is going on 
