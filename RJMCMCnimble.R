# Demoing tools for RJMCMC in NIMBLE
# AEB
# August 27, 2020
# Converse Gardner Lab Meeting
# Model Averaging exercises from Dormann et al 2018

# Fake data story (normal): 
# We measured N = 100 nudibranchs to determine how much 8 covariates could help explain their size (in cm). 
# TRUTH (normal)
# beta0 <- 5	# intercept
# Beta <- c(0.6, -0.5, -0.4, 0.3, -0.2, 0.1, -0.05, 0.0)	# covariate effects
# sig <- 0.5

# Fake data story (logit): 
# We tested N = 100 nudibranchs to see if they failed (0) or succeeded (1) in solving a maze. 
# There are 8 covariates associated to varying degrees with nudibranch maze-solving ability
# TRUTH (logit)
# beta0 <- qlogis(.6)	# intercept
# Beta <- c(0.6, -0.5, -0.4, 0.3, -0.2, 0.1, -0.05, 0.0)	# covariate effects

# Reference
# https://r-nimble.org/variable-selection-in-nimble-using-reversible-jump-mcmc

##### LIBRARIES ##### 
library(nimble)
library(coda)

##### DATA ##### 
logitdat <- read.csv("logitSimulatedData1April2020.csv")[,-1]
normaldat <- read.csv("normalSimulatedData1April2020.csv")[,-1]

# there are two approaches
# (1) do you need to estimate inclusion nodes? use indicator variables
  # covered here
# (2) do you know a prior probability for variable inclusion?
  # covered at reference above

##### REGRESSION IN NIMBLE #####

# NORMAL MODEL
normal <- nimbleCode({
  # PRIORS
  sigma ~ dunif(0, 20) # prior on variance
  beta0 ~ dunif(-10, 10) # prior on intercept
  for(i in 1:numVars) {
    beta[i] ~ dunif(-1, 1)
  }
  # LIKELIHOOD
  for(i in 1:N) {
    pred.y[i] <- beta0 + inprod(X[i, 1:numVars], beta[1:numVars])
    y[i] ~ dnorm(pred.y[i], sd = sigma)
  }
})

normalConstants <- list(N = 100, numVars = 8)
normalInits <- list(sigma = 1, 
                    beta0 = 0.5,
                    beta = runif(normalConstants$numVars, -1, 1))

normalData  <- list(y = normaldat[, 1], X = normaldat[, 2:9])
normalModel <- nimbleModel(code = normal, constants = normalConstants,
                                    inits = normalInits, data = normalData)

normalConf <- configureMCMC(normalModel)
mcmcnormal <- buildMCMC(normalConf)
cnormalModel <- compileNimble(normalModel)
CMCMCnormal <- compileNimble(mcmcnormal, project = normalModel)
system.time(samplesnormal <- runMCMC(CMCMCnormal, niter = 100000, nburnin = 10000))

summary(samplesnormal)
# These estimates look pretty good!
plot(as.mcmc(samplesnormal))

# LOGIT MODEL
logit <- nimbleCode({
  # PRIORS
  beta0 ~ dunif(-1, 1)
  for(i in 1:numVars) {
    beta[i] ~ dunif(-1,1)
  }
  # LIKELIHOOD
  for(i in 1:N) {
    logit(p[i]) <- beta0 + inprod(X[i, 1:numVars], beta[1:numVars])
    y[i] ~ dbern(p[i])
  }
})

logitConstants <- list(N = 100, numVars = 8)
logitInits <- list(beta0 = 0.5, beta = runif(logitConstants$numVars, -1, 1))

logitData  <- list(y = logitdat[, 1], X = logitdat[, 2:9])
logitModel <- nimbleModel(code = logit, constants = logitConstants,
                                   inits = logitInits, data = logitData)

logitConf <- configureMCMC(logitModel)
mcmclogit <- buildMCMC(logitConf)
clogitModel <- compileNimble(logitModel)
CMCMClogit <- compileNimble(mcmclogit, project = logitModel)

system.time(sampleslogit <- runMCMC(CMCMClogit, niter = 100000, nburnin = 10000))

summary(sampleslogit)
plot(as.mcmc(sampleslogit))

##### MODEL SELECTION WITH INDICATOR VARIABLES IN NIMBLE #####

# NORMAL MODEL
# TODO
# why not psi ~ dbern(0.5), or unique psi for each variable
normalindicator <- nimbleCode({
  sigma ~ dunif(0, 20)  ## uniform prior per Gelman (2006)
  psi ~ dunif(0,1)    ## prior on inclusion probability
  beta0 ~ dunif(-10, 10)
  z0 ~ dbern(psi)
  zbeta0 <- z0 * beta0
  for(i in 1:numVars) {
    z[i] ~ dbern(psi) ## indicator variable for each coefficient
    beta[i] ~ dunif(-1, 1) 
    zbeta[i] <- z[i] * beta[i]  ## indicator * beta
  }
  for(i in 1:N) {
    pred.y[i] <- zbeta0 + inprod(X[i, 1:numVars], zbeta[1:numVars])
    y[i] ~ dnorm(pred.y[i], sd = sigma)
  }
})

normalindicatorConstants <- list(N = 100, numVars = 8)
normalindicatorInits <- list(sigma = 1, psi = 0.5,
                             beta0 = runif(1, -10, 10), 
                             z0 = sample(0:1, 1, 0.5),
                             beta = runif(normalindicatorConstants$numVars, -1, 1),
                             z = sample(0:1, normalindicatorConstants$numVars, 0.5))

normalindicatorData  <- list(y = normaldat[, 1], X = normaldat[, 2:9])
normalindicatorModel <- nimbleModel(code = normalindicator, constants = normalindicatorConstants,
                                    inits = normalindicatorInits, data = normalindicatorData)

normalindicatorConf <- configureMCMC(normalindicatorModel)
normalindicatorConf$addMonitors('z', 'z0')
normalindicatorConf$printSamplers(c("z[1]", "beta[1]"))
mcmcnormalindicator <- buildMCMC(normalindicatorConf)
cnormalindicatorModel <- compileNimble(normalindicatorModel)
CMCMCnormalindicator <- compileNimble(mcmcnormalindicator, project = normalindicatorModel)
system.time(samplesnormalindicator <- runMCMC(CMCMCnormalindicator, niter = 100000, nburnin = 10000))

summary(samplesnormalindicator)
colSums(samplesnormalindicator)/dim(samplesnormalindicator)[1]
# works but does not mix very well
plot(as.mcmc(samplesnormalindicator))

# LOGIT MODEL
logitindicator <- nimbleCode({
  psi ~ dunif(0,1)    ## prior on inclusion probability
  beta0 ~ dunif(-1, 1)
  z0 ~ dbern(psi)
  zbeta0 <- z0 * beta0
  for(i in 1:numVars) {
    z[i] ~ dbern(psi) ## indicator variable for each coefficient
    beta[i] ~ dunif(-1, 1)
    zbeta[i] <- z[i] * beta[i]  ## indicator * beta
  }
  for(i in 1:N) {
    logit(p[i]) <- zbeta0 + inprod(X[i, 1:numVars], zbeta[1:numVars])
    y[i] ~ dbern(p[i])
  }
})

logitindicatorConstants <- list(N = 100, numVars = 8)
logitindicatorInits <- list(psi = 0.5, beta0 = runif(1, -1, 1), z0 = sample(0:1, 1, 0.5),
                            beta = runif(logitindicatorConstants$numVars, -1, 1),
                            z = sample(0:1, logitindicatorConstants$numVars, 0.5))

logitindicatorData  <- list(y = logitdat[, 1], X = logitdat[, 2:9])
logitindicatorModel <- nimbleModel(code = logitindicator, constants = logitindicatorConstants,
                                   inits = logitindicatorInits, data = logitindicatorData)

logitindicatorConf <- configureMCMC(logitindicatorModel)
logitindicatorConf$addMonitors('z', 'z0')
logitindicatorConf$printSamplers(c("z[1]", "beta[1]"))
mcmclogitindicator <- buildMCMC(logitindicatorConf)
clogitindicatorModel <- compileNimble(logitindicatorModel)
CMCMClogitindicator <- compileNimble(mcmclogitindicator, project = logitindicatorModel)
system.time(sampleslogitindicator <- runMCMC(CMCMClogitindicator, niter = 100000, nburnin = 10000))

summary(sampleslogitindicator)
colSums(sampleslogitindicator)/dim(sampleslogitindicator)[1]
# works but doesn't mix very well
plot(as.mcmc(sampleslogitindicator))


##### RJMCMC IN NIMBLE ##### 
# we demonstrate the case using indicator variables, but this can also be done
# without. see reference link above. 

# NORMAL MODEL
normalindicator <- nimbleCode({
  sigma ~ dunif(0, 20)  ## uniform prior per Gelman (2006)
  psi ~ dunif(0,1)    ## prior on inclusion probability
  beta0 ~ dunif(-10, 10)
  z0 ~ dbern(psi)
  zbeta0 <- z0 * beta0
  for(i in 1:numVars) {
    z[i] ~ dbern(psi) ## indicator variable for each coefficient
    beta[i] ~ dunif(-1, 1) 
    zbeta[i] <- z[i] * beta[i]  ## indicator * beta
  }
  for(i in 1:N) {
    pred.y[i] <- zbeta0 + inprod(X[i, 1:numVars], zbeta[1:numVars])
    y[i] ~ dnorm(pred.y[i], sd = sigma)
  }
})

normalindicatorConstants <- list(N = 100, numVars = 8)
normalindicatorInits <- list(sigma = 1, psi = 0.5,
                         beta0 = runif(1, -10, 10), 
                         z0 = sample(0:1, 1, 0.5),
                         beta = runif(normalindicatorConstants$numVars, -1, 1),
                         z = sample(0:1, normalindicatorConstants$numVars, 0.5))

normalindicatorData  <- list(y = normaldat[, 1], X = normaldat[, 2:9])
normalindicatorModel <- nimbleModel(code = normalindicator, constants = normalindicatorConstants,
                                inits = normalindicatorInits, data = normalindicatorData)

normalindicatorConf <- configureMCMC(normalindicatorModel)
normalindicatorConf$addMonitors('z', 'z0')

# NOTE z0 is not one of the indicator nodes
configureRJ(normalindicatorConf,
            targetNodes = 'beta',
            indicatorNodes = 'z',
            control = list(mean = 0, scale = .2)) # mean and scale for proposal distrib (should coeff go back in model)

normalindicatorConf$printSamplers(c("z[1]", "beta[1]"))

mcmcnormalindicatorRJ <- buildMCMC(normalindicatorConf)

cnormalindicatorModel <- compileNimble(normalindicatorModel)
CMCMCnormalindicatorRJ <- compileNimble(mcmcnormalindicatorRJ, project = normalindicatorModel)

system.time(samplesnormalindicator <- runMCMC(CMCMCnormalindicatorRJ, niter = 100000, nburnin = 10000))

summary(samplesnormalindicator)
colSums(samplesnormalindicator)/dim(samplesnormalindicator)[1]
# Looks good!
plot(as.mcmc(samplesnormalindicator))

# LOGIT MODEL
logitindicator <- nimbleCode({
  psi ~ dunif(0,1)    ## prior on inclusion probability
  beta0 ~ dunif(-1, 1)
  z0 ~ dbern(psi)
  zbeta0 <- z0 * beta0
  for(i in 1:numVars) {
    z[i] ~ dbern(psi) ## indicator variable for each coefficient
    beta[i] ~ dunif(-1, 1)
    zbeta[i] <- z[i] * beta[i]  ## indicator * beta
  }
  for(i in 1:N) {
    logit(p[i]) <- zbeta0 + inprod(X[i, 1:numVars], zbeta[1:numVars])
    y[i] ~ dbern(p[i])
  }
})

logitindicatorConstants <- list(N = 100, numVars = 8)
logitindicatorInits <- list(psi = 0.5, beta0 = runif(1, -1, 1), z0 = sample(0:1, 1, 0.5),
                             beta = runif(logitindicatorConstants$numVars, -1, 1),
                             z = sample(0:1, logitindicatorConstants$numVars, 0.5))

logitindicatorData  <- list(y = logitdat[, 1], X = logitdat[, 2:9])
logitindicatorModel <- nimbleModel(code = logitindicator, constants = logitindicatorConstants,
                                    inits = logitindicatorInits, data = logitindicatorData)

logitindicatorConf <- configureMCMC(logitindicatorModel)
logitindicatorConf$addMonitors('z', 'z0')

# NOTE z0 is not one of the indicator nodes
# TODO
# perhaps play with configureRJ control settings
configureRJ(logitindicatorConf,
            targetNodes = 'beta',
            indicatorNodes = 'z',
            control = list(mean = 0, scale = .2)) # mean and scale for proposal distrib (should coeff go back in model)

logitindicatorConf$printSamplers(c("z[1]", "beta[1]"))

mcmclogitindicatorRJ <- buildMCMC(logitindicatorConf)

clogitindicatorModel <- compileNimble(logitindicatorModel)
CMCMClogitindicatorRJ <- compileNimble(mcmclogitindicatorRJ, project = logitindicatorModel)

system.time(sampleslogitindicator <- runMCMC(CMCMClogitindicatorRJ, niter = 100000, nburnin = 10000))

summary(sampleslogitindicator)
colSums(sampleslogitindicator)/dim(sampleslogitindicator)[1]
# These don't look so good...any ideas?
plot(as.mcmc(sampleslogitindicator))

# TODO
# how to compute model weights from this

