#conducts model averaging using 9 different methods on two data sets and returns RMSE, MA predictions, and weights for each MA method

library(MuMIn)
source("mod_ave_funcs.R")

#### Normal data ####
dat_normal<- read.csv("normalSimulatedData1April2020.csv")[,-1]

methods<-c("AIC","BIC","Cp","loo","minimum_variance","bootstrap","stacking","jackknife","cos_squared")

all_MA_norman<-all_methods_MA(data_type = "gaussian", dat = dat_normal, methods =  methods[1:9], covs =  2:5)

#### Logit data ####
dat.logit <-  read.csv("logitSimulatedData1April2020.csv")[,-1]
dat.logit$y<-as.numeric(dat.logit$y) #loo function doesn't tale integer response

all_MA_logit<-all_methods_MA(data_type = "binomial", dat = dat.logit, methods =  methods[1:9], covs =  2:5)


#plot results...