#conducts model averaging using 9 different methods on two data sets and returns RMSE, MA predictions, and weights for each MA method

library(MuMIn)
source("mod_ave_funcs.R")

#### Normal data ####
#load data
dat_normal<- read.csv("normalSimulatedData1April2020.csv")[,-1]

#vector of methods to use in function below
methods<-c("AIC","BIC","Cp","loo","minimum_variance","bootstrap","stacking","jackknife","cos_squared")

#conduct model averaging (this takes a few minutes)
all_MA_normal<-all_methods_MA(data_type = "gaussian", dat = dat_normal, methods =  methods[1:9], covs =  2:5)

#### Logit data ####
#load data
dat.logit <-  read.csv("logitSimulatedData1April2020.csv")[,-1]
dat.logit$y<-as.numeric(dat.logit$y) #loo function doesn't take integer response

#conduct model averaging (this takes a few minutes)
all_MA_logit<-all_methods_MA(data_type = "binomial", dat = dat.logit, methods =  methods[1:9], covs =  2:5)


#plot results... this is just a quick start. feel free to expand or beautify

##  Normal  ##
#RMSE
all_MA_normal$rmse

par(mfrow=c(1,2))
#RMSE barplot
barplot(all_MA_normal$rmse,las=2)

#Weights boxplot
boxplot(all_MA_normal$weights,las=2)

#true vs predicted
par(mfrow=c(3,3),mar=rep(0,4),oma=c(4,4,2,2))
for ( i in 1:9){
plot(dat_normal$y[51:100],all_MA_normal$ma_preds[,i],pch=19,yaxt="n",xaxt="n")
text(x=1.1*min(dat_normal$y[51:100]),y=.99*max(all_MA_normal$ma_preds[,i]),labels =names(all_MA_normal$rmse[i]),pos=4 )
abline(0,1,lwd=2)
}
mtext("true value",1,1,outer = TRUE)
mtext("MA prediction",2,1,outer = TRUE)

#boxplot predictions
boxplot(all_MA_normal$ma_preds,las=2)



#RMSE logit
all_MA_logit$rmse
barplot(all_MA_logit$rmse,las=2)
