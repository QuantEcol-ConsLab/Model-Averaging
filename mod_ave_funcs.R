# this script contians functions for model averaging and comparison of methods

#function. Takes data, data type, vector of model averaging methods, and covariate columns to use. Returns a list of RMSE, predictions,weights, and run time for each method.

all_methods_MA<-function(data_type,dat,methods,covs){
  n_methods<-length(methods)      # number of methods     
  methods_vec<-numeric(n_methods) # vector that will go in list below
  names(methods_vec)=methods
  
  n_test_data<-nrow(dat)/2        #number of data points in the test data
  n_mods<-2^length(covs)         # number of models being averaged
  
  
  #list to hold outputs from model averaging methods: weights, predictions, rmse, and run time
  compiled_outputs<-list(
    rmse=methods_vec,
    weights=matrix(NA,nrow=n_mods,ncol=n_methods,dimnames = list(1:n_mods,methods)),
    ma_preds=matrix(NA,nrow=n_test_data,ncol=n_methods,dimnames = list(1:n_test_data,methods)),
    time=methods_vec
  )
  
  for(i in 1:n_methods){
    ma_out<-model_average_function(method= methods[i],data_in = dat,covs=covs,data_type = data_type)
    
    compiled_outputs$rmse[i] <- ma_out$rmse
    compiled_outputs$weights[,i] <- ma_out$weights
    compiled_outputs$ma_preds[,i] <- ma_out$wtpreds
    compiled_outputs$time[i] <- ma_out$run_time
  } 
  
  return(compiled_outputs)  
  
}





# The functions take as inputs training and testing data, and return 1) model weights, 2) model predictions to test data, 3) RMSE for test data

# methods include: CP,  minimum variance, bottstrapping, AIC, BIC, leave one out cross validation, jackknife, and stacking

model_average_function<-function(
  method,   #model averaing method
  data_in,  #data
  covs,     # vector of columns of covariates to include
  data_type # for flm function "family" argument. e.g., "binomial" or "gaussian"
){

 #training data
  train<-data_in[((nrow(data_in)/2)+1):(nrow(data_in)),c(1,covs)]
 
  #testing response
  test_y<-data_in[1:(nrow(data_in)/2),1]
  
  #testing covariate values
  test_covs<-data_in[1:(nrow(data_in)/2),covs]
  
  start_time<-Sys.time()
  
  ## fit global model
  fm1 <- glm(y~., data=train, na.action = "na.fail", family=data_type,x=TRUE)
  
  #get all possible models
  dfm1<-dredge(fm1)
  
  #get all of the model output
  fm1.mods.all<-get.models(dfm1, subset=TRUE)
  
#calculate weights  

  #information criterion
  if(method%in%c("AIC","BIC","Cp","loo")){
    #Get the information criterion value
    if(method=="loo"){
    mod_criterion_value<-sapply(fm1.mods.all, method,type="rmse")      
    }else{
    mod_criterion_value<-sapply(fm1.mods.all, method)
    }
    # Calculate the weights for each model based on Cp
    weights<-exp(-0.5*(mod_criterion_value-min(mod_criterion_value)))/sum(exp(-0.5*(mod_criterion_value-min(mod_criterion_value))))
  }else
      #Bates-granger/minimal variance weights
      if(method=="minimum_variance"){
        weights<-(BGWeights(fm1.mods.all, data=train))#some are negative
      }else
        
        if(method=="bootstrap"){
          set.seed(2020)
          weights<-(bootWeights(fm1.mods.all, data=train,R=1000))
        }else
          
          if(method=="stacking"){
            weights<-(stackingWeights(fm1.mods.all, data=train, R=1000))[1,]
          }else
            if(method=="jackknife"){
              weights<-(jackknifeWeights(fm1.mods.all, data=train))
            }else
              if(method=="cos_squared"){
                weights<-(cos2Weights(fm1.mods.all, data=train))
              }
  


  
  
#predict for all models
preds_all_mods <-sapply(fm1.mods.all, predict, newdata = test_covs)
#weigh the models by the BG/MV weights
wtpreds<-preds_all_mods%*%weights
#RMSE 
rmse<-sqrt(mean((wtpreds-test_y )^2))

#time it took to run
total_time<-Sys.time()-start_time

#return

return(list(
  weights=weights,
  wtpreds=wtpreds,
  rmse=rmse,
  run_time=total_time,
  method=method,
  true_test=test_y
))
}

