#Statistika Deskriptif
library(readxl)
IPM<-read_excel("IPM Skripsi.xlsx")
IPM <- data.frame(IPM)
colnames(IPM) <- c("KAB_KOTA","LONG","LAT","Y","X1","X2","X3","X4","X5","X6","X7")
summary(IPM)
var(IPM$Y);var(IPM$X1);var(IPM$X2);var(IPM$X3);var(IPM$X4)
var(IPM$X5);var(IPM$X6);var(IPM$X7)

regresi<-lm(formula=Y~X1+X2+X3+X4+X5+X6+X7,data=IPM)
#Heterogenitas (Breusch-Pagan)
library(lmtest)
bptest(regresi)

# Multikolinearitas
library(car)
vif(regresi)

#' GWRLASSO: a hybrid model that uses the LASSO model for important variable selection and GWR model with exponential kernel for prediction at an unknown location based on the selected variables.
library(qpdf)
library(numbers)
library(glmnet)
library(Matrix)
library(spdep)

min_max_scale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

set.seed(123)
GWRLASSO_exponential<- function(data_sp,bw,split_value,exponential_kernel,nfolds) {
  
  # Step1: Generation of simulated spatial population with spatial coordinates
  data_sp<-as.data.frame(data_sp)
 
  # Step2: Split data into training and testing sets
  train_idx <- sample(nrow(data_sp), split_value * nrow(data_sp))
  train_data <- data_sp[train_idx,]
  nrow(train_data)
  test_data <- data_sp[-train_idx, ]
  nrow(test_data)
  x_train <-subset(train_data[,-c(1,ncol(train_data)-1,ncol(train_data))])
  y_train <-train_data[,1]
  x_test <- subset(test_data,select = -c(1,ncol(test_data)-1,ncol(test_data)))
  y_test <- test_data[,1]
  
  # Step3: Important variable selection using LASSO
  # fit LASSO model for variable selection
  
  train_data_LASSO<-subset(train_data,select = -c(ncol(train_data)-1,ncol(train_data)))
  
  x_lasso <- model.matrix(y ~ ., train_data_LASSO)
  y_lasso <- train_data_LASSO[, "y"]
  y_lasso <- y_lasso[1:nrow(x_lasso)]
  
  ### select important variables based on LASSO model
  cvfit <- cv.glmnet(x_lasso, y_lasso, alpha = 1, standardize = TRUE)
  coef_path <- predict(cvfit, type="coefficients", s=cvfit$lambda.min)
  selected_vars <- row.names(coef_path)[which(coef_path != 0)[-1]] # exclude intercept
  
  # Extracting important variables
  imp_var_lasso <- x_train[, selected_vars]
  
  # extract optimal lambda value
  opt_lambda <- cvfit$lambda.min
  
  # Step5 : Fitting of the GWR model on training data
  
  coords<-cbind(data_sp[,ncol(data_sp)-1],data_sp[,ncol(data_sp)])
  dists<- as.matrix(dist(coords))
  
  # Define the Exponential kernel function
  exponential_kernel <- function(dists, bw) {
    exp(-(dists/bw))
  }
  weights <- exponential_kernel(dists, bw)
  dim(weights)
  Xt <- as.matrix(cbind(rep(1, nrow(x_train)),imp_var_lasso))
  dim(Xt)
  y_tr<- matrix(y_train)
  dim(y_tr)
  W.train<-list()
  for (i in 1:ncol(weights)){
    t <- diag(weights[,i],nrow=nrow(x_train), ncol=nrow(x_train))
    W.train[[i]]<-t
  }
  W.train
  
  Beta.train<-list()
  for (i in 1:length(W.train)){
    lm<- solve((t(Xt)%*%W.train[[i]]%*%Xt))%*%t(Xt)%*%W.train[[i]]%*%y_tr
    Beta.train[[i]]<-lm
  }
  Beta.train
  
  X.tr_row<-list()
  for(i in 1:nrow(Xt)){
    
    X.tr_row[[i]]<-(Xt[i,])
  }
  
  
  y_hat.train<-mapply("%*%", X.tr_row,Beta.train)
  y_hat.train
  
  
  # Step6 : Make predictions at the test locations
  
  lasso_imp_X_test<-x_test[, selected_vars]
  
  coords_test<-cbind(test_data[,ncol(test_data)-1],test_data[,ncol(train_data)])
  dists_test <- as.matrix(dist(coords_test))
  weights_test <- exponential_kernel(dists_test, bw)
  dim(weights_test)
  
  Xtest <- as.matrix(cbind(rep(1, nrow(x_test)),lasso_imp_X_test))
  dim(Xtest)
  ytest<- matrix(y_test)
  dim(ytest)
  
  W.test<-list()
  for (i in 1:ncol(weights_test)){
    test <- diag(weights_test[,i],nrow=nrow(x_test), ncol=nrow(x_test))
    W.test[[i]]<-test
  }
  W.test
  
  Beta.test<-list()
  for (i in 1:length(W.test)){
    lm_test<- solve((t(Xtest)%*%W.test[[i]]%*%Xtest))%*%t(Xtest)%*%W.test[[i]]%*%ytest
    Beta.test[[i]]<-lm_test
  }
  Beta.test
  
  X.test_row<-list()
  for(i in 1:nrow(Xtest)){
    
    X.test_row[[i]]<-(Xtest[i,])
  }
  
  
  y_hat.test<-mapply("%*%", X.test_row,Beta.test)
  y_hat.test
  
  # Step7 : Create summary output
  # Compute accuracy measures
  
  rmse_test <- sqrt(mean((y_hat.test - y_test)^2)) / mean(y_test)
  r2_test <- cor(y_hat.test, y_test)^2
  
  # Return results
  results_GWRLASSO <- list(Important_vars = selected_vars,
                           Optimum_lamda = opt_lambda,
                           GWR_y_pred_train = y_hat.train,
                           GWR_y_pred_test = y_hat.test,
                           R_square = r2_test,
                           rmse=rmse_test,
                           Beta.train = Beta.train,
                           Beta.test = Beta.test,
                           Weight.train = W.train,
                           Weight.test = W.test)
  
  return(results_GWRLASSO)
}

library(readxl)
Dataset <- read_excel("IPM Skripsi.xlsx")
Latitude<-Dataset$longitude
Longitude<-Dataset$longitude
y<-Dataset$Y
library(dplyr)
x <- select(Dataset,-Y,-longitude,-latitude,-`Kabupaten/Kota`)
data_df<-data.frame(y,x,Latitude,Longitude)
data_sp<-data.frame(y,x,Latitude,Longitude)
coordinates(data_sp) <- c("Longitude","Latitude")
#data_sp <- as(data_sp,"Spatial")
# Tentukan kernel bandwidth
library(GWmodel)
bw <- bw.gwr(y~X1+X2+X3+X4+X5+X6+X7,
             approach = "CV",
             adaptive = F, 
             data=data_sp)
bw

GWRLASSO_exp<-GWRLASSO_exponential(data_df,0.8,bw=bw,exponential_kernel,10)
GWRLASSO_exp$Beta

n <- length(IPM$`Kabupaten/Kota`)

GWRLASSO_exp$Beta

Beta <- data.frame(GWRLASSO_exp$Beta)
# yhat
y_hat <- mapply("%*%", x,GWRLASSO_exp$Beta)
# residual
residual <- IPM$Y - GWRLASSO_exp