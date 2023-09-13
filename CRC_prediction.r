library(ggplot2)
library(tidyr)
library(dplyr)
library(FactoMineR)
library(data.table)
library(caret)
library(reshape2)
library(vegan)
library(patchwork)
library(ecodist)
library(MASS)
library(glmnet)
library(pROC)


# 70 -30 train test ml function (logistic regression) 
train_test = function(data, samples){
  
  data_ml = merge(data, samples, by = "sample")
  data_ml$Group = as.factor(data_ml$Group)
  data_ml$sample = NULL
  data_ml$population = NULL

  aucs = list()
  accs = list()
  n = 0
  for (j in 1:20) {
    
    
      set.seed(n)
      
      # get train and test data
      validation_index = createDataPartition(data_ml$Group, p=0.30, list=FALSE)
      train <- data_ml[-validation_index,]
      test <- data_ml[validation_index,]
      
      yt <- as.factor(train$Group)
      train$Group = NULL
      xt <- train
      
      ytest <- as.factor(test$Group)
      test$Group = NULL
      xtest <- test
      
      # fit model
      cvglmfit  <- cv.glmnet(as.matrix(xt), yt,family = "binomial")
      s0 = cvglmfit$lambda.min
      fit = glmnet(as.matrix(xt), yt,family = "binomial")
      pred = predict(fit,newx = as.matrix(xtest),type = "class", s = s0)
      
      # get AUC
      auc = assess.glmnet(fit, newx = as.matrix(xtest), newy = ytest,s = s0)
      
      # get accuracy
      cm = confusionMatrix(table(pred,ytest))  
      acc = cm$overall[1]
      
      # store data
      n = n + 1
      aucs[n] = auc$auc
      accs[n] = acc
      
      
  }
  r = as.data.frame(cbind(aucs,accs))
  names(r) = c("auc","acc")
  r$auc = as.numeric(r$auc)
  r$acc = as.numeric(r$acc)
  return(r)
}

# train test 70 - 30%, with random forest
train_test_rf = function(data,samples){

  data_ml = merge(data, samples, by = "sample")
  data_ml$Group = as.factor(data_ml$Group)
  data_ml$sample = NULL
  data_ml$population = NULL
  aucs = list()
  accs = list()
  n = 0
  for (j in 1:20) {
  
    set.seed(n)
    
    # random sample 1000 kos / novel fams 
    data_ml_sub = data_ml[,sample(ncol(data_ml), 1000)]
    data_ml_sub$Group = data_ml$Group
    data_ml = data_ml_sub

     # get train and test data
    validation_index = createDataPartition(data_ml$Group, p=0.30, list=FALSE)
    train <- data_ml[validation_index,]
    test <- data_ml[-validation_index,]
    
    xt <- train
    xtest <- test
    
    # train model
    control <- trainControl(method="cv", number=10)
    metric <- "Accuracy"
    fit.rf <- train(Group~., data=xt, method="rf", metric=metric, trControl=control)
    
    # calculate accuracy
    predictions <- predict(fit.rf, xtest)
    results = confusionMatrix(predictions, as.factor(xtest$Group))
    acc = results$overall[1]
    
    # calculate AUC
    predictions <- predict(fit.rf, xtest,  type="prob")
    auc = roc(xtest$Group, predictions$CTR)
    
    # combine in table
    n = n + 1
    aucs[n] = auc$auc
    accs[n] = acc
     
  }
  r = as.data.frame(cbind(aucs,accs))
  names(r) = c("auc","acc")
  r$auc = as.numeric(r$auc)
  r$acc = as.numeric(r$acc)
  return(r)
}

####
# load data
####

# nfams data
data = read.csv("abs_per_fam_all.f_fams.h.tab",sep = ',',header = T)
names = data[,1]
data = data[,-1]
colnames= names(data)
data = data.frame(transpose(data))
names(data) = names
data$sample = colnames
nfam_data = data

# kos data 
data = read.csv("kos_abs.tab",header = T)
names = data[,1]
data = data[,-1]
colnames= names(data)
data = data.frame(transpose(data))
names(data) = names
data$sample = colnames
ko_data = data

# kos + nfams data
data_b = merge(ko_data,nfam_data,by = "sample")

# load metadata
samples = read.table("metadata.tsv.f",header = T)
samples = samples[,c(1,6,11)]
names(samples) = c("sample","Group","population")


#####
# ml logistic regression 70% train - 30% test
#####


# run model
ttko = train_test(ko_data,samples)
ttko$origin = "KOs"
ttnfam = train_test(nfam_data,samples)
ttnfam$origin = "nfams"
ttboth = train_test(data_b,samples)
ttboth$origin = "KOs + nfams"
tt = rbind(ttko,ttnfam,ttboth)

write.csv(tt,"prediction_results/linear_regression.30Train_70Test.lambda.min.tab")

#####
# ml random forest 30% train - 70% test, 1000 random nfams / kos / combination
#####


# run model
ttko = train_test_rf(ko_data,samples)
ttko$origin = "KOs"
ttnfam = train_test_rf(nfam_data,samples)
ttnfam$origin = "nfams"
ttboth = train_test_rf(data_b,samples)
ttboth$origin = "KOs + nfams"
tt = rbind(ttko,ttnfam,ttboth)

write.csv(tt,"prediction_results/random_forest.30Train_70Test.tab")


