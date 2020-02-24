### Libraries

library(caret)
library(fda)
library(fda.usc)
library(foreach)
library(doParallel)
library(pROC)
library(ggROC)
library(glmnet)
library(dplyr)

### Data
### Establish a clear set of inputs that can be used generally for statistical estimation.
### We can also start working to establish .Rdata files that have the data frame stored correctly.

data.raw <- read.table('lupustherm.txt')
temperatures <- seq(45, 90, 0.1)
thermogram.length <- length(temperatures)
thermograms.raw <- data.raw %>% select(1:thermogram.length)
class.raw <- data.raw %>% select(thermogram.length+1)


### thermogram.fda
### A function that takes thermogram data and produces a 'smoothed' version based on nbasis or lambda-pen
###
### x : a matrix (for thermograms this should be ...), accepta a data.frame also
### argvals : the temperatures
### basis.size : the n.basis used to produce fda representation
### 
### NEEDS: Derivative.

thermogram.fda <- function(x, argvals, nbasis, nderiv=NULL)
{
  if(typeof(x)=='list') x<-as.matrix(x)
  thermogram.length <- length(argvals)
  samples.n <- dim(x)[1]
  
  fdata.raw <- fdata(x, argvals=argvals)
  reduced.basis.spline <- create.bspline.basis(range=range(argvals), nbasis = nbasis)
  tryCatch(Smoothing.Matrix <- S.basis(argvals, basis=reduced.basis.spline), 
           error=function(x) {Smoothing.Matrix<<-diag(thermogram.length)})
  
  ### Tryign to decide how to handle if the 'smoothing' is singular
  ### Need to just return an unsmoothed matrix, but its global and I can't find another way yet.
  
  ### Need to also evaluate using penalized smoothing
  
  fdata.smooth <- fdata.raw
  fdata.smooth$data <- fdata.raw$data%*%Smoothing.Matrix
  
  if(!is.null(nderiv))
  {
    fdata.smooth <- fdata.deriv(fdata.smooth, nderiv=nderiv)
  }
  
  return(fdata.smooth)
}

fdata.working <- thermogram.fda(thermograms.raw, temperatures, 300)
plot(fdata.working)
# rm(Smoothing.Matrix)


### Fold Creator

fold.creator<-function(classes, folds=10, trials=10)
{
  folds.list<-list()
  for(j in  1:trials)
    folds.list[[j]]<-caret::createFolds(factor(classes), k=folds, list=FALSE)
  return(folds.list)
}

### Logistic Regression
### Function for performing LR and cross-validation.
### Returns the resulting predictions of test.x

lr.fit<-function(train.x, test.x, train.classes)
{
  glm.step<-glm(train.classes~., family="binomial", data=as.data.frame(train.x))
  return(glm.step)
}

lasso.fit<-function(train.x, test.x, train.classes)
{
  cv.step<-cv.glmnet(x=train.x, y=train.classes, family="binomial", alpha=1)
  return(cv.step)
}

enet.fit<-function(train.x, test.x, train.classes)
{
  cv.step<-cv.glmnet(x=train.x, y=train.classes, family="binomial", alpha=0.5)
  return(cv.step)
}

ridge.fit<-function(train.x, test.x, train.classes)
{
  cv.step<-cv.glmnet(x=train.x, y=train.classes, family="binomial", alpha=0)
  return(cv.step)
}

# pred.step<-predict(cv.step, test.x, type="response", s=cv.step$lambda.min)
# 
# coef(cv.step, s=cv.step$lambda.min)

# if(model.chosen=='lr.fit') func<-lr.fit
# func<-lr.fit
# ### loop over iterations
# i=1
# 
# ### loop over folds
# fld=1
# 
# train.index <- which(folds.list[[i]]!=fld)
# train.x <- predictor.set[train.index,]
# test.x <- predictor.set[-train.index,]
# train.classes <- classes[train.index]
# n.test <- dim(test.x)[1]
# 
# if(model.chosen=='lr.fit') 
# {
#   func<-lr.fit
#   model.result<-func(train.x, test.x, train.classes)
#   coefficients <- model.result$coefficients[-1]
#   coefficients[which(is.na(output$Coef))]<-0
#   predictions <- predict.glm(model.result, as.data.frame(test.x), type="response")
#   output <- list(Coef=coefficients, Pred.Test=predictions, Test.Info=c(i, fld, n.test))
# }
# 
# plot(temperatures, coefficients)
# 
# 
# which(is.na(output$Coef))
# 
# 
# plot(predictions)
# 
# model.chosen='lr.fit'
# predictor.set <- fdata.working$data
# classes <- class.raw$V452
# folds <- 10
# repeats <- 10
# folds.list <- fold.creator(classes, folds, repeats)

validator.step<-function(predictor.set, classes, model.chosen, folds.list, index, fld)
{
  # index=1
  # fld=1
  train.index <- which(folds.list[[index]]!=fld)
  train.x <- predictor.set[train.index,]
  test.x <- predictor.set[-train.index,]
  train.classes <- classes[train.index]
  n.test <- dim(test.x)[1]
  
  if(model.chosen=='lr.fit') 
  {
    func<-lr.fit
    model.result<-func(train.x, test.x, train.classes)
    coefficients <- model.result$coefficients[-1]
    coefficients[which(is.na(model.result$Coef))]<-0
    predictions <- predict.glm(model.result, as.data.frame(test.x), type="response")
    output <- list(Coef=coefficients, Pred.Test=predictions, Test.Info=c(index, fld, n.test))
  }
  
  if(model.chosen=='lasso.fit') 
  {
    func<-lasso.fit
    model.result<-func(train.x, test.x, train.classes)
    coefficients <- coef(model.result, s=model.result$lambda.min)[-1]
    coefficients[which(is.na(model.result$Coef))]<-0
    predictions <- predict(model.result, test.x, type="response", s=model.result$lambda.min)
    output <- list(Coef=coefficients, Pred.Test=predictions, Test.Info=c(index, fld, n.test), Penalty.Lambda=c(model.result$lambda.min))
  }
  
  if(model.chosen=='enet.fit') 
  {
    func<-enet.fit
    model.result<-func(train.x, test.x, train.classes)
    coefficients <- coef(model.result, s=model.result$lambda.min)[-1]
    coefficients[which(is.na(model.result$Coef))]<-0
    predictions <- predict(model.result, test.x, type="response", s=model.result$lambda.min)
    output <- list(Coef=coefficients, Pred.Test=predictions, Test.Info=c(index, fld, n.test), Penalty.Lambda=c(model.result$lambda.min))
  }
  
  if(model.chosen=='ridge.fit') 
  {
    func<-ridge.fit
    model.result<-func(train.x, test.x, train.classes)
    coefficients <- coef(model.result, s=model.result$lambda.min)[-1]
    coefficients[which(is.na(model.result$Coef))]<-0
    predictions <- predict(model.result, test.x, type="response", s=model.result$lambda.min)
    output <- list(Coef=coefficients, Pred.Test=predictions, Test.Info=c(index, fld, n.test), 
                   Penalty.Lambda=c(model.result$lambda.min))
  }
  
  ### Adaptive LASSO
  ### Adaptive ENET
  ### LDA
  ### Penalized LDA
  ### Boosting
  ### Bagging
  ### Random Forest
  ### SVM?
  
  return(output)
}



# model.chosen='lr.fit'
# predictor.set <- fdata.working$data
# classes <- class.raw$V452
# folds <- 10
# repeats <- 10
# folds.list <- fold.creator(classes, folds, repeats)
# 
# validator.step(predictor.set, classes, model.chosen, folds.list, 1, 1)

# length(folds.list)
# max(folds.list[[1]])

validation.kcv<-function(predictor.set, classes, model.chosen, folds.list, parallel.on=0)
{
  iterations <- length(folds.list)
  folds <- max(folds.list[[1]])
  if(parallel.on)
  {
    cat('No Parallel Yet')
  }
  if(!parallel.on)
  {
    all.trials <- foreach(j=1:iterations, .combine=append, .export=c('validator.step', 'model.chosen', 'folds.list')) %do%
    {
      cat(j)
      trial.temp<-foreach(k=1:folds, .combine=append, .export=c('validator.step', 'model.chosen', 'folds.list')) %do%
      {
        output <- validator.step(predictor.set, classes, model.chosen, folds.list, j, k)
        return(list(output))
      }
      return(list(trial.temp))
    }
  }
  
  return(all.trials)
}




# # all.trials[[1]]
# library(foreach)
# library(doParallel)

model.chosen='lr.fit'
# model.chosen='lasso.fit'
predictor.set <- fdata.working$data
classes <- class.raw$V452
folds <- 4
repeats <- 2
folds.list <- fold.creator(classes, folds, repeats)

test.1 <- validation.kcv(predictor.set, classes, model.chosen, folds.list, parallel.on=0)

test.1

### Always [[repeats]][[folds]]

validation.COEF <- function(validation.kcv.output)
{
  # validation.kcv.output <- test.1
  repeats <- length(validation.kcv.output)
  folds <- length(validation.kcv.output[[1]]) # should always be at least 1
  n.coef <- length(validation.kcv.output[[1]][[1]]$Coef)
  coef.results <- matrix(nrow=repeats*folds, ncol=n.coef)
  for(j in 1:repeats)
  {
    for(k in 1:folds) coef.results[folds*(j-1)+k,]<- validation.kcv.output[[j]][[k]]$Coef
  }
  return(coef.results)
}

coefficient.final <- validation.COEF(test.1)
boxplot(coefficient.final)

validation.penalty <- function(validation.kcv.output)
{
  repeats <- length(validation.kcv.output)
  folds <- length(validation.kcv.output[[1]]) # should always be at least 1
  n.coef <- length(validation.kcv.output[[1]][[1]]$Coef)
  penalty.results <- numeric()
  for(j in 1:repeats)
  {
    for(k in 1:folds) penalty.results[folds*(j-1)+k]<- validation.kcv.output[[j]][[k]]$Penalty.Lambda
  }
  return(penalty.results)
}


penalty.final <- validation.penalty(test.1)
hist(penalty.final)

validation.ROC <- function(validation.kcv.output, classes, folds.list)
{
  require(ROCR)
  require(pROC)
  
  repeats <- length(validation.kcv.output)
  folds <- length(validation.kcv.output[[1]]) # should always be at least 1
  n.coef <- length(validation.kcv.output[[1]][[1]]$Coef)
  ROC.results <- list()
  for(j in 1:repeats)
  {
    for(k in 1:folds) 
      {
        preds.temp <- validation.kcv.output[[j]][[k]]$Pred.Test
        
        train.index <- which(folds.list[[j]]!=k)
        train.x <- predictor.set[train.index,]
        test.x <- predictor.set[-train.index,]
        train.classes <- classes[train.index]
        test.classes<-classes[-train.index]
        
        roc.temp <- roc(as.numeric(test.classes), as.numeric(preds.temp))
        rocr.temp <- prediction(preds.temp, test.classes)
        
        ROC.results[[folds*(j-1)+k]]<-list(roc=roc.temp, rocr=rocr.temp)
      }
  }
  
  return(ROC.results)
  
}

ROC.output <- validation.ROC(test.1, classes, folds.list)

plot(ROC.output[[1]][[1]])
for(j in 1:(repeats*folds)) lines(ROC.output[[j]][[1]])





acc.temp <- performance(ROC.output[[1]][[2]], 'acc')
plot(acc.temp)

for(j in 2:(repeats*folds)) 
  {
  acc.temp <- performance(ROC.output[[j]][[2]], 'acc')
  lines(acc.temp)
  }


roc.temp <- roc(as.numeric(test.classes), as.numeric(preds.temp))

auc(roc.temp)

roc.temp$thresholds

test.classes
preds.temp

library(ggROC)

plot(roc.temp)







library(ROCR)


pred.rocr <- prediction(preds.temp, test.classes)
cutoffs.temp <- pred.rocr@cutoffs
acc.temp <- performance(pred.rocr, 'acc')
sens.temp <- performance(pred.rocr, 'sens')
spec.temp <- performance(pred.rocr, 'spec')

pred.rocr@predictions

acc.temp <- performance(pred.rocr, 'acc')

acc.temp@y.values
acc.temp@x.values

performance(pred.rocr, 'sens', 'spec')

perf <- performance(pred.rocr,"tpr","fpr")
plot(perf,colorize=TRUE)





